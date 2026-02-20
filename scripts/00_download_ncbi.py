import os
import requests
import time
import sys
import re
from requests.exceptions import ChunkedEncodingError, ConnectionError, ReadTimeout

# ================= CONFIGURATION =================
BOLD_BASE_URL = "https://v4.boldsystems.org/index.php/API_Public/sequence"
NCBI_API_KEY = "f8ffdbb3a603606188fe140fb934f68d8b08"  # User API Key
EMAIL = "roshguadiana@gmail.com" 

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/91.0.4472.124 Safari/537.36"
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_REF_DIR = os.path.join(BASE_DIR, "../data/fasta_ref")
taxid_cache = {}

PRESETS = {
    "FISH_UK": {
        "taxa": ["Actinopterygii", "Chondrichthyes", "Petromyzontida", "Myxini"],
        "geo": "United Kingdom"
    },
    "PLANTS_UK": {
        "taxa": ["Magnoliopsida", "Liliopsida", "Pinopsida", "Polypodiopsida"],
        "geo": "United Kingdom"
    }
}

VALID_SUBTREES = {
    # --- PECES ---
    "Actinopterygii": "7898",   # Peces óseos (atunes, bacalaos, etc.)
    "Chondrichthyes": "7777",   # Tiburones, rayas, quimeras
    "Petromyzontida": "7745",   # Lampreas
    "Myxini": "7762",           # Mixinos
    
    # --- PLANTAS (De tus presets) ---
    "Magnoliopsida": "3398",    # Dicotiledóneas
    "Liliopsida": "4081",       # Monocotiledóneas
    "Pinopsida": "3313",        # Coníferas
    "Polypodiopsida": "241806"  # Helechos
}
# =================================================================
# =================================================

def check_api_key_status():
    if not NCBI_API_KEY:
        print("\n⚠️  WARNING: NCBI API Key not detected (Slow mode).")
    else:
        print("\n⚡ TURBO MODE: NCBI API Key active.")

def get_ncbi_taxid(name, expected_lineage):
    """
    Queries NCBI Taxonomy using strict Subtree mathematical graph matching.
    """
    name = name.strip()
    name = re.sub(r'[|_]', ' ', name)
    if not name: return None
    
    cache_key = f"{name}_{expected_lineage}"
    if cache_key in taxid_cache: 
        return taxid_cache[cache_key]

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    # 1. Buscamos el nodo padre en nuestro diccionario
    parent_taxid = VALID_SUBTREES.get(expected_lineage)
    
    if parent_taxid:
        # FIREWALL MATEMÁTICO: Filtramos por descendencia de nodo (Infalible)
        strict_query = f"{name}[Organism] AND txid{parent_taxid}[Subtree]"
    else:
        # Fallback de texto por si introduces un taxón nuevo que no está en el diccionario
        strict_query = f"{name}[Organism] AND {expected_lineage}[Lineage]"
        
    params = {'db': 'taxonomy', 'term': strict_query, 'retmode': 'json', 'email': EMAIL}
    if NCBI_API_KEY: params['api_key'] = NCBI_API_KEY
    
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=10)
        if r.status_code == 200:
            data = r.json()
            id_list = data.get('esearchresult', {}).get('idlist', [])
            taxid = id_list[0] if id_list else None
            
            taxid_cache[cache_key] = taxid
            time.sleep(0.1 if NCBI_API_KEY else 0.34)
            return taxid
    except:
        pass
    return None

def resolve_taxid_from_fasta_header(header_line, expected_lineage):
    """
    Extracts the organism name from the BOLD header and passes it to the firewall.
    """
    clean_line = header_line.lstrip('>').strip()
    parts = clean_line.split('|')
    if len(parts) < 2: return None, None, None

    raw_name = parts[1].strip()
    
    # First barrier: Try matching the full Species name
    tid = get_ncbi_taxid(raw_name, expected_lineage)
    if tid: return tid, raw_name, "Species"

    genus = raw_name.split(' ')[0]
    if len(genus) > 2 and genus != raw_name:
        # Second barrier: Try matching the Genus if Species fails
        tid = get_ncbi_taxid(genus, expected_lineage)
        if tid: return tid, genus, "Genus"

    return None, None, None

def download_with_retries(taxon, geo, max_retries=9):
    """
    Downloads data from BOLD and processes it on the fly, passing the taxon as the filter.
    """
    params = {'taxon': taxon, 'format': 'fasta'}
    if geo: params['geo'] = geo

    for attempt in range(1, max_retries + 1):
        try:
            print(f"   ⬇️  Attempt {attempt}/{max_retries}...")
            
            r = requests.get(BOLD_BASE_URL, params=params, headers=HEADERS, stream=True, timeout=300)
            
            if r.status_code != 200:
                print(f"   ❌ HTTP Error {r.status_code}. Retrying...")
                time.sleep(5)
                continue

            # Pass the 'taxon' as the expected lineage to the processor
            return process_stream_with_progress(r, taxon)

        except (ChunkedEncodingError, ConnectionError, ReadTimeout) as e:
            print(f"\n   ⚠️  Network issue: {e}")
            wait_time = attempt * 5 
            if attempt < max_retries:
                print(f"   ⏳ Waiting {wait_time}s before retrying...")
                time.sleep(wait_time)
            else:
                print("   ❌ Retries exhausted.")
                return None
    return None

def process_stream_with_progress(response, expected_lineage):
    """
    Processes the BOLD stream, filtering out contaminants using the taxonomic firewall.
    """
    valid_sequences = []
    stats = {'species': 0, 'genus': 0, 'blocked_contaminants': 0}
    
    current_header = None
    current_seq = []
    downloaded_bytes = 0
    start_time = time.time()
    
    print("   📊 Starting download and strict filtering...")

    for line in response.iter_lines():
        if not line: continue
        
        downloaded_bytes += len(line) + 1 
        
        # Terminal UI update
        if len(valid_sequences) % 50 == 0:
            mb = downloaded_bytes / (1024 * 1024)
            elapsed = time.time() - start_time
            speed = mb / elapsed if elapsed > 0 else 0
            
            msg = f"\r   🚀 {mb:.2f} MB received | {speed:.2f} MB/s | Valid Seqs: {len(valid_sequences)} | Blocked: {stats['blocked_contaminants']}"
            sys.stdout.write(msg)
            sys.stdout.flush()

        line_str = line.decode('utf-8').strip()
        
        if line_str.startswith('>'):
            if current_header and current_seq:
                taxid, name, rank = resolve_taxid_from_fasta_header(current_header, expected_lineage)
                
                # THE FIREWALL ACTS HERE
                if taxid:
                    orig_id = current_header.split('|')[0].lstrip('>')
                    new_header = f">{orig_id}|kraken:taxid|{taxid} {name}"
                    valid_sequences.append(f"{new_header}\n{''.join(current_seq)}")
                    if rank == "Species": stats['species'] += 1
                    else: stats['genus'] += 1
                else:
                    # Discarded: eDNA contaminants (Humans, Bacteria, Insects, etc.)
                    stats['blocked_contaminants'] += 1
            
            current_header = line_str
            current_seq = []
        else:
            current_seq.append(line_str)

        # Process the final block in the stream
    if current_header and current_seq:
        taxid, name, rank = resolve_taxid_from_fasta_header(current_header, expected_lineage)
        if taxid:
            orig_id = current_header.split('|')[0].lstrip('>')
            new_header = f">{orig_id}|kraken:taxid|{taxid} {name}"
            valid_sequences.append(f"{new_header}\n{''.join(current_seq)}")
            if rank == "Species": stats['species'] += 1
            else: stats['genus'] += 1
        else:
            stats['blocked_contaminants'] += 1

    sys.stdout.write("\n") # Clean terminal line
    
    if valid_sequences:
        return valid_sequences, stats
    return None, None

def main():
    print("=== BOLD DOWNLOADER: SECURE MODE ===")
    check_api_key_status()

    db_name = input("1. DB Name (e.g., FISH): ").strip().upper()
    if not db_name: sys.exit(1)
    
    target_dir = os.path.join(DATA_REF_DIR, db_name)
    os.makedirs(target_dir, exist_ok=True)

    mode = input("\n[1] Manual\n[2] Preset: UK Fish\n[3] Preset: UK Plants\nSel: ").strip() or "1"
    
    if mode == "2": taxa, geo = PRESETS["FISH_UK"]["taxa"], PRESETS["FISH_UK"]["geo"]
    elif mode == "3": taxa, geo = PRESETS["PLANTS_UK"]["taxa"], PRESETS["PLANTS_UK"]["geo"]
    else:
        taxa = [t.strip() for t in input("Taxa (comma separated): ").split(",")]
        geo = input("Geo: ").strip()
    
    print(f"\n📂 Directory: {target_dir}")
    
    total_files = 0
    for taxon in taxa:
        filename = f"{taxon}_{'UK' if 'United Kingdom' in str(geo) else 'Global'}_Kraken.fasta"
        filepath = os.path.join(target_dir, filename)

        if os.path.exists(filepath) and os.path.getsize(filepath) > 1000:
            print(f"\n✅ {taxon}: Already exists ({os.path.getsize(filepath)/1024/1024:.2f} MB). Skipping.")
            total_files += 1
            continue

        print(f"\n🧬 Processing: {taxon} (with strict Lineage filtering)")
        result = download_with_retries(taxon, geo)
        
        if result:
            seqs, stats = result
            if seqs:
                with open(filepath, 'w') as f:
                    f.write("\n".join(seqs))
                print(f"   💾 Saved: {filename}")
                print(f"      📊 Final: {stats['species']} Species | {stats['genus']} Genera | {stats['blocked_contaminants']} Blocked")
                total_files += 1
        else:
            print(f"   ❌ Failed to download or process {taxon}.")

    if total_files > 0:
        print(f"\n🎉 DONE! Next step: ./01_build_db.sh {db_name}")

if __name__ == "__main__":
    main()