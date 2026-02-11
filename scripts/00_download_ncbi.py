import os
import requests
import time
import sys
import re
from requests.exceptions import ChunkedEncodingError, ConnectionError, ReadTimeout

# ================= CONFIGURACIÓN =================
BOLD_BASE_URL = "https://v4.boldsystems.org/index.php/API_Public/sequence"
NCBI_API_KEY = "f8ffdbb3a603606188fe140fb934f68d8b08"  # e.g. "a1b2c3d4e5..."
EMAIL = "roshguadiana@gmail.com" 
#
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
# =================================================

def check_api_key_status():
    if not NCBI_API_KEY:
        print("\n⚠️  AVISO: No se detectó NCBI API Key (Modo lento).")
    else:
        print("\n⚡ TURBO MODE: NCBI API Key activa.")

def get_ncbi_taxid(name):
    name = name.strip()
    name = re.sub(r'[|_]', ' ', name)
    if not name: return None
    if name in taxid_cache: return taxid_cache[name]

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {'db': 'taxonomy', 'term': name, 'retmode': 'json', 'email': EMAIL}
    if NCBI_API_KEY: params['api_key'] = NCBI_API_KEY
    
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=10)
        if r.status_code == 200:
            data = r.json()
            id_list = data.get('esearchresult', {}).get('idlist', [])
            taxid = id_list[0] if id_list else None
            taxid_cache[name] = taxid
            time.sleep(0.1 if NCBI_API_KEY else 0.34)
            return taxid
    except:
        pass
    return None

def resolve_taxid_from_fasta_header(header_line):
    clean_line = header_line.lstrip('>').strip()
    parts = clean_line.split('|')
    if len(parts) < 2: return None, None, None

    raw_name = parts[1].strip()
    tid = get_ncbi_taxid(raw_name)
    if tid: return tid, raw_name, "Species"

    genus = raw_name.split(' ')[0]
    if len(genus) > 2 and genus != raw_name:
        tid = get_ncbi_taxid(genus)
        if tid: return tid, genus, "Genus"

    return None, None, None

def download_with_retries(taxon, geo, max_retries=9): # Subimos a 5 reintentos
    params = {'taxon': taxon, 'format': 'fasta'}
    if geo: params['geo'] = geo

    for attempt in range(1, max_retries + 1):
        try:
            print(f"   ⬇️  Intento {attempt}/{max_retries}...")
            
            # Timeout muy alto para conexiones lentas
            r = requests.get(BOLD_BASE_URL, params=params, headers=HEADERS, stream=True, timeout=300)
            
            if r.status_code != 200:
                print(f"   ❌ Error HTTP {r.status_code}. Reintentando...")
                time.sleep(5)
                continue

            return process_stream_with_progress(r)

        except (ChunkedEncodingError, ConnectionError, ReadTimeout) as e:
            print(f"\n   ⚠️  Corte de red: {e}")
            wait_time = attempt * 5 # Backoff progresivo: 5s, 10s, 15s...
            if attempt < max_retries:
                print(f"   ⏳ Esperando {wait_time}s antes de reintentar...")
                time.sleep(wait_time)
            else:
                print("   ❌ Se agotaron los reintentos.")
                return None
    return None

def process_stream_with_progress(response):
    """
    Procesa el stream mostrando MB descargados en tiempo real.
    """
    valid_sequences = []
    stats = {'species': 0, 'genus': 0, 'failed': 0}
    
    current_header = None
    current_seq = []
    
    # Variables para la barra de progreso
    downloaded_bytes = 0
    start_time = time.time()
    
    # Intentamos obtener tamaño total (usualmente será None en BOLD)
    total_size = response.headers.get('content-length')
    total_size = int(total_size) if total_size else None

    print("   📊 Iniciando descarga...")

    # Iterar por líneas pero contando bytes
    for line in response.iter_lines():
        if not line: continue
        
        # Actualizar progreso
        downloaded_bytes += len(line) + 1 # +1 por el salto de línea
        
        # Calcular velocidad cada 100 líneas para no saturar la terminal
        if len(valid_sequences) % 50 == 0:
            mb = downloaded_bytes / (1024 * 1024)
            elapsed = time.time() - start_time
            speed = mb / elapsed if elapsed > 0 else 0
            
            # Barra visual personalizada para tu Hyprland ;)
            msg = f"\r   🚀 {mb:.2f} MB recibidos | {speed:.2f} MB/s | Seqs: {len(valid_sequences)}"
            sys.stdout.write(msg)
            sys.stdout.flush()

        line_str = line.decode('utf-8').strip()
        
        if line_str.startswith('>'):
            if current_header and current_seq:
                taxid, name, rank = resolve_taxid_from_fasta_header(current_header)
                if taxid:
                    orig_id = current_header.split('|')[0].lstrip('>')
                    new_header = f">{orig_id}|kraken:taxid|{taxid} {name}"
                    valid_sequences.append(f"{new_header}\n{''.join(current_seq)}")
                    if rank == "Species": stats['species'] += 1
                    else: stats['genus'] += 1
                else:
                    stats['failed'] += 1
            
            current_header = line_str
            current_seq = []
        else:
            current_seq.append(line_str)

    # Procesar último bloque
    if current_header and current_seq:
        taxid, name, rank = resolve_taxid_from_fasta_header(current_header)
        if taxid:
            orig_id = current_header.split('|')[0].lstrip('>')
            new_header = f">{orig_id}|kraken:taxid|{taxid} {name}"
            valid_sequences.append(f"{new_header}\n{''.join(current_seq)}")
            if rank == "Species": stats['species'] += 1
            else: stats['genus'] += 1

    sys.stdout.write("\n") # Limpiar línea final
    
    if valid_sequences:
        return valid_sequences, stats
    return None, None

def main():
    print("=== BOLD DOWNLOADER: MONITOR MODE ===")
    check_api_key_status()

    db_name = input("1. Nombre DB (ej. FISH): ").strip().upper()
    if not db_name: sys.exit(1)
    
    target_dir = os.path.join(DATA_REF_DIR, db_name)
    os.makedirs(target_dir, exist_ok=True)

    mode = input("\n[1] Manual\n[2] Preset: UK Fish\n[3] Preset: UK Plants\nSel: ").strip() or "1"
    
    if mode == "2": taxa, geo = PRESETS["FISH_UK"]["taxa"], PRESETS["FISH_UK"]["geo"]
    elif mode == "3": taxa, geo = PRESETS["PLANTS_UK"]["taxa"], PRESETS["PLANTS_UK"]["geo"]
    else:
        taxa = [t.strip() for t in input("Taxones: ").split(",")]
        geo = input("Geo: ").strip()
    
    print(f"\n📂 Directorio: {target_dir}")
    
    total_files = 0
    for taxon in taxa:
        filename = f"{taxon}_{'UK' if 'United Kingdom' in str(geo) else 'Global'}_Kraken.fasta"
        filepath = os.path.join(target_dir, filename)

        if os.path.exists(filepath) and os.path.getsize(filepath) > 1000:
            print(f"\n✅ {taxon}: Ya existe ({os.path.getsize(filepath)/1024/1024:.2f} MB). Saltando.")
            total_files += 1
            continue

        print(f"\n🐟 Procesando: {taxon}")
        result = download_with_retries(taxon, geo)
        
        if result:
            seqs, stats = result
            if seqs:
                with open(filepath, 'w') as f:
                    f.write("\n".join(seqs))
                print(f"   💾 Guardado: {filename}")
                print(f"      📊 Final: {stats['species']} Especies | {stats['genus']} Géneros")
                total_files += 1
        else:
            print(f"   ❌ Falló descarga de {taxon}.")

    if total_files > 0:
        print(f"\n🎉 ¡LISTO! Siguiente: ./01_build_db.sh {db_name}")

if __name__ == "__main__":
    main()
