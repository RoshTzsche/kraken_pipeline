#!/bin/bash
# Script: 01_build_db.sh
# Description: Builds a Kraken2 DB using a CENTRALIZED taxonomy and MODULAR references.
# Usage: ./01_build_db.sh [DB_NAME]
# Example: ./01_build_db.sh PLANTS (Looks for fasta files in data/fasta_ref/PLANTS)

# --- CONFIGURATION ---
DB_NAME=${1:-"CUSTOM_DB"}

# UPDATE: References are now separated by subfolder matching the DB name
FASTA_SOURCE="../data/fasta_ref/$DB_NAME"

DBS_ROOT="../data/dbs"
DB_DIR="$DBS_ROOT/$DB_NAME"
CENTRAL_TAXO_DIR="../data/taxonomy" 
THREADS=12

echo "=========================================="
echo "🏗️  KRAKEN2 DB BUILDER: $DB_NAME"
echo "=========================================="

# 0. Check Input References first
if [ ! -d "$FASTA_SOURCE" ]; then
    echo "❌ ERROR: Reference directory not found."
    echo "   Expected path: $FASTA_SOURCE"
    echo "   Please create the folder 'data/fasta_ref/$DB_NAME' and put your .fasta files there."
    echo "   If you already have your folder, please put the name in the running script"
    echo "   example: ./01_build_db.sh PLANTS - to use the data/fasta_ref/PLANTS folder"
    echo "   Check the README for more information"
    exit 1
fi

# Check if folder is empty
shopt -s nullglob
files=("$FASTA_SOURCE"/*.fasta "$FASTA_SOURCE"/*.fa)
if [ ${#files[@]} -eq 0 ]; then
    echo "❌ ERROR: No .fasta or .fa files found in $FASTA_SOURCE"
    exit 1
fi

# 1. Prepare Directories
mkdir -p "$DB_DIR"
mkdir -p "$CENTRAL_TAXO_DIR"

# 2. CENTRALIZED TAXONOMY LOGIC (Symlinks)
if [ -f "$CENTRAL_TAXO_DIR/names.dmp" ] && [ -f "$CENTRAL_TAXO_DIR/nodes.dmp" ]; then
    echo "✅ Central taxonomy found."
else
    echo "🔍 Central taxonomy is empty. Checking if we can migrate from an old DB..."
    
    # Check if user has taxonomy in a previous DB to move it (Save space)
    MOVED=false
    for other_db in "$DBS_ROOT"/*; do
        # Look for a folder that has taxonomy BUT is not a symlink (! -L)
        if [ -d "$other_db/taxonomy" ] && [ -f "$other_db/taxonomy/names.dmp" ] && [ ! -L "$other_db/taxonomy" ]; then
            echo "⚠️  Found real taxonomy data in: $other_db"
            echo "   To save space, we should MOVE this to the central folder."
            read -p "❓ Move taxonomy from $other_db to central storage? (y/n): " answer
            
            if [[ "${answer,,}" == "y" ]]; then
                echo "📦 Moving files... (This saves 60GB of duplication)"
                mv "$other_db/taxonomy/"* "$CENTRAL_TAXO_DIR/"
                # Remove empty folder and replace with link
                rmdir "$other_db/taxonomy"
                ln -s "../../taxonomy" "$other_db/taxonomy"
                echo "✅ Migration complete. Old DB is now linked too."
                MOVED=true
                break
            fi
        fi
    done

    # If still not found after checking neighbors, offer download
    if [ "$MOVED" = false ]; then
        # Check again in case the migration failed or wasn't approved
        if [ ! -f "$CENTRAL_TAXO_DIR/names.dmp" ]; then
             echo "⚠️  No taxonomy found."
             read -p "❓ Download NCBI taxonomy to Central Storage now? (y/n): " answer
             if [[ "${answer,,}" == "y" ]]; then
                 echo "⬇️  Downloading taxonomy (this may take a while)..."
                 # Download to temp folder then move to keep structure clean
                 kraken2-build --download-taxonomy --db "$CENTRAL_TAXO_DIR/temp_dl" --threads "$THREADS"
                 mv "$CENTRAL_TAXO_DIR/temp_dl/taxonomy/"* "$CENTRAL_TAXO_DIR/"
                 rm -rf "$CENTRAL_TAXO_DIR/temp_dl"
                 echo "✅ Taxonomy downloaded to central storage."
             else
                 echo "🛑 Aborted. You need taxonomy in $CENTRAL_TAXO_DIR to proceed."
                 exit 1
             fi
        fi
    fi
fi

# 3. LINKING TAXONOMY
echo "🔗 Linking taxonomy..."
if [ ! -L "$DB_DIR/taxonomy" ]; then
    # Remove empty folder if exists to allow linking
    [ -d "$DB_DIR/taxonomy" ] && rmdir "$DB_DIR/taxonomy" 2>/dev/null
    ln -s "../../taxonomy" "$DB_DIR/taxonomy"
    echo "✅ Symlink created."
else
    echo "✅ Link already exists."
fi

# 4. Add FASTA files
echo "➕ Adding FASTA files from: $FASTA_SOURCE"
for fasta_file in "${files[@]}"; do
    echo "   📄 $(basename "$fasta_file")"
    kraken2-build --add-to-library "$fasta_file" --db "$DB_DIR" --threads "$THREADS" --no-masking
done

# 5. Build
echo "🔨 Building final index..."
kraken2-build --build --db "$DB_DIR" --threads "$THREADS"

if [ -f "$DB_DIR/taxo.k2d" ]; then
    echo "🎉 SUCCESS! Database '$DB_NAME' ready."
else
    echo "❌ ERROR: Build failed."
fi