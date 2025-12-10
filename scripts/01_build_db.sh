#!/bin/bash
# Script: 01_build_db.sh
# Description: Builds a Kraken2 DB using a CENTRALIZED taxonomy to save disk space.
#              Uses symbolic links instead of copying massive files.
# Usage: ./01_build_db.sh [DB_NAME]

# --- CONFIGURATION ---
DB_NAME=${1:-"CUSTOM_DB"}

FASTA_DIR="../data/fasta_ref"
DBS_ROOT="../data/dbs"
DB_DIR="$DBS_ROOT/$DB_NAME"

#Central location for the 60GB taxonomy files
CENTRAL_TAXO_DIR="../data/taxonomy" 

THREADS=12

echo "=========================================="
echo "🏗️  KRAKEN2 DB BUILDER (Linked Taxonomy)"
echo "=========================================="
echo "Target DB: $DB_NAME"

# 1. Prepare Directories
mkdir -p "$DB_DIR"
mkdir -p "$CENTRAL_TAXO_DIR"

# 2. CENTRALIZED TAXONOMY LOGIC
# We want names.dmp and nodes.dmp to live ONLY in ../data/taxonomy/
# and symlink them to specific DBs.

if [ -f "$CENTRAL_TAXO_DIR/names.dmp" ] && [ -f "$CENTRAL_TAXO_DIR/nodes.dmp" ]; then
    echo "✅ Central taxonomy found."
else
    echo "🔍 Central taxonomy is empty. Checking if we can migrate from an old DB..."
    
    # Check if user has taxonomy in a previous DB to move it (Save space)
    MOVED=false
    for other_db in "$DBS_ROOT"/*; do
        if [ -d "$other_db/taxonomy" ] && [ -f "$other_db/taxonomy/names.dmp" ] && [ ! -L "$other_db/taxonomy" ]; then
            echo "⚠️  Found real taxonomy data in: $other_db"
            echo "   To save space, we should MOVE this to the central folder."
            read -p "❓ Move taxonomy from $other_db to central storage? (y/n): " answer
            
            if [[ "${answer,,}" == "y" ]]; then
                echo "📦 Moving files... (This saves 60GB duplicating)"
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

    # If still not found, offer download
    if [ "$MOVED" = false ]; then
        echo "⚠️  No taxonomy found."
        read -p "❓ Download NCBI taxonomy to Central Storage now? (y/n): " answer
        if [[ "${answer,,}" == "y" ]]; then
            # We use kraken2-build to download into the central dir
            # Note: kraken2-build expects the folder to be the DB root, so we trick it slightly
            # or just download, but let's use the tool's native download command targeting central dir context
            # Actually, standard download puts it in /taxonomy. 
            echo "⬇️  Downloading..."
            kraken2-build --download-taxonomy --db "$CENTRAL_TAXO_DIR/temp_dl" --threads "$THREADS"
            # Move out of temp structure to our clean central structure
            mv "$CENTRAL_TAXO_DIR/temp_dl/taxonomy/"* "$CENTRAL_TAXO_DIR/"
            rm -rf "$CENTRAL_TAXO_DIR/temp_dl"
        else
            echo "🛑 Aborted. You need taxonomy in $CENTRAL_TAXO_DIR to proceed."
            exit 1
        fi
    fi
fi

# 3. LINKING 
# Create a symbolic link: DB/taxonomy -> ../../taxonomy
echo "🔗 Linking taxonomy..."

# Remove existing folder if it exists (to replace with link)
if [ -d "$DB_DIR/taxonomy" ] && [ ! -L "$DB_DIR/taxonomy" ]; then
    echo "⚠️  Warning: This DB has a physical copy of taxonomy."
    echo "   Skipping link creation to avoid accidental data loss."
    echo "   (If you want to save space, delete $DB_DIR/taxonomy manually and re-run)."
elif [ ! -L "$DB_DIR/taxonomy" ]; then
    # Create the link relative path: from inside data/dbs/NAME/ -> ../../taxonomy
    ln -s "../../taxonomy" "$DB_DIR/taxonomy"
    echo "✅ Symlink created: $DB_DIR/taxonomy -> $CENTRAL_TAXO_DIR"
else
    echo "✅ Symlink already exists."
fi

# 4. Add FASTA files
shopt -s nullglob
files=("$FASTA_DIR"/*.fasta "$FASTA_DIR"/*.fa)

if [ ${#files[@]} -eq 0 ]; then
    echo "❌ ERROR: No .fasta files found in $FASTA_DIR"
    exit 1
fi

echo "➕ Adding FASTA files..."
for fasta_file in "${files[@]}"; do
    echo "   📄 $(basename "$fasta_file")"
    kraken2-build --add-to-library "$fasta_file" --db "$DB_DIR" --threads "$THREADS" --no-masking
done

# 5. Build
echo "🔨 Building final index..."
kraken2-build --build --db "$DB_DIR" --threads "$THREADS"

if [ -f "$DB_DIR/taxo.k2d" ]; then
    echo "🎉 SUCCESS! Database '$DB_NAME' ready (using shared taxonomy)."
else
    echo "❌ ERROR: Build failed."
fi