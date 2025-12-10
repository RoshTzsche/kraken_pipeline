#!/bin/bash
# Script: 01_build_db.sh
# Description: Builds a custom Kraken2 database using local FASTA files.
#              Offers to download taxonomy if missing.
# Usage: Run from the 'scripts/' directory: ./01_build_db.sh

# --- CONFIGURATION ---
# Paths relative to the scripts/ directory
FASTA_DIR="../data/fasta_ref"
DB_DIR="../data/dbs/CUSTOM_DB"
THREADS=12

echo "=========================================="
echo "🏗️  KRAKEN2 DATABASE BUILDER"
echo "=========================================="

# 1. Create DB directory
mkdir -p "$DB_DIR"

# 2. Check Taxonomy (Interactive)
# Kraken requires the 'taxonomy' folder containing names.dmp and nodes.dmp
if [ -f "$DB_DIR/taxonomy/names.dmp" ] && [ -f "$DB_DIR/taxonomy/nodes.dmp" ]; then
    echo "✅ Taxonomy detected successfully in $DB_DIR/taxonomy/"
else
    echo "⚠️  Taxonomy not found in $DB_DIR/taxonomy/"
    echo "---------------------------------------------------"
    echo "   You have two options:"
    echo "   1. Download it now from NCBI (Requires Internet)."
    echo "   2. Quit and manually copy an existing 'taxonomy' folder."
    echo "---------------------------------------------------"
    
    # Interactive prompt
    read -p "❓ Do you want to download the taxonomy now? (y/n): " answer

    # Convert answer to lowercase for easier comparison
    if [[ "${answer,,}" == "y" ]]; then
        echo "⬇️  Starting taxonomy download (this may take a while)..."
        # Download command
        kraken2-build --download-taxonomy --db "$DB_DIR" --threads "$THREADS"
        
        # Verify download success
        if [ -f "$DB_DIR/taxonomy/names.dmp" ]; then
            echo "✅ Taxonomy downloaded successfully."
        else
            echo "❌ Download failed or incomplete."
            exit 1
        fi
    else
        echo "🛑 Process aborted by user."
        echo "   Please manually copy your 'taxonomy' folder into: $DB_DIR"
        exit 0
    fi
fi

# 3. Add FASTA files to the library
echo "➕ Adding FASTA files from $FASTA_DIR..."
count=0

# Enable nullglob to handle empty directories gracefully
shopt -s nullglob
files=("$FASTA_DIR"/*.fasta "$FASTA_DIR"/*.fa)

if [ ${#files[@]} -eq 0 ]; then
    echo "⚠️  No .fasta or .fa files found in $FASTA_DIR"
    exit 1
fi

for fasta_file in "${files[@]}"; do
    echo "   📄 Processing: $(basename "$fasta_file")"
    
    # --no-masking speeds up the process by skipping DustMasker (optional)
    kraken2-build --add-to-library "$fasta_file" --db "$DB_DIR" --threads "$THREADS" --no-masking
    ((count++))
done

echo "📦 Added $count files to the library."

# 4. Build the final index
echo "🔨 Building final index..."
kraken2-build --build --db "$DB_DIR" --threads "$THREADS"

# 5. Final Verification
if [ -f "$DB_DIR/taxo.k2d" ]; then
    echo "🎉 SUCCESS! Database ready at: $DB_DIR"
else
    echo "❌ ERROR: Build failed. 'taxo.k2d' was not generated."
fi