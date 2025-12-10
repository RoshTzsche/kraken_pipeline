#!/bin/bash
# Script: 02_run_kraken.sh
# Description: Classifies Paired-End (PE) reads against a SPECIFIC DB.
# Usage: ./02_run_kraken.sh [DB_NAME]
# Example: ./02_run_kraken.sh FISH

# --- CONFIGURATION ---
# 1. Capture the first argument as the DB Name (default: CUSTOM_DB)
DB_NAME=${1:-"CUSTOM_DB"}

# 2. Set dynamic paths based on that name
DB_PATH="../data/dbs/$DB_NAME"
INPUT_DIR="../data/raw_fastq"

# 3. Create a specific subfolder for results to keep organized
OUTPUT_DIR="../results/reports/$DB_NAME"
THREADS=12

echo "=========================================="
echo "🐙 KRAKEN2 CLASSIFIER -> DB: $DB_NAME"
echo "=========================================="

# Verify DB existence
if [ ! -f "$DB_PATH/taxo.k2d" ]; then
    echo "❌ ERROR: Database '$DB_NAME' not found in $DB_PATH"
    echo "   Available databases:"
    ls -1 ../data/dbs/
    echo "   Did you run '01_build_db.sh $DB_NAME' first?"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "📂 Reports will be saved to: $OUTPUT_DIR"
echo "▶️  Starting sample processing..."

# Loop to find R1 files
for r1 in "$INPUT_DIR"/*_1.fastq.gz; do
    
    if [ ! -e "$r1" ]; then
        echo "⚠️  No FASTQ files found in $INPUT_DIR"
        break
    fi

    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    if [ ! -f "$r2" ]; then
        echo "⚠️  Skipping $(basename "$r1"): Missing R2 pair."
        continue
    fi

    filename=$(basename "$r1")
    # Clean sample name
    sample=$(echo "$filename" | sed 's/cleaned//' | sed 's/_EKDN.*//' | sed 's/^_//')

    report_file="$OUTPUT_DIR/${sample}_report.txt"

    if [ -s "$report_file" ]; then
        echo "   ⏭️  $sample already processed. Skipping."
        continue
    fi

    echo "   🦠 Analyzing $sample against $DB_NAME..."

    # Run Kraken2
    kraken2 --db "$DB_PATH" \
            --threads "$THREADS" \
            --paired "$r1" "$r2" \
            --report "$report_file" \
            --use-names \
            --output /dev/null

    if [ $? -ne 0 ]; then
        echo "      ❌ Error processing $sample"
        rm -f "$report_file"
    fi

done

echo "=========================================="
echo "🎉 Done! Results in: $OUTPUT_DIR"