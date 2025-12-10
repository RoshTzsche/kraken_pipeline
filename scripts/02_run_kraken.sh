#!/bin/bash
# Script: 02_run_kraken.sh
# Description: Classifies Paired-End (PE) reads against the selected DB.
# Usage: Run from the 'scripts/' directory: ./02_run_kraken.sh

# --- CONFIGURATION ---
# Adjust your DB path here if necessary
DB_PATH="../data/dbs/CUSTOM_DB" 

INPUT_DIR="../data/raw_fastq"
OUTPUT_DIR="../results/reports"
THREADS=12

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "🐙 KRAKEN2 CLASSIFIER PIPELINE"
echo "=========================================="

# Verify DB existence
if [ ! -f "$DB_PATH/taxo.k2d" ]; then
    echo "❌ ERROR: Database not found at $DB_PATH"
    echo "   Did you run '01_build_db.sh' first?"
    exit 1
fi

echo "▶️  Starting sample processing..."

# Loop to find R1 files (assuming _1.fastq.gz extension)
for r1 in "$INPUT_DIR"/*_1.fastq.gz; do
    
    # Check if files actually exist
    if [ ! -e "$r1" ]; then
        echo "⚠️  No FASTQ files found in $INPUT_DIR"
        break
    fi

    # 1. Define pairs and sample names
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    
    # Verify R2 exists
    if [ ! -f "$r2" ]; then
        echo "⚠️  Skipping $(basename "$r1"): Missing R2 pair."
        continue
    fi

    # Clean sample name (Adjust sed commands to match your filename structure)
    filename=$(basename "$r1")
    # Example: cleaned_SAMPLE1_EKDN_1.fastq.gz -> SAMPLE1
    sample=$(echo "$filename" | sed 's/cleaned//' | sed 's/_EKDN.*//' | sed 's/^_//')

    report_file="$OUTPUT_DIR/${sample}_report.txt"

    # 2. Check if already processed
    if [ -s "$report_file" ]; then
        echo "   ⏭️  $sample already processed. Skipping."
        continue
    fi

    echo "   🦠 Analyzing: $sample"

    # 3. Run Kraken2
    # --output /dev/null discards the heavy read-by-read output (we only need the report)
    kraken2 --db "$DB_PATH" \
            --threads "$THREADS" \
            --paired "$r1" "$r2" \
            --report "$report_file" \
            --use-names \
            --output /dev/null

    if [ $? -eq 0 ]; then
        echo "      ✅ Report generated."
    else
        echo "      ❌ Error processing $sample"
        rm -f "$report_file" # Remove corrupt/empty report
    fi

done

echo "=========================================="
echo "🎉 Classification complete! Reports saved in: $OUTPUT_DIR"