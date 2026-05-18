#!/bin/bash

INPUT_DIR="../results/final_tables"
METADATA="../data/metadata.xlsx"

# 1. Function to parse the filename and extract the Organism name
#    This ensures our output files don't overwrite each other.
get_org_name() {
    local filepath="$1"
    local filename=$(basename "$filepath")
    
    if [[ "$filename" == "taxonomic_classification_clean.xlsx" ]]; then
        echo "Microbiome"
    elif [[ "$filename" == Taxonomy_*_Cumulative_Reads.xlsx ]]; then
        # Native Bash string manipulation to extract the organism
        local org="${filename#Taxonomy_}"
        org="${org%_Cumulative_Reads.xlsx}"
        echo "$org"
    else
        echo "Unknown_Organism"
    fi
}

# 2. Export the function so GNU parallel subshells can access it
export -f get_org_name

echo "Starting parallel generation of PCoA plots across multiple ranks..."

# 3. Use find to grab all relevant .xlsx files and pipe them to parallel
#    {1} represents the filepath piped from find.
#    {2} represents the taxonomic rank from the list at the bottom.
find "$INPUT_DIR" -type f -name "*.xlsx" ! -name "Taxonomy_ALL_Cumulative_Reads.xlsx" | \
parallel --verbose --jobs 4 \
  python 05_generate_PCoA_PieChart.py \
    -d {1} \
    -r {2} \
    -m "$METADATA" \
    -c Treatment \
    -id SampleID \
    --mode pcoa \
    -fmt png \
    --unknown drop_all \
    -org '$(get_org_name "{1}")' \
    :::: - \
    ::: species genus family phylum

echo "Parallel processing complete!"