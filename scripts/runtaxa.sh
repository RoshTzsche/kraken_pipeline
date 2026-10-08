#!/bin/bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
INPUT_DIR="../results/final_tables"
METADATA="../data/metadata.xlsx"
get_org_name() {
    local filename
    filename=$(basename "$1")
    if [[ "$filename" == "taxonomic_classification_clean.xlsx" ]]; then
        echo "Microbiome"
    else
        filename="${filename#Taxonomy_}"
        echo "${filename%_Cumulative_Reads.xlsx}"
    fi
}
export -f get_org_name
find "$INPUT_DIR" -type f -name "*.xlsx" ! -name "Taxonomy_ALL_Cumulative_Reads.xlsx" | \
parallel --verbose --jobs 4 \
  python 06_generate_Violin_ANOVA.py \
    -d {1} -m "$METADATA" -c Time -id SampleID \
    -r {2} -t {3} -org '$(get_org_name "{1}")' \
    -fmt png --mode taxa \
    :::: - ::: genus species ::: 0.005 0.01 0.02 0.007
