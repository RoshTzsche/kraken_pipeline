#!/bin/bash

INPUT_DIR="../results/final_tables"
METADATA="../data/metadata.xlsx"

# 1. Función para extraer el nombre del Organismo (idéntica a la anterior)
get_org_name() {
    local filepath="$1"
    local filename=$(basename "$filepath")
    
    if [[ "$filename" == "taxonomic_classification_clean.xlsx" ]]; then
        echo "Microbiome"
    elif [[ "$filename" == Taxonomy_*_Cumulative_Reads.xlsx ]]; then
        local org="${filename#Taxonomy_}"
        org="${org%_Cumulative_Reads.xlsx}"
        echo "$org"
    else
        echo "Unknown_Organism"
    fi
}
export -f get_org_name

echo "============================================================"
echo " Iniciando generación de Alpha Diversity (Violin Plots) "
echo "============================================================"
# Modo Alpha: Cruza Archivos x Rangos (genus, species). No usa -t.
find "$INPUT_DIR" -type f -name "*.xlsx" ! -name "Taxonomy_ALL_Cumulative_Reads.xlsx" | \
parallel --verbose --jobs 4 \
  python 06_generate_Violin_ANOVA.py \
    -d {1} \
    -m "$METADATA" \
    -c Treatment \
    -id SampleID \
    -r {2} \
    -org '$(get_org_name "{1}")' \
    -fmt png \
    --mode alpha \
    :::: - \
    ::: genus species



echo "============================================================"
echo " ¡Procesamiento paralelo completado!                        "
echo "============================================================"
