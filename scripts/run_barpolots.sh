#!/bin/bash

INPUT_DIR="../results/final_tables"
METADATA="../data/metadata.xlsx"

# 1. Definimos una función que parsea el nombre del archivo y extrae el Organismo
get_org_name() {
    local filepath="$1"
    local filename=$(basename "$filepath")
    
    if [[ "$filename" == "taxonomic_classification_clean.xlsx" ]]; then
        echo "Microbiome"
    elif [[ "$filename" == Taxonomy_*_Cumulative_Reads.xlsx ]]; then
        # Manipulación de strings nativa de Bash
        local org="${filename#Taxonomy_}"
        org="${org%_Cumulative_Reads.xlsx}"
        echo "$org"
    else
        echo "Unknown_Organism"
    fi
}

# 2. Exportamos la función
export -f get_org_name

echo "Iniciando generación de Barplots en paralelo..."

# 3. Cambio crucial: usamos ':::: -' para leer los archivos del find
find "$INPUT_DIR" -type f -name "*.xlsx" ! -name "Taxonomy_ALL_Cumulative_Reads.xlsx" | \
parallel --verbose --jobs 4 \
  python 04_generate_Barplots.py \
    -d {1} \
    -m "$METADATA" \
    -c Treatment \
    -r {2} \
    -t {3} \
    -org '$(get_org_name "{1}")' \
    -fmt png \
    -ord "Control" "Partial" "Untreated" \
    :::: - \
    ::: phylum \
    ::: 0.005 0.01 0.02 0.007

echo "¡Procesamiento paralelo completado!"
