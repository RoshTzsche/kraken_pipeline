echo "============================================================"
echo " Iniciando generación de Taxa Abundance (Violin Plots)  "
echo "============================================================"
# Modo Taxa: Cruza Archivos x Rangos (genus, species) x Umbrales
find "$INPUT_DIR" -type f -name "*.xlsx" ! -name "Taxonomy_ALL_Cumulative_Reads.xlsx" | \
parallel --verbose --jobs 4 \
  python 06_generate_Violin_ANOVA.py \
    -d {1} \
    -m "$METADATA" \
    -c Treatment \
    -id SampleID \
    -r {2} \
    -t {3} \
    -org '$(get_org_name "{1}")' \
    -fmt png \
    --mode taxa \
    :::: - \
    ::: genus species \
    ::: 0.005 0.01 0.02 0.007
