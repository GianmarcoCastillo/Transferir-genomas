#!/bin/bash
tsv="mags_filtrados_rest-PAS.tsv"
destino="/media/mordor/gianmarco/tesis_mag/bakta_analisis_MAGs/analisis_mags_100k"
mkdir -p "$destino"

tail -n +2 "$tsv" | cut -f1 | tr -d '\r' | while read acc; do
    for d in /media/mordor/gianmarco/tesis_mag/gtdbtk_analisis_MAGs_completo/mags_lote_{1..8}; do
        f=$(find "$d" -maxdepth 1 -type f -name "${acc}*.fna" -print -quit 2>/dev/null)
        [ -n "$f" ] && cp "$f" "$destino/" && break
    done
done
