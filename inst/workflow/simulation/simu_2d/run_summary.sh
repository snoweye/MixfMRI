#!/bin/sh

mkdir -p new
mkdir -p new/summary

# mkdir -p new/summary_plotpv
# mkdir -p new/summary_plotdensity
# mkdir -p new/summary_plotfclust
mkdir -p new/summary_table

mkdir -p new/summary_alternative

Rscript script_summary/summary.r \
        > new/summary/summary.txt

# Rscript script_summary/summary_plotpv.r \
#         > new/summary_plotpv/summary_plotpv.txt
# Rscript script_summary/summary_plotdensity.r \
#         > new/summary_plotdensity/summary_plotdensity.txt
# Rscript script_summary/summary_plotfclust.r \
#         > new/summary_plotfclust/summary_plotfclust.txt
Rscript script_summary/summary_table.r \
        > new/summary_table/summary_table.txt
Rscript script_summary/summary_table_new4.r \
        > new/summary_table/summary_table_new4.txt
Rscript script_summary/summary_table_new5.r \
        > new/summary_table/summary_table_new5.txt

Rscript script_summary/summary_jaccard_table.r \
        > new/summary_table/summary_jaccard_table.txt

Rscript script_summary/summary_alternative.r \
        > new/summary_alternative/summary_alternative.txt

Rscript script_summary/summary_jaccard_alternative.r \
        > new/summary_alternative/summary_jaccard_alternative.txt

