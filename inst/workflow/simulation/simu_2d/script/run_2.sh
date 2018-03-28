#!/bin/sh

# PREFIX=../2d_common
PREFIX=`Rscript -e "cat(MixfMRI:::get.workflow())"`/2d_common
NP=8

# Rscript ${PREFIX}/01-plot_pv.r
# mpiexec -np ${NP} \
#         Rscript ${PREFIX}/02-spmd_fclust.r > ./output/log.out
Rscript ${PREFIX}/04-get_ic_lrt.r > ./summary.txt

### new
# mpiexec -np ${NP} \
#         Rscript ${PREFIX}/13-spmd_merge_lrt2.r > ./output/log.merge.out
# Rscript ${PREFIX}/14-get_ic_lrt.r > ./summary.new.txt

### new2
# mpiexec -np ${NP} \
#         Rscript ${PREFIX}/23-spmd_merge_lrt2.r > ./output/log.merge2.out
# Rscript ${PREFIX}/24-get_ic_lrt.r > ./summary.new2.txt

### new3, merging P1
Rscript ${PREFIX}/33-spmd_merge_lrt2.r > ./output/log.merge3.out

### new4, merging P2
Rscript ${PREFIX}/33-spmd_merge_lrt3.r > ./output/log.merge4.out

### new5, merging P3
Rscript ${PREFIX}/33-spmd_merge_lrt4.r > ./output/log.merge5.out

