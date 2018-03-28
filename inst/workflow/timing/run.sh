#!/bin/sh

NP=4
mpiexec -np ${NP} \
        Rscript ./02-spmd_fclust.r > log_${NP}.txt

NP=8
mpiexec -np ${NP} \
        Rscript ./02-spmd_fclust.r > log_${NP}.txt

NP=16
mpiexec -np ${NP} \
        Rscript ./02-spmd_fclust.r > log_${NP}.txt

NP=32
mpiexec -np ${NP} \
        Rscript ./02-spmd_fclust.r > log_${NP}.txt

NP=64
mpiexec -np ${NP} \
        Rscript ./02-spmd_fclust.r > log_${NP}.txt
