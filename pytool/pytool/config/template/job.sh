#!/bin/bash

#PJM -g {{ user_id }}
#PJM -L rscgrp=small
#PJM --mpi max-proc-per-node={{ proc_per_node }}
#PJM -L node={{ node }}
#PJM -L elapse={{ elapse }}
#PJM -x PJM_LLIO_GFSCACHE=/vol0004

if [ "$1" == "--local" ]; then
  export OMP_NUM_THREADS={{ thread }}
  mpirun -np {{ proc }} spdyn inp/{{ mode }}{{ index }}.inp | tee out/{{ mode }}{{ index }}.log
else
  . /vol0004/apps/oss/spack/share/spack/setup-env.sh
  spack load genesis@2.1.1~mixed~single
  export OMP_NUM_THREADS={{ thread }}
  mpiexec -std-proc out/{{ mode }}{{ index }}.log -n ${PJM_MPI_PROC} spdyn inp/{{ mode }}{{ index }}.inp
fi
