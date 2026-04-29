#!/bin/bash

set -eo pipefail
mkdir -p ./log

# Initialize conda in the current shell
eval "$(conda shell.bash hook)"
# Load modules and activate environment
conda activate cpdb_env

MIN_PERCENT=0.05
REGION_LIST=(Sinushorn-SAnode-AtriumR) # (LA RA SAN PV)
# YOUNG_OR_OLD=(young) # (young old)

for REGION in "${REGION_LIST[@]}"; do
    sbatch \
      -p icelake-himem \
      -A POLONIUS-TEICHMANN-SL3-CPU \
      -N 1 \
      -n 1 \
      -c 4 \
      -t 1:00:00 \
      -J "cellphone_${REGION}" \
      -o ./log/cellphone_%j_${REGION}.out \
      -e ./log/cellphone_%j_${REGION}.err \
      --wrap="python ./run_cpdb_statistical.py --min-percent $MIN_PERCENT --region $REGION"
done