#!/bin/bash

condatainer scheduler | head -7

# Slurm test
condatainer create grcm39/salmon/1.10.2/gencodeM6
condatainer run src/run_chain_dep.sh --auto-install # chain dep submit (Slurm)
condatainer run src/run_chain_external.sh --auto-install # should fail (external img cannot auto-install)

# LSF test (Cross-scheduler test)
condatainer run --debug src/run_chain_lsf.sh
# PBS test (Cross-scheduler test)
condatainer run --debug src/run_chain_pbs.sh

# GPU jobs test
condatainer create -p src/pytorch -f src/pytorch.def
condatainer run src/gpu_test_slurm.sh
