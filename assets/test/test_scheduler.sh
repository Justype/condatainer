#!/bin/bash

condatainer scheduler | head -7

condatainer create grcm39/salmon/1.10.2/gencodeM6
condatainer run src/run_chain_dep.sh --auto-install # chain dep submit
condatainer run src/run_chain_external.sh --auto-install # should fail (external img cannot auto-install)
