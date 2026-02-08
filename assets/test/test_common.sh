#!/bin/bash
# Basic tests for condatainer commands
condatainer avail
condatainer avail --remote
condatainer list

# Create; execute; and remove conda overlays
condatainer create samtools/1.16
condatainer exec -o samtools/1.16 samtools --version | head -3
condatainer e -n samtools/1.16 -- samtools --version | head -3
condatainer remove samtools/1.16 -y

# Create; execute; and remove apptainer containers
condatainer create ubuntu22/code-server
condatainer exec -o ubuntu22/code-server code-server --version
condatainer exec -b ubuntu22/code-server code-server --version
condatainer e -n ubuntu22/code-server -- code-server --version
condatainer exec -o ubuntu22/code-server cat /etc/os-release
condatainer remove ubuntu22/code-server -y

# Create; execute; remote docker containers
condatainer create -n r4.2.2 -s docker://posit/r-base:4.2.2-noble-amd64
condatainer exec -o r4.2.2 R --version
condatainer e r4.2.2 -- bash -c "type R"

# Create; execute; and remove script built overlays
condatainer create rpvg/1.0
condatainer exec -o rpvg/1.0 rpvg --help 2>&1 | head
condatainer rm rpvg/1.0 -y

# Create; execute; and remove script built overlays which requires link input
condatainer create cellranger/8.0.1
condatainer exec -o cellranger/8.0.1 cellranger --help 2>&1 | head

# Create; execute; and remove script built reference data overlays
condatainer create grcm39/transcript-gencode/M6
condatainer exec -o grcm39/transcript-gencode/M6 bash -c 'echo $TRANSCRIPT_FASTA'
condatainer e -n grcm39/transcript-gencode/M6 -- bash -c 'head -2 $TRANSCRIPT_FASTA'
condatainer rm grcm39/transcript-gencode/M6 -y

# Check out the script dependencies
condatainer check src/check_module.sh
condatainer check src/check_ml.sh
condatainer check src/check_dep.sh

# Auto install
condatainer check src/check_module.sh --auto-install

# Local run
condatainer run src/check_module.sh --local | head -5
condatainer run src/check_ml.sh --local 2>&1 | head -5
condatainer run src/check_dep.sh | head -5

## Overlay create info check chown
condatainer o --sparse test
condatainer overlay info test.img
condatainer overlay create -s 1g src/test
condatainer overlay info src/test.img
condatainer overlay check src/test.img
condatainer overlay chown src/test.img
condatainer e src/test.img r4.2.2 -- R -e "options(repos = c(CRAN = \"https://packagemanager.posit.co/cran/__linux__/noble/latest\")); install.packages('pak'); .libPaths(); sessionInfo()"
condatainer overlay chown -u 1000 -g 1000 src/test.img
debugfs -R 'stat /upper' src/test.img | head
debugfs -R 'stat /upper/opt' src/test.img | head
condatainer overlay chown --root src/test.img
debugfs -R 'stat /upper' src/test.img | head

## Using img overlay (current uid)
condatainer exec -o test.img ls /ext3/env
condatainer e test.img # Should open bash shell
condatainer e test.img -- mm-install python=3.9 -y
condatainer e -r test.img -- mm-install r-base=4.2 -y # Should fail (not writable)
condatainer e -r test.img -- touch /ext3/env/testfile.txt # Should fail (not writable)
condatainer exec -o test.img mm-install r-base=4.2 -y # Should fail (not writable)
condatainer e test.img -r -- mm-list # CNT_CONDA_PREFIX is not set! cannot work
condatainer exec -o test.img micromamba list | head -5 # should work

## Using img overlay (root uid)
condatainer exec -o src/test.img --fakeroot id # should show uid 0 (auto detect root)
condatainer e src/test.img -- apt update
condatainer e -r src/test.img -- apt update # Should fail (not writable)

## img overlay resize
condatainer overlay resize -s 5g test.img
condatainer overlay info test.img | grep Used
condatainer overlay resize -s 8192 test.img
condatainer overlay info test.img | grep Used

## instance
condatainer instance -h
condatainer instance start -o test.img -w test_instance
condatainer instance list
condatainer instance exec test_instance python --version # should show 3.9.*
condatainer instance exec test_instance mm-update python=3.10 -y
condatainer instance exec test_instance python --version # should show 3.10.*
condatainer instance stats test_instance
condatainer instance stop test_instance
