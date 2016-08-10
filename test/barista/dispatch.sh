#!/bin/bash -x
# This script dispatches a number of uge jobs which output analysis data
# of barista runs on different size data sets

if [[ -e /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load modules modules-init modules-gs modules-noble
fi

set -o nounset
set -o errexit 
set -o xtrace

#Global UGE options
SIM_ROOT=$HOME
QUEUE="testing.q"
# Which shell to use
SHELL="/bin/bash"
# Send an email notification when the jobs stops because of: end (e) or abort (a)
NOTIFY="ea"
# Location to write UGE diagnostic files
OUTPUT="$SIM_ROOT"
# -w: remove the job if an error occurs, -j: merge stderr and stdout.
MISC_OPTS="-w e -j yes"
UGE_OPTS="-q $QUEUE -S $SHELL -m $NOTIFY -o $OUTPUT $MISC_OPTS"
prefix="/net/noble/vol1/data/crux-datasets/2013wu-variation"
FASTA="$prefix/ipi.HUMAN.v3.74.fasta.gz"
MS2="$prefix/ms2-centroid"

# TODO: what kind of step size should we use?
STEP_SIZE=10
for i in {0...NUM_SPECTRA...STEP_SIZE}; do
    # Name of the job
    NAME="barista_$i"
    # TODO: is this how to specificy args?
    ARGS="-v $FASTA -v $MS2 -v $i"
    # TODO: these requirements should be dynamic on input size. How?
    # TODO: do these jobs really need exclusive use of the nodes?
    REQS="-l h_rt=24:0:0 exculsive=true disk_free=50G"
    # run jobs in uge
    qsub $UGE_OPTS $ARGS -N "$NAME_time" $REQS time_job.sh
    qsub $UGE_OPTS $ARGS -N "$NAME_memory" $REQS memory_job.sh
fi
