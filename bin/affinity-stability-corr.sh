#!/bin/sh
########################################################################
# Usage: affinity-stability-corr.sh ALLELE PEPTIDE_FILE OUTPUT_FILE [SAMPLE_SIZE]
########################################################################

if [ $# -lt 3 ]
then
    echo "Usage:" `basename $0` "ALLELE PEPTIDE_FILE OUTPUT_FILE [SAMPLE_SIZE]"
    exit 1
fi

if [ -z "${JAM_HOME}" ]
then
    echo "Environment variable JAM_HOME is not set; exiting."
    exit 1
fi

if [ -z "${PEPMHC_HOME}" ]
then
    echo "Environment variable PEPMHC_HOME is not set; exiting."
    exit 1
fi

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.app.AffinityStabilityCorr "$@"
