#!/bin/sh
########################################################################
# Usage: affinity-proxy-builder.sh \
#        ALLELE_FILE PRED_METHOD PEPTIDE_FILE SAMPLE_SIZE
########################################################################

if [ $# -ne 4 ]
then
    echo "Usage: `basename $0` ALLELE_FILE PRED_METHOD PEPTIDE_FILE SAMPLE_SIZE"
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

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.stab.AffinityProxyBuilder "$@"
