#!/bin/sh
########################################################################
# Usage: process-peptide-source.sh \
#        MISSENSE_DIR SELF_PEP_FILE BARCODE_FILE PEP_SOURCE_DIR
########################################################################

if [ $# -lt 4 ]
then
    echo "Usage: `basename $0` MISSENSE_DIR SELF_PEP_FILE BARCODE_FILE PEP_SOURCE_DIR [AGPRO_PROP_FILE]"
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

JVM_FLAGS="-Xmx96g -verbose:gc"

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME $JVM_FLAGS pepmhc.neo.PeptideSourceProcessor "$@"
