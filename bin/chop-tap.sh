#!/bin/sh
########################################################################
# Usage: chop-tap.sh [JVM OPTIONS] ENSEMBL_FILE
########################################################################

if [ $# -ne 2 ]
then
    echo "Usage:" `basename $0` "[JVM OPTIONS] ENSEMBL_FILE OUTPUT_FILE"
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

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.app.ChopTAP "$@"
