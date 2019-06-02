#!/bin/sh
########################################################################
# Usage: gene-processor.sh [JVM OPTIONS] INPUT_FILE OUTPUT_FILE
########################################################################

if [ $# -lt 2 ]
then
    echo "Usage:" `basename $0` "[JVM OPTIONS] INPUT_FILE OUTPUT_FILE"
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

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.proc.GeneProcessor "$@"
