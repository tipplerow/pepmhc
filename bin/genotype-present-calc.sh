#!/bin/sh
########################################################################
# Usage: genotype-present-calc.sh [JVM OPTIONS] FILE1 [FILE2 ...]
########################################################################

if [ $# -lt 1 ]
then
    echo "Usage:" `basename $0` "[JVM OPTIONS] FILE1 [FILE2 ...]"
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

${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.app.GenotypePresentCalc "$@"
