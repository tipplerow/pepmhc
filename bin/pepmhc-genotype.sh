#!/bin/sh
########################################################################
# Usage: pepmhc-genotype.sh [JVM OPTIONS] FILE1 [FILE2 ...]
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

WORKDIR=`pwd`

cd `dirname $0`
cd ..

export PEPMHC_HOME=`pwd`
echo PEPMHC_HOME=$PEPMHC_HOME

cd $WORKDIR
${JAM_HOME}/bin/jam-run.sh $PEPMHC_HOME pepmhc.app.GenotypePresentCalc "$@"
