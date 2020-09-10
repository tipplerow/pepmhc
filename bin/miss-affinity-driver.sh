#!/bin/sh
########################################################################
# Usage: miss-driver-driver.sh [JVM OPTIONS] PROP_FILE1 [PROP_FILE2 ...]
########################################################################

if [ $# -lt 1 ]
then
    echo "Usage:" `basename $0` "[JVM OPTIONS] PROP_FILE1 [PROP_FILE2 ...]"
    exit 1
fi

if [ -z "${PEPMHC_HOME}" ]
then
    echo "Environment variable PEPMHC_HOME is not set; exiting."
    exit 1
fi

${PEPMHC_HOME}/bin/pepmhc-run.sh pepmhc.miss.MissDriverDriver "$@"
