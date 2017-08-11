#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "  "
    echo "  ./purge.sh ../rootfiles/<systematic>/<analysis>"
    echo "  ./purge.sh ../minitrees/<systematic>/<analysis>"
    echo "  "
    exit -1
fi

FOLDER="$1"

pushd $FOLDER

rm S* Z* W* V* *03Feb* g* G* H* TTZ* TTW*  

popd
