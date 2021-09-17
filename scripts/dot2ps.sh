#!/bin/sh

if [ $# -ne 1 ] ; then
    echo "Usage: > $0 dot_file"
    exit
fi

# Requires graphviz (dot executable)

ofn=`echo $1 | sed -e s/\.dot$/\.ps/`
dot -Tps $1 -o ${ofn}
