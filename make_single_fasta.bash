#!/bin/bash

echo "$0 $1 $2 $3"
echo
echo

if [[ "$1" != "" ]]; then
    BASEDIR=$1
else
    BASEDIR=.
fi

if [[ "$2" != "" ]]; then
    OUTFILENAME=$2
else
    echo "Usage: $0 parent_dir output_file_name"
    exit
fi

FILEEXT=""

if [[ "$3" != "" ]]; then
    FILEEXT=$3 # this should be escaped by the user if necessary, i.e. specify "\.fna" not ".fna"
    NAMECLAUSE="-name \*$FILEEXT"
fi

if [[ -e $OUTFILENAME ]]; then
    echo "FATAL: file $OUTFILENAME already exists; delete it first"
    exit
fi

FIND="find $BASEDIR -type f $NAMECLAUSE -exec cat {} >> '$OUTFILENAME' \;"

echo $FIND
eval $FIND




