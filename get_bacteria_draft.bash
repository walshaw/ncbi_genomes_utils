#!/bin/bash

DIRFILE=Bacteria_DRAFT__dir_20150429.lst # should be pathless dir names

#URL_PREFIX=ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT/
REMOTE_HOST=ftp.ncbi.nlm.nih.gov
REMOTE_PATH=genomes/Bacteria_DRAFT/
REMOTE_SUFFIX=.fna.tgz

if [[ "$1" != "" ]]; then
    DIRFILE=$1
fi

if [[ "$2" != "" ]]; then
    REMOTE_HOST=$2
fi

if [[ "$3" != "" ]]; then
    REMOTE_PATH=$3
fi

if [[ "$4" != "" ]]; then
    REMOTE_SUFFIX=$4
fi



for DIR in $( cat $DIRFILE ); do

#    echo $DIR
    MKDIR="mkdir $DIR"
    NCFTPGET="ncftpget $REMOTE_HOST $DIR ${REMOTE_PATH}$DIR/*$REMOTE_SUFFIX"

    if [[ -d $MKDIR ]]; then
        echo "dir $DIR already exists"
    else
        echo $MKDIR
        eval $MKDIR
    fi

    echo $NCFTPGET
    eval $NCFTPGET

done


