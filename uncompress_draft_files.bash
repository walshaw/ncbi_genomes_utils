#!/bin/bash

echo "$0 $1 $2"
echo
echo

if [[ "$1" != "" ]]; then
    BASEDIR=$1
else
    BASEDIR=.
fi

if [[ "$2" != "" ]]; then
    DIRROOTNAME=$2
else
    DIRROOTNAME=bacteria_draft
fi

DRAFT_DIR_COMPRESSED=$BASEDIR/compressed/$DIRROOTNAME
DRAFT_DIR_UNCOMPRESSED=$BASEDIR/$DIRROOTNAME

for GENOME_DIR in $(ls $DRAFT_DIR_COMPRESSED); do

    COMPRESSED_FILES=$(ls $DRAFT_DIR_COMPRESSED/$GENOME_DIR)
    for FILE in $COMPRESSED_FILES; do
        FILE_FULLPATH=$DRAFT_DIR_COMPRESSED/$GENOME_DIR/$FILE

        if [[ "$FILE" != *.tgz ]]; then
             echo "unexpected file name: $FILE"
        fi

        UNCDIR="$DRAFT_DIR_UNCOMPRESSED/$GENOME_DIR"

        UNC_FILE_FULLPATH=$UNCDIR/$UNCOMPRESSED_FILE
        echo "$FILE_FULLPATH -> $UNC_FILE_FULLPATH"

        if [[ ! -d $UNCDIR ]]; then
            MKDIR1="mkdir $UNCDIR"
            echo $MKDIR1
            eval $MKDIR1
        fi

#        CD1="cd $UNCDIR"
        TAR="tar -C $UNCDIR -zvxf $FILE_FULLPATH"
#        CD2="cd -"
#        echo $CD1
        #eval $CD1
        echo $TAR
        eval $TAR
#        echo $CD2
        #eval $CD2

    done

done


