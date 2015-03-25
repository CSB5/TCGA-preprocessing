#!/bin/bash

CLINICAL=$1
TISSUE=$2
PREFIX=${TISSUE// /_}

echo "CLINICAL:$CLINICAL"
echo "TISSUE:$TISSUE"
echo "PREFIX:$PREFIX"

grep "$TISSUE" $CLINICAL | cut -f 1 | sort | uniq > ${PREFIX}-TCGA_IDs.txt
