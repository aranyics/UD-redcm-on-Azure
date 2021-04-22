#!/bin/bash

SUBJECT=s01
PROJ=vdmodel_4node_${SUBJECT}
ROOTPROJDIR="$( cd $( dirname ${BASH_SOURCE[0]} ) >/dev/null 2>&1 && pwd )"
#PROJDIR=${ROOTPROJDIR}/${PROJ}
DATADIR=${ROOTPROJDIR}/../../../data
RESULTSDIR=${ROOTPROJDIR}/../../../results/modelspace/${PROJ}
OUTDIR=${RESULTSDIR}/dcm_ms

SCRIPTDIR=${ROOTPROJDIR}
RDCMEXEC=${SCRIPTDIR}/rdcm_wrapper.R
MERGEEXEC=${SCRIPTDIR}/rdcm_merge_estimates.R

AMX=${DATADIR}/A.csv
BMX=${DATADIR}/B.csv

START=1
MODELS=`cat $AMX | wc -l`
BATCH=1000

echo "Models: $MODELS"
echo "Batch size: $BATCH"

mkdir -p $OUTDIR

k=1


for j in `seq $START $BATCH $MODELS`
do

    DCMDIR=${RESULTSDIR}/dcm/out_${j}

    SUBSPACE=$(( j + BATCH - 1 ))

    SIZE=${BATCH}


    for i in `seq 1 10 $SIZE`
    do

	SUBSUBSPACE=$(( i + 10 - 1 ))

	echo "$MERGEEXEC $DCMDIR -x $i -y $SUBSUBSPACE -o $OUTDIR"
	$MERGEEXEC $DCMDIR -x $i -y $SUBSUBSPACE -o $OUTDIR 2>&1 </dev/null

	exit

    done

done

exit
