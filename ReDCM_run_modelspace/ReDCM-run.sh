#!/bin/bash

SUBJECT=${1:-}
PROJ=vdmodel_4node_${SUBJECT}
ROOTPROJDIR="$( cd $( dirname ${BASH_SOURCE[0]} ) >/dev/null 2>&1 && pwd )"
#PROJDIR=${ROOTPROJDIR}/${PROJ}
DATADIR=${ROOTPROJDIR}/../../../data
RESULTSDIR=${ROOTPROJDIR}/../../../results/modelspace

SCRIPTDIR=${ROOTPROJDIR}
RDCMEXEC=${SCRIPTDIR}/rdcm_wrapper.R
MERGEEXEC=${SCRIPTDIR}/rdcm_merge_estimates.R

DCMTEMP=${DATADIR}/DCM/DCM_${SUBJECT}.mat
MXDIR=${RESULTSDIR}/mx
AMX=${DATADIR}/A.csv
BMX=${DATADIR}/B.csv

START=25001
MODELS=26000
#MODELS=`cat $AMX | wc -l`
BATCH=1000

echo "Models: $MODELS"
echo "Batch size: $BATCH"

mkdir -p $MXDIR

k=1


for j in `seq $START $BATCH $MODELS`
do

    DCMDIR=${RESULTSDIR}/${PROJ}/dcm/out_${j}
    mkdir -p $DCMDIR
    chmod 777 -R $DCMDIR

    SUBSPACE=$(( j + BATCH - 1 ))
    AMXCURR=${MXDIR}/A_${j}.csv
    BMXCURR=${MXDIR}/B_${j}.csv
    cat ${AMX} | awk 'NR>='${j}' && NR<='${SUBSPACE} > ${AMXCURR}
    cat ${BMX} | awk 'NR>='${j}' && NR<='${SUBSPACE} > ${BMXCURR}

    SIZE=`cat $AMXCURR | wc -l`

    echo "Start batch: jobs $j - $(( j + SIZE - 1 )) "
    #echo "Start batch: jobs $j - $(( j + SIZE - 1 )) " | mail -s "rDCM - Start new batch (automatic mail)" -c "emri@pet.dote.hu" aranyics11@gmail.com


    for i in `seq 1 10 $SIZE`
    do

	SUBSUBSPACE=$(( i + 10 - 1 ))
	CURR=$(( i + j - 1 ))
	CURR2=$(( i + j + 10 - 2 ))
	OFFSET=$(( j - 1 ))

	echo "$RDCMEXEC $DCMTEMP -a $AMXCURR -b $BMXCURR -x $i -y $SUBSUBSPACE -n $OFFSET -o $DCMDIR"
	nohup $RDCMEXEC $DCMTEMP -a $AMXCURR -b $BMXCURR -x $i -y $SUBSUBSPACE -n $OFFSET -o $DCMDIR 2>&1 </dev/null &

        sleep 5
	#exit

	while [ `ps -eaf | grep "/usr/lib/R/bin/exec/R" | wc -l` -gt 1 ]; do
	    echo "sleep 1m"
	    #sleep 60
	    sleep 5
	done

    done

done

exit
