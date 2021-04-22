#!/bin/bash


PROJ=vdmodel_4node_s01
ROOTPROJDIR=/mnt/raid6_data/user/aranyics/rdcm2/ReDCM-nosge
PROJDIR=${ROOTPROJDIR}/${PROJ}

SCRIPTDIR=${ROOTPROJDIR}

DCMDIR=${PROJDIR}/dcm

#echo "Watcher for rdcm jobs (refreshes every 8 hours)"
#echo "$PROJDIR"

MODELS=`cat ${ROOTPROJDIR}/A.csv | wc -l`

#echo "Models: $MODELS"


CURRENT_TIME=`echo $SECONDS`
COMPUTED=`find ${DCMDIR} -maxdepth 2 -type f | wc -l`

while [ 1 ]; do

    DONE=$(($COMPUTED - 1))
    START_TIME=${CURRENT_TIME}

    #sleep for 2 hours
    #sleep 28800
    sleep 7200

    CURRENT_TIME=`echo $SECONDS`
    COMPUTED=`find ${PROJDIR}/dcm -maxdepth 2 -type f | wc -l`

    ELAPSED_TIME=$(($CURRENT_TIME - $START_TIME))
    COMP_RATE=$(($COMPUTED - $DONE))
    REMAINED=$(($MODELS - $COMPUTED))
    ESTIMATED_RUNTIME=$(($REMAINED / $COMP_RATE * $ELAPSED_TIME))

    #printf "Estimated runtime (of $REMAINED models): $(($ESTIMATED_RUNTIME/60)) min $(($ESTIMATED_RUNTIME%60)) sec ($COMPUTED / $MODELS done)\tJobs: `qstat -f | grep "aranyics     r" | wc -l`\n"
    printf "Estimated runtime (of $REMAINED models): $(($ESTIMATED_RUNTIME/60)) min $(($ESTIMATED_RUNTIME%60)) sec ($COMPUTED / $MODELS done)\tJobs: `qstat -f | grep "aranyics     r" | wc -l`\n" | mail -s "rDCM - Estimated runtime (automatic mail)" aranyics11@gmail.com

done

exit
