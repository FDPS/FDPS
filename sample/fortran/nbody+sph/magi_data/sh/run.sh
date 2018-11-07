#!/bin/sh
###############################################################


###############################################################
# global configurations
###############################################################
MAGI_INSTALL_DIR=/work4a/local/MAGI/ymiki-magi-0944b759e2dc/build/bin
EXEC=${MAGI_INSTALL_DIR}/magi
###############################################################
###############################################################
# number of N-body particles
if [ -z "$NTOT" ]; then
    NTOT=2097152
    # NTOT=4194304
    # NTOT=8388608
    # NTOT=16777216
fi
###############################################################
# dump file generation interval in the N-body simulation in units of minute (required for GOTHIC)
if [ -z "$SAVE" ]; then
    SAVE=140.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# FILE: the name of the simulation
# CONFIG: cfg/$CONFIG specifies the physical properties of the initial condition
# EPS: Plummer softening length in units of astrophysical unit system
# ETA: safety parameter to control time stepping (required for GOTHIC)
# FINISH: final time of the N-body simulation in units of astrophysical unit system
# INTERVAL: snapshot interval in the N-body simulation in units of astrophysical unit system (required for GOTHIC)
###############################################################
FILE=Galaxy
CONFIG=Galaxy.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1175.0
INTERVAL=25.0
###############################################################
# set input arguments
OPTION="-file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################


###############################################################
# job execution via command line
###############################################################
# make necessary directories
DOC_DIR=doc # some documents are output here
DAT_DIR=dat # particle data is output here
LOG_DIR=log # runtime log is output here
mkdir -p $DOC_DIR $DAT_DIR $LOG_DIR
# set stdout and stderr
JOB_ID=$$
LOG=${LOG_DIR}/$FILE.l
STDOUT=${LOG_DIR}/$FILE.$JOB_ID.out
STDERR=${LOG_DIR}/$FILE.$JOB_ID.err
###############################################################
# start logging
TIME=`date`
echo "start: $TIME" >> $LOG
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" >> $LOG
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" >> $LOG
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" >> $LOG
###############################################################
