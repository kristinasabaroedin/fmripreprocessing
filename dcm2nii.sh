#!/bin/bash

timestamp() {
  date +"%T"
}

# Static directories
# Dropbox
SCRIPTDIR=~/Dropbox/scripts/projects/OCDPG
SUBJIDS=$(<$SCRIPTDIR/sublists/scratch.txt)
# HDD
DISK=/media/lindenmp/SG8_4TB
PROJDIR=$DISK/Research_Projects/OCDPG
DATADIR=$PROJDIR/data


# Module toggles (on/off)
	MODULE1=1 #mrconvert

for SUBJ in $SUBJIDS; do
	# Dynamic directories
	DICOMDIR=$DATADIR/$SUBJ/dicoms
	RESTDIR=$DATADIR/$SUBJ/rfMRI
	T1DIR=$DATADIR/$SUBJ/t1
	# clean outputs and reinitialise
	if [ ! -d $RESTDIR ]; then mkdir $RESTDIR; else echo "Cleaning rfMRI outputs"; rm -rf $RESTDIR; mkdir $RESTDIR; fi
	if [ ! -d $T1DIR ]; then mkdir $T1DIR; else echo "Cleaning t1 outputs"; rm -rf $T1DIR; mkdir $T1DIR; fi

	################################ MODULE 1: MRtrix mrconvert #######################################
	if [ $MODULE1 = "1" ]; then
		echo -e "\nRunning MODULE 1: MRtrix mrconvert: $SUBJ \n"
		
		# fMRI epi
		mrconvert $DICOMDIR/rfMRI $RESTDIR/epi.nii

		# t1
		mrconvert $DICOMDIR/t1 $T1DIR/t1.nii
		
		echo -e "\nFinished MODULE 1: MRtrix mrconvert: $SUBJ \n"
	else
		echo -e "\nSkipping MODULE 1: MRtrix mrconvert: $SUBJ \n"
	fi
	###################################################################################################

done