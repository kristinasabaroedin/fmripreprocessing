#!/bin/bash

timestamp() {
  date +"%T"
}

# Static directories
# Dropbox
SCRIPTDIR=~/Dropbox/scripts/projects/OCDPG
SUBJIDS=$(<$SCRIPTDIR/sublists/SubjectIDs.txt)

INDIR=/Volumes/SG8_4TB/Research_Projects/OCDPG/data
OUTDIR=/Volumes/SG8_4TB/Research_Projects/OCDPG/motion_reports 

if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; else echo "Cleaning outputs"; rm -rf $OUTDIR; mkdir $OUTDIR; fi

for SUBJ in $SUBJIDS; do
	echo "copying $SUBJ"
	cp $INDIR/$SUBJ/rfMRI/prepro_report_motion*.pdf $OUTDIR/$SUBJ.pdf
done