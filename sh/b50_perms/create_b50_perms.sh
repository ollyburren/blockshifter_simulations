#!/bin/sh

TXT_FILE='sub_b50_bs.txt'
RSCRIPT='/scratch/ob219/bs_sim/blockshifter_simulations/R/blockshifter_sim.R'
OUT_DIR_STUB='/scratch/ob219/bs_sim/'
BLOCK_NUMBER=50
OUT_DIR="${OUT_DIR_STUB}${BLOCK_NUMBER}_pir_nopir_bs/"
DATA_DIR="${OUT_DIR_STUB}${BLOCK_NUMBER}_pir_nopir/"
if [ ! -e $OUT_DIR ]
then
	echo "CREATING OUT DIR $OUT_DIR";
	`mkdir $OUT_DIR`;
fi
[ -e $TXT_FILE ] && rm $TXT_FILE;
for i in `\ls ${DATA_DIR}*.RData` 
do
	if [ ! -e "${DATA_DIR}${i}" ]
	then
		echo "/home/ob219/git/slurmer/qRscript.rb -r $RSCRIPT -p $i -o $OUT_DIR"
	fi
done
