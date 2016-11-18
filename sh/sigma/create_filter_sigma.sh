#!/bin/sh

TXT_FILE='create_sigma.txt'
DATA_DIR='/scratch/ob219/bs_sim/support/chr1_EUR_LD/'
RSCRIPT='/scratch/ob219/bs_sim/blockshifter_simulations/R/create_filtered_sigma.R'
OUT_DIR='/scratch/ob219/bs_sim/sigma_pir/'
[ -e $TXT_FILE ] && rm $TXT_FILE;
for i in `\ls $DATA_DIR`;
do
	ext="${i##*.}"
	if [ $ext = 'RData' ]  && [ ! -e "${OUT_DIR}$i" ]
	 then
		echo "Rscript --vanilla $RSCRIPT  -f $DATA_DIR/$i -o $OUT_DIR" >> create_sigma.txt;
	fi
done
