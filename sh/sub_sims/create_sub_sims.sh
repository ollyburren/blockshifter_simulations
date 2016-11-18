#!/bin/sh

TXT_FILE='sub_sims.txt'
DATA_DIR='/scratch/ob219/bs_sim/sigma_pir/'
RSCRIPT='/scratch/ob219/bs_sim/blockshifter_simulations/R/sim_gwas.R'
OUT_DIR='/scratch/ob219/bs_sim/sims/'
[ -e $TXT_FILE ] && rm $TXT_FILE;
for i in `\ls $DATA_DIR`;
do
	ext="${i##*.}"
	if [ $ext = 'RData' ]
	 then
		for j in {1..10}
		do
			if [ ! -e "${OUT_DIR}${j}_${i}" ]
			then
				echo "Rscript --vanilla $RSCRIPT  -f $DATA_DIR/$i -o $OUT_DIR -i $j" >> $TXT_FILE;
			fi
		done
	fi
done
