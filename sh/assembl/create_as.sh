#!/bin/sh

TXT_FILE='sub_a.txt'
RSCRIPT='/scratch/ob219/bs_sim/blockshifter_simulations/R/assemble_sim_gwas.R'
OUT_DIR_STUB='/scratch/ob219/bs_sim/'
BLOCK_NUMBER=50
OUT_DIR="${OUT_DIR_STUB}${BLOCK_NUMBER}_pir_nopir/"
if [ ! -e $OUT_DIR ]
then
	echo "CREATING OUT DIR $OUT_DIR";
	`mkdir $OUT_DIR`;
fi
[ -e $TXT_FILE ] && rm $TXT_FILE;
for j in {1..10}
do
	if [ ! -e "${OUT_DIR}${j}-${i}.RData" ]
	then
		echo "Rscript --vanilla $RSCRIPT  -b $BLOCK_NUMBER -o $OUT_DIR -i $j" >> $TXT_FILE;
	fi
done
