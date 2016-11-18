#!/bin/sh

TXT_FILE='missing_files.txt'
DATA_DIR='/scratch/ob219/bs_sim/sigma_pir/'
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
				echo "${OUT_DIR}${j}_${i}" >> $TXT_FILE;
			fi
		done
	fi
done
