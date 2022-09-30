#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 2
#$ -P quanah


source ~/conda/etc/profile.d/conda.sh
conda activate selection


#This script is to calculate the log likeligood ratio change from the null and alternative model of CODEML 
mkdir likelihood_values/

#Extract likelihood values from all the genes

for d in OG0*/; do PX=$(basename $d | awk -F '_' '{print $1}'); (cd $d && grep "lnL" */*model | sed 's/(.[^)]*)/\t/g' | sed 's/.*\///g' | sed 's/:[^:]*://g' | sed 's/_model/\t/g' | awk '{print $1, $2}' | awk 'ORS=NR%2?"\t":"\n"' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $PX'_lnl' && mv $PX'_lnl' ../likelihood_values/); done

#Concatenate values
cd likelihood_values/
ln -s ../lrt_construction.py

cat *_lnl >> all_lnl_values

rm *_lnl

python3 lrt_construction.py

#Calculate LRT
sed -i 's/_alternative//g' LRT_results

#Calculate P value with Chi2

awk '{print $1","$5}' LRT_results | sed '1d' > lrt_temp; less lrt_temp | while IFS=, read OG LRT; do CHI2=$(chi2 1 $LRT); echo "$OG,$LRT,$CHI2" >> chi2_temp; done

#Create final table with OG, LRT ahd CHI2 value:

sed ':a;N;$!ba;s/,\n/,/g' chi2_temp > probability_values_temp
sed -i -e 's/,/\t/g;s/df//g;s/prob//g;s/ch2//g;s/=/\t/g' probability_values_temp

sed '1 i\Orthogroup\tLRT\tdf\tprob\tprob exp' probability_values_temp > codeml_prob_results.txt

rm *_temp
