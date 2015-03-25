#!/bin/bash

mkdir -p $4
cd $4

table=$1
QSUB="qsub -pe OpenMP 1 -l h_rt=2:0:0,h_vmem=10G -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH"

#for i in {2..375}
tumor_s=`expr $2 + 1`
echo "$tumor_s"
#for i in {$tumor_s..$3}
#for i in {9..426}
for (( i=$tumor_s; i<=$3; i++ ))
do
	outPrefix=`head -n 1 $table | cut -f $i`
	echo "head -n 1 $table | cut -f $i"
#	echo "PREFIX $outPrefix"
	head -n 1 $table | cut -f 1-$2,$i > ${outPrefix}.pool.dat
	echo "head -n 1 $table | cut -f 1-$2,$i > ${outPrefix}.pool.dat"
	j=`expr $i + 1`
	echo "tail -n+2 $table | cut -f 1-$tumor_s,$j >> ${outPrefix}.pool.dat"
	tail -n+2 $table | cut -f 1-$tumor_s,$j >> ${outPrefix}.pool.dat
	
	echo "$QSUB /home/chiakhb/SCRIPTS/TCGA_PREPROCESSING/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_DEseq DESeq"
	$QSUB /home/chiakhb/SCRIPTS/TCGA_PREPROCESSING/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_DEseq DESeq

	echo "$QSUB /home/chiakhb/SCRIPTS/TCGA_PREPROCESSING/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_foldChange foldChange"
	$QSUB /home/chiakhb/SCRIPTS/TCGA_PREPROCESSING/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_foldChange foldChange


done

