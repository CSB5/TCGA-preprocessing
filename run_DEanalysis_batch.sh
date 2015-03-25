#!/bin/bash
table=$1
#table="../RNA_SEQ_COUNT_MATRIX_8normals.dat"
QSUB="qsub -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH"

#for i in {2..375}
tumor_s=`expr $2 + 1`
echo "$tumor_s"

ana_type=$4
out_dir=$5

#for i in {$tumor_s..$3}
#for i in {9..426}
for (( i=$tumor_s; i<=$3; i++ ))
do
	outPrefix=`head -n 1 $table | cut -f $i`
	file_name=$out_dir${outPrefix}.${ana_type}.dat
	echo " **** $file_name"
	echo "head -n 1 $table | cut -f $i"
#	echo "PREFIX $outPrefix"
	head -n 1 $table | cut -f 1-$2,$i > $file_name
	echo "head -n 1 $table | cut -f 1-$2,$i > $file_name"
	j=`expr $i + 1`
	echo "tail -n+2 $table | cut -f 1-$tumor_s,$j >> $file_name"
	tail -n+2 $table | cut -f 1-$tumor_s,$j >> $file_name
	
	if [ "$ana_type" = "DE_pool" ]; then
	    echo "$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R $file_name $2 res DESeq"
	    $QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R $file_name $2 res DESeq
	else
	    echo "$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE.R $file_name $2 res DESeq"
	    $QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE.R $file_name $2 res DESeq
	fi

	#echo "$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_foldChange foldChange"
	#$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_foldChange foldChange
	
	#echo "$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_DEseq DESeq"
	#$QSUB /home/chiakhb/SCRIPTS/oncoIMPACT/call_DE_pool.R ${outPrefix}.pool.dat $2 DE_result_cqn_DEseq DESeq

	
done

