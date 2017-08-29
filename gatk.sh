#!/bin/bash

# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

while [[ $# -gt 1 ]]
do 
key="$1"

case $key in
	-o|--output-dir)
	OUTPUT_DIR="$2"
	shift # past argument
	;;
	-f|--fwd)
	FWD="$2"
	shift # past argument
	;;
	-r|--rev)
	REV="$2"
	shift # past argument
	;;
	-c|--condition)
	CONDITION="$2"
	shift # past argument
	;;
	*)
		# unknown option
	;;
esac
shift # past argument or value
done

# https://gatkforums.broadinstitute.org/gatk/discussion/2799
# (howto) Map and mark duplicates

# Map reads with bwa
bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' $GENOME_FASTA $FWD $REV > ${OUTPUT_DIR}/${CONDITION}/aligned_reads.sam

# Convert to BAM, sort and mark duplicates
java -jar $picard SortSam INPUT=${OUTPUT_DIR}/${CONDITION}/aligned_reads.sam OUTPUT=${OUTPUT_DIR}/${CONDITION}/sorted_reads.bam SORT_ORDER=coordinate
java -jar $picard MarkDuplicates INPUT=${OUTPUT_DIR}/${CONDITION}/sorted_reads.bam OUTPUT=${OUTPUT_DIR}/${CONDITION}/dedup_reads.bam METRICS_FILE=${OUTPUT_DIR}/${CONDITION}/metrics.txt
java -jar $picard BuildBamIndex INPUT=${OUTPUT_DIR}/${CONDITION}/dedup_reads.bam

# https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#2.2
# Call variants with HaplotypeCaller
java -jar $GATK -T HaplotypeCaller -R $GENOME_FASTA -I ${OUTPUT_DIR}/${CONDITION}/dedup_reads.bam -o ${OUTPUT_DIR}/${CONDITION}/raw_variants.vcf

# https://software.broadinstitute.org/gatk/documentation/article?id=2806
# (howto) Apply hard filters to a call set
java -jar $GATK -T SelectVariants -R $GENOME_FASTA -V ${OUTPUT_DIR}/${CONDITION}/raw_variants.vcf -selectType SNP -o ${OUTPUT_DIR}/${CONDITION}/raw_snps.vcf
java -jar $GATK -T VariantFiltration -R $GENOME_FASTA -V ${OUTPUT_DIR}/${CONDITION}/raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o ${OUTPUT_DIR}/${CONDITION}/filtered_snps.vcf
java -jar $GATK -T SelectVariants -R $GENOME_FASTA -V ${OUTPUT_DIR}/${CONDITION}/raw_variants.vcf -selectType INDEL -o ${OUTPUT_DIR}/${CONDITION}/raw_indels.vcf
java -jar $GATK -T VariantFiltration -R $GENOME_FASTA -V ${OUTPUT_DIR}/${CONDITION}/raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o ${OUTPUT_DIR}/${CONDITION}/filtered_indels.vcf

# http://snpeff.sourceforge.net/SnpEff_manual.html#run
# Annotate using snpEff
java -Xmx4g -jar $snpEff Arabidopsis_thaliana ${OUTPUT_DIR}/${CONDITION}/filtered_snps.vcf -stats ${OUTPUT_DIR}/${CONDITION}/snpEff_summary.snps.html > ${OUTPUT_DIR}/${CONDITION}/filtered_snps.ann.vcf
java -Xmx4g -jar $snpEff Arabidopsis_thaliana ${OUTPUT_DIR}/${CONDITION}/filtered_indels.vcf -stats ${OUTPUT_DIR}/${CONDITION}/snpEff_summary.indels.html > ${OUTPUT_DIR}/${CONDITION}/filtered_indels.ann.vcf
