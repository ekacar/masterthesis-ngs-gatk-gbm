#!/bin/bash


qc="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/qc"
ref="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta"
known_sites="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/dbsnp_146.hg38.vcf.gz"
aligned_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/aligned_reads"
reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/fastq"
merged_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/merged_reads"
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"
genomic_regions="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED"



#lanes=("1" "2" "3" "4")
samples=("COMP5_1_S126518" "COMP5_2_S126517" "COMP5_5_S126520" "COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")

#Using Fastqc, we check the quality of individual bam files belonging to different samples from our patients in the conda environment qc
#for lane in "${lanes[@]}"; do
#  for sample in "${samples[@]}"; do
#      fastqc "/data/Unit_LMM/selberherr-group/kacar/rawdata/L${lane}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.bam" -o ${qc} 
#  done
#done


#To use mosdepth, we need to index our bam files
#for sample in "${samples[@]}"; do
# gatk BuildBamIndex \
#      I="${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam"
#      O="${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bai"
#done

#Calculating the coverage in case of WES of genomic regions 
for sample in "${samples[@]}"; do
       mosdepth --by ${genomic_regions} -t 4 ${qc}/coverage/${sample}  ${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam
done


