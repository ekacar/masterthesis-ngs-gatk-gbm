#!/bin/bash

#Shotcuts for my data and other resources
home="/data/Unit_LMM/selberherr-group/kacar"
qc="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/qc"
ref="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta"
resources="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources"
aligned_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/aligned_reads"
reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/fastq"
merged_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/merged_reads"
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"

#NGS Sequencing was undertaken for each subregion and patient sample on 4 lanes
lanes=("1" "2" "3" "4")
samples=("COMP5_1_S126518" "COMP5_2_S126517" "COMP5_5_S126520" "COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")

#Index the reference genome using bwa index
bwa index ${ref}

#Turn ubam to Fastq file as preparation for bwa alignment, then use the output for bwa alignment and use MergeBamAlignment to make sam to bam file
for lane in "${lanes[@]}"; do
  for sample in "${samples[@]}"; do
    gatk SamToFastq \
       I="${home}/rawdata/L${lane}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.bam" \
       F="${reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}_r1.fastq" \
       F2="${reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}_r2.fastq"
  done
done

#Alignment with reference using bwa mem
for lane in "${lanes[@]}"; do
  for sample in "${samples[@]}"; do
     bwa mem -t 8 "${ref}" "${reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}_r1.fastq" "${reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}_r2.fastq" > "${aligned_reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.sam"
  done
done

#MergeBamAlignment needs the reference file as a dictionary
gatk CreateSequenceDictionary \
      R=/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta \
      O=/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.dict


#From SAM to BAM, applying the info from unaligned BAM to aligned BAM at the end!
for lane in "${lanes[@]}"; do
  for sample in "${samples[@]}"; do
    gatk MergeBamAlignment \
      ALIGNED="${aligned_reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.sam" \
      UNMAPPED="${home}/rawdata/L${lane}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.bam" \
      O="${aligned_reads}/BSF_1274_H3CTKDSX5_${lane}%23${sample}.bam" \
      R="${ref}"
  done
done

#Merge all bam files from one sample
for sample in "${samples[@]}"; do
  gatk MergeSamFiles \
    -I "${aligned_reads}/BSF_1274_H3CTKDSX5_1%23${sample}.bam" \
    -I "${aligned_reads}/BSF_1274_H3CTKDSX5_2%23${sample}.bam" \
    -I "${aligned_reads}/BSF_1274_H3CTKDSX5_3%23${sample}.bam" \
    -I "${aligned_reads}/BSF_1274_H3CTKDSX5_4%23${sample}.bam" \
    -O "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_merged.bam"
done

#Sort aligned bam files before marking the duplicates
for sample in "${samples[@]}"; do
    gatk SortSam \
       -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_merged.bam" \
       -O "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_merged_sorted.bam" \
       --SORT_ORDER coordinate
done

#Mark Duplicates in my aligned bam files
for sample in "${samples[@]}"; do
    gatk MarkDuplicates \
      -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_merged_sorted.bam" \
      -O "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup.bam" \
      -M "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_marked_dup_metrics.txt"
done


#Get indexed feature files for Vcf files from databases
gatk IndexFeatureFile \
    -I "${resources}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" 

gatk IndexFeatureFile\
    -I "${resources}/dbsnp_146.hg38.vcf.gz"

gatk IndexFeatureFile\
    -I "${resources}/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz" 



#Recalibrate base quality scores using GATK:
for lane in "${lanes[@]}"; do
  for sample in "${samples[@]}"; do
    gatk BaseRecalibrator \
      -R "${ref}" \
      -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup.bam" \
      --known-sites "${resources}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
      --known-sites "${resources}/dbsnp_146.hg38.vcf.gz" \
      --known-sites "${resources}/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz"\
      -O "${merged_reads}/BSF_1274_H3CTKDSX5#${sample}_dedup_recal_data.table"
  done
done

#Applying recalibration on the bam files
for sample in "${samples[@]}"; do
    gatk ApplyBQSR \
      -R "${ref}" \
      -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup.bam" \
      -bqsr "${merged_reads}/BSF_1274_H3CTKDSX5#${sample}_dedup_recal_data.table" \
      -O "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam"
done

