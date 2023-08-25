#!/bin/bash

home="/data/Unit_LMM/selberherr-group/kacar"
qc="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/qc"
ref="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta"
resources="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources"
aligned_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/aligned_reads"
reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/fastq"
merged_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/merged_reads"
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"


samples=("COMP5_1_S126518" "COMP5_2_S126517" "COMP5_5_S126520" "COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")

#gatk IndexFeatureFile \
#     -I "${resources}/panel_of_normals/1000g_pon.hg38.vcf.gz"

#Somatic Variant Calling with Mutect2
#for sample in "${samples[@]}"; do
#    gatk Mutect2 \
#        -R ${ref} \
#        -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam" \
#        --germline-resource "${resources}/gnomAD/af-only-gnomad.hg38.vcf.gz" \
#        --panel-of-normals "${resources}/panel_of_normals/1000g_pon.hg38.vcf.gz" \
#        --f1r2-tar-gz "${results}/BSF_1274_H3CTKDSX5#${sample}_f1r2.tar.gz"\
#        -O "${results}/BSF_1274_H3CTKDSX5#${sample}_vcf.gz" \
#        -bamout "${merged_reads}/BSF_1274_H3CTKDSX5#${sample}_assembledhaplotypes.bam"
#done

#Get the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter
#for sample in "${samples[@]}"; do
#   gatk LearnReadOrientationModel\
#      -I "${results}/BSF_1274_H3CTKDSX5#${sample}_f1r2.tar.gz" \
#      -O "${results}/BSF_1274_H3CTKDSX5#${sample}_read_orientation_model.tar.gz"
#done

##Estimate cross-sample contamination

#GetPileupSummaries to summarize read support for a set number of known variant sites.
#for sample in "${samples[@]}"; do
# gatk GetPileupSummaries \
#     -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam" \
#     -V ${resources}/small_exac_common_3.hg38.vcf.gz \
#     -L ${resources}/small_exac_common_3.hg38.vcf.gz \
#     -O "${results}/BSF_1274_H3CTKDSX5#${sample}_pileupsummaries.table"
#done

##Estimate contamination with CalculateContamination.The tool gives the fraction contamination. This estimation informs downstream filtering by FilterMutectCalls.
#for sample in "${samples[@]}"; do
# gatk CalculateContamination \
#    -I "${results}/BSF_1274_H3CTKDSX5#${sample}_pileupsummaries.table" \
#    -O "${results}/BSF_1274_H3CTKDSX5#${sample}_contamination.table"
#done

#Filter somatic SNVs and indels called by Mutect2
#for sample in "${samples[@]}"; do
#  gatk FilterMutectCalls \
#     -R ${ref} \
#     -V "${results}/BSF_1274_H3CTKDSX5#${sample}_vcf.gz" \
#     --contamination-table "${results}/BSF_1274_H3CTKDSX5#${sample}_contamination.table"\
#     --stats "${results}/BSF_1274_H3CTKDSX5#${sample}_vcf.gz.stats" \
#     --orientation-bias-artifact-priors "${results}/BSF_1274_H3CTKDSX5#${sample}_read_orientation_model.tar.gz" \
#    -O "${results}/BSF_1274_H3CTKDSX5#${sample}_filtered_vcf.gz" 
#done

#Select a subset of variants from a VCF file by excluding all the filtered vcfs
#for sample in "${samples[@]}"; do
# gatk SelectVariants \
#    -R ${ref} \
#    -V "${results}/BSF_1274_H3CTKDSX5#${sample}_filtered_vcf.gz" \
#    --exclude-filtered \
#    -OVI \
#    -O "${results}/BSF_1274_H3CTKDSX5#${sample}_final_vcf.gz"
#done 


#Functional Annotation with Funcotator
for sample in "${samples[@]}"; do
  gatk Funcotator\
    -V "${results}/BSF_1274_H3CTKDSX5#${sample}_final_vcf.gz" \
    -R ${ref} \
    --ref-version hg38 \
    --data-sources-path "${resources}/funcotator_dataSources.v1.6.20190124s" \
    --output "${results}/BSF_1274_H3CTKDSX5#${sample}_final1.maf" \
    --output-file-format MAF \
    --annotation-default tumor_barcode:${sample}
done
