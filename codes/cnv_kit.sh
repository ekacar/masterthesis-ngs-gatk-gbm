
#!/bin/bash

cnv_kit="/home/kacar/apps/cnvkit-master/cnvlib"
tumor="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/merged_reads"
genomic_regions="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED"
ref="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta"
resources="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources"
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"


wget https://github.com/etal/cnvkit/blob/master/data/access-10kb.hg38.bed

samples=("COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")

for sample in "${samples[@]}"; do
  cnvkit.py batch ${tumor}/*dedup_recal.bam -n -t ${genomic_regions} -f ${ref} \
     --access ${resources}/access-10kb.hg38.bed \
     --output-reference ${sample}.cnn -d ${results}/cna_calling/cnvkit
done
