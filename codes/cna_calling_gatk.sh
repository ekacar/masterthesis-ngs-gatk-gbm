#!/bin/bash
ref="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/Homo_sapiens_assembly38.fasta"
resources="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources"
merged_reads="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/merged_reads"
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"
genomic_regions="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-ssv-tumoronly/resources/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED"

samples=("COMP5_1_S126518" "COMP5_2_S126517" "COMP5_5_S126520" "COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")

#Converting the whole exome BED file to a Picard Interval List. 
#gatk BedToIntervalList \
#      I="${genomic_regions}" \
#      O="${resources}/hg38_Twist_ILMN_Exome_2.0_Plus_Panel.interval_list" \
#      SD="${ref}"

#Preprocess the whole exome interval list within the reference.fasta
#gatk PreprocessIntervals \
#    -L "${resources}/hg38_Twist_ILMN_Exome_2.0_Plus_Panel.interval_list"\
#    -R "${ref}" \
#    --bin-length 0 \
#    --interval-merging-rule OVERLAPPING_ONLY \
#    -O "${resources}/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_preprocessed.interval_list"

#Create a fasta index file for CollectReadCounts
#samtools faidx ${ref}

#Collecting read counts at specified intervals
#for sample in "${samples[@]}"; do
#   gatk CollectReadCounts \
#     -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam"\
#     -L "${resources}/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_preprocessed.interval_list"\
#     -R "${ref}"\
#     --format HDF5\
#     -imr OVERLAPPING_ONLY\
#     -O "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}.read_counts.hdf5"
#done
########This "self-made read counts file didn't match well with the Pon read counts from Broad Institute.
########Therefore, I will continue with the read counts file from Broad Institute Server.

##At first, we preprocess our whole exome interval list in our reference.fasta
#gatk PreprocessIntervals \
#    -L ${resources}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list\
#    -R ${ref} \
#    --bin-length 0 \
#    --interval-merging-rule OVERLAPPING_ONLY \
#    -O ${resources}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_preprocessed_targets.interval_list##

#Collecting read counts at specified intervals
#for sample in "${samples[@]}"; do
#   gatk CollectReadCounts \
#     -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam"\
#     -L ${resources}/panel_of_normals/targets_C.preprocessed.interval_list \
#     -R "${ref}"\
#     --format HDF5\
#     -imr OVERLAPPING_ONLY\
#     -O "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}.read_counts.hdf5"
#done

#Download PoN for CNV Analysis (hg38)
#wget storage.googleapis.com/gatk-test-data/cnv/somatic/wes-do-gc.pon.hdf5

#Remove noise from sample coverage using the PoN
#for sample in "${samples[@]}"; do
#  gatk --java-options "-Xmx7g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" DenoiseReadCounts \
#     -I "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}.read_counts.hdf5" \
#     --count-panel-of-normals "${resources}/panel_of_normals/cnvponC.pon.hdf5"\
#     --standardized-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.standardizedCR.tsv" \
#     --denoised-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.denoisedCR.tsv"
#done

#Annotating intervals with GC content, mappability, and segmental-duplication content
#gatk AnnotateIntervals \
#          -R ${ref} \
#          -L "${resources}/panel_of_normals/targets_C.preprocessed.interval_list"\
#          --interval-merging-rule OVERLAPPING_ONLY \
#          -O "${resources}/panel_of_normals/targets_C.annotated.interval_list" 

##Remove noise from sample coverage this time with annotated_intervals instead of PoNs
#for sample in "${samples[@]}"; do
#  gatk DenoiseReadCounts \
#          -I "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}.read_counts.hdf5" \
#          --annotated-intervals "${resources}/panel_of_normals/targets_C.annotated.interval_list"\
#          --standardized-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.standardizedCR.tsv" \
#          --denoised-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.denoisedCR.tsv"
#done

#Plotting Denoised Copy Ratios
#for sample in "${samples[@]}"; do
#  gatk PlotDenoisedCopyRatios \
#     --standardized-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.standardizedCR.tsv" \
#     --denoised-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.denoisedCR.tsv" \
#     --sequence-dictionary "${resources}/Homo_sapiens_assembly38.dict" \
#     --minimum-contig-length 46709983 \
#     --output "${results}/cna_calling/plots" \
#     --output-prefix "BSF_1274_H3CTKDSX5#${sample}_denoised"
#done

##Bring the plot to my local comp
#scp kacar@i121srv02.vu-wien.ac.at:/home/kacar/masterproject/Project_1/Dataflow/cna_plots/T_clean.denoised.png emrah@emrah-VirtualBox
#rsync -avz --partial --progress -e "ssh -i /home/emrah/.ssh/id_rsa -p 12121" kacar@i121srv02.vu-wien.ac.at:/home/kacar/masterproject/Project_1/Dataflow/cna_plots/T_clean.denoised.png .

#for sample in "${samples[@]}"; do
#   gatk --java-options "-Xmx3g" CollectAllelicCounts \
#     -L "${resources}/small_exac_common_3.hg38.vcf.gz"\
#     -I "${merged_reads}/BSF_1274_H3CTKDSX5%23${sample}_dedup_recal.bam" \
#     -R ${ref} \
#     -O "${results}/cna_calling/${sample}.allelicCounts.tsv"
#done

#Modelling segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts
#for sample in "${samples[@]}"; do
# gatk --java-options "-Xmx7g" ModelSegments \
#    --denoised-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.denoisedCR.tsv" \
#    --allelic-counts "${results}/cna_calling/${sample}.allelicCounts.tsv"\
#    --number-of-changepoints-penalty-factor 5.0 \
#    --output "${results}/cna_calling" \
#    --output-prefix "BSF_1274_H3CTKDSX5#${sample}_segments_new"
#done

#Plot Modeled Segments
#for sample in "${samples[@]}"; do
# gatk PlotModeledSegments \
#    --denoised-copy-ratios "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_T_clean.denoisedCR.tsv" \
#    --segments "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_new.modelFinal.seg" \
#    --sequence-dictionary "${resources}/Homo_sapiens_assembly38.dict" \
#    --minimum-contig-length 46709983 \
#    --output "${results}/cna_calling/plots" \
#    --output-prefix "BSF_1274_H3CTKDSX5#${sample}_modelled_new"
#done

#Calling copy-ratio segments as amplified, deleted, or copy-number neutral
for sample in "${samples[@]}"; do
 gatk CallCopyRatioSegments \
    -I "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_new.cr.seg" \
    -O "${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_new.called.seg"
done
