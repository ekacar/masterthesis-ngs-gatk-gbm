#!/bin/bash
vcf="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"
titan="/home/kacar/apps/TitanCNA-master/scripts/R_scripts"
samples=("COMP5_1_S126518" "COMP5_2_S126517" "COMP5_5_S126520" "COMP9_5_S126523" "COMP9_6_S126526" "COMP9_7_S126525" "GA11_1_S126532" "GA11_2_S126531" "GA11_3_S126534")
results="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results"
titanresults="/data/Unit_LMM/selberherr-group/kacar/ngs-gatk-trial/results/cna_calling/titancna"

#Translate CNA results from GATK4 to TITAN 
#conda activate  titancna

#Convert GATK hets from ModelSegments to Titan input:
#for sample in "${samples[@]}"; do
#   cat ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_new.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_titan.hets.tsv
#done

#for sample in "${samples[@]}"; do
#   cat ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_new.cr.seg | awk 'BEGIN{print"chr\tstart\tend\tlog2_TNratio_corrected"} /^[^@C]/ { print $1,$2,$3,$5 }' | tr ' ' '\t' > ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_titan.cr.tsv
#done

#Running Titan for CNA Analysis
#for sample in "${samples[@]}"; do
#    Rscript ${titan}/titanCNA.R --id ${sample} --hetFile ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_titan.hets.tsv --cnFile ${results}/cna_calling/BSF_1274_H3CTKDSX5#${sample}_segments_titan.cr.tsv --numClusters 1 --numCores 5 --normal_0 0.5 --ploidy_0 2 --estimatePloidy TRUE --outDir ${titanresults}
#done

#Creating ssm_txt for PhyloWGS
#for sample in "${samples[@]}"; do
#    create_phylowgs_inputs.py -s 5000 --vcf-type sample=mutect_smchet sample=${vcf}/ssv_calling/BSF_1274_H3CTKDSX5#${sample}_final_vcf.gz --regions=all  --output-variants ${vcf}/phylo/ssm_data_${sample}.txt
#done

#Creating cnv.txt for PhyloWGS from Titan CNA output
#for sample in "${samples[@]}"; do
#   python /home/kacar/miniconda3/envs/phylowgs/share/phylowgs/parser/parse_cnvs_new.py -f titan -c 0.81 --cnv-output ${results}/phylo/${sample}.cnv.txt  ${titanresults}/${sample}_cluster01.segs.txt 
#done

#Erase CNAs with the same start and end point
#for sample in "${samples[@]}"; do
#     awk -F'\t' '$2 != $3' ${results}/phylo/${sample}.cnv.txt > ${results}/phylo/${sample}_corr.cnv.txt
#    python /home/kacar/miniconda3/envs/phylowgs/share/phylowgs/parser/create_phylowgs_inputs.py -s 5000 --cnvs ${sample}=${results}/phylo/${sample}_corr.cnv.txt --vcf-type ${sample}=vardict ${sample}=${vcf}/ssv_calling/BSF_1274_H3CTKDSX5#${sample}_final_vcf --output-cnvs ${vcf}/phylo/${sample}_cnv_data.txt
#done

for sample in "${samples[@]}"; do
    python2 /home/kacar/miniconda3/envs/phylowgs/share/phylowgs/multievolve.py -n 7 -O ${results}/phy$ --ssms ${vcf}/phylo/ssm_data_${sample}.txt --cnvs ${vcf}/phylo/${sample}_cnv_data.txt 
done

#write_results.py example_data trees.zip example_data.summ.json.gz example_data.muts.json.gz example_data.mutass.zip
