#BSUB -W 30:00
#BSUB -n 15
#BSUB -M 100000

echo Running Souporcell 

module load souporcell
module load bedtools

source /users/katy2h/.bashrc
set location = /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor1/2/Souporcell_In
set clusters = 2
set threads = 15
set bam_file = atac_possorted_bam.bam
set barcodes_file = barcodes.tsv
set reference_genome = /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor1/15/souporcell_in
set Skip_Remap = True
set No_Umi = True
set variants_file = /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Control_Genome/common_variants_with_chr_2_005.vcf
set output = souporcell_out_commonvariants_atac_skip

conda activate souporcell
cd ${location}
souporcell_pipeline.py -i ${bam_file} -b ${barcodes_file} -f ${reference_genome} --common_variants ${variants_file} --skip_remap ${Skip_Remap} --no_umi ${No_Umi} -o ${output} -t ${threads} -k ${clusters}