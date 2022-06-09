#BSUB -W 30:00
#BSUB -n 10
#BSUB -M 100000

module load souporcell
module load bedtools

source /users/katy2h/.bashrc

conda activate souporcell

cd /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor4/15/Souporcell_In

souporcell_pipeline.py -i atac_possorted_bam.bam -b barcodes.tsv -f /data/miraldiNB/Katko/Projects/Genomes/refdata-GRCh38-1.0.0/fasta/genome.fa --common_variants /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Control_Genome/common_variants_with_chr_2_005.vcf --skip_remap True --no_umi True --ignore True -o souporcell_out_commonvariants_atac -t 10 -k 2
