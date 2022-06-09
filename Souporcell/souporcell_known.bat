#BSUB -W 20:00
#BSUB -n 10
#BSUB -M 100000

module load souporcell
module load bedtools

source /users/katy2h/.bashrc

conda activate souporcell

cd /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor4/15/Souporcell_In

souporcell_pipeline.py -i gex_possorted_bam.bam -b barcodes.tsv -f /data/miraldiNB/Katko/Projects/Genomes/refdata-GRCh38-1.0.0/fasta/genome.fa --known_genotypes /data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Control_Genome/Control_Donor4_15_allchr.vcf --known_genotypes_sample_names Control sample -o souporcell_out_commonvariants_known_allchr -t 10 -k 2
