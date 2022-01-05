#BSUB -W 10:00
#BSUB -n 10

module load souporcell

cd /data/weirauchlab/collab/miraldi/kat6ti/Barski_Multiome/data/resting_combined_set2

souporcell_pipeline.py -i gex_possorted_bam.bam -b /raw_feature_bc_matrix/barcodes.tsv -f /data/weirauchlab/databank/genomes/hg38/hg38.fa -o souporcell_out -k 2
