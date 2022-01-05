# CD4-multiome
Contains code for the analysis of the CD4 multiome data produced from the Barski lab. 

# SouporCell

Souporcell module on the cluster already exists with every all the parts it needs. You will just need to activate the correct conda envrionment (already created). Inputs are the postsorted bam (cellranger output), the barcodes tsv file (cellranger output) and the reference fasta file. -k is the number of clusters (genotypes). -t is the number of threads (= number of cores).

1. Module load souporcell
2. conda activate souporcell
3. souporcell_pipeline.py -i gex_possorted_bam.bam -b barcodes.tsv -f hg38.fna -o souporcell_out -t 10 -k 2

Problems found and how I fixed them

1. No package named pystan
    1. This happened simply because I didnt activate the souporcell envrionment. To fix this, see step 2
2. fa.fai file may not be writable
    1. This happened because the .fa reference file I was using wasn't writeable because it was in the weirauchlab data bank. There may be a better fix, but I just downloaded my own Hg38 reference file from the internet and placed it in my writable folder.
