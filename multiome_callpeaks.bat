#BSUB -n 2
#BSUB -M 100000
#BSUB -W 10:00 

module load R/4.1.1

Rscript multiome_callpeaks.R
