#BSUB -n 1
#BSUB -M 200000
#BSUB -W 30:00 

module purge
module load R/4.1.1

Rscript mergeall.r
