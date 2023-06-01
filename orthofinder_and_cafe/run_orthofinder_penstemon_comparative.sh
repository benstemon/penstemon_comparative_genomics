#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH --out=slurm-OrthoFinder.%j.out
#SBATCH -p wessinger-48core
#SBATCH --job-name=orthofinder


cd $SLURM_SUBMIT_DIR


#load conda env with modules numpy and scipy
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#load necessary modules
module load blast

#set important variables
orthofinder="/work/bs66/software/OrthoFinder_source"
indir="/work/bs66/testing_comparative_genomics/orthogroup_testing/single_isoform_CDS_fastas"
NThreads=12


#run OrthoFinder basic full analysis
#include option -d to indicate sequences are nucleotide sequences
python $orthofinder/orthofinder.py -t $NThreads -a $NThreads -d -f $indir

