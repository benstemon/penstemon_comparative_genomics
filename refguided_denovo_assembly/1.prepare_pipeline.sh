#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=prepare_genome

cd $SLURM_SUBMIT_DIR

#script to prepare the reference genome

source /home/bs66/.bashrc
source /home/bs66/.bash_profile



#name of the working directory (directory with programsm including sort_and_filter_fasta.py)
workPathFiles="/work/bs66/testing_comparative_genomics/refguided_denovo_testing"

#path to reference genome, and name of reduced reference genome
ref="/work/bs66/project_compare_genomes/barbatus/Pbar.2022.LG.fa"
refRed="${workPathFiles}/reduced_reference.fa"




  # PREPARE REFERENCE -- no scaffolds < 10kb
  conda activate mapping_etc
  
  python3 sort_and_filter_fasta.py $ref $refRed 10000
  samtools faidx $refRed
  bwa index $refRed

