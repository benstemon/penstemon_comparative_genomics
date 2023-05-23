#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=run_gemoma

#########################################################


cd $SLURM_SUBMIT_DIR

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#specify n threads
NThreads=8


#specify locations of reference genome fastas and gff files
ref1fa="/work/bs66/project_compare_genomes/barbatus/Pbar.2022.LG.fa"
ref1gff="/work/bs66/project_compare_genomes/barbatus/M4_annotation_putative_function_domain_added_blast_tomato.genemodels.noseq.gff"

ref2fa="/work/bs66/project_compare_genomes/davidsonii_unfiltered/annot_Pdavidsonii_genome.fasta"
ref2gff="/work/bs66/project_compare_genomes/davidsonii_unfiltered/annot_Pdavidsonii_genome.gff"

ref3fa="/work/bs66/project_compare_genomes/petiolatus/petiolatus_genome.fasta"
ref3gff="/work/bs66/project_compare_genomes/petiolatus/Rnd1.all.maker.snapdragon.noseq.gff"


#######################################################
  # ANNOTATE GENOME WITH GEMOMA
  
  #load conda environment
  conda activate gemoma
  mkdir outfiles_gemoma

  #run GeMoMa
  #make new names so that they can be filtered appropriately later on (or not)
  GeMoMa GeMoMaPipeline threads=$NThreads AnnotationFinalizer.r=NO p=true o=false\
   t=ragtag.scaffold.fasta outdir=outfiles_gemoma/ \
   s=own i=barbatus a=$ref1gff g=$ref1fa \
   s=own i=davidsonii a=$ref2gff g=$ref2fa \
   s=own i=petiolatus a=$ref3gff g=$ref3fa


