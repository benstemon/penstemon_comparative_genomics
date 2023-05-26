#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=refdenovo_assembly
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bs66@mailbox.sc.edu

#########################################################
#  Reference-guided de novo assembly - SPAdes
# ====================================================
# by Heidi Lischer, 2015/2016
# Updated by Ben Stone, 2023
#########################################################

# set variables #########################################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#raw PE reads. 
rawreads="/work/bs66/testing_comparative_genomics/refguided_denovo_testing/rawreads_smallii"

#other important directories and files
workPathFiles="/work/bs66/testing_comparative_genomics/refguided_denovo_testing"
refRed="${workPathFiles}/reduced_reference.fa"



NThreads=16      # set the number of threads of every parallelizable step


# paired-end libraries -------------------
name="smallii"           # set name of your species


# set work path ---------------------------
workPath=${workPathFiles}/${name}_spades
# log file
log=${workPath}/log_${name}_spades.txt




# Programs --------------------------------
progPrepareReference="${workPathFiles}/sort_and_filter_fasta.py"
progSPAdes="/work/bs66/software/SPAdes-3.15.4-Linux/bin/spades.py"
progGetBlocks="${workPathFiles}/GetBlocks.jar"
progRemovShortSeq="${workPathFiles}/RemoveShortSeq.jar"
progFastaToAmos="${workPathFiles}/FastaToAmos.jar"
progFastaStats="${workPathFiles}/FastaStats.jar"
progSplitSeqLowCov="${workPathFiles}/SplitSeqLowCov.jar"
progWriteSoapConfig="${workPathFiles}/WriteSoapConfig.jar"
prognucmer="/work/bs66/software/mummer-4.0.0rc1/nucmer"
progPicard="/share/apps/gcc/4.8.5/picard2018/picard.jar"
progQuast="/work/bs66/software/quast-5.2.0/quast.py"


#load modules on HPC
module load java
module load bedtools

#necessary conda installs
#samtools
#bamtools
#bwa
#fastqc
#multiqc
#fastp
#soapdenovo2
#soapdenovo2-prepare


#########################################################



# run pipeline ##########################################
mkdir ${workPath}


# 1. Step: quality/adapter trimming and quality check:
#######################################################
  # INITIAL QUALITY CHECK ----------
  echo "1a. quality check of raw reads..."
  echo "1a. quality check of raw reads..." > $log
  
  #activate conda environment with QC packages installed
  conda activate QC
  
  #run fastqc and summarize with multiqc
  mkdir $workPath/fastqc-out_rawreads
  cd $rawreads
  files=(*.fastq.gz)
  fastqc "${files[@]}" -t $NThreads -o $workPath/fastqc-out_rawreads
  multiqc $workPath/fastqc-out_rawreads -o $workPath/fastqc-out_rawreads


  # MERGE ILLUMINA LANES (if reads across >1 lane)
  mkdir $workPath/merged_reads
  for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
  
      do echo "Merging R1"
  cat "$i"_L00*_R1_001.fastq.gz > "$i"_merged_L001_R1_001.fastq.gz
         echo "Merging R2"
  cat "$i"_L00*_R2_001.fastq.gz > "$i"_merged_L001_R2_001.fastq.gz
  done;
  mv *merged* $workPath/merged_reads


  # TRIM ADAPTERS AND PERFORM QC ON READS
  echo "performing fastp step..."
  echo "performing fastp step..." >> $log
  mkdir $workPath/filtered_reads
  cd $workPath/merged_reads
  for r1in in $workPath/merged_reads/*_R1_001.fastq.gz; 
  do
      r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
      r1out="${r1in##*/}"
      r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
      fastp -i "$r1in" -I "$r2in" --detect_adapter_for_pe -l 30 --out1 $workPath/filtered_reads/"${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R1_001.fastq.gz}" --out2 $workPath/filtered_reads/"${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R2_001.fastq.gz}" -x -c -w $NThreads -h $workPath/filtered_reads/"${r1out/merged_L001_R1_001.fastq.gz/html}" -j $workPath/filtered_reads/"${r1out/merged_L001_R1_001.fastq.gz/json}"
  done


  # POST-QC QUALITY CHECK
  echo "1b. quality check of filtered reads..."
  echo "1b. quality check of filtered reads..." >> $log
  mkdir $workPath/fastqc-out_filtered
  cd $workPath/filtered_reads
  
  # perform fastqc and multiqc on quality-filtered reads
  files=(*.fastq.gz)
  fastqc "${files[@]}" -t $NThreads -o $workPath/fastqc-out_filtered
  multiqc $workPath/fastqc-out_filtered -o $workPath/fastqc-out_filtered
  
  #return to base conda
  conda deactivate


#originally this created mateRead1TrimUnPair and mateRead2TrimUnPair
#originally, fastQC -> Trimmomatic -> fastqc. Now, replace Trimmomatic with fastp




# 2. Step: map reads against reference 
#          and define blocks and superblocks
#######################################################
  
  #load conda env with mapping tools installed
  conda activate mapping_etc
  
  # MAP READS AGAINST REFERENCE ----------
  echo "2a. run reference mapping..."
  echo "2a. run reference mapping..." >> $log
  
  cd $workPath/filtered_reads
  
  #for all R1s in the filtered reads dir
  for i in *R1_001.fastq.gz
  do
      #label reads and naming scheme
      read1=$i
      read2="${read1/L001_R1/L001_R2}"
      outname="${read1/_trimmed_L001_R1_001.fastq.gz/}"
      
      
      #map reads and send EVERYTHING to all.sorted.bam
      bwa mem -t $NThreads -M $refRed $read1 $read2 | \
       samtools fixmate -@ $NThreads -m -u -O bam - - | \
       samtools sort -@ $NThreads -u - | \
       samtools markdup -@ $NThreads - ${outname}_all.sorted.bam
      samtools index ${outname}_all.sorted.bam
      
      
      #send only mapped reads to mapped.sorted.bam
      #bams should retain their sort property from before
      samtools view -@ $NThreads -b -F 4 ${outname}_all.sorted.bam -o ${outname}.mapped.sorted.bam
      samtools index ${outname}.mapped.sorted.bam
      
      #send unmapped reads to unmapped.sorted.bam
      samtools view -@ $NThreads -b -f 4 ${outname}_all.sorted.bam > ${outname}.unmapped.sorted.bam
      
      #from unmapped reads, send PAIRED reads unmapped.failPair fastqs
      samtools fastq -@ $NThreads -f 9 -N ${outname}.unmapped.sorted.bam -1 ${outname}.unmapped.failPair.1.fastq -2 ${outname}.unmapped.failPair.2.fastq
      
      #from unmapped reads, send UNPAIRED reads to new failUnPair fastqs
      samtools fastq -@ $NThreads -F 8 -f 64 -N ${outname}.unmapped.sorted.bam -1 ${outname}.unmapped.failUnPair.1.fastq
      samtools fastq -@ $NThreads -F 8 -f 128 -N ${outname}.unmapped.sorted.bam -2 ${outname}.unmapped.failUnPair.2.fastq
      
      
      
      #get summary statistics from bamtools
      bamtools stats -in ${outname}_all.sorted.bam >> $log
      echo "--> ${outname}_all.sorted.bam" >> $log
      
      #filter for mapping quality >=10 and read stats to log
      samtools view -b -q 10 ${outname}.mapped.sorted.bam > ${outname}.mapped.sorted.filtered.bam
      bamtools stats -in ${outname}.mapped.sorted.filtered.bam >> $log
      echo "--> ${outname}.mapped.sorted.filtered.bam" >> $log
      
      
      
      # GET BLOCKS AND SUPERBLOCKS ----------
      echo "2b. get blocks and superblocks..."
      echo "2b. get blocks and superblocks..." >> $log
      
      #generate coverage file with the mapped.sorted reads
      covFile=${workPath}/${outname}_coverage.txt
      bedtools genomecov -ibam ${outname}.mapped.sorted.bam -bga > $covFile
      
      
      #generate coverage file with only properly paired mapped reads (also sort by name?)
      samtools view -@ $NThreads -b -f 2 ${outname}.mapped.sorted.bam | \
       samtools sort -@ $NThreads - -n | \
       bedtools bamtobed -i - -bedpe | \
       awk '$1 == $4' | \
       cut -f 1,2,6 | \
       sort -k 1,1 | \
       bedtools genomecov -i - -bga -g ${refRed}.fai > ${covFile%.txt}Paired.txt
#this gives lots of warnings about reads marked as paired, but not occurring next to it in bam file
#I am unsure if this is a legitimately serious problem. Going to push forward for now.

      
      #prepare outfiles for blocks and superblocks
      blocks=${workPath}/${outname}_blocks.txt
      superblocks=${workPath}/${outname}_superblocks.txt
      
      #get blocks -- minimum 10 (PE) reads
      #get superblocks -- >= 12kb, minimal overlap of 300bp (max overlap = 3*300 bp)
      java -jar ${progGetBlocks} -i ${covFile} -paired ${covFile%.txt}Paired.txt -o ${blocks} -oSuper ${superblocks} -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 100000
      
  done


#Step 2 -- mapping and defining blocks -- should be considered complete,
#barring determining that the warnings about paired reads require fixing.
#After consideration, I think that is OK, because we are filtering out those reads.
#notes: originally, step 2 used bowtie to map andtrimmomatic for QC and trimming.
#now, I am using bwa mem to map, and fastp for QC/trimming
#I also simplified some of the samtools steps
#block and superblock construction remains the same


# 3. Step: do deNovo assembly within superblocks (SPAdes)
#######################################################
  echo "3a. deNovo assembly within superblocks..."
  echo "3a. deNovo assembly within superblocks..." >> $log
  
  cd ${workPath}
  
  #create dir for infiles and move block files to here
  mkdir ${workPath}/assembly_infiles
  mv *blocks.txt ${workPath}/assembly_infiles
  
  
  # SEQUENCE GENERATION
  
  #find sequence names overlapping regions in the superblocks
  #these include all mapped reads -- both paired and unpaired 
  samtools view -@ $NThreads -b filtered_reads/${outname}.mapped.sorted.bam $(tr '\n' ' ' < assembly_infiles/${outname}_superblocks.txt) | \
   samtools sort -@ $NThreads - -n -o assembly_infiles/superblock.sequences_${outname}.bam
  
  
  # SUBSEQUENCE GENERATION
  #generate fastq files for these superblock sequences
  #NOTE that this does not work with samtools fastq -- unsure, but related to skipping reads
  #perhaps because of name-sort... in any case, bedtools should work
  cd ${workPath}/assembly_infiles
  bedtools bamtofastq -i superblock.sequences_${outname}.bam -fq 1_subseq_${outname}_R1.fastq -fq2 1_subseq_${outname}_R2.fastq
  
  
  #extract paired reads with one pair unmapped
  #first, do this for the first in the pair (-f 72)
  samtools view -b -f 72 superblock.sequences_${outname}.bam | \
   bamtools convert -format fastq | \
   paste - - - - | \
   sort -k1,1 -t ' ' | \
   tr '\t' '\n' > 2_subseq_${outname}_R1.fastq
  
  samtools view -f 72 superblock.sequences_${outname}.bam | \
   cut -f 1 | \
   awk "{print \$0\"/2\"}" > ${outname}_R2.txt
  
  seqtk subseq ${workPath}/filtered_reads/${outname}.unmapped.failUnPair.2.fastq ${outname}_R2.txt | \
   paste - - - - | \
   sort -k1,1 -t ' ' | \
   tr '\t' '\n' > 2_subseq_${outname}_R2.fastq
  
  
  #then extract reads for second in the pair (-f 136)
  #NOTE: I am pretty sure there is an error here in original pipeline
  #should be writing to R2, R1, R1
  samtools view -b -f 136 superblock.sequences_${outname}.bam | \
   bamtools convert -format fastq | \
   paste - - - - | \
   sort -k1,1 -t ' ' | \
   tr '\t' '\n' > 3_subseq_${outname}_R2.fastq
  
  samtools view -f 136 superblock.sequences_${outname}.bam | \
   cut -f 1 | \
   awk "{print \$0\"/1\"}" > ${outname}_R1.txt
  
  seqtk subseq ${workPath}/filtered_reads/${outname}.unmapped.failUnPair.1.fastq ${outname}_R1.txt | \
   paste - - - - | \
   sort -k1,1 -t ' ' | \
   tr '\t' '\n' > 3_subseq_${outname}_R1.fastq
   
   
  #combine all of the unpaired sequences into a single file.
  cat 1_subseq_${outname}_R1.fastq 2_subseq_${outname}_R1.fastq 3_subseq_${outname}_R1.fastq > superblock.subseq_${outname}_R1.fastq
  cat 1_subseq_${outname}_R2.fastq 2_subseq_${outname}_R2.fastq 3_subseq_${outname}_R2.fastq > superblock.subseq_${outname}_R2.fastq
  
  #perform a bit of cleanup
  rm 1_* 2_* 3_*
  
  
  
  # DE NOVO ASSEMBLY OF MAPPED READS
  echo "3b. de novo assembly -- mapped reads..."
  echo "3b. de novo assembly -- mapped reads..." >> $log
  
  mkdir ${workPath}/assembly_outfiles
  cd ${workPath}/assembly_infiles
  
  #gzip the fastq files for input into SPAdes
#  gzip -c superblock.subseq_${outname}_R1.fastq > superblock.subseq_${outname}_R1.fq.gz
#  gzip -c superblock.subseq_${outname}_R2.fastq > superblock.subseq_${outname}_R2.fq.gz
  
  python $progSPAdes \
   -1 superblock.subseq_${outname}_R1.fastq -2 superblock.subseq_${outname}_R2.fastq \
   -t $NThreads -o ${workPath}/assembly_outfiles


  # DE NOVO ASSEMBLY OF UNASSEMBLED (UNPAIRED) READS ----------
  echo "3c. deNovo assembly of unassembled reads..."
  echo "3c. deNovo assembly of unassembled reads..." >> $log
  
  cd ${workPath}/filtered_reads
  
  #merge unpaired files
  cat ${outname}.unmapped.failUnPair.1.fastq ${outname}.unmapped.failUnPair.1.fastq > ${outname}.FailUnpairMerged.fastq
  
  #run SPAdes assembly on all unmapped reads
  mkdir ${workPath}/assembly_outfiles_unmapped
  cd ${workPath}/filtered_reads
  
  python $progSPAdes \
   -1 ${outname}.unmapped.failPair.1.fastq -2 ${outname}.unmapped.failPair.2.fastq \
   -s ${outname}.FailUnpairMerged.fastq \
   -t $NThreads -o ${workPath}/assembly_outfiles_unmapped
  
  
  
# 4. Step: get non-redundant supercontigs
####################################################### 
  echo "4a. get supercontigs..."
  echo "4a. get supercontigs..." >> $log
  cd ${workPath}
  
  #Remove short sequences (<500 bp) from the assembly with unpaired reads
  java -jar ${progRemovShortSeq} -i assembly_outfiles_unmapped/contigs.fasta \
   -o assembly_outfiles_unmapped/contigs_500.fasta -length 500 >> $log
  
  
  #merge the contigs from both assemblies into new file (AKA superblockSeq)
  cat assembly_outfiles/contigs.fasta assembly_outfiles_unmapped/contigs_500.fasta > deNovo_Superblocks.fa
  
  #Again, remove short seqs (this time <200) -- NOTE: none removed -- commenting out currently
  java -jar ${progRemovShortSeq} -i deNovo_Superblocks.fa -o deNovo_Superblocks_200.fa \
   -length 200 -u >> $log
   rm deNovo_Superblocks.fa
  
  
  # REMOVE REDUNDANCY WITH AMOSCMP
  
  #load the refdenovo conda environment and set up new AMOS folder
  conda deactivate
  conda activate refdenovo
  mkdir ${workPath}/AMOScmp
  cd ${workPath}/AMOScmp
  
  #assemble all assembled superblocks with AMOScmp to supercontigs
  #changed parameters in AMOScmp: (casm-layout -t 1000 (maximum ignorable trim length), make-consensus -o 10 (minimum overlap base)) 
  java -jar ${progFastaToAmos} -i ${workPath}/deNovo_Superblocks_200.fa -o superblockseq_Amos.afg
  
  #in case this is relevant later ... 
  supercontigs=Amos_supercontigs
  amosSupercontigs=${amosFolder}/Amos_supercontigs.fasta
  
  #can in theory run the whole thing in this step
  #but the nucmer loaded with AMOS on bioconda isn't functioning correctly
  #echo "run AMPScmp..." >> $log
  #AMOScmp -D TGT=superblockseq_Amos.afg -D REF=$refRed Amos_supercontigs
  
  
  # running AMPScmp step by step and use multithread nucmer to spead it up
  ## Building AMOS bank
  echo "4b.  build AMPS bank..." >> $log
  bank-transact -c -z -b Amos_supercontigs.bnk -m superblockseq_Amos.afg

  ## Collecting clear range sequences
  echo "4b.  clear range sequences..." >> $log
  dumpreads Amos_supercontigs.bnk > Amos_supercontigs.seq

  ## Running nucmer
  echo "4b.  run nucmer..." >> $log
  #note -- the conda install version of this will not work.
  #so the software must be installed with the mummer package on github
  $prognucmer --maxmatch -t $NThreads --prefix=Amos_supercontigs $refRed Amos_supercontigs.seq

  ## Running layout -- max ignorable trim length = 1000
  echo "4b.  run layout..." >> $log
  casm-layout -t 1000 -U Amos_supercontigs.layout -C Amos_supercontigs.conflict -b Amos_supercontigs.bnk Amos_supercontigs.delta

  ## Running consensus
  echo "4b.  run consensus..." >> $log
  make-consensus -o 10 -B -b Amos_supercontigs.bnk

  ## Outputting contigs
  echo "4b.  output contigs..." >> $log
  bank2contig Amos_supercontigs.bnk > Amos_supercontigs.contig

  ## Outputting fasta
  echo "4b.  output fasta..." >> $log
  bank2fasta -b Amos_supercontigs.bnk > Amos_supercontigs.fasta
  
  
# 5. Step: map reads on supercontigs
#          and de novo assemble unmapped reads
####################################################### 
  echo "5a. map reads on supercontigs and correct them..."
  echo "5a. map reads on supercontigs and correct them..." >> $log
  
  #make seqnames unique
  java -jar ${progRemovShortSeq} -i Amos_supercontigs.fasta -o Amos_supercontigs_unique.fa -length 1 -u
  
  #get statistics
  echo Amos_supercontigs_unique.fa >> $log
  java -jar ${progFastaStats} -i Amos_supercontigs_unique.fa -min 200 >> $log
  
  

  # A. MAPPPING TO THE NEW SUPERCONTIGS
  #index the new reference -- the supercontigs
  # deactivate refdenovo conda environment and load the mapping environment again
  conda deactivate
  conda activate mapping_etc
  cd ${workPath}/AMOScmp
  bwa index Amos_supercontigs_unique.fa
  
  #map trimmed reads back to the supercontigs reference
  cd ${workPath}
  bwa mem -t $NThreads -M AMOScmp/Amos_supercontigs_unique.fa \
   filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
   filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz | \
   samtools fixmate -@ $NThreads -m -u -O bam - - | \
   samtools sort -@ $NThreads -u - | \
   samtools markdup -@ $NThreads - AMOScmp/${outname}_supercont_all.sorted.bam
  samtools index AMOScmp/${outname}_supercont_all.sorted.bam
  
  
  #create bam for unmapped reads
  cd AMOScmp
  samtools view -@ $NThreads -b -f 4 ${outname}_supercont_all.sorted.bam > ${outname}_supercont_unmapped.sorted.bam
  
  #from unmapped reads, send PAIRED to separate fastqs
  samtools view -@ $NThreads -b -f 9 ${outname}_supercont_unmapped.sorted.bam > ${outname}_supercont_unmapped_pair.sorted.bam
  samtools fastq -@ $NThreads -N ${outname}_supercont_unmapped_pair.sorted.bam \
   -1 ${outname}_supercont_FailPair.1.fastq \
   -2 ${outname}_supercont_FailPair.2.fastq
  rm ${outname}_supercont_unmapped_pair.sorted.bam
  
  #From unmapped reads, send UNPAIRED reads to a single failUnPair fastq
  samtools view -b -F 8 ${outname}_supercont_unmapped.sorted.bam | \
   bamtools convert -format fastq -out ${outname}_supercont_failUnPair.fastq
  
  
  #get statistics for these reads
  bamtools stats -in ${outname}_supercont_all.sorted.bam >> $log
  echo "--> ${outname}_supercont_all.sorted.bam" >> $log
  
  
  #filter for mapping quality >=10 and get stats for these as well
  samtools view -b -F 4 -q 10 ${outname}_supercont_all.sorted.bam > ${outname}_supercont_all.sorted.filtered.bam
  bamtools stats -in ${outname}_supercont_all.sorted.filtered.bam >> $log
  echo "--> ${outname}_supercont_all.sorted.filtered.bam" >> $log
  
  

  # B. DENOVO ASSEMBLY OF THE UNMAPPED/UNSASSEMBLED READS
  echo "5b. deNovo assemble unassembled reads..."
  echo "5b. deNovo assemble unassembled reads..." >> $log
  
  
  #assemble with spades
  mkdir ${workPath}/AMOScmp/unassembled-assembly_outfiles
  python $progSPAdes \
   -1 ${outname}_supercont_FailPair.1.fastq -2 ${outname}_supercont_FailPair.2.fastq \
   -s ${outname}_supercont_failUnPair.fastq \
   -t $NThreads -o ${workPath}/AMOScmp/unassembled-assembly_outfiles
  
  
  #after assembly, remove contigs shorter than 200 bp -- check names prior to running
  cd unassembled-assembly_outfiles
  java -jar ${progRemovShortSeq} -i contigs.fasta -o unass_contigs_200.fa -length 100 >> $log
  echo unass_contigs_200.fa >> $log
  java -jar ${progFastaStats} -i unass_contigs_200.fa -min 200 >> $log


# 6. Step: map reads to all supercontigs and correct them 
########################################################## 
  echo "6a. merge contigs..." 
  echo "6a. merge contigs..." >> $log
  
  mkdir ${workPath}/merged_corr
  cd ${workPath}/merged_corr

  
  #create file for merged unique and unmapped supercontigs
  cat ${workPath}/AMOScmp/Amos_supercontigs_unique.fa ${workPath}/AMOScmp/unassembled-assembly_outfiles/unass_contigs_200.fa > ${outname}_supercontSeq_Unass.fa
  
  echo "${outname}_supercontSeq_Unass.fa" >> $log
  java -jar ${progFastaStats} -i ${outname}_supercontSeq_Unass.fa -min 200 >> $log
  
  
  #create new reference for this concatenated file
  bwa index ${outname}_supercontSeq_Unass.fa 
  
  
  
  #map reads to the merged supercontigs
  #input are original read1 and read2, with reference as the new merged supercontigs
  echo "6b. map reads to contigs..." 
  echo "6b. map reads to contigs..." >> $log

  bwa mem -t $NThreads -M ${outname}_supercontSeq_Unass.fa \
   ${workPath}/filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
   ${workPath}/filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz | \
   samtools fixmate -@ $NThreads -m -u -O bam - - | \
   samtools sort -@ $NThreads -u - | \
   samtools markdup -@ $NThreads - ${outname}_all.sorted.bam
  samtools index ${outname}_all.sorted.bam
  
  #read stats to log
  bamtools stats -in ${outname}_all.sorted.bam >> $log
  echo "--> ${outname}_all.sorted.bam" >> $log
  
  #filter for mapping quality >=10
  samtools view -b -F 4 -q 10 ${outname}_all.sorted.bam > ${outname}.filtered.sorted.bam
  bamtools stats -in ${outname}.filtered.sorted.bam >> $log
  echo "--> ${outname}.filtered.sorted.bam" >> $log
  
  
  
  # ERROR CORRECTION
  #add RG header of the second file (lost while merging)
  #because I don't have different libraries, I think this should be fine. I don't need to merge...
  #... as far as I know... Hopefully this doesn't break!
  
  #index the filtered sorted bam and generate a fasta index for the supercontseq_Unass.fa
  samtools index ${outname}.filtered.sorted.bam
  samtools faidx ${outname}_supercontSeq_Unass.fa
  
  #realign reads
  #original verion uses picard and gatk3 to align reads around indels
  #this is no longer functional, and gatk3 is deprecated. So I will use bcftools
  echo "6c. realign reads..." 
  echo "6c. realign reads..." >> $log
  bcftools mpileup --threads $NThreads -Ou -f ${outname}_supercontSeq_Unass.fa ${outname}.filtered.sorted.bam -Ou | \
   bcftools call --threads $NThreads -c -Ou | \
   bcftools norm --threads $NThreads -f ${outname}_supercontSeq_Unass.fa -Oz -o tmp.vcf.gz

  #index vcf and create consensus genome
  bcftools index tmp.vcf.gz
  cat ${outname}_supercontSeq_Unass.fa | bcftools consensus tmp.vcf.gz > ${outname}_supercontSeq_Unass_corr.fa
  rm tmp.vcf.gz
  rm tmp.vcf.gz.csi
  

  #remove start and end N, and a min of at least 100 bp
  echo ${outname}_supercontSeq_Unass_corrWN.fa >> $log
  java -jar ${progRemovShortSeq} -i ${outname}_supercontSeq_Unass_corr.fa -o ${outname}_supercontSeq_Unass_corrWN.fa -length 100 -n >> $log
  
  #get statistics
  echo ${outname}_supercontSeq_Unass_corrWN.fa >> $log
  java -jar ${progFastaStats} -i ${outname}_supercontSeq_Unass_corrWN.fa -min 200 >> $log
  
  
  
  # SPLIT SEQUENCES AT PLACES WITH NO COVERAGE
  #index WN.fa file
  bwa index ${outname}_supercontSeq_Unass_corrWN.fa
  
  
  #map the original trimmed reads, yet again, to the WN fasta file
  echo "6d. map reads to WN fasta..." 
  echo "6d. map reads to WN fasta..." >> $log

  bwa mem -t $NThreads -M ${outname}_supercontSeq_Unass_corrWN.fa \
   ${workPath}/filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
   ${workPath}/filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz | \
   samtools fixmate -@ $NThreads -m -u -O bam - - | \
   samtools sort -@ $NThreads -u - | \
   samtools markdup -@ $NThreads - ${outname}_corrWN_all.sorted.bam
  samtools index ${outname}_corrWN_all.sorted.bam

  #get stats
  bamtools stats -in ${outname}_corrWN_all.sorted.bam >> $log
  echo "--> ${outname}_corrWN_all.sorted.bam" >> $log
  
  #filter for MQ >=10
  samtools view -b -F 4 -q 10 ${outname}_corrWN_all.sorted.bam > ${outname}_corrWN_filtered.sorted.bam
  bamtools stats -in ${outname}_corrWN_filtered.sorted.bam >> $log
  echo "--> ${outname}_corrWN_filtered.sorted.bam" >> $log

  #get genome coverage information
  bedtools genomecov -ibam ${outname}_corrWN_filtered.sorted.bam -bga > ${outname}_supercontSeq_Unass_corrWN_filteredCov.txt
  
  #and only with properly paired reads...
  samtools faidx ${outname}_supercontSeq_Unass_corrWN.fa
  samtools view -bf 0x2 ${outname}_corrWN_filtered.sorted.bam | \
   samtools sort - -n | \
   bedtools bamtobed -i - -bedpe | \
   awk '$1 == $4' | \
   cut -f 1,2,6 | \
   sort -k 1,1 | \
   bedtools genomecov -i - -bga -g ${outname}_supercontSeq_Unass_corrWN.fa.fai > ${outname}_supercontSeq_Unass_corrWN_filteredPairedCov.txt

  
  #split sequences with low coverage
  echo "6e. split sequences with low coverage..." 
  echo "6e. split sequences with low coverage..." >> $log
  java -jar ${progSplitSeqLowCov} -i ${outname}_supercontSeq_Unass_corrWN_filteredCov.txt \
   -paired ${outname}_supercontSeq_Unass_corrWN_filteredPairedCov.txt \
   -o ${outname}_supercontSeq_Unass_corrWN_filteredNotCov.txt -mCov 1 -fasta ${outname}_supercontSeq_Unass_corrWN.fa \
   -fastaOut ${outname}_supercontSeq_Unass_corrWN_splitFiltered.fa >> $log
  echo ${outname}_supercontSeq_Unass_corrWN_splitFiltered.fa >> $log
  java -jar ${progFastaStats} -i ${outname}_supercontSeq_Unass_corrWN_splitFiltered.fa -min 200 >> $log
  
  


# 7. Step: scaffolding and gap closing
#######################################################
  echo "7a. scaffolding with SOAPdenovo..." 
  echo "7a. scaffolding with SOAPdenovo..." >> $log
  
  #make new directory for scaffolding
  cd ${workPath}/merged_corr
  mkdir scaffold_gapClosed
  cd scaffold_gapClosed
  
  
  #write SOAPdenovo config file
  #I am just going to use the average insert size estimated from bwa, which is 316
  java -jar ${progWriteSoapConfig} -insLength 316 \
   -r1 ${workPath}/filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
   -r2 ${workPath}/filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz \
   -max 151 -ru 2 -o soap.config
  
  # USE SOAPDENOVO2 TO SCAFFOLD AND CLOSE GAPS
  # unload the mapping environment again, and reload the ref denovo conda environment
  conda deactivate
  conda activate refdenovo
  
  #generate more config files
  #using kmer 61 because insert size is similar to guide and that's what they used...
  finalFusion -D -c ${workPath}/merged_corr/${outname}_supercontSeq_Unass_corrWN_splitFiltered.fa \
   -K 61 -g ${outname}_61 -p $NThreads
  
  #map and scaffold
  SOAPdenovo-127mer map -s soap.config -g ${outname}_61 -p $NThreads
  SOAPdenovo-127mer scaff -g ${outname}_61 -p $NThreads -F
  
  # remove scaffolds < 200 bp, 500bp, and 1kb
  echo ${outname}_61.scafSeq >> $log
  java -jar ${progRemovShortSeq} -i ${outname}_61.scafSeq -o ${outname}_scafSeq.fa -length 200 >> $log
  java -jar ${progRemovShortSeq} -i ${outname}_scafSeq.fa -o ${outname}_scafSeq_500.fa -length 500 >> $log
  java -jar ${progRemovShortSeq} -i ${outname}_scafSeq.fa -o ${outname}_scafSeq_1000.fa -length 1000 >> $log
  
  
  #get statistics
  echo ${outname}_scafSeq.fa >> $log
  java -jar ${progFastaStats} -i ${outname}_scafSeq.fa -min 200 >> $log
  java -jar ${progFastaStats} -i ${outname}_scafSeq.fa -min 500 >> $log
  java -jar ${progFastaStats} -i ${outname}_scafSeq.fa -min 1000 >> $log
  
  
  
  
#  # FINAL STEP -- IN ORIGINAL PIPELINE
#  # MAP READS AGAINST THESE SCAFFOLDS
#  conda deactivate
#  conda activate mapping_etc
#  
#  bwa index ${outname}_scafSeq.fa
#  bwa mem -t $NThreads -M ${outname}_scafSeq.fa \
#   ${workPath}/filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
#   ${workPath}/filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz | \
#   samtools fixmate -@ $NThreads -m -u -O bam - - | \
#   samtools sort -@ $NThreads -u - | \
#   samtools markdup -@ $NThreads - ${outname}_all.sorted.bam
#  samtools index ${outname}_all.sorted.bam
#  
#  
#  #get bamtools stats
#  bamtools stats -in ${outname}_all.sorted.bam >> $log
#  echo "--> ${outname}_all.sorted.bam" >> $log
#  
#  #filter for mapping quality >=10
#  samtools view -b -F 4 -q 10 ${outname}_all.sorted.bam > ${outname}_filtered.sorted.bam
#  
#  #get bamtools stats on filtered reads
#  bamtools stats -in ${outname}_filtered.sorted.bam >> $log
#  echo "--> ${outname}_filtered.sorted.bam" >> $log
  
  
  
  #IMPROVEMENT -- ASSEMBLE NEW SCAFFOLDS IN SPADES WITH THESE AS UNTRUSTED CONTIGS
  #original merged reads
  echo "7b. perform de novo assembly using skeleton scaffolds as untrusted contigs..." 
  echo "7b. perform de novo assembly using skeleton scaffolds as untrusted contigs..." >> $log

  mkdir ${workPath}/final_spades_assembly
  python $progSPAdes \
   -1 ${workPath}/filtered_reads/${outname}_trimmed_L001_R1_001.fastq.gz \
   -2 ${workPath}/filtered_reads/${outname}_trimmed_L001_R2_001.fastq.gz \
   --untrusted-contigs ${outname}_scafSeq.fa \
   -t $NThreads -o ${workPath}/final_spades_assembly

  
  
  # final assembly scaffolding with RagTag
  # this only works with python 3.5, so need separate conda environment to function
  conda deactivate
  conda activate ragtag
  cd  ${workPath}/final_spades_assembly
  mkdir ${workPath}/final_genome_assembly
  ragtag.py scaffold -t $NThreads $refRed scaffolds.fasta -o ${workPath}/final_genome_assembly
  
  
  
# 8. Step: Find BUSCO scores from the genome
  echo "8. Find BUSCO scores..." 
  echo "8. Find BUSCO scores..." >> $log
  
  #load conda environment with busco installed
  conda deactivate
  conda activate busco
  cd ${workPath}/final_genome_assembly
  
  #run busco
  busco -c $NThreads -i ragtag.scaffold.fasta -l eudicots_odb10 --contig_break 100 \
   --download_path ${workPathFiles} -o busco_output -m genome 
  

# 9. File consolidation and cleanup
  echo "9. file consolidation and cleanup..." 
  echo "9. file consolidation and cleanup..." >> $log
  
  #make final output directory
  mkdir ${workPathFiles}/${name}_final_output
  
  #move ragtag and busco info to final output
  mv ragtag.scaffold.fasta ${workPathFiles}/${name}_final_output/${outname}_refdenovo-genome.fasta
  mv ragtag.scaffold.stats ${workPathFiles}/${name}_final_output/${outname}_refdenovo-genome.stats
  mv ragtag.scaffold.confidence.txt ${workPathFiles}/${name}_final_output/${outname}_refdenovo-genome.confidence.txt
  mv busco_output/short_summary.specific.eudicots_odb10.busco_output.txt ${workPathFiles}/${name}_final_output/${outname}_busco_output.txt
  
  #move fastqc reports to final output
  mv ${workPath}/fastqc-out_filtered ${workPathFiles}/${name}_final_output
  mv ${workPath}/fastqc-out_rawreads ${workPathFiles}/${name}_final_output
  
  #move filtered reads to final output, in case they are wanted for future mapping/variant calling
  cd ${workPath}/filtered_reads
  mv *.fastq.gz ${workPathFiles}/${name}_final_output


# 10. QUAST final scaffold stats
  echo "10. QUAST final scaffold stats..." 
  echo "10. QUAST final scaffold stats..." >> $log
  
  #run quast
  cd ${workPathFiles}/${name}_final_output
  python $progQuast -o quast_scaffold_stats -t $NThreads ${outname}_refdenovo-genome.fasta
  
  mv $log ${workPathFiles}/${name}_final_output

