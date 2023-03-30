# 1. Install jcvi

I ended up using installation of jcvi on my laptop.
steps for installation:
```
conda env create -n jcvi
conda activate jcvi
conda install pip
pip install jcvi
```

# 2. Install LAST
Conda install will not work for some reason. Last can be found [here](https://gitlab.com/mcfrith/last). After using makefile, copy files in last/bin to software bin and add to PATH variable
```
cp * /Users/benstone/software
sudo echo /Users/benstone/software >> /etc/paths
```

# 3. Install Latex
LaTex can be found [here](https://www.latex-project.org/get/).


# 4. Preparing input files
We need genome.fa files and .gff3 files for each of the reference genomes of interest.
a. Generate mRNA-only bed file with gffread.  Input: genome.gff3
* #NOTE: smallii annotations subject to change, because they are liftovers. Current processing doesn't remove anything
```shell
#for barbatus:
gffread M4_annotation_putative_function_domain_added_blast_tomato.genemodels.noseq.gff --bed --keep-genes --sort-alpha -o M4_annotation_barbatus-genes_CMKEVH.bed

#for davidsonii:
gffread annot_Pdavidsonii_genome_FUNCTIONAL-INCLUDED.gff --bed --keep-genes --sort-alpha -o annot_Pdavidsonii_FUNCTIONAL-INCLUDED-genes_CMKEVH.bed

#for petiolatus:
gffread Rnd1.all.maker.snapdragon.noseq.gff --bed --keep-genes --sort-alpha -o petiolatus-genes_CMKEVH.bed

#for smallii:
gffread smallii_NAMECHANGE_final_annotation.gff --bed --keep-genes --sort-alpha -o smallii_NAMECHANGE_PGA-genes_CMKEVH.bed
```

b. filter the bed file to include only a single isoform (longest) with `filter_isoforms_bedfile-*species*.py`. These are species-specific scripts dependent on the gene model names given to species in their annotation files.
```
#for barbatus:
python filter_isoforms_bedfile-barbatus.py M4_annotation_barbatus-genes_CMKEVH.bed M4_annotation_barbatus-genes_CMKEVH-single-isoform.bed

#for davidsonii:
python filter_isoforms_bedfile-davidsonii.py annot_Pdavidsonii_FUNCTIONAL-INCLUDED-genes_CMKEVH.bed annot_Pdavidsonii_FUNCTIONAL-INCLUDED-genes_CMKEVH-single-isoform.bed

#for petiolatus:
python filter_isoforms_bedfile-petiolatus.py petiolatus-genes_CMKEVH.bed petiolatus-genes_CMKEVH-single-isoform.bed

#for smallii:
python filter_isoforms_bedfile-smallii.py smallii_NAMECHANGE_PGA-genes_CMKEVH.bed smallii_NAMECHANGE_PGA-genes_CMKEVH-single-isoform.bed
```


c. Generate a CDS protein fasta file and .cds nucleotide file from the filtered bed file with gffread. Input: genome.gff3, and single-isoform.bed. Additionally:
	Discard CDS with in-frame stop codons
	Paramters: -C -M -K -E -V -H.
```
#for barbatus:
gffread -y M4_annotation_barbatus_CMKEVH-single-isoform-protein.fasta -C -M -K -E -V -H --sort-alpha -g Pbar.2022.LG.fa M4_annotation_barbatus-genes_CMKEVH-single-isoform.bed

gffread -x M4_annotation_barbatus_CMKEVH-single-isoform-nucleotide.cds -C -M -K -E -V -H --sort-alpha -g Pbar.2022.LG.fa M4_annotation_barbatus-genes_CMKEVH-single-isoform.bed


#for davidsonii:
gffread -y annot_Pdavidsonii_FUNCTIONAL-INCLUDED_CMKEVH-single-isoform-protein.fasta -C -M -K -E -V -H --sort-alpha -g annot_Pdavidsonii_genome.fasta annot_Pdavidsonii_FUNCTIONAL-INCLUDED-genes_CMKEVH-single-isoform.bed

gffread -x annot_Pdavidsonii_FUNCTIONAL-INCLUDED_CMKEVH-single-isoform-nucleotide.cds -C -M -K -E -V -H --sort-alpha -g annot_Pdavidsonii_genome.fasta annot_Pdavidsonii_FUNCTIONAL-INCLUDED-genes_CMKEVH-single-isoform.bed


#for petiolatus:
gffread -y petiolatus_CMKEVH-single-isoform-protein.fasta -C -M -K -E -V -H --sort-alpha -g petiolatus_genome.fasta petiolatus-genes_CMKEVH-single-isoform.bed

gffread -x petiolatus_CMKEVH-single-isoform-nucleotide.cds -C -M -K -E -V -H --sort-alpha -g petiolatus_genome.fasta petiolatus-genes_CMKEVH-single-isoform.bed


#for smallii:
gffread -y smallii_NAMECHANGE_PGA_CMKEVH-single-isoform-protein.fasta -C -M -K -E -V -H --sort-alpha -g smallii_NAMECHANGE_PGA_assembly.fasta smallii_NAMECHANGE_PGA-genes_CMKEVH-single-isoform.bed

gffread -x smallii_NAMECHANGE_PGA_CMKEVH-single-isoform-nucleotide.cds -C -M -K -E -V -H --sort-alpha -g smallii_NAMECHANGE_PGA_assembly.fasta smallii_NAMECHANGE_PGA-genes_CMKEVH-single-isoform.bed

```
We should now have single-isoform bedfiles and the respective fasta sequence files for all of the genomes of interest. *Note that I subsequently copied these output bed files into a new directory, and changed file names to be more concise*


# 5. Running MCscan (Python version) implementation of jcvi
These steps will mirror closely those found in the jcvi [manual](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)). I am currently setting up barbatus and davidsonii as the "central" genomes, as they are probably the two most complete genomes +annotated assemblies.
```
conda activate jcvi

#create anchors files and dotplots
python -m jcvi.compara.catalog ortholog barbatus davidsonii --cscore=.99 --dbtype=nucl
python -m jcvi.compara.catalog ortholog barbatus smallii --cscore=.99 --dbtype=nucl
python -m jcvi.compara.catalog ortholog davidsonii petiolatus --cscore=.99 --dbtype=nucl


#extract subsets of blocks from anchorfile --minspan=30 --
python -m jcvi.compara.synteny screen --simple barbatus.davidsonii.anchors barbatus.davidsonii.anchors.new
python -m jcvi.compara.synteny screen --simple barbatus.smallii.anchors barbatus.smallii.anchors.new
python -m jcvi.compara.synteny screen --simple davidsonii.petiolatus.anchors davidsonii.petiolatus.anchors.new


#after organizing the seqids.txt and layout.txt files...
python -m jcvi.graphics.karyotype seqids.txt layout.txt -o karyotype_barbatus.davidsonii.petiolatus.smallii.pdf


```




