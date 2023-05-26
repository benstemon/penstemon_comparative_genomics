## Overview
This readme is a work-in-progress. Java scripts and, more generally, pipeline conceptualization is forked from (https://bitbucket.org/HeidiLischer/refguideddenovoassembly_pipelines/src/master/)

The scripts perform the reference-guided de novo assembly, with modifications, using SPAdes as the main assembler. All scripts needed to perform the assembly are included within this repo.

Brief overview:
1. Prepare the reference genome with [`1.prepare_pipeline.sh`](1.prepare_pipeline.sh). This removes all scaffolds  from the reference genome < specified amount (default = 10kb)

2. The main script: [`2.refGuidedDeNovoAssembly_SPAdes.sh`](2.refGuidedDeNovoAssembly_SPAdes.sh). There are several parameters that need adjusting, including the paths to programs, raw reads, etc. Ensure that all of the necessary programs are installed in conda environments or as executables, as directed.

3. Annotate genomes with GeMoMa: [`3.annotate_genome.sh`](3.annotate_genome.sh). This generates liftover annotations from desired set of reference annotations onto the DeNovo genome.