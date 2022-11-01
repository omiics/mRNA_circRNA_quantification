# mRNA and circRNA quantification

## Dependencies

The script is implemented in python2.7, but can be run as python3 if the import 'from pathlib2 import Path' is changed to 'from pathlib import Path'.

The reason it is implemented to work for python2.7 is that the find_circ scripts are execlusively in python2.7.

Included is a conda_environment.yml which can be used to see environment this software was tested with and what versions everything was run with.

To get it to work the following python2 packages should be installed:
* For the python script:
    * argparse
    * plumbum
    * pathlib2
* For find_circ:
    * pysam
    * numpy

For CIRI2.pl you need to have atleast perl v5.8

External software:
* STAR
* samtools
* subread (For featureCounts)
* bwa
* bowtie2

## Installation

To install the conda environment we suggest to use the yml file:
```
conda create -f conda_environment.yml
conda activate mRNA_circRNA_quant
```

Then you need to install the bioconda version of argparse (which matches the python2.7 verison).
```
conda install -c bioconda argparse
```
or pip install
```
pip install argparse
```

## Help message

```
mRNA and circRNA quantification (v2.1.0)
usage: python run_quantification.py [<args>] <fastq1> <fastq2>

positional arguments:
  {full,STAR,circRNA,CIRI2,find_circ}
                        Choose which pipelines should be run.
  fastq1                First gzip fastq file from paired end sequencing
  fastq2                Second gzip fastq file from paired end sequencing

optional arguments:
  -h, --help            show this help message and exit
  --output-prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Output prefix path default is current directory and
                        output: './output'
  --STARindex STARINDEX, --star STARINDEX
                        STAR genome index directory
  --mRNA-annotation MRNA_ANNOTATION, -a MRNA_ANNOTATION
                        Annotation .gtf file used to quantify mRNAs
  --BWAindex BWAINDEX, --bwa BWAINDEX
                        BWAindex prefix
  --Bowtie2Index BOWTIE2INDEX, --bt2 BOWTIE2INDEX
                        Bowtie2Index prefix
  --reference-fasta REFERENCE_FASTA, -r REFERENCE_FASTA
                        Reference fasta file
  --Chromosome-dir CHROMOSOME_DIR
                        Provide a directory with chromosomes
  --script-path SCRIPT_PATH, -s SCRIPT_PATH
                        Path to where the scripts are installed
  --create-path         Should the output prefix path be automatically created
  --threads THREADS, -t THREADS
                        How many threads should be used in all of the external
                        programs
  --keep-temp           Specify that temporary files should be kept.

```

## Examples on how to run it

The script can be used to run 3 different tools.

* STAR + featureCounts
* CIRI2
* find_circ

The first input to the script is which "run mode" should be used by the script. The possible options are as follows:

* full
  - Run the full script with all of the tools
* circRNA
  - Run all circRNA tools CIRI2 and find_circ
* STAR
  - Only run star
* CIRI2
  - Only run CIRI2
* find_circ
  - Only run find_circ

For each of the modes you only need to provide what is required for that script to run. Below are examples with the required inputs for each mode.

The input fastq files can be labelled with .fq.gz or .fastq.gz

### Full

For the full script you need to provide all of the inputs.

The script only runs with sample at a time, where it needs specified all of the different index files, and annotations for running the processes. To speed up the run you can use -t or --threads to provide a number of threads that the external programs should use.

--create-path can be used to automatically make the path provided in the output-prefix to create data divided into subdirectories.

```
python2 run_quantification.py full test_sample_1.fq.gz test_sample_2.fq.gz \
    -o output/test_sample/test_sample -s path/to/scripts \ 
    --star path/to/STARIndex \
    --bwa path/to/BWAIndex/GRCh37.p13.genome \
    --bt2 path/to/Bowtie2Index/GRCh37.p13.genome \
    -a path/to/annotations/gencode.v37lift37.annotation.gtf \
    --Chromosome-dir path/to/human/Chromosomes \
    -r /mnt/storage/genomes/human/GRCh37.p13.genome.fa \
    --create-path -t 10
```

### STAR

Using STAR, you are only required to provide the STAR index and an annotation file.

```
python2 run_quantification.py STAR input1.fq.gz input2.fq.gz \
    -o outputs/sample_name/sample_name \
    --star path/to/STARIndex \
    -a path/to/annotations/gencode.v37lift37.annotation.gtf \
    --create-path -t 10
```

### circRNA

For running both of the circRNA tools, you need nearly all of the inputs except the STAR index.

```
python2 run_quantification.py circRNA input1.fq.gz input2.fq.gz \
    -o output/sample/test_sample -s path/to/scripts \
    --bwa path/to/BWAIndex/GRCh37.p13.genome \
    --bt2 path/to/Bowtie2Index/GRCh37.p13.genome \
    -a path/to/annotations/gencode.v37lift37.annotation.gtf \
    --Chromosome-dir path/to/human/Chromosomes \
    -r /mnt/storage/genomes/human/GRCh37.p13.genome.fa \
    --create-path -t 10
```

### CIRI2

For running CIRI2 you need the script_path, bwa index, annotations file and the reference file.

```
python2 run_quantification.py circRNA input1.fq.gz input2.fq.gz \
    -o output/sample/test_sample -s path/to/scripts \
    --bwa path/to/BWAIndex/GRCh37.p13.genome \
    -a path/to/annotations/gencode.v37lift37.annotation.gtf \
    -r /mnt/storage/genomes/human/GRCh37.p13.genome.fa \
    --create-path -t 10
```

### find_circ

For running find circ you need the script_path, bowtie2 index and a Chromosome directory with each individual chromosome in its own fasta file.

```
python2 run_quantification.py circRNA input1.fq.gz input2.fq.gz \
    -o output/sample/test_sample -s path/to/scripts \
    --bt2 path/to/Bowtie2Index/GRCh37.p13.genome \
    --Chromosome-dir path/to/human/Chromosomes \
    --create-path -t 10
```