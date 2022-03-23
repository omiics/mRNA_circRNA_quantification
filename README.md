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

## Help message

```
mRNA and circRNA quantification (v1.0.0)
usage: python run_quantification.py [<args>] <fastq1> <fastq2>

positional arguments:
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
```

## Example on how to run it

The script only runs with sample at a time, where it needs specified all of the different index files, and annotations for running the processes. To speed up the run you can use -t or --threads to provide a number of threads that the external programs should use.

--create-path can be used to automatically make the path provided in the output-prefix to create data devided into subdirectories.

```
python2 run_quantification.py test_sample_1.fq.gz test_sample_2.fq.gz \
    -o output/test_sample -s path/to/scripts \ 
    --star path/to/STARIndex \
    --bwa path/to/BWAIndex/GRCh37.p13.genome \
    --bt2 path/to/Bowtie2Index/GRCh37.p13.genome \
    -a path/to/annotations/gencode.v37lift37.annotation.gtf \
    --Chromosome-dir path/to/human/Chromosomes \
    -r /mnt/storage/genomes/human/GRCh37.p13.genome.fa \
    --create-path -t 10
```