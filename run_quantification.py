from __future__ import annotations
import argparse
import sys
import os
from pathlib import Path
from plumbum.cmd import STAR, samtools, featureCounts

def star_mapping(fq1, fq2, bam_prefix, genomeDir, annotation, threads=1):
    
    bam_file = "{}.Aligned.sortedByCoord.out.bam".format(bam_prefix)
    
    STAR['--runThreadN', threads, '--genomeDir', genomeDir, 
         '--sjdbGTFfile', annotation, '--sjdbOverhang', 79,
         '--readFilesIn', fq1, fq2, '--readFileCommand', 'zcat',
         '--outSAMtype', 'BAM', 'SortedByCoordinate',
         '--outFileNamePrefix', bam_prefix].run()

    samtools['index', bam_file].run()

    return bam_file

def featureCounts(bamfile, output_prefix, annotation, threads=1):
    
    output_files = []

    for i in range(0,2):
        output_files.append("{}.s{}_overlap".format(output_prefix, i))
        featureCounts['-T', threads, '-p', '-C', '-O', '-s', i,
                      '-o', output_files[-1], bamfile]

    return output_files

def CIRI2():
    pass

def find_circ():
    pass


class ArgumentCaller():
    
    __version__ = "1.0.0"

    def __init__(self):
        print("mRNA and circRNA quantification (v{})".format(self.__version__))
        parser = argparse.ArgumentParser(
            usage = "python run_quantification.py [<args>] <fastq1> <fastq2>"
        )

        parser.add_argument("fastq1", help="First gzip fastq file from paired end sequencing")
        parser.add_argument("fastq2", help="Second gzip fastq file from paired end sequencing")
        parser.add_argument("--output_prefix", default="./output", help="")
        parser.add_argument("--STARindex", help="STAR genome index directory")
        parser.add_argument("--mRNA-annotation", help="Annotation .gtf file used to quanit")
        parser.add_argument("--threads", default=1, type=int, help="How many threads should be used in all of the external programs")

        args = parser.parse_args(sys.argv[1:])

        # Check that the fastq1 and fastq2 look like .fq.gz files
        # Check that the path to the output_prefix exists
        # Check that the STARindex provided is a directory and contains the standard filenames

        print("running STAR mapping")
        bam_file = star_mapping(args.fastq1, args.fastq2, args.output_prefix, 
                                genomeDir=args.STARindex, annotation=args.mRNA_annotations, threads=args.threads)

        print("running featureCounts")
        feature_count_files = featureCounts(bam_file, args.output_prefix,
                                            annotation=args.mRNA_annotations, threads=args.threads)
        


if __name__=="__main__":
    ArgumentCaller()
