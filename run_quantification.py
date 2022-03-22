from __future__ import annotations
import argparse
import sys
import os
from pathlib import Path
from plumbum.cmd import STAR, samtools, featureCounts, perl, bwa, bowtie2
from plumbum.cmd import python2, grep, gzip

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


def CIRI2(fq1, fq2, output_prefix, scripts_path, annotation, reference_fasta, BWAIndex, threads=1):
    
    bwa_sam = "{}_ciri2_temp.sam".format(output_prefix)
    
    print("\tRunning bwa")
    
    (bwa['mem', '-t', threads, BWAIndex, fq1, fq2] > bwa_sam).run()
    
    output_circRNA = "{}.CIRI2.circRNAs.txt".format(output_prefix)
    log_file = "{}.CIRI2.log".format(output_prefix)
    
    print("\tRunning CIRI2.pl")
    
    perl[os.path.join(scripts_path, "CIRI2.pl"), '-I', bwa_sam, '-O', output_circRNA,
         '-F', reference_fasta, '-A', annotation, '-G', log_file].run()
    
    return output_circRNA
    

def find_circ(fq1, fq2, output_prefix, scripts_path, BowtieIndex, chromoDir, threads=1):
    
    temp_bam = "{}.mapped.bam".format(output_prefix)
    log_file = "{}.bt2_firstpass.log".format(output_prefix)
    
    print("\tRunning bowtie2")
    
    bowtie_call = (bowtie2['-p', threads, '--very-sensitive', '--mm', '-M20', 
                          '--score-min=C,-15,0', '-x', BowtieIndex, '-q', 
                          '-1', fq1, '-2', fq2] >= log_file) | samtools['sort', '-'] > temp_bam
    bowtie_call.run()
    
    print("\tExtracting all unmapped reads")
    
    unmapped_bam = "{}.unmapped.bam".format(output_prefix)
    
    (samtools['view', '-hf', 4, temp_bam] | samtools['view', '-Sb', '-'] > unmapped_bam).run()
    
    print("\tProducing anchors")
    
    anchors = "{}.anchors.fq.gz".format(output_prefix)
    
    (python2[os.path.join(scripts_path, "find_circ", "unmapped2anchors.py"), unmapped_bam] | gzip > anchors).run()
    
    print("\tRunning the final bowtie")
    
    log_sites = "{}.sites.log".format(output_prefix)
    sites_bed = "{}.sample.sites.bed".format(output_prefix)
    sites_reads = "{}.sample.sites.reads".format(output_prefix)
    
    bowtie_call = bowtie2['-p', threads, '--reorder', '--mm', '-M20', 
                          '--scores-min=C,-15,0', '-x', BowtieIndex, 
                          '-U', anchors] | python2[os.path.join(scripts_path, "find_circ", "unmapped2anchors.py"),
                                                   '-G', chromoDir, '-p', output_prefix, 
                                                   '-s', log_sites] > sites_bed >= sites_reads
    bowtie_call.run()
    
    circ_candidates = "{}.circ_candidates.bed".format(output_prefix)
    
    print("\t Applying thresholds")
    threshold_call = grep['circ_', sites_bed] | grep['-v', 'chrM']
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "sum.py"), '-2,3']
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-20', 1]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-19', 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-18', 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), 7, 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '9,10', 35]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-21', 10000]
    threshold_call = threshold_call > circ_candidates
    
    circ_candidates_35x2 = "{}.circ_candidates_map35x2.bed".format(output_prefix)
    
    print("\t Applying thresholds")
    threshold_call = grep['circ_', sites_bed] | grep['-v', 'chrM']
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "sum.py"), '-2,3']
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-20', 1]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-19', 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-18', 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), 7, 2]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '9,10', 35]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '10', 35]
    threshold_call = threshold_call | python2[os.path.join(scripts_path, "find_circ", "scorethresh.py"), '-21', 10000]
    threshold_call = threshold_call > circ_candidates_35x2
    
    return circ_candidates, circ_candidates_35x2
    
    
class ArgumentCaller():
    
    __version__ = "1.0.0"

    def __init__(self):
        print("mRNA and circRNA quantification (v{})".format(self.__version__))
        parser = argparse.ArgumentParser(
            usage = "python run_quantification.py [<args>] <fastq1> <fastq2>"
        )

        parser.add_argument("fastq1", help="First gzip fastq file from paired end sequencing")
        parser.add_argument("fastq2", help="Second gzip fastq file from paired end sequencing")
        parser.add_argument("--output_prefix", '-o', default="./output", help="Output prefix path default is current directory and output: './output'")
        parser.add_argument("--STARindex", '--star', help="STAR genome index directory")
        parser.add_argument("--mRNA-annotation", '-a', help="Annotation .gtf file used to quantify mRNAs")
        parser.add_argument("--BWAindex", '--bwa', help="BWAindex prefix")
        parser.add_argument("--Bowtie2Index", '--bt2', help="Bowtie2Index prefix")
        parser.add_argument("--reference-fasta", '-r', help="Reference fasta file")
        parser.add_argument("--Chromosome-dir", help="Provide a directory with chromosomes")
        parser.add_argument("--script-path", '-s', default="~/mRNA_circRNA_quantification", help="Path to where the scripts are installed")
        parser.add_argument("--threads", '-t', default=1, type=int, help="How many threads should be used in all of the external programs")

        args = parser.parse_args(sys.argv[1:])

        # Check that the fastq1 and fastq2 look like .fq.gz files
        # Check that the path to the output_prefix exists
        # Check that the STARindex provided is a directory and contains the standard filenames
        # Check that the BWAindex follows the standard naming convention
        # Check that the Bowtie2Index follows the standard naming convention
        
        # Check all of the software is available

        print("Running STAR mapping")
        bam_file = star_mapping(args.fastq1, args.fastq2, args.output_prefix, 
                                genomeDir=args.STARindex, annotation=args.mRNA_annotations, threads=args.threads)

        print("Running featureCounts")
        feature_count_files = featureCounts(bam_file, args.output_prefix,
                                            annotation=args.mRNA_annotations, threads=args.threads)
        
        print("Running CIRI2")
        output_circRNA = CIRI2(args.fastq1, args.fastq2, args.output_prefix, args.scripts_path,
                               annotation=args.mRNA_annotations, reference_fasta=args.reference_fasta,
                               BWAIndex=args.BWAindex, threads=args.threads)
        
        print("Running find_circ")
        circ_candidates, circ_candidates_35x2 = find_circ(args.fastq1, args.fastq2, 
                                                          args.output_prefix, args.scripts_path,
                                                          BowtieIndex=args.Bowtie2Index, chromoDir=args.Chromosome_dir,
                                                          threads=args.threads)
        
        
        # Print the location all of the output files
        print("Finished")
        
        print("Output featureCounts:")
        for fc_file in feature_count_files:
            print("\t{}".format(fc_file))
        
        print("Output CIRI2 circRNA: {}".format(output_circRNA))
        
        print("Output from find_circ:")
        print("\t{}".format(circ_candidates))
        print("\t{}".format(circ_candidates_35x2))
        

if __name__=="__main__":
    ArgumentCaller()
