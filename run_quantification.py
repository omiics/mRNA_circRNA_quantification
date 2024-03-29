import argparse
import sys
import os
import warnings
from pathlib2 import Path
from plumbum.cmd import STAR, samtools, featureCounts, perl, bwa, bowtie2
from plumbum.cmd import python2, grep, gzip

def star_mapping(fq1, fq2, bam_prefix, genomeDir, annotation, threads=1):
    
    bam_file = "{}.Aligned.sortedByCoord.out.bam".format(bam_prefix)
    
    STAR['--runThreadN', threads, '--genomeDir', genomeDir, 
         '--sjdbGTFfile', annotation, '--sjdbOverhang', 79,
         '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
         '--outSAMtype', 'BAM', 'SortedByCoordinate',
         '--outFileNamePrefix', bam_prefix+"."].run()

    samtools['index', bam_file].run()

    return bam_file


def run_featureCounts(bamfile, output_prefix, annotation, threads=1):
    
    output_files = []

    for i in range(0,2):
        output_files.append("{}.s{}_overlap".format(output_prefix, i))
        featureCounts['-T', threads, '-p', '-C', '-O', '-s', i,
                      '-a', annotation, '-o', output_files[-1], bamfile].run()

    return output_files


def CIRI2(fq1, fq2, output_prefix, scripts_path, annotation, reference_fasta, BWAIndex, keep_temp, threads=1):
    
    bwa_sam = "{}_ciri2_temp.sam".format(output_prefix)
    
    print("\tRunning bwa")
    
    (bwa['mem', '-t', threads, BWAIndex, fq1, fq2] > bwa_sam).run()
    
    output_circRNA = "{}.CIRI2.circRNAs.txt".format(output_prefix)
    log_file = "{}.CIRI2.log".format(output_prefix)
    
    print("\tRunning CIRI2.pl")
    
    perl[os.path.join(scripts_path, "CIRI2.pl"), '-I', bwa_sam, '-O', output_circRNA,
         '-F', reference_fasta, '-A', annotation, '-G', log_file].run()
    
    if not keep_temp:
        # Delete temp files
        os.remove(bwa_sam)
    
    return output_circRNA
    

def find_circ(fq1, fq2, output_prefix, scripts_path, BowtieIndex, chromoDir, keep_temp, threads=1):
    
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
    
    if not keep_temp:
        # Delete temp file
        os.remove(temp_bam)
    
    print("\tProducing anchors")
    
    anchors = "{}.anchors.fq.gz".format(output_prefix)
    
    (python2[os.path.join(scripts_path, "find_circ", "unmapped2anchors.py"), unmapped_bam] | gzip > anchors).run()
    
    if not keep_temp:
        # Delete temp file
        os.remove(unmapped_bam)
    
    print("\tRunning the final bowtie")
    
    bowtie_log = "{}.bt2_finalpass.log".format(output_prefix)
    log_sites = "{}.sites.log".format(output_prefix)
    sites_bed = "{}.sample.sites.bed".format(output_prefix)
    sites_reads = "{}.sample.sites.reads".format(output_prefix)
    
    bowtie_call = (bowtie2['-p', threads, '--reorder', '--mm', '-M20', 
                          '--score-min=C,-15,0', '-q', '-x', BowtieIndex, 
                          '-U', anchors] >= bowtie_log) | (python2[os.path.join(scripts_path, "find_circ", "find_circ_v2.py"),
                                                   '-G', chromoDir, '-p', output_prefix, 
                                                   '-s', log_sites] >= sites_reads) > sites_bed
    bowtie_call.run()
    
    if not keep_temp:
        # Delete temp file
        os.remove(anchors)
    
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
    
    threshold_call.run()
    
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
    
    threshold_call.run()
    
    return circ_candidates, circ_candidates_35x2
    
    
class ArgumentCaller():
    
    __version__ = "2.1.0"

    def __init__(self):
        print("mRNA and circRNA quantification (v{})".format(self.__version__))
        parser = argparse.ArgumentParser(
            usage = "python run_quantification.py [<args>] <fastq1> <fastq2>"
        )

        parser.add_argument("run_mode", choices=["full", "STAR", "circRNA", "CIRI2", "find_circ"], help="Choose which pipelines should be run.")
        parser.add_argument("fastq1", help="First gzip fastq file from paired end sequencing")
        parser.add_argument("fastq2", help="Second gzip fastq file from paired end sequencing")
        parser.add_argument("--output-prefix", '-o', default="./output", help="Output prefix path default is current directory and output: './output'")
        parser.add_argument("--STARindex", '--star', help="STAR genome index directory")
        parser.add_argument("--mRNA-annotation", '-a', help="Annotation .gtf file used to quantify mRNAs")
        parser.add_argument("--BWAindex", '--bwa', help="BWAindex prefix")
        parser.add_argument("--Bowtie2Index", '--bt2', help="Bowtie2Index prefix")
        parser.add_argument("--reference-fasta", '-r', help="Reference fasta file")
        parser.add_argument("--Chromosome-dir", help="Provide a directory with chromosomes")
        parser.add_argument("--script-path", '-s', default="~/mRNA_circRNA_quantification", help="Path to where the scripts are installed")
        parser.add_argument("--create-path", action="store_true", help="Should the output prefix path be automatically created")
        parser.add_argument("--threads", '-t', default=1, type=int, help="How many threads should be used in all of the external programs")
        parser.add_argument("--keep-temp", action="store_true", help="Specify that temporary files should be kept.")

        args = parser.parse_args(sys.argv[1:])

        # Check that the fastq1 and fastq2 look like .fq.gz files
        self.check_fastq(args.fastq1)
        self.check_fastq(args.fastq2)
        # Check that the path to the output_prefix exists
        self.check_output_prefix(args.output_prefix, args.create_path)
        # Check that the STARindex provided is a directory and contains the standard filenames
        if args.run_mode in ["full", "STAR"]:
            self.check_STARindex(args.STARindex)
        # Check that the BWAindex follows the standard naming convention
        if args.run_mode in ["full", "circRNA", "CIRI2"]:
            self.check_BWAindex(args.BWAindex)
        # Check that the Bowtie2Index follows the standard naming convention
        if args.run_mode in ["full", "circRNA", "find_circ"]:
            self.check_Bowtie2Index(args.Bowtie2Index)
        # Check the script-path
        if args.run_mode in ["full", "circRNA", "CIRI2", "find_circ"]:
            self.check_script_path(args.script_path)
        
        # Check the annotation
        if args.run_mode in ["full", "STAR", "circRNA", "CIRI2"]:
            self.check_annotation(args.mRNA_annotation)
        
        # Check the Chromosome-dir
        if args.run_mode in ["full", "circRNA", "find_circ"]:
            self.check_chromosome_dir(args.Chromosome_dir)
        
        # Check the reference fasta
        if args.run_mode in ["full", "circRNA", "CIRI2"]:
            self.check_reference_fasta(args.reference_fasta)

        if args.run_mode in ["full", "STAR"]:
            print("Running STAR mapping")
            bam_file = star_mapping(args.fastq1, args.fastq2, args.output_prefix, 
                                    genomeDir=args.STARindex, annotation=args.mRNA_annotation, threads=args.threads)

        if args.run_mode in ["full", "STAR"]:
            print("Running featureCounts")
            feature_count_files = run_featureCounts(bam_file, args.output_prefix,
                                                annotation=args.mRNA_annotation, threads=args.threads)
        
        if args.run_mode in ["full", "circRNA", "CIRI2"]:
            print("Running CIRI2")
            output_circRNA = CIRI2(args.fastq1, args.fastq2, args.output_prefix, args.script_path,
                                   annotation=args.mRNA_annotation, reference_fasta=args.reference_fasta,
                                   BWAIndex=args.BWAindex, keep_temp=args.keep_temp, threads=args.threads)
        
        if args.run_mode in ["full", "circRNA", "find_circ"]:
            print("Running find_circ")
            circ_candidates, circ_candidates_35x2 = find_circ(args.fastq1, args.fastq2, 
                                                             args.output_prefix, args.script_path,
                                                             BowtieIndex=args.Bowtie2Index, chromoDir=args.Chromosome_dir,
                                                             keep_temp=args.keep_temp,
                                                             threads=args.threads)
        
        # Print the location all of the output files
        print("Finished")
        
        if args.run_mode in ["full", "STAR"]:
            print("Output featureCounts:")
            for fc_file in feature_count_files:
                print("\t{}".format(fc_file))
        
        if args.run_mode in ["full", "circRNA", "CIRI2"]:
            print("Output CIRI2 circRNA: {}".format(output_circRNA))
        
        if args.run_mode in ["full", "circRNA", "find_circ"]:
            print("Output from find_circ:")
            print("\t{}".format(circ_candidates))
            print("\t{}".format(circ_candidates_35x2))
        
    def check_fastq(self, fastq):
        if not (fastq.endswith(".fq.gz") or fastq.endswith(".fastq.gz")):
            raise Exception("Provided fastq file '{}' does not endwith .fq.gz or .fastq.gz, is this a gzip fastq file?".format(fastq))
        if not os.path.exists(fastq):
            raise Exception("Provided fastq file '{}' does not exist!".format(fastq))
            
    def check_output_prefix(self, output_prefix, create_path=False):
        
        dirname = os.path.dirname(output_prefix)
        if create_path:
            path = Path(dirname)
            path.mkdir(parents=True, exist_ok=True)
        else:
            if not os.path.exists(dirname):
                raise Exception("Path to output_prefix {} does not exist. If it should please add the --create-path option".format(dirname))
    
    def check_STARindex(self, STARindex):
        
        if STARindex is None:
            print("No STARindex has been provided. Please provide a directory with STARindex to be mapped against.")
            print("Use --STARindex or --star and provide the directory where STARindex is stored.")
            sys.exit(1)
        
        if not os.path.exists(STARindex):
            raise Exception("Provided STARindex {} does not exists!".format(STARindex))
            
        if not os.path.isdir(STARindex):
            raise Exception("Provided STARindex {} is not a directory!".format(STARindex))
        
        # Check if specific files exist in the directory
        files = ["Genome", "SA", "SAindex"]
        for filename in files:
            if not os.path.exists(os.path.join(STARindex, filename)):
                print("The provided STARindex {} does not look like a valid STARindex.".format(STARindex))
                print("It is missing the following essential file: {}".format(filename))
                
    
    def check_BWAindex(self, BWAindex):
        
        if BWAindex is None:
            print("No BWAindex has been provided. Please provide a prefix to BWAindex without the extensions.")
            print("Use --BWAindex or --bwa and provide the prefix BWAindex  files.")
            sys.exit(1)
        
        if os.path.exists(BWAindex):
            if os.path.isdir(BWAindex):
                raise Exception("Provided BWAindex {} is a directory! Please provide a prefix to the ".format(BWAindex))
            # Does it end with one of the extensions
            extensions = [".amb", ".ann", ".bwt", ".fa", ".pac", ".sa"]
            for extension in extensions:
                if BWAindex.endswith(extension):
                    print("Looks like you provided one of the files in BWAindex {}".format(BWAindex))
                    print("Please only provide the prefix to the index files: {}".format(BWAindex[:-len(extension)]))
                    sys.exit(1)
            raise Exception("Provided a file {} that exists and does not look like it belongs in a BWAindex".format(BWAindex))
        
        # Check if specific files exist in the directory
        extensions = [".bwt", ".sa", ".amb"]
        for extension in extensions:
            if not os.path.exists(BWAindex+extension):
                print("The provided BWAindex {} does not look like a valid BWAindex.".format(BWAindex))
                print("It is missing the following essential file: {}".format(BWAindex+extension))
                sys.exit(1)
    
    def check_Bowtie2Index(self, Bowtie2Index):
        
        if Bowtie2Index is None:
            print("No Bowtie2Index has been provided. Please provide a prefix to Bowtie2Index without the extensions.")
            print("Use --Bowtie2Index or --bt2 and provide the prefix to Bowtie2Index files.")
            sys.exit(1)
        
        if os.path.exists(Bowtie2Index):
            if os.path.isdir(Bowtie2Index):
                raise Exception("Provided Bowtie2Index {} is a directory! Please provide a prefix to the ".format(BWAindex))
            # Does it end with one of the extensions
            extensions = [".fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.1.bt2"]
            for extension in extensions:
                if Bowtie2Index.endswith(extension):
                    print("Looks like you provided one of the files in Bowtie2Index {}".format(Bowtie2Index))
                    print("Please only provide the prefix to the index files: {}".format(Bowtie2Index[:-len(extension)]))
                    sys.exit(1)
            raise Exception("Provided a file {} that exists but does not look like it belongs in a Bowtie2Index".format(Bowtie2Index))
        
        # Check if specific files exist in the directory
        extensions = [".1.bt2", ".rev.1.bt2"]
        for extension in extensions:
            if not os.path.exists(Bowtie2Index+extension):
                print("The provided Bowtie2Index {} does not look like a valid Bowtie2Index.".format(Bowtie2Index))
                print("It is missing the following essential file: {}".format(Bowtie2Index+extension))
                sys.exit(1)
    
    def check_script_path(self, script_path):
        
        if not os.path.exists(script_path):
            raise("Script-path {} does not exists! Unable to run CIRI2 and find_circ".format(script_path))
        
        if not os.path.isdir(script_path):
            raise("Script-path {} is not a directory, please provide a directory with the CIRI2.pl script and the find_circ scripts".format(script_path))
        
        if not os.path.exists(os.path.join(script_path, "CIRI2.pl")):
            raise("CIRI2.pl does not exist in the Script-path {}".format(script_path))
        
        if os.path.exists(os.path.join(script_path, "find_circ")):
            if not os.path.isdir(os.path.join(script_path, "find_circ")):
                raise("find_circ is not a directory in the Script-path")
            filenames = ['find_circ_v2.py', 'sum.py', 'scorethresh.py', 'unmapped2anchors.py']
            for filename in filenames:
                if not os.path.exists(os.path.join(script_path, "find_circ", filename)):
                    raise("Script {} does not exist in {}/find_circ!".format(filename, script_path))
        else:
            raise("find_circ does not exist in the Script-path {}".format(script_path))
        
    def check_annotation(self, annotation_file):
        
        if annotation_file is None:
            print("No mRNA annotation file has been provided. Please provide a gtf file.")
            print("Use --mRNA-annotation or -a to provide a gtf file with annotations")
            sys.exit(1)
            
        if not os.path.exists(annotation_file):
            raise("mRNA annotation file {} does not exist!".format(annotation_file))
        
    def check_chromosome_dir(self, chromosome_dir):
        
        if chromosome_dir is None:
            print("No chromosome dir has been provided please provide a directory with all of the chromosomes in individual fasta files.")
            print("Or provide a reference fasta file")
            print("Use --Chromosome-dir to provide such a directory.")
            sys.exit(1)
            
        if not os.path.exists(chromosome_dir):
            raise("Provided chromosome dir {} does not exist!".format(chromosome_dir))
        
        if not (chromosome_dir.endswith(".fasta") or chromosome_dir.endswith(".fa")):
            # Check if it is a directory
            if not os.path.isdir(chromosome_dir):
                raise("Provided chromosome dir {} is not a directory".format(chromosome_dir))
        
    def check_reference_fasta(self, reference_fasta):
        
        if reference_fasta is None:
            print("No reference fasta file has been provided!")
            print("Use --reference-fasta or -r to provide a reference fasta file")
            sys.exit(1)
        
        if not os.path.exists(reference_fasta):
            raise("Provided reference fasta file {} does not exist!".format(reference_fasta))
        

if __name__=="__main__":
    ArgumentCaller()
