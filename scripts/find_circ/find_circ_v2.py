#!/usr/bin/env python
import os
import sys
import re
import mmap
import pysam
import numpy
import resource
from logging import debug,warning,error,getLogger
from optparse import OptionParser
from datetime import datetime


COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

class mmap_fasta(object):
    def __init__(self,fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline-self.ldata
        self.skip_char = row[self.ldata:]
        #print "SKIP",self.skip,self.skip_char
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def __getslice__(self,start,end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs
        
        #print "ofs",ofs_start,ofs_end
        
        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char,"")
        L = end-start
        if len(s) == L:
            return s
        else:
            return s+"N"*(L-len(s))
            
        #self.mmap.close ()
        #self.f.close ()
        return 

class Accessor(object):
    supports_write = False

    def get_data(self,start,end,sense):
        return []

    def get_oriented(self,start,end,sense):
        data = self.get_data(start,end,sense)
        if sense == "-": #  self.sense_specific and
            return data[::-1]
        else:
            return data

    def get_sum(self,start,end,sense):
        return self.get_data(start,end,sense).sum()

    def flush(self):
        pass

def to_bool(obj):
    if obj == "False":
        return False
    else:
        return bool(obj)
        
class Track(object):
    """
    Abstraction of chromosome-wide adressable data like sequences, coverage, scores etc.
    Actual access to the data is delegated to accessor objects which are instantiated on-the-fly for
    each chromosome (strand) upon first access and then cached.
    Use of mmap for the accessors is recommended and implemented for sequences and numpy (C-type)
    arrays.

    See io/track_accessors.py for more examples.
    """

#    def __init__(self,path,accessor,sense_specific=True,description="unlabeled track",system="hg19",dim=1,auto_flush=False,mode="r",**kwargs):
    def __init__(self,path,accessor,auto_flush,sense_specific=True,description="unlabeled track",system="hg19",dim=1,mode="r",**kwargs):
        self.path = path
        self.mode = mode
        self.acc_cache = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = to_bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = getLogger("Track('%s')" % path)
        self.system = system

        self.logger.debug("Track(auto_flush=%s)" % (str(auto_flush)))
        kwargs['sense_specific'] = self.sense_specific
        kwargs['mode'] = self.mode
        kwargs['system'] = self.system
        kwargs['description'] = self.description
        kwargs['system'] = self.system
        kwargs['dim'] = self.dim

        if "w" in mode:
            # make sure the path exists right away so that the accessors 
            # can flush the actual data there!
            from sequence_data.io import ensure_path
            ensure_path(self.path)

    def load(self,chrom,sense):

        # automatically flush buffers whenever a new chromosome is seen. reduces memory-footprint for sorted input data
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches." % chrom)
            self.flush_all()
        self.last_chrom = chrom
         
        ID = self.get_identifier(chrom,sense)
        if not ID in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor" % (chrom,sense))
            self.acc_cache[ID] = self.accessor(self.path,chrom,sense,**(self.kwargs))
            #print len(self.acc_cache)

        return self.acc_cache[ID]

    def save(self):
        """
        If a track is opened in "rw" or "w" mode this will save the track-definition config files and flush all accessors.
        Saving of the actual data is performed by the accessors that support writing upon a call to flush.
        """
        if not "w" in self.mode:
            self.logger.warning("save() called on a read-only opened track. Ignored!")
            return

        if not self.accessor.supports_write:
            self.logger.warning("save() called on a track with only read-access supporting accessors. Ignored!")
            return
      
        self.logger.debug("save(): writing '%s'" % self.path)

        def to_str(obj):
            # convert simple data-types to their string representation
            # but classes and more complex types to their names.
            return getattr(obj,"__name__",str(obj))

        kwarg_str = "\n".join(["%s=%s" % (k,to_str(self.kwargs[k])) for k in sorted(self.kwargs.keys()) if k != "mode"])
        file(os.path.join(self.path,"track.rc"),"w+").write(trackrc % dict(accessor=self.accessor.__name__,kwargs=kwarg_str))
        self.flush_all()

    def __del__(self):
        if "w" in self.mode:
            self.save()

    def flush(self,chrom,sense):
        ID = self.get_identifier(chrom,sense)
        if ID in self.acc_cache:
            self.logger.warning("Flushing %s%s" % (chrom,sense))
            del self.acc_cache[ID]

    def flush_all(self):
        for a in self.acc_cache.values():
            a.flush()
        self.acc_cache = {}

    def get(self,chrom,start,end,sense):
        #print "get", chrom, start, end, sense
        acc = self.load(chrom,sense)
        return acc.get_data(start,end,sense)

    def get_oriented(self,chrom,start,end,sense):
        acc = self.load(chrom,sense)
        return acc.get_oriented(start,end,sense)

    def get_sum(self,chrom,start,end,sense):
        acc = self.load(chrom,sense)
        return acc.get_sum(start,end,sense)
        
    def get_identifier(self,chrom,sense):
        if self.sense_specific:
            return chrom+sense
        else:
            return chrom

class GenomeAccessor(Accessor):
    def __init__(self,path,chrom,sense,**kwargs):
        debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))
        #print "init genomeAccessor",chrom,sense        
        fname = os.path.join(path,chrom+".fa")
        try:
            self.data = mmap_fasta(fname)
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only Ns)" % fname)
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented

    def get_data(self,start,end,sense):
        if start < 0 or end < 0:
            return self.get_dummy(start,end,sense)
        #UCSC convention: start with 1, end is inclusive
        if sense == "+":
            return self.data[start:end]
        else:
            return complement(self.data[start:end])

    def get_oriented(self,start,end,sense):
        if end < 0:
            return self.get_dummy(start,end,sense)
        elif start < 0:
            return self.get_dummy(start,0,sense) + self.get_oriented(0,end,sense)
        if sense == "+":
            return self.data[start:end]
        else:
            try:
                return rev_comp(self.data[start:end])
            except KeyError:
                print start,end,sense

    def get_dummy(self,start,end,sense):
        return "N"*(end-start)



usage = """

  bowtie2 anchors.qfa.gz | %prog > candidates.bed 2> candidates.reads

"""

parser = OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (alternatively use -G <path_to_genome_folder>)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome. Point to folder with one fasta file for each chromosome.")
parser.add_option("-a","--anchor",dest="asize",type=int,default=20,help="anchor size (default=20)")
parser.add_option("-w","--wiggle",dest="wiggle",type=int,default=2,help="maximum nts a linear splice site may be away from circ to be counted as competitor (default=2)")
parser.add_option("-m","--margin",dest="margin",type=int,default=2,help="maximum nts the BP is allowed to reside inside an anchor (default=2)")
parser.add_option("-d","--maxdist",dest="maxdist",type=int,default=2,help="maximum mismatches (no indels) allowed in anchor extensions (default=2)")
parser.add_option("-p","--prefix",dest="prefix",default="",help="prefix to prepend to each name")
parser.add_option("-q","--min_uniq_qual",dest="min_uniq_qual",type=int,default=2,help="minimal uniqness (alignment score margin to next-best hit) for anchor alignments (default=2)")
parser.add_option("-s","--stats",dest="stats",default="runstats.log",help="where to put stats (default='runstats.log'")
parser.add_option("-r","--reads2samples",dest="reads2samples",default="",help="path to tab-separated two-column file with read-name prefix -> sample ID mapping")
parser.add_option("-f","--autoflush",dest="autoflush",action="store_true",default=False,help="Use autoflush to reduce memory usage (default=FALSE)")

options,args = parser.parse_args()

if options.system:
    # you have the rest of the library installed and set up, great
    import sequence_data.systems
    genome = getattr(sequence_data.systems,options.system).genome
else:
    genome = Track(options.genome,GenomeAccessor,options.autoflush)
    
    
from itertools import izip_longest
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

if options.reads2samples:
    samples = [line.rstrip().split('\t') for line in file(options.reads2samples)]
   # loci_tissue  = [line.rstrip().split('\t') for line in file(options.reads2samples)]
else:
    samples = []
    
samples.append(('','unknown'))

        
minmapscore = options.asize * (-2)
from collections import defaultdict

class Hit(object):
    def __init__(self):
        self.reads = []
        self.uniq = set()
        self.uniq2 = set()
        self.mapquals = []
        self.tissues = defaultdict(int)
        self.rel_exp = defaultdict(int)
        self.edits = []
        self.overlaps = []
        self.n_hits = []
        self.shs = []
        
    def add(self,read,A,B,dist,ov,n_hits,dshs=-1):
        self.reads.append(read)
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)
        self.shs.append(dshs)

        # Alignment Score - Secondbest hit score
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS') - aopt.get('XS',minmapscore)
        qB = bopt.get('AS') - bopt.get('XS',minmapscore)
       
        self.mapquals.append((qA+qB,qA,qB))
        
        for (prefix,tiss) in samples:            
            if A.qname.startswith(prefix):
                self.tissues[tiss] += 1                
                break
    
        self.uniq.add((read,tiss))
        self.uniq.add((rev_comp(read),tiss))
        
        self.uniq2.add ((A.pos, tiss))
        self.uniq2.add ((B.pos, tiss))


    def scores(self,chrom,start,end,sense):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) / 2
        n_uniq2 = len(self.uniq2) / 2
        
        total_mq,best_qual_A,best_qual_B = sorted(self.mapquals,reverse=True)[0] # Retrieve best mapping pair (based on total_mq)

        total_A = 0;
        total_B = 0;
        
        for x in self.mapquals:
           #total_mq,best_qual_A,best_qual_B = self.mapquals[x]
           total_A += x[1]
           total_B += x[2]
           
           
        avr_qual_A = total_A / len(self.mapquals)
        avr_qual_B = total_B / len(self.mapquals)       
        
        wiggle = numpy.arange(-options.wiggle,options.wiggle+1)   
        
        spliced_at_begin = 0
        for x in wiggle:
            begin = (chrom,start+x,sense)
            if begin in loci2:
                spliced_at_begin += loci2[begin]
                #for l in loci[begin]:
                #    #spliced_at_begin += len(l.reads)
                #    spliced_at_begin += 1

        spliced_at_end = 0
        for x in wiggle:
            ending = (chrom,end+x,sense)
            if ending in loci2:
                spliced_at_end += loci2[ending]
                #for l in loci[ending]:
                #    #spliced_at_end += len(l.reads)
                #    spliced_at_end += 1

        tissues = sorted(self.tissues.keys())
        
       
                         
        #print self.edits,self.overlaps,self.n_hits
        
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        #tiss_rel_exp = [str(spliced_tissue[k]) for k in tissues]
        tiss_rel_exp = tiss_counts
        
        
        #return (n_reads,n_uniq,best_qual_A,best_qual_B,spliced_at_begin,spliced_at_end,tissues,tiss_counts,min(self.edits),min(self.overlaps),min(self.n_hits))
        return (n_reads,n_uniq,n_uniq2,best_qual_A,best_qual_B,avr_qual_A,avr_qual_B,spliced_at_begin,spliced_at_end,tissues,tiss_counts,tiss_rel_exp,min(self.edits),min(self.overlaps),min(self.n_hits))
                
        
loci = defaultdict(list)
loci2 = defaultdict(int)
#loci_tissue = defaultdict(list)
circs = defaultdict(Hit)
#circs_ol = defaultdict(Hit)
splices = defaultdict(Hit)
#splices_ol = defaultdict(Hit)
#trans = defaultdict(Hit)
tinyc = defaultdict(Hit)

N = defaultdict(int)

from numpy import chararray as carray
from numpy import fromstring,byte

def trunc_divmod(a, b):
    q = a / b
    q = -int(-q) if q<0 else int(q)
    r = a - b * q
    return q, r

def mismatches(a,b):
    if len(a) != len(b):
       return 100
    a,b = fromstring(a,dtype=byte), fromstring(b,dtype=byte)
    return (a != b).sum()
    
def find_breakpoints(A,B,read,chrom,margin=options.margin,maxdist=options.maxdist):

    
        
    L = len(read) #100
    hits = []
    #print "readlen",L
    #print " "*2+read
    eff_a = options.asize-margin  #20-2 = 18
    internal = read[eff_a:-eff_a].upper()   #from 18 to 92
        
    #margin: maximum nts the BP is allowed to reside inside an anchor (default=2)
    
    flank = L - 2*eff_a + 2  #=100 - 2*18 + 2 = 100-36+2=66
    
    #aend one past last aligned residue

    A_flank = genome.get(chrom, A.aend - margin, A.aend-margin + flank, '+').upper()
    B_flank = genome.get(chrom, B.pos - flank+margin, B.pos+margin, '+').upper()

    
    #print " "*2+A.seq[:-margin]+A_flank.lower()
    #print " "*(eff_a)+B_flank.lower()+B.seq[margin:]
    #print "testing..."
    #print " "*(eff_a+2)+internal
    l = L - 2*eff_a
    for x in range(l+1):
        spliced = A_flank[:x] + B_flank[x+2:]
        dist = mismatches(spliced,internal)        
        
        bla = A_flank[:x].lower() + B_flank[x+2:]
        #print " "*(eff_a+2)+bla,dist

        ov = 0
        if x < margin:
            ov = margin-x
        if l-x < margin:
            ov = margin-(l-x)
        
        if dist <= maxdist:
            gt = A_flank[x:x+2]
            ag = B_flank[x:x+2]
            #print x,gt,ag
            #print "MATCH", A_flank[x:x+2].lower(),B_flank[-(l-x)-2:-(l-x)].upper()
            
            start,end = B.pos+margin-l+x,A.aend-margin+x+1
            start,end = min(start,end),max(start,end)
            if gt == 'GT' and ag == 'AG':
                #print "FOUND ONE PLUS STRAND"
                #print x,L
                hits.append((dist,ov,chrom,start,end,'+'))
            elif gt == 'CT' and ag == 'AC':
                #print "FOUND ONE MINUS STRAND"
                hits.append((dist,ov,chrom,start,end,'-'))

    if len(hits) < 2:
        return hits

    # return only hits that are tied with the best candidate by edit distance and anchor overlap. 
    # Hits are still sorted, with low edit distance beating low anchor overlap
    hits = sorted(hits)
    best = hits[0]
    return [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]
    
def find_breakpoints_tiny(A,B,read,chrom,margin=options.margin,maxdist=options.maxdist):

#    def mismatches(a,b):
#        a,b = fromstring(a,dtype=byte), fromstring(b,dtype=byte)
#        return (a != b).sum()
        
    BA_dist = B.pos - A.pos
    
    #if A.is_reverse:
    #   BA_dist = A.pos - B.pos
    
        
    L = len(read) #100
    hits = []
    #print "readlen",L
    #print " "*2+read
    eff_a = options.asize-margin  #20-2 = 18
    internal = read[eff_a:-eff_a].upper()   #from 18 to 92
    
    exon_size = [0] * 5
    
    exon_size[0] = (L - BA_dist - options.asize)
    if A.is_reverse:
       exon_size[0] = exon_size[0] * -1
    
    for x in range (1,4):
    
       exon_size[x] = -1
       
       if exon_size[0]%(x+1) == 0:
          exon_size[x] = (L - BA_dist - options.asize) / (x+1)
           
       
               
    #margin: maximum nts the BP is allowed to reside inside an anchor (default=2)
    
    flank = L - 2*eff_a + 2  #=100 - 2*18 + 2 = 100-36+2=66
    
    #aend one past last aligned residue

    A_flank = genome.get(chrom, A.aend - margin, A.aend-margin + flank, '+').upper()
    B_flank = genome.get(chrom, B.pos - flank+margin, B.pos+margin, '+').upper()

    #if exon_size3 == -1 and exon_size2 == -1 and exon_size1 == -1:
       #print " "*2+A.seq[:-margin]+A_flank.lower()
       #print " "*(eff_a)+B_flank.lower()+B.seq[margin:]
       #print "testing..."
       #print " "*(eff_a+2)+internal
       
    l = L - 2*eff_a
    
    
    for e in range(0,4):
       
       if exon_size[e] == -1:
          continue
       

       #if e >= 2:
          #print exon_size,e,chrom,A.pos,B.pos,L,BA_dist,read
          #print " "*2+A.seq[:-margin]+A_flank.lower()
          #print " "*(eff_a)+B_flank.lower()+B.seq[margin:]
          #print "testing..."
          #print " "*(eff_a+2)+internal   
          
       for x in range(l+1-(exon_size[e]*e)):
       #for x in range(l+1):
           repeat_unit = ""
           repeat_seq = ""
           
           if e > 0:
              repeat_unit = internal[x:x+exon_size[e]]              
              for r in range (e):
                 repeat_seq += repeat_unit
           
           spliced = A_flank[:x] + repeat_seq + B_flank[x+2+len(repeat_seq):]
           #spliced = A_flank[:x] + B_flank[x+2:]
           
           #print spliced,internal           
           
           dist = mismatches(spliced,internal)  

           #print dist           
           bla = A_flank[:x].lower() + repeat_seq.upper () + B_flank[x+2+len(repeat_seq):].lower ()
            
           #if e >= 2:         
              #print " "*(eff_a+2)+bla,dist,repeat_seq,repeat_unit
           
           #if exon_size3 == -1 and exon_size2 == -1 and exon_size1 == -1:
              #bla = A_flank[:x].lower() + B_flank[x+2:]
              #print " "*(eff_a+2)+bla,dist

           ov = 0
           if x < margin:
               ov = margin-x
           if l-x < margin:
               ov = margin-(l-x)
           
           if dist <= maxdist:
               gt = A_flank[x:x+2]
               ag = B_flank[x+len(repeat_seq):x+2+len(repeat_seq)]
               
               #if e >= 2:
                  #print x,gt,ag
                  #print "MATCH", A_flank[x:x+2].lower(),B_flank[-(l-x)-2:-(l-x)].upper()
               
               start,end = B.pos+margin-l+x+len(repeat_seq),A.aend-margin+x+1
               start,end = min(start,end),max(start,end)
               if gt == 'GT' and ag == 'AG':
                   #print "FOUND ONE PLUS STRAND"
                   #print x,L
                   hits.append((dist,ov,chrom,start,end,'+'))
               elif gt == 'CT' and ag == 'AC':
                   #print "FOUND ONE MINUS STRAND"
                   hits.append((dist,ov,chrom,start,end,'-'))


       if len(hits) == 1:
          return hits

       if len(hits) > 1:
          # return only hits that are tied with the best candidate by edit distance and anchor overlap. 
          # Hits are still sorted, with low edit distance beating low anchor overlap
          hits = sorted(hits)
          best = hits[0]
          return [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]



    return hits 


def timedelta_total_seconds(td):
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / float(10**6)      


nlines = 100


div = 1000
line = 0
d1 = datetime.now()

sam = pysam.Samfile('-','r')
try:
    for A,B in grouper(2,sam):
        #print A
        #print B
        
        line = line+1
                                
        q, r = trunc_divmod(line, div)
        
        if r == 0:
            with open("progress.txt", "w") as text_file:
               text_file.write("Done: {0} of {1}".format(line, nlines))
               ts_pr_line = (datetime.now () - d1) / line
               time_to_go = (nlines-line) * ts_pr_line
               eta = datetime.now () + time_to_go
               text_file.write("ETA: {0}".format (str (eta)))
                   
        #if line > div*2:
        #   div = div*2
        
        N['total'] += 1
        if A.is_unmapped or B.is_unmapped:
            N['unmapped'] += 1
            continue
        if A.tid != B.tid:
            N['other_chrom'] += 1
            continue
        if A.is_reverse != B.is_reverse:
            N['other_strand'] += 1
            continue

      

        dist = B.pos - A.pos
        

        if numpy.abs(dist) < options.asize*2 and A.tid == B.tid and A.is_reverse == B.is_reverse:
            N['overlapping_anchors1'] += 1
            #if (A.is_reverse and dist >= -options.asize*2) or (not A.is_reverse and dist <= options.asize*2): #circ?
               
            read = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read)
                            
            bp = find_breakpoints_tiny(A,B,read,chrom)
            if not bp:
                N['tinyc_no_bp'] += 1
            else:
                N['tinyc_reads'] += 1

            n_hits = len(bp) 
            for h in bp:
                #print h
                # for some weird reason for circ we need a correction here
                dist,ov,chrom,start,end,sense = h
                h = (chrom,start+1,end-1,sense)
                tinyc[h].add(read,A,B,dist,ov,n_hits)
            
       
            
                                
        elif (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
            # get original read sequence
            read = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read)
                            
            bp = find_breakpoints(A,B,read,chrom)
            if not bp:
                N['circ_no_bp'] += 1
            elif numpy.abs(dist) < options.asize:
                N['circ_reads_overlap'] += 1
            else:
                N['circ_reads'] += 1

            n_hits = len(bp) 
            for h in bp:
                #print h
                # for some weird reason for circ we need a correction here
                dist,ov,chrom,start,end,sense = h
                h = (chrom,start+1,end-1,sense)
                circs[h].add(read,A,B,dist,ov,n_hits)

        elif (A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0):
            read = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            name = A.qname
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read)
                            
            bp = find_breakpoints(A,B,read,chrom)
            if not bp:
                N['splice_no_bp'] += 1
            else:
                N['spliced_reads'] += 1
            n_hits = len(bp)                
            for h in bp:
                #print h
                dist,ov,chrom,start,end,sense = h
                h = (chrom,start,end,sense)
                
                #No need to remember read
                #splices[h].add(read,A,B,dist,ov,n_hits)
                
                                
                # remember the spliced reads at these sites
                #loci[(chrom,start,sense)].append(splices[h])
                #loci[(chrom,end,sense)].append(splices[h])
                
                loci2[(chrom,start,sense)] = loci2[(chrom,start,sense)] + 1
                loci2[(chrom,end,sense)]   = loci2[(chrom,end,sense)] + 1
                
                               
                #tissue specific loci ### Removed again to reduce mem usage
                #for (prefix,tiss) in samples:
                #   if name.startswith(prefix):
                #       loci_tissue[(chrom,start,sense,tiss)].append(splices[h])
                #       loci_tissue[(chrom,end,sense,tiss)].append(splices[h])
                       
                
            
    #break
except KeyboardInterrupt:
    pass

def output(cand,prefix,output_reads):
    #                      1        2      3     4       5         6        7         8          9             10            11           12          13                 14                  15     16            17                    18            19            20
    print "#","\t".join(['chrom','start','end','name','n_reads','strand','n_uniq','n_uniq2','best_qual_A','best_qual_B','avr_qual_A','avr_qual_B','spliced_at_begin','spliced_at_end','tissues','tiss_counts','tiss_counts_spliced','edits','anchor_overlap','breakpoints'])
    n = 1
    for c,hit in cand.items():
        #print c
        chrom,start,end,sense = c
        n_reads,n_uniq,n_uniq2,best_qual,best_other,avr_qual_A,avr_qual_B,spliced_at_begin,spliced_at_end,tissues,tiss_counts,tiss_rel_exp,min_edit,min_anchor_ov,n_hits = hit.scores(chrom,start,end,sense)
        
        
        if best_other < options.min_uniq_qual:
            N['anchor_not_uniq'] += 1
            continue
        
        name = "%s%s_%06d" % (options.prefix,prefix,n)
        n += 1
        
        
        #sys.stderr.write("%s\t%s\n" % (name,"\t".join(sorted(reads))))
        
        if output_reads == True:
           for r in hit.reads:
               sys.stderr.write("%s\t%s\n" % (name,r))

               
        bed = [chrom,start-1,end,name,n_reads,sense,n_uniq,n_uniq2,best_qual,best_other,avr_qual_A,avr_qual_B,spliced_at_begin,spliced_at_end,",".join(tissues),",".join(tiss_counts),",".join(tiss_rel_exp),min_edit,min_anchor_ov,n_hits]
        print "\t".join([str(b) for b in bed])

stats = file(options.stats,"w")
stats.write(str(N)+"\n")


#print (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
#print (memory_usage())

output(circs,"circ", True)
#output(circs_ol,"c_overlap")
output(splices,"norm", False)
#output(trans,"trans", True)
output(tinyc,"tinyc", True)
#output(splices_ol,"n_overlap")
