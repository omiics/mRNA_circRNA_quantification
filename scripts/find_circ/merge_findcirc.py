#!/usr/bin/env python
import sys
import glob
from collections import defaultdict
dReads = defaultdict (list)
dPos = defaultdict (int)

files = []



for x in range (1, len(sys.argv)):
   files.append (sys.argv[x])


files.sort()
   
for file in files:
   #print file
   f = open(file, 'r')
   
   for line in f:
     
      cols = line.split ('\t')
      
      if (cols[0][0] == "#"):
         continue

      
      junction = cols[4]
      
            
      #non_junction
      SA = cols[12]
      SD = cols[13]
      name = cols[3]
      
      #if name.find ("circ") == -1:
      #   continue
      
      chrom = cols[0]
      start = cols[1]
      end = cols[2]
      strand = cols[5]
      
      #if int(end)-int(start) > 100000:
      #   continue
            
      dReads[(chrom, start, end, strand, file)] = (junction, SA, SD, name)
      dPos[(chrom, start, end, strand)] = 1
      
output = "\t".join (["#chrom", "start", "end", "names", "total_junction", "strand", "total_SA", "total_SD"])
      
for f in files:
   #output = "\t".join ([output, "\t".join ([f + "_junction",f+"_SA",f+"_SD"])])
   output = "\t".join ([output, "\t".join ([f + "_junction"])])


print output
      

      

#print dReads
for (chrom, start, end, strand) in dPos:
   if dPos[(chrom, start, end, strand)] == 1:
     
      output = "\t".join ([chrom, start, end])
     
      output_data = []
      
      name_cat = ""
      total_junction = 0
      
      total_SA = 0
      total_SD = 0
          
      for f in files:
         #print f
         if dReads[(chrom, start, end, strand, f)]:
            (junction, SA, SD, name) = dReads[(chrom, start, end, strand, f)]         
            #output_data.append ("\t".join ([junction, SA, SD]))
            output_data.append ("\t".join ([junction]))

            total_junction = total_junction + int (junction)
            total_SA = total_SA + int (SA)
            total_SD = total_SD + int (SD)
            name_cat = ",".join ([name_cat, name])
         else:
            #output_data.append ("\t".join (['0', '0', '0']))
            output_data.append ("\t".join (['0']))

           
      
      if total_junction < 2:
         continue
      
      print output + "\t" + name_cat + "\t" + "\t".join ([str(total_junction), strand, str(total_SA), str(total_SD)]) + "\t" + "\t".join (output_data)
       
    
    
