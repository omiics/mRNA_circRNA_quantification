#!/bin/bash


a1=cx115_Index12_1
a2=cx115_Index12_2
data_ref=~/primary_data/ssc_cortex_RNA-seq/cortex_scriptSeq/

genome=~/pig_genome/Sscrofa10.2/genome-fastas/
bowtie2=Sscrofa10.2


for sample in $a1 $a2
	do
	echo "starting for sample "$sample
	bowtie2 -p16 --very-sensitive --phred64 --mm -M20 --score-min=C,-15,0 -x $bowtie2 -q -U $data_ref$a1.fq.gz 2> bowtie2.$sample.log | samtools view -hbuS - | samtools sort - $sample
	samtools view -hf 4 $sample.bam | samtools view -Sb - > unmapped_$sample.bam
	./unmapped2anchors.py unmapped_$sample.bam | gzip > $sample_anchors.fq.gz

	mkdir circ.$sample
	bowtie2 -p 16 --reorder --mm -M20 --score-min=C,-15,0 -q -x $bowtie2 -U $sample_anchors.fq.gz | ./find_circ.py -G $genome -p $sample.junction_ -s circ.$sample/sites.log > circ.$sample/sites.bed 2> circ.$sample/sites.reads
	done

echo Done

