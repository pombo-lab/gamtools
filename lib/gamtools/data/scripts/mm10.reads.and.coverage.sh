#!/bin/bash
#$ -N mm10.reads&coverage
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G

export PATH=/data/pombo/software/bedtools_recent:$PATH

FILES=/path/to/bams/*.rmdup.bam
windowsize=50000
name=GAM_pool_bla_bla
genome_file=/data/pombo/Sasha/Common_files/Genomic_windows/mm10.chrom.sizes

#=========================================================================================

bedtools makewindows -g $genome_file -w $windowsize > ${name}.at.${windowsize}.windows.bed
cp ${name}.at.${windowsize}.windows.bed ${name}.at.${windowsize}.reads.windows.bed
cp ${name}.at.${windowsize}.windows.bed ${name}.at.${windowsize}.coverage.windows.bed

for f in $FILES
do
	bedtools coverage -a ${name}.at.${windowsize}.windows.bed -b $f | tee >(awk '{print $4}' > ${name}.at.${windowsize}.reads.bed) >(awk '{print $5}' > ${name}.at.${windowsize}.coverage.bed)
	paste ${name}.at.${windowsize}.reads.windows.bed ${name}.at.${windowsize}.reads.bed > ${name}.at.${windowsize}.merged.bed && mv ${name}.at.${windowsize}.merged.bed ${name}.at.${windowsize}.reads.windows.bed
	paste ${name}.at.${windowsize}.coverage.windows.bed ${name}.at.${windowsize}.coverage.bed > ${name}.at.${windowsize}.merged.bed && mv ${name}.at.${windowsize}.merged.bed ${name}.at.${windowsize}.coverage.windows.bed
done

echo -e "chrom\tstart\tstop" $FILES > ${name}.at.${windowsize}.list_of_files.txt
tr ' ' \\t < ${name}.at.${windowsize}.list_of_files.txt > ${name}.at.${windowsize}.list_of_files.tab.txt && mv ${name}.at.${windowsize}.list_of_files.tab.txt ${name}.at.${windowsize}.list_of_files.txt
cat ${name}.at.${windowsize}.list_of_files.txt ${name}.at.${windowsize}.reads.windows.bed > ${name}.at.${windowsize}.reads.windows.header.bed && mv ${name}.at.${windowsize}.reads.windows.header.bed ${name}.reads.at.${windowsize}.table
cat ${name}.at.${windowsize}.list_of_files.txt ${name}.at.${windowsize}.coverage.windows.bed > ${name}.at.${windowsize}.coverage.windows.header.bed && mv ${name}.at.${windowsize}.coverage.windows.header.bed ${name}.coverage.at.${windowsize}.table
rm ${name}.at.${windowsize}.list_of_files.txt ${name}.at.${windowsize}.reads.bed ${name}.at.${windowsize}.coverage.bed ${name}.at.${windowsize}.windows.bed ${name}.at.${windowsize}.reads.windows.bed ${name}.at.${windowsize}.coverage.windows.bed 


