#! /bin/bash

echo -e "Sample\tReads_sequenced\tReads_mapped\tUnique_reads_mapped"

for fq in $*; 
 do 
    if [ ${fq: -10} == '.rmdup.bam' ]; then continue; fi;
    sample=$(basename $fq);
    nlines=$(zcat ${fq} | wc -l);
    nreads=$(echo "${nlines} / 4" | bc);
    nmap=$(samtools view -c ${fq%.*}.bam);
    if [ -z $nmap ]; then $nmap=0; fi;
    nnodup=$(samtools view -c ${fq%.*}.rmdup.bam);
    if [ -z $nnodup ]; then $nnodup=0; fi;
    echo -e "${sample%%.*}\t${nreads}\t${nmap}\t${nnodup}"; 
 done
