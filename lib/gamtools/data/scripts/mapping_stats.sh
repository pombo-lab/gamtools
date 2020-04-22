#! /bin/bash

echo -e "Sample\tReads_sequenced\tReads_mapped\tUnique_reads_mapped\tPercent_mapped_reads"

for fq in $*; 
 do 
    if [ ${fq: -10} == '.rmdup.bam' ]; then continue; fi;
    cat_cmd=cat;
    if [ ${fq##*.} '==' 'gz' ]; then cat_cmd=zcat; fi;
    sample=$(basename $fq);
    name_no_ext=$(echo -n $fq | tr "." "\n" | sed '/^gz$/ d' | sed '$ d' | perl -pe 'chomp if eof' | tr "\n" ".");
    nlines=$(${cat_cmd} ${fq} | wc -l);
    nreads=$(echo "${nlines} / 4" | bc);
    nmap=$(samtools view -c ${name_no_ext}.bam);
    if [ -z $nmap ]; then $nmap=0; fi;
    nnodup=$(samtools view -c ${name_no_ext}.rmdup.bam);
    if [ -z $nnodup ]; then $nnodup=0; fi;
    pct_mapped=$(echo "scale=2; 100*$nmap/$nreads" | bc);
    echo -e "${sample%%.*}\t${nreads}\t${nmap}\t${nnodup}\t${pct_mapped}"; 
 done
