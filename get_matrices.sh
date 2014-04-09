#! /bin/bash
set -u
set -e

segmentation=$1

python -u ~/repos/gam/gam_matrix.py create -o -a -z lzf $segmentation
chroms=`tail -n +2 $segmentation | cut -f 1,1 | sort | uniq | grep -v _random | tr "\n" " "`
for chrom in $chroms
do
    echo $chrom
    python -u ~/repos/gam/gam_matrix.py chrom -p 3 -s ${segmentation}.matrix.hdf5 -c $chrom
done
