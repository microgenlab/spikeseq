#!/bin/bash
#this script should be run using source

echo "Activating environment using fastq-tools through conda"

conda activate fastq-tools

echo "Sampling with replacement"
echo "Read sampling for each replicate will start with 500 reads up to 10000"

for f in {1..5}
do
echo "$f"
for d in *.fastq
do
echo $d
mkdir dir_$d_$f;
i=0
end=9500
while [ $i -le $end ]; do
    echo $i
    i=$(($i+500))
    fastq-sample -n $i $d -s $f -o samp_$i;
    mv samp_* dir_$d_$f;
done
done
done

