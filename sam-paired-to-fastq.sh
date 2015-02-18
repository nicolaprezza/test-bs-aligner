#converts a sam file containing paired-end alignments in 2 fastq files

#ASSUMPTION: all pairs have the reads on different strands: if read1 is on forward (resp. reverse), then read2 must be on reverse (resp. forward)!
#NOTE: the script does not preserve the order of sam entries

#usage: sam-paired-to-fastq.sh <input.sam> <prefix_name>
#behaviour: creates files prefix_name_1.fq and prefix_name_2.fq

#requires: gawk, awk
#requires: fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)

samfile=$1
prefix=$2

#separate read 1 and read 2 in 2 different files (read 1 and read 2). We remove header lines with grep.
cat $samfile | grep -v ^@ | awk 'NR % 2 == 1' > ${prefix}_1_temp
cat $samfile | grep -v ^@ | awk 'NR % 2 == 0' > ${prefix}_2_temp

#separate read 1: first all read 1 on fwd strand, then all read 1 on rev strand

cat ${prefix}_1_temp | gawk '{if (!and($2,0x10)) {print}}' > ${prefix}_read1_fwd_sam
cat ${prefix}_1_temp | gawk '{if (and($2,0x10)) {print}}' > ${prefix}_read1_rev_sam

rm ${prefix}_1_temp

#separate read 2: first all read 2 on rev strand, then all read 2 on fwd strand

cat ${prefix}_2_temp | gawk '{if (and($2,0x10)) {print}}' > ${prefix}_read2_rev_sam
cat ${prefix}_2_temp | gawk '{if (!and($2,0x10)) {print}}' > ${prefix}_read2_fwd_sam

rm ${prefix}_2_temp

#read 1 fwd can be converted as is
cat ${prefix}_read1_fwd_sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${prefix}_1.fq
rm ${prefix}_read1_fwd_sam

#read 1 rev must be reversed complemented
cat ${prefix}_read1_rev_sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${prefix}_1_rev.fq
rm ${prefix}_read1_rev_sam

cat ${prefix}_1_rev.fq | fastx_reverse_complement -Q33 >> ${prefix}_1.fq
rm ${prefix}_1_rev.fq

#read 2 rev must be reversed complemented
cat ${prefix}_read2_rev_sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${prefix}_2_rev.fq
rm ${prefix}_read2_rev_sam

cat ${prefix}_2_rev.fq | fastx_reverse_complement -Q33 > ${prefix}_2.fq
rm ${prefix}_2_rev.fq

#read 2 fwd can be converted as is
cat ${prefix}_read2_fwd_sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> ${prefix}_2.fq
rm ${prefix}_read2_fwd_sam
