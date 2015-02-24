#simulates BS conversion in the genome + adds SNPs and indels and generates :
#1) true methylation annotations
#2) 2 fastq files with the simulated reads (paired-end)

#requires: SimSeq (https://github.com/jstjohn/SimSeq) installed in ~/workspace/SimSeq/
#requires: samtools (http://www.htslib.org/)
#requires: fastx-mutate-tools (https://github.com/nicolaprezza/fastx-mutate-tools) installed in ~/workspace/fastx-mutate-tools/
#requires: realpath
#requires: fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
#requires: gawk, awk
#requires: erne (http://erne.sourceforge.net/)

genome=$1
output_folder=`realpath $2`
simseq_folder=`realpath $3`
fastx_mutate_tools_folder=`realpath $4`
snp_frequency=$5
indel_open_freq=$6
indel_ext_freq=$7
number_of_pairs=$8
read_length=$9
bs_failure_rate=${10}
bs_conversion_rate=${11}

simseq_jar=$simseq_folder/SimSeqNBProject/store/SimSeq.jar
fastx_mutate_tools_jar=$fastx_mutate_tools_folder/FastxMutateTools.jar

#error_profile_1=$simseq_folder/examples/illumina_hg19_erne_map.txt
#error_profile_2=$error_profile_1

error_profile_1=$simseq_folder/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt
error_profile_2=$simseq_folder/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt

shopt -s expand_aliases
alias simseq='java -jar -Xmx4000m $simseq_jar'
alias fastx-mutate-tools='java -jar $fastx_mutate_tools_jar'

#1) insert SNPs in the original genome

genome_with_snps=$output_folder/genome_with_snps.fa

fastx-mutate-tools snp --input $genome --output $genome_with_snps --snp $snp_frequency

sam_header=${genome_with_snps}.fai

samtools faidx $genome_with_snps

#2) simulate the reads that failed to be BS-converted

number_bs_failed_pairs=$(((bs_failure_rate*number_of_pairs)/100))
failed_bs_sam=$output_folder/failed_bs.sam

reads1=$output_folder/reads_1.fq
reads2=$output_folder/reads_2.fq

simseq -1 $read_length -2 $read_length --error $error_profile_1 --error2 $error_profile_2 -n $number_bs_failed_pairs -o $failed_bs_sam -p failed_bs -r $genome_with_snps 

#sam to fastq

./sam-paired-to-fastq.sh $failed_bs_sam $output_folder/reads

rm $failed_bs_sam

#  3) generate simulated BS sam files

#3.1) simulate methylation on fw strand

genome_with_snps_bs_fw=$output_folder/genome_with_snps_bs_fw.fa

annotations_fwd=$output_folder/meth_annotations_fwd.bed
annotations_rev=$output_folder/meth_annotations_rev.bed

fastx-mutate-tools snp --input $genome_with_snps --output $genome_with_snps_bs_fw --bs-fwd $bs_conversion_rate --output-methylation ${annotations_fwd}_temp --1-based

#3.2) simulate sam fw strand

number_bs_pairs=$((((100-bs_failure_rate)*number_of_pairs)/100))
fwd_bs_sam=$output_folder/fwd_bs.sam

simseq -1 $read_length -2 $read_length --error $error_profile_1 --error2 $error_profile_2 -n $number_bs_pairs -o $fwd_bs_sam -p fwd_bs -r $genome_with_snps_bs_fw

rm $genome_with_snps_bs_fw

#3.3) simulate methylation on rev strand

genome_with_snps_bs_rev=$output_folder/genome_with_snps_bs_rev.fa

fastx-mutate-tools snp --input $genome_with_snps --output $genome_with_snps_bs_rev --bs-rev $bs_conversion_rate --output-methylation ${annotations_rev}_temp --1-based

rm $genome_with_snps

cat ${annotations_fwd}_temp | grep '+' > $annotations_fwd
rm ${annotations_fwd}_temp
cat ${annotations_rev}_temp | grep '-' > $annotations_rev
rm ${annotations_rev}_temp

#3.4) simulate sam rev strand

rev_bs_sam=$output_folder/rev_bs.sam

simseq -1 $read_length -2 $read_length --error $error_profile_1 --error2 $error_profile_2 -n $number_bs_pairs -o $rev_bs_sam -p rev_bs -r $genome_with_snps_bs_rev

rm $genome_with_snps_bs_rev

#4) split read 1 and read 2 in sam files

fwd_bs_sam_read1_temp=$output_folder/fwd_bs_1_temp.sam
fwd_bs_sam_read2_temp=$output_folder/fwd_bs_2_temp.sam
rev_bs_sam_read1_temp=$output_folder/rev_bs_1_temp.sam
rev_bs_sam_read2_temp=$output_folder/rev_bs_2_temp.sam

bs_sam_read1=$output_folder/bs_1.sam
bs_sam_read2=$output_folder/bs_2.sam

cat $fwd_bs_sam | awk 'NR % 2 == 1' > $fwd_bs_sam_read1_temp
cat $fwd_bs_sam | awk 'NR % 2 == 0' > $fwd_bs_sam_read2_temp
cat $rev_bs_sam | awk 'NR % 2 == 1' > $rev_bs_sam_read1_temp
cat $rev_bs_sam | awk 'NR % 2 == 0' > $rev_bs_sam_read2_temp

rm $fwd_bs_sam $rev_bs_sam

#5) extract from fwd only pairs that have the first read on fwd and 2nd read on reverse AND extract from rev only pairs that have the first read on rev and 2nd read on forward

samtools view -t $sam_header -S -F 0x10 $fwd_bs_sam_read1_temp > $bs_sam_read1
rm $fwd_bs_sam_read1_temp
samtools view -t $sam_header -S -f 0x10 $fwd_bs_sam_read2_temp > $bs_sam_read2
rm $fwd_bs_sam_read2_temp

#sam to fastq
cat $bs_sam_read1 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> $reads1
rm $bs_sam_read1
cat $bs_sam_read2 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${reads2}_temp
rm $bs_sam_read2

#reverse complement query2 (because in SAM reads on reverse are reversed complemented)
cat ${reads2}_temp | fastx_reverse_complement -Q33 >> ${reads2}
rm ${reads2}_temp

samtools view -t $sam_header -S -F 0x10 $rev_bs_sam_read1_temp > $bs_sam_read2 #forward reads that are the first of their pair. They become read2
rm $rev_bs_sam_read1_temp
samtools view -t $sam_header -S -f 0x10 $rev_bs_sam_read2_temp > $bs_sam_read1 #reverse reads that are the second of their pair. They become read1
rm $rev_bs_sam_read2_temp

rm $sam_header

#sam to fastq
cat $bs_sam_read1 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${reads1}_temp
rm $bs_sam_read1
cat $bs_sam_read2 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> ${reads2}
rm $bs_sam_read2

#reverse complement query2 (because in SAM reads on reverse are reversed complemented)
cat ${reads1}_temp | fastx_reverse_complement -Q33 >> ${reads1}
rm ${reads1}_temp

#6) insert indels and create final reads

reads_indel1=$output_folder/query1_ind.fq
reads_indel2=$output_folder/query2_ind.fq

fastx-mutate-tools indel --input ${reads1} --output $reads_indel1 --open $indel_open_freq --extend $indel_ext_freq
rm ${reads1}
fastx-mutate-tools indel --input ${reads2} --output $reads_indel2 --open $indel_open_freq --extend $indel_ext_freq
rm ${reads2}

#7) filter

output_prefix=$output_folder/query
erne-filter --query1 $reads_indel1 --query2 $reads_indel2 --output-prefix $output_prefix

rm $reads_indel1 $reads_indel2


