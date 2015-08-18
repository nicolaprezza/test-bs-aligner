#generate dataset with default parameters
#usage: testcase-bs-aligner-default.sh <input_genome> <output_folder>
#behaviour: generates in output_folder a set of reads with SNPs, BS conversions, indels, and sequencing errors. BS conversions are consistent with a simulated methylome which is output as a .bed file. Moreover, there is a certain % of bisulfite error (i.e. some reads are not BS converted) 

#requires: SimSeq (https://github.com/jstjohn/SimSeq) installed in ~/workspace/SimSeq/
#requires: samtools (http://www.htslib.org/)
#requires: fastx-mutate-tools (https://github.com/nicolaprezza/fastx-mutate-tools) installed in ~/workspace/fastx-mutate-tools/
#requires: realpath
#requires: fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
#requires: gawk, awk

input_genome=$1
output_folder=`realpath $2`
simseq_folder=~/workspace/SimSeq
fastx_mutate_tools_folder=~/workspace/fastx-mutate-tools
snp_frequency=0.005 #on average, 1 snp every 200 bases
indel_open_freq=0.0003 #equivalent to inserting 1M indels in the human genome
indel_ext_freq=0.8 #average indel length is 5
number_of_pairs=1500
#number_of_pairs=6000000
read_length=100
bs_failure_rate=2 #percentage of read pairs that fail to be bs-converted
bs_conversion_rate=0.5 #all Cs on fw/rev strands are converted to Ts with this probability
create_bed=true #create .bed files with true methylation annotations?

./testcase-bs-aligner.sh $input_genome $output_folder $simseq_folder $fastx_mutate_tools_folder $snp_frequency $indel_open_freq $indel_ext_freq $number_of_pairs $read_length $bs_failure_rate $bs_conversion_rate $create_bed
