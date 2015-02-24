#compares meth annotations in bismark cov format with those produced by the script testcase-bs-aligner.sh
#usage: compare_meth_annotations_cov_format.sh <annotations in bismark cov format> <input_output_folder> <max_diff> <min coverage>
#where input_output_folder is the folder containing the files meth_annotations_fwd.bed and meth_annotations_rev.bed. max_diff is a value in [0,1]. A position is considered only if covered with >= <min coverage> bases. A methylation call is considered as correct if the absolute difference between the call and the correct beta value does not exceed max_diff (included: |beta-correct_beta|<=max_diff)

annot_bismark=$1
io_folder=`realpath $2`
max_diff=$3
min_covg=$4

#merge true meth annotations in one single bed file

meth_annotations_fwd_bed=$io_folder/meth_annotations_fwd.bed 
meth_annotations_rev_bed=$io_folder/meth_annotations_rev.bed

total=`bedtools intersect -a $meth_annotations_fwd_bed -b $1 -wa -wb | awk -v min_covg=$min_covg '($10+$11>=min_covg && $2==$7){print}'|wc -l`
total=$((total+`bedtools intersect -a $meth_annotations_rev_bed -b $1 -wa -wb | awk -v min_covg=$min_covg '($10+$11>=min_covg && $2==$7){print}'|wc -l`))

correct=`bedtools intersect -a $meth_annotations_fwd_bed -b $1 -wa -wb | awk -v min_covg=$min_covg -v maxdiff=$max_diff '(($10+$11>=min_covg) && ($2==$7) && (sqrt( ($5-($9/100))*($5-($9/100))) )<=maxdiff ){print}'|wc -l` 

correct=$((correct+`bedtools intersect -a $meth_annotations_rev_bed -b $1 -wa -wb | awk -v min_covg=$min_covg -v maxdiff=$max_diff '(($10+$11>=min_covg) && ($2==$7) && (sqrt( ($5-($9/100))*($5-($9/100))) )<=maxdiff ){print}'|wc -l`))

printf "\nnumber of correct calls/total number of calls: "$correct"/"$total"\n"
