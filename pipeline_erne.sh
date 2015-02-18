#pipeline with erne (http://erne.sourceforge.net/)

#clear working diretory and simulate dataset
rm wd/*; ./testcase-bs-aligner-default.sh ~/workspace/datasets/genomes/50000/50000.fasta wd/

#align with the BS-aligner erne-bs5
erne-bs5 --reference ~/workspace/datasets/genomes/50000/50000.ebm --query1 wd/query1.fq --query2 wd/query2.fq --output wd/erne_bs5_out.bam

#call methylation with the caller erne-meth, generating methylation annotations in bismark format 
erne-meth --fasta ~/workspace/datasets/genomes/50000/50000.fasta --input wd/erne_bs5_out.bam --output-prefix wd/erne_meth_out --annotations-bismark

#count number of calls that are 20% from the correct value
./compare_meth_annotations.sh wd/erne_meth_out_bismark.bed wd/ 0.2