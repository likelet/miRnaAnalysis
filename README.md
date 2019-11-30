# miRnaAnalysis
A collection of miRna Analysis pipes, current function could be used for predicting the miRNA based on the sequences from input dataset.

# Installation
Before use this pipe users should Install the following softare and the configure the `config` file in this folder when use this pipeline at the first time 

    Nextflow
    miranda 
    targetScan
    java 
    R 
    vendiagram package 



# Predicting miRNA target by several tools 

    nextflow run miRnaAnalysis/predictMiRNA.nf -profile c2 --fasta candidates_seq.fa 

`--fasta` the mRNA sequences in fasta  format for predicting.

