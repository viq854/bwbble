###Compilation

Run ```make all``` in this folder to compile the program.  

###Command Line Arguments  

The following commands and options are available with the BWBBLE aligner.  

Usage: ```bwbble command [options]```   

Command:  
  ```index``` index the reference multi-genome in the FASTA format   
  ```bwbble index <seq_fasta>```   

  ```align```   align reads  
  ```bwbble align [options] <seq_fasta> <reads_fastq> [output_aln]```  
 
Options:  
   M mismatch penalty 
   O  gap open penalty 
   E  gap extend penalty 
   n  maximum number of differences in the alignment (gaps and mismatches) 
   l  length of the seed (seed := first l chars of the read) 
   k  maximum number of differences in the seed 
   o  maximum number of gap opens 
   e  maximum number of gap extends 
   t  run multi-threaded with t threads 
   S align with a single-genome reference 

  aln2sam convert alignment results to SAM file format for single-end mapping
  bwbble aln2sam [options] seq_fasta reads_fastq alns_aln 
