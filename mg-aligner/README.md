### Compilation

Run `make all` in this folder to compile the program.

### Command Line Arguments

The following commands and options are available with the BWBBLE aligner.

Usage: `bwbble command [options]`

#### Commands:

1.`index` index the reference multi-genome in the FASTA format  
`bwbble index <seq_fasta>`

2.`align` align reads  
`bwbble align [options] <seq_fasta> <reads_fastq> [output_aln]`

##### Options:

`-M <arg>` mismatch penalty (default: 3)  
`-O <arg>` gap open penalty (default: 11)  
`-E <arg>` gap extend penalty (default: 4)  
`-n <arg>` maximum number of differences in the alignment (gaps and mismatches) (default: 0)  
`-l <arg>` length of the seed (seed := first `l` chars of the read) (default: 32)  
`-k <arg>` maximum number of differences in the seed (default: 2)  
`-o <arg>` maximum number of gap opens (default: 1)  
`-e <arg>` maximum number of gap extends (default: 6)  
`-t <arg>` run multi-threaded with t threads (default: 1)  
`-S` align with a single-genome reference

3.`aln2sam` convert alignment results to SAM file format for single-end mapping  
`bwbble aln2sam [options] <seq_fasta> <reads_fastq> <alns_aln> <output_sam>`

### Example

`make clean; make`  
`./bwbble index ../test_data/chr21.fa`  
`./bwbble align ../test_data/chr21.fa ../test_data/sim_chr21_N100.fastq ../test_data/sim_chr21_N100.aln`  
`./bwbble aln2sam ../test_data/chr21.fa ../test_data/sim_chr21_N100.fastq ../test_data/sim_chr21_N100.aln ../test_data/sim_chr21_N100.sam`
