BWBBLE: Short Read Alignment with Populations of Genomes
===

###Description

BWBBLE is a BWT-based aligner that maps short reads to a compressed linear representation of a collection of genomes (the reference 'multi-genome').  

The reference multi-genome is created from a single genome build (e.g. GRCh37) augmented with a set of genomic variants obtained from a given collection of genomes (e.g. the 1000 Genomes Project). The SNPs present in the population are incorporated by extending the 4-letter A/C/T/G nucleotide code of the single reference to the 16-letter IUPAC code, which allows us to encode which nucleotides have been observed at a given position across all the available sequenced genomes (for example, the IUPAC character R can be used to represent that both A and G were observed). Indels are incorporated into the multi-genome by appending them at the end of the reference padded at both ends with the bases surrounding the indel at its position in the reference.  

The reference multi-genome is stored using the standard FASTA file format. We provide a script to generate this FASTA file from a single reference FASTA file (e.g. GRCh37) and a set of VCF files containing the SNPs and indels to be incorporated into the reference multi-genome.  

Given the reference multi-genome, the BWBBLE aligner operates in 3 steps:  
1. reference indexing (this creates the BWT of the reference multi-genome and auxiliary data structures to enable read mapping) 
2. read mapping  
3. alignment results evaluation and reporting  

The reads are accepted in the FASTQ file format and the alignment results are reported using the SAM file format.  

The BWBBLE read alignment alrogithm is based on the BWT backwards search procedure, with the difference that a read nucleotide can now match multiple bases in the reference multi-genome (namely, all the IUPAC codes that include this nucleotide; e.g. read base A can match the IUPAC characters A, D, H, M, N, R, V and W). In order to tolerate sequencing errors and other variations of the reads from the reference multi-genome, similar to other aligners, BWBBLE allows up to a certain number of mismatches and gaps. Our inexact read alignment to the multi-genome reference algorithm is an extension of the inexact search algorithm used by BWA (Li and Durbin, 2009).  

In order to reduce the search space and improve the performance of the alignment algorithm we have also adapted many heuristics from the existing BWA aligner. In particular, no sub-optimal alignments are explored if the number of best hits exceeds a given threshold. Furthermore, if the number of differences in the best hit, ```n_best```, is less than ```n```, then ```n``` is reset to ```n_best + 1```. It is also possible to limit how many differences are allowed in the seed sequence (the first ```k``` base pairs at the beginning of the read, where ```k``` is the seed length and can be adjusted) and disallow insertions and deletions at the ends of the read. Other search parameters available to the users of BWA (e.g. mismatch, gap open and gap extension penalties) have also been incorporated to provide a similar user interface and improve the efficiency of the BWBBLE implementation.

###Program Structure

The program consists of the following two subdirectories:  
```mg-aligner```: BWBBLE aligner functionality (reference indexing, read-mapping, and alignment results evaluation)  
```mg-ref```: scripts to construct the reference multi-genome FASTA file  

###Command Line Arguments

Please check the subdirectories ```mg-aligner``` and ```mg-ref``` for details on how to run each program component.

###License

MIT License  

###System Requirements

OpenMP

###Authors

BWBBLE aligner functionality (```bwbble/mg-aligner```): Victoria Popic (viq@stanford.edu)  
Reference multi-genome FASTA file construction (```bwbble/mg-ref```): Lin Huang (linhuang@cs.stanford.edu)  

###Publication
Lin Huang\*, Victoria Popic\*, and Serafim Batzoglou. "Short read alignment with populations of genomes." Bioinformatics 29.13 (2013): i361-i370. (\*joint first authors)  

###Support

Please contact Victoria Popic (viq@stanford.edu) for support related to the alignment functionality and Lin Huang (linhuang@cs.stanford.edu) for help with the multi-genome construction.
