bwbble mg-ref
======

Questions or comments may be directed to Lin Huang at linhuang@cs.stanford.edu

## compile

In this directory compile mg-ref by running

make clean; make

## data_prep

data_prep [option] input1.vcf input2.vcf ...

This command extracts SNPs and INDELs from a list of VCF files, write the SNPs into mg-ref-output/SNP.extract.chrxx.data, and write the INDELs into mg-ref-output/INDEL.extract.chrxx.data, where xx is a chromosome identified from the reference genome. A VCF file can include variants from multiple chromosomes, and variants from the same chromosome can be recorded in multiple VCF files.

Option -c clears all mg-ref-output/SNP.extract.chrxx.data and mg-ref-output/INDEL.extract.chrxx.data files before usage.

## comb

comb [options] ref.fasta ref_w_snp.fasta ref_w_snp_and_bubble.fasta bubble.data

This command reads SNPs from mg-ref-output/SNP.extract.chrxx.data and INDELs from mg-ref-output/INDEL.extract.chrxx.data, and combine them with the reference genome ref.fasta. The output ref_w_snp.fasta combines the reference genome with the SNPs, ref_w_snp_and_bubble.fasta combines the reference genome with the SNPs and the INDELs.

An INDEL is padded at both ends with w base pairs. This parameter can be specified with option -w.

Options -i and -a specify minimum and maximum SNP occurrence in the VCF files. If the number of the occurrence of a SNP in the VCF files is lower than the specified minimum value or higher than the specified maximum value, this SNP will not be incorporated into ref_w_snp.fasta and ref_w_snp_and_bubble.fasta.


## sam_pad

sam_pad bubble.data sam.input sam.output

This command parses the output sam file of mg-aligner. For the reads aligned to bubble branches, this command computes their positions in the reference genome and appends this information to the end of the correspondings lines in the sam file. The chromsome identifier is specified with TAG bC, the 1-based leftmost mapping position of the first matching base in the reference genome is specified with the TAG bP. If the leftmost mapping position of a read locates within the INDEL of a bubble branch, the bP field specifies a position range at the reference genome. Consider the following example (a). If a read is aligned to position 3 at the bubble branch, its corresponding position at the reference genome is 3-7. Consider example (b). If a read is aligned to position 6-8 at the bubble branch, its corresponding position at the reference genome is 12.

(a)

                12345 67890 12345 6
        ref     ACGTA ATATA TAGGC T
        bubble  ACG-- --ATA TA
                123     456 78

(b)

                12345 67890 12--345 6
        ref     ACGTA ATATA TA--GGC T
        bubble         TATA TACTGGC T
                       1234 5678901 2
