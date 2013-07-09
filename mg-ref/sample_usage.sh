#!/bin/bash

#######################################
## compile mg-ref
#######################################
make clean; make


#######################################
## clean up the output directory
#######################################
if [ -d "mg-ref-output" ]; then
  echo "directory mg-ref-output found"
  rm mg-ref-output/SNP.extract.chr*.data
  rm mg-ref-output/INDEL.extract.chr*.data
else
  echo "create directory mg-ref-output"
  mkdir mg-ref-output
fi


#######################################
## download and extract variants from vcf files
#######################################
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
do
  echo "start download $chr"
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
  echo "start unzip $chr"
  gunzip ALL.$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
  echo "start extract $chr"
  ./data_prep ALL.$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf
  echo "remove $chr"
  rm ALL.$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf
done


#######################################
## download GRCh37 if not exist
#######################################
if [ -f "human_g1k_v37.fasta" ]
then
  echo "file human_g1k_v37.fasta found"
else
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
  gunzip human_g1k_v37.fasta.gz
fi


#######################################
## create reference multi-genome:
## reference_w_snp_and_bubble.fasta in this sample
#######################################
echo "start comb"
output_dir=mg-ref-output
./comb human_g1k_v37.fasta $output_dir/reference_w_snp.fasta $output_dir/reference_w_snp_and_bubble.fasta $output_dir/meta.data
