/*
 * Program BWBBLE mg-aligner for Short-Read Alignment to a Multi-Genome
 * by Victoria Popic (viq@stanford.edu) 2012
 *
 * MIT License
 *
 * Copyright (c) 2012 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bwt.h"
#include "align.h"
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>

static int usage() {
	printf("Usage:   bwbble command [options] \n");
	printf("Command: index    index sequences in the FASTA format\n");
	printf("         align    exact or inexact read alignment\n");
	printf("         fasta2ref    constructs a single linear reference from the input file \n");
	printf("         aln2sam  convert alignment results to SAM file format for single-end mapping\n");
	printf("\n");
	return 1;
}

static int index_usage() {
	printf("Usage: bwbble index [options] <seq_fasta> \n");
	printf("Options: e    file with the SA precomputed by the external memory eSAIS algorithm.\n");
	printf("\n");
	return 1;
}

static int align_usage() {
	printf("Usage: bwbble align [options] <seq_fasta> <reads_fastq> <output_aln> \n");
	printf("Options: M    mismatch penalty (default: 3)\n");
	printf("         O    gap open penalty (default: 11) \n");
	printf("         E    gap extend penalty (default: 4) \n");
	printf("         n    maximum number of differences in the alignment (gaps and mismatches) (default: 0)\n");
	printf("         l    length of the seed (seed := first seed_length chars of the read) (default: 32)\n");
	printf("         k    maximum number of differences in the seed (default: 2)\n");
	printf("         o    maximum number of gap opens (default: 1)\n");
	printf("         e    maximum number of gap extends (default: 6) \n");
	printf("         t    run multi-threaded with t threads (default: 1)\n");
	printf("         S    align with a single-genome reference\n");
	printf("         P    use pre-calculated partial alignment results\n");
	printf("\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0) {
		if(argc < 3) {
			index_usage();
			exit(1);
		}

		char* use_esa_file = NULL;
		int c;
		while ((c = getopt(argc-1, argv+1, "e:")) >= 0) {
			switch (c) {
				case 'e': use_esa_file = optarg; break;
				case '?': index_usage(); return 1;
				default: return 1;
			}
		}
		index_bwt(argv[optind+1], use_esa_file);
	}
	else if (strcmp(argv[1], "align") == 0) {
		if(argc < 5) {
			align_usage();
			exit(1);
		}
		aln_params_t* params = (aln_params_t*) calloc(1, sizeof(aln_params_t));
		set_default_aln_params(params);

		int c;
		while ((c = getopt(argc-1, argv+1, "M:O:E:n:k:o:e:l:m:t:SP")) >= 0) {
				switch (c) {
					case 'M': params->mm_score = atoi(optarg); break;
					case 'O': params->gapo_score = atoi(optarg); break;
					case 'E': params->gape_score = atoi(optarg); break;
					case 'n': params->max_diff = atoi(optarg); break;
					case 'k': params->max_diff_seed = atoi(optarg); break;
					case 'o': params->max_gapo = atoi(optarg); break;
					case 'e': params->max_gape = atoi(optarg); break;
					case 'l': params->seed_length = atoi(optarg); break;
					case 'm': params->max_entries = atoi(optarg); break;
					case 't': params->n_threads = atoi(optarg); break;
					case 'S': params->is_multiref = 0; break;
					case 'P': params->use_precalc = 1; break;
					case '?': align_usage(); return 1;
					default: return 1;
				}
		}
		
		align_reads(argv[optind+1], argv[optind+2], argv[optind+3], params);
		free(params);
	} else if (strcmp(argv[1], "fasta2ref") == 0) {
                if(argc < 3) {
                        printf("Usage: bwbble fasta2ref <seq_fasta> \n");
                        exit(1);
                }
		char* refFname  = (char*) malloc(strlen(argv[optind+1]) + 5);
		char* annFname  = (char*) malloc(strlen(argv[optind+1]) + 5);
        	sprintf(refFname, "%s.ref", argv[optind+1]);
        	sprintf(annFname, "%s.ann", argv[optind+1]);
		unsigned char *seq;
        	bwtint_t seqLen;
                fasta2ref(argv[optind+1], refFname, annFname, &seq, &seqLen);
        	free(annFname);
        	free(refFname);
		free(seq);
	} else if (strcmp(argv[1], "aln2sam") == 0) {
		if(argc < 6) {
			printf("Usage: bwbble aln2sam [-S, -n] <seq_fasta> <reads_fastq> <alns_aln> <out_sam> \n");
			exit(1);
		}
		int is_multiref = 1;
		int max_diff = 6;
		//int n_occ = 3;
		int c;
		while ((c = getopt(argc-1, argv+1, "n:S:o")) >= 0) {
			switch (c) {
				case 'S': is_multiref = 0; break;
				case 'n': max_diff = atoi(optarg); break;
				//case 'o': n_occ = atoi(optarg); break;
				case '?': printf("Unknown option \n"); break;
				default: return 1;
			}
		}
		alns2sam(argv[optind+1], argv[optind+2], argv[optind+3], argv[optind+4], is_multiref, max_diff);
	} else {
		printf("Error: Unknown command '%s'\n", argv[1]);
		usage();
	}
	return 0;
}
