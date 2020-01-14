//-------------------------------------------------------
// Sequence IO and Compression Functionality
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#ifndef BWTSNP_IO_H
#define BWTSNP_IO_H

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "io.h"

void fasta_error(const char *fastaFname)
{
	printf("Error: File %s does not comply with the FASTA file format \n", fastaFname);
	exit(1);
}

void fastq_error(const char *fastqFname)
{
	printf("Error: File %s does not comply with the FASTQ file format \n", fastqFname);
	exit(1);
}

/* Reference I/O */

// reads the sequence data from the FASTA file, encodes and packs it into a .pac file,
// sequence annotations are stored into an .ann file
void fasta2pac(const char *fastaFname, const char *pacFname, const char *annFname)
{
	FILE *fastaFile = (FILE *)fopen(fastaFname, "r");
	if (fastaFile == NULL)
	{
		printf("fasta2pac: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}
	FILE *pacFile = (FILE *)fopen(pacFname, "wb");
	if (pacFile == NULL)
	{
		printf("fasta2pac: Cannot open PAC file: %s!\n", pacFname);
		exit(1);
	}
	FILE *annFile = (FILE *)fopen(annFname, "wb");
	if (annFile == NULL)
	{
		printf("fasta2pac: Cannot open ANN file: %s!\n", annFname);
		exit(1);
	}

	unsigned char seqBuffer[SEQ_BATCH_ALLOC_LEN];
	unsigned char packedSeqBuffer[SEQ_BATCH_ALLOC_LEN / CHARS_PER_BYTE];
	bwtint_t allocatedSeqLen = SEQ_BATCH_ALLOC_LEN;
	char *seq = (char *)malloc(allocatedSeqLen * sizeof(char));
	bwtint_t totalSeqLen = 0;
	bwtint_t seqBufferCount = 0;

	fasta_annotations_t *annotations = (fasta_annotations_t *)calloc(1, sizeof(fasta_annotations_t));
	bwtint_t allocatedAnnNum = ANN_ALLOC_LEN;
	annotations->seq_anns = (seq_annotation_t *)malloc(allocatedAnnNum * sizeof(seq_annotation_t));

	char c = (char)getc(fastaFile);
	if (c != '>')
		fasta_error(fastaFname);

	while (!feof(fastaFile))
	{
		if (allocatedAnnNum == annotations->num_seq)
		{
			allocatedAnnNum <<= 1;
			annotations->seq_anns = (seq_annotation_t *)realloc(annotations->seq_anns, allocatedAnnNum * sizeof(seq_annotation_t));
		}
		seq_annotation_t *seqAnnotation = &(annotations->seq_anns[annotations->num_seq]);

		c = (char)getc(fastaFile);

		// sequence description line (> ...)
		bwtint_t seqNameLen = 0;
		while (c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(fastaFile))
		{
			seqAnnotation->name[seqNameLen] = c;
			seqNameLen++;
			c = (char)getc(fastaFile);
		}
		seqAnnotation->name[seqNameLen] = '\0';
		while (c != '\n' && !feof(fastaFile))
		{
			c = (char)getc(fastaFile);
		}
		if (feof(fastaFile))
			fasta_error(fastaFname);

		// sequence data
		bwtint_t seqLen = 0;
		while (c != '>' && !feof(fastaFile))
		{
			if (c != '\n')
			{
				if (c >= 'a' && c <= 'z')
				{
					c += 'A' - 'a';
				}
				// reallocate twice as much memory
				if (seqLen >= allocatedSeqLen)
				{
					allocatedSeqLen <<= 1;
					seq = (char *)realloc(seq, sizeof(char) * allocatedSeqLen);
				}

				*(seq + seqLen) = c;
				seqLen++;
			}
			c = (char)getc(fastaFile);
		}
		// add $ as a separator between chromosomes and bubbles
		// s.t. no match is found spanning multiple sequences
		*(seq + seqLen) = '$';
		seqLen++;
		printf("Done reading a sequence of size %" PRIbwtint_t " from FASTA\n", seqLen);

		// pack and store the sequence
		for (bwtint_t i = 0; i < seqLen; i++)
		{
			if (seqBufferCount >= SEQ_BATCH_ALLOC_LEN)
			{
				pack_byte(seqBuffer, packedSeqBuffer, SEQ_BATCH_ALLOC_LEN);
				fwrite(packedSeqBuffer, 1, SEQ_BATCH_ALLOC_LEN / CHARS_PER_BYTE, pacFile);
				seqBufferCount = 0;
			}
			seqBuffer[seqBufferCount] = nt16_table[(unsigned int)seq[i]];
			//if(seqBuffer[seqBufferCount] == 10) { // ambiguous base
			// replace by a random ACGT base to improve performance
			//seqBuffer[seqBufferCount] = nt4_gray[lrand48()&3];
			//}
			seqBufferCount++;
		}

		// record the sequence range in the concatenated genome
		seqAnnotation->start_index = totalSeqLen;
		seqAnnotation->end_index = totalSeqLen + seqLen - 1;

		totalSeqLen += seqLen;
		annotations->num_seq++;
	}

	// pack any remaining chars
	if (seqBufferCount > 0)
	{
		pack_byte(seqBuffer, packedSeqBuffer, seqBufferCount);
		fwrite(packedSeqBuffer, 1, ceil((float)seqBufferCount / CHARS_PER_BYTE), pacFile);
		seqBufferCount = 0;
	}
	c = (char)(totalSeqLen % CHARS_PER_BYTE);
	fwrite(&c, 1, 1, pacFile);

	// store annotations info
	fprintf(annFile, "%" PRIbwtint_t "\t%d\n", totalSeqLen, annotations->num_seq);
	for (int i = 0; i < annotations->num_seq; i++)
	{
		seq_annotation_t seqAnnotation = annotations->seq_anns[i];
		fprintf(annFile, "%s\t%" PRIbwtint_t "\t%" PRIbwtint_t "\n", seqAnnotation.name, seqAnnotation.start_index, seqAnnotation.end_index);
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", totalSeqLen);

	fclose(fastaFile);
	fclose(pacFile);
	fclose(annFile);

	free(annotations->seq_anns);
	free(annotations);
	free(seq);
}

void ref2seq(const char *refFname, unsigned char **seq, bwtint_t *totalSeqLen)
{
	// input sequence
	FILE *refFile = (FILE *)fopen(refFname, "r");
	if (refFile == NULL)
	{
		printf("ref2seq: Cannot open .ref file: %s!\n", refFname);
		exit(1);
	}

	bwtint_t allocatedSeqLen = 50 * 1024 * 1024;
	*seq = (unsigned char *)malloc(allocatedSeqLen * sizeof(unsigned char));
	bwtint_t seqLen = 0;
	unsigned char c = (unsigned char)getc(refFile);
	while (!feof(refFile))
	{
		if (seqLen >= allocatedSeqLen)
		{
			allocatedSeqLen <<= 1;
			*seq = (unsigned char *)realloc(*seq, sizeof(unsigned char) * allocatedSeqLen);
			if (*seq == NULL)
			{
				printf("ref2seq: Could not allocate memory for the input sequence, size alloc'd = %" PRIbwtint_t "\n", allocatedSeqLen);
				exit(1);
			}
		}
		(*seq)[seqLen] = c;
		seqLen++;
		c = (unsigned char)getc(refFile);
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", seqLen);
	*totalSeqLen = seqLen;
}

// reads the input sequence data from the FASTA file, concatenates the subsequences,
// appends the reverse complement, encodes the resulting sequence and writes it to a .ref file;
// sequence annotations are stored into an .ann file
void fasta2ref(const char *fastaFname, const char *refFname, const char *annFname, unsigned char **seq, bwtint_t *totalSeqLen)
{
	// input sequence
	FILE *fastaFile = (FILE *)fopen(fastaFname, "r");
	if (fastaFile == NULL)
	{
		printf("fasta2ref: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}

	// annotations
	FILE *annFile = (FILE *)fopen(annFname, "wb");
	if (annFile == NULL)
	{
		printf("fasta2ref: Cannot open .ann file: %s!\n", annFname);
		exit(1);
	}

	// Write the initial header for the ANN file
	int readAnnNum = 0;
	fprintf(annFile, "%50" PRIbwtint_t "\t%50d\n", 0l, readAnnNum);

	// output combined ref sequence
	FILE *refFile;
	int save_ref = 0;
	if (refFname != NULL)
	{
		save_ref = 1;
		refFile = (FILE *)fopen(refFname, "wt");
		if (refFile == NULL)
		{
			printf("fasta2ref: Cannot open .ref file: %s!\n", refFname);
			exit(1);
		}
	}

	bwtint_t allocatedSeqLen = 50 * 1024 * 1024;
	*seq = (unsigned char *)malloc(allocatedSeqLen * sizeof(unsigned char));
	bwtint_t seqLen = 0;

	char c = (char)getc(fastaFile);
	if (c != '>')
		fasta_error(fastaFname);

	while (!feof(fastaFile))
	{
		seq_annotation_t *seqAnnotation = malloc(sizeof(seq_annotation_t));

		// sequence description line (> ...)
		c = (char)getc(fastaFile);
		bwtint_t seqNameLen = 0;
		while (c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(fastaFile))
		{
			seqAnnotation->name[seqNameLen] = c;
			seqNameLen++;
			c = (char)getc(fastaFile);
		}
		seqAnnotation->name[seqNameLen] = '\0';
		while (c != '\n' && !feof(fastaFile))
		{
			c = (char)getc(fastaFile);
		}
		if (feof(fastaFile))
			fasta_error(fastaFname);

		// sequence data
		bwtint_t subseqLen = 0;
		while (c != '>' && !feof(fastaFile))
		{
			if (c != '\n')
			{
				if (c >= 'a' && c <= 'z')
				{
					c += 'A' - 'a';
				}
				// reallocate twice as much memory
				if (seqLen >= allocatedSeqLen)
				{
					allocatedSeqLen <<= 1;
					*seq = (unsigned char *)realloc(*seq, sizeof(unsigned char) * allocatedSeqLen);
					if (*seq == NULL)
					{
						printf("fasta2ref: Could not allocate memory for the input sequence, size alloc'd = %" PRIbwtint_t "\n", allocatedSeqLen);
						exit(1);
					}
				}

				(*seq)[seqLen] = nt16_table[(unsigned int)c];
				seqLen++;
				subseqLen++;
				if (save_ref)
				{
					putc((int)(*seq)[seqLen - 1], refFile);
				}
			}
			c = (char)getc(fastaFile);
		}
		// add $ as a separator between chromosomes and bubbles
		// s.t. no match is found spanning multiple sequences
		(*seq)[seqLen] = nt16_table[(unsigned int)'$'];
		seqLen++;
		subseqLen++;
		if (save_ref)
		{
			putc((*seq)[seqLen - 1], refFile);
		}
		printf("Done reading a subsequence of size %" PRIbwtint_t " from FASTA\n", subseqLen);

		// record the sequence range in the concatenated genome
		seqAnnotation->start_index = seqLen - subseqLen;
		seqAnnotation->end_index = seqLen - 1;
		fprintf(annFile, "%s\t%" PRIbwtint_t "\t%" PRIbwtint_t "\n", seqAnnotation->name, seqAnnotation->start_index, seqAnnotation->end_index);
		free(seqAnnotation);

		readAnnNum++;
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", seqLen);

	// store annotations info
	fseek(annFile, 0, SEEK_SET);
	fprintf(annFile, "%50" PRIbwtint_t "\t%50d\n", seqLen, readAnnNum);

	// add the reverse complement
	*seq = (unsigned char *)realloc(*seq, sizeof(unsigned char) * 2 * seqLen);
	if (*seq == NULL)
	{
		printf("fasta2ref: Could not allocate memory for the reference sequence (including complement), of length = %" PRIbwtint_t "\n", 2 * seqLen);
		exit(1);
	}
	for (bwtint_t i = 0; i < seqLen; i++)
	{
		(*seq)[2 * seqLen - i - 1] = iupacCompl[(int)(*seq)[i]];
	}
	(*totalSeqLen) = 2 * seqLen;

	if (save_ref)
	{
		for (bwtint_t i = 0; i < seqLen; i++)
		{
			putc((*seq)[seqLen + i], refFile);
		}
		printf("Wrote %" PRIbwtint_t " chars to ref file\n", 2 * seqLen);
	}

	fclose(fastaFile);
	fclose(annFile);
	if (save_ref)
		fclose(refFile);
}

// load annotations from file
fasta_annotations_t *annf2ann(const char *annFname)
{
	FILE *annFile = (FILE *)fopen(annFname, "r");
	if (annFile == NULL)
	{
		printf("annf2ann: Cannot open ANN file: %s!\n", annFname);
		perror(annFname);
		exit(1);
	}

	fasta_annotations_t *annotations = (fasta_annotations_t *)calloc(1, sizeof(fasta_annotations_t));
	bwtint_t totalSeqLen;
	if (fscanf(annFile, "%" SCNbwtint_t "\t%d\n", &totalSeqLen, &(annotations->num_seq)) != 2)
	{
		printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
		exit(1);
	}
	annotations->seq_anns = (seq_annotation_t *)malloc(annotations->num_seq * sizeof(seq_annotation_t));

	for (int i = 0; i < annotations->num_seq; i++)
	{
		seq_annotation_t *seqAnnotation = &(annotations->seq_anns[i]);
		if (fscanf(annFile, "%[^\n\t]\t%" SCNbwtint_t "\t%" SCNbwtint_t "\n", (char *)&(seqAnnotation->name), &(seqAnnotation->start_index), &(seqAnnotation->end_index)) < 3)
		{
			printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
			exit(1);
		}
	}
	fclose(annFile);
	return annotations;
}

void free_ann(fasta_annotations_t *annotations)
{
	free(annotations->seq_anns);
	free(annotations);
}

// the packed reference sequence is read from the PAC file, uncompressed,
// concatenated with its reverse complement, and stored into seq
void pac2seq(const char *pacFname, unsigned char **seq, bwtint_t *totalSeqLen)
{
	FILE *pacFile = (FILE *)fopen(pacFname, "rb");
	if (pacFile == NULL)
	{
		printf("pac2seq: Cannot open PAC file %s!\n", pacFname);
		exit(1);
	}
	fseek(pacFile, -1, SEEK_END);
	// position in the file (# bytes from the beginning of file)
	bwtint_t pacFileLen = ftell(pacFile);
	if (pacFileLen < 0)
	{
		printf("pac2seq: Cannot determine the length of the PAC file!\n");
		exit(1);
	}
	unsigned char endByte;
	if (fread(&endByte, sizeof(unsigned char), 1, pacFile) < 1)
	{
		printf("pac2seq: Cannot read the PAC file!\n");
		exit(1);
	}
	bwtint_t seqLength = pacFileLen * CHARS_PER_BYTE - endByte;

	// read the packed sequence
	fseek(pacFile, 0, SEEK_SET);
	unsigned char *packedSeq = (unsigned char *)malloc(sizeof(char) * pacFileLen);
	if (fread(packedSeq, 1, pacFileLen, pacFile) < pacFileLen)
	{
		printf("pac2seq: Could not read the expected length of the PAC file!\n");
		exit(1);
	}
	fclose(pacFile);

	// unpack the sequence
	*seq = (unsigned char *)malloc(sizeof(char) * (2 * seqLength + 1));
	unpack_byte(packedSeq, *seq, seqLength);

	// add the reverse complement
	for (bwtint_t i = 0; i < seqLength; i++)
	{
		(*seq)[2 * seqLength - i - 1] = iupacCompl[(int)(*seq)[i]];
	}
	(*totalSeqLen) = 2 * seqLength;

	free(packedSeq);
}

void seq2rev_compl(unsigned char *seq, bwtint_t seqLen, unsigned char **rcSeq)
{
	*rcSeq = (unsigned char *)malloc(sizeof(char) * (seqLen + 1));
	for (bwtint_t i = 0; i < seqLen; i++)
	{
		(*rcSeq)[seqLen - i - 1] = iupacCompl[seq[i]];
	}
}

/* Reads I/O */

// loads the read sequences from the FASTQ file
reads_t *fastq2reads(const char *readsFname, int skip, int limit)
{
	FILE *readsFile = (FILE *)fopen(readsFname, "r");
	if (readsFile == NULL)
	{
		printf("load_reads_fastq: Cannot open reads file: %s !\n", readsFname);
		exit(1);
	}
	reads_t *reads = (reads_t *)calloc(1, sizeof(reads_t));
	int allocatedReads = NUM_READS_ALLOC;
	reads->reads = (read_t *)malloc(allocatedReads * sizeof(read_t));
	reads->count = 0;

	skip *= 4;
	while (skip-- > 0 && !feof(readsFile))
	{
		while (!feof(readsFile) && '\n' != (char)getc(readsFile))
		{
			// Do nothing with this character, it is being skipped while we look
			// for a new-line.
		}
	}

	char c;
	while (!feof(readsFile))
	{
		if (reads->count >= allocatedReads)
		{
			allocatedReads += NUM_READS_ALLOC;
			reads->reads = (read_t *)realloc(reads->reads, allocatedReads * sizeof(read_t));
		}

		read_t *read = &(reads->reads[reads->count]);

		c = (char)getc(readsFile);
		while (c != '@' && !feof(readsFile))
		{
			c = (char)getc(readsFile);
		}
		if (feof(readsFile))
			break;

		// line 1 (@ ...)
		int seqNameLen = 0;
		c = (char)getc(readsFile);
		while (c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(readsFile))
		{
			read->name[seqNameLen] = c;
			seqNameLen++;
			c = (char)getc(readsFile);
		}
		read->name[seqNameLen] = '\0';
		if (feof(readsFile))
			fastq_error(readsFname);

		while (c != '\n' && !feof(readsFile))
		{
			c = (char)getc(readsFile);
		}
		if (feof(readsFile))
			fastq_error(readsFname);

		// line 2 (sequence letters)
		int readLen = 0;
		int allocatedReadLen = READ_LENGTH_ALLOC;
		read->seq = (char *)malloc(allocatedReadLen * sizeof(char));
		read->rc = (char *)malloc(allocatedReadLen * sizeof(char));
		read->qual = (char *)malloc(allocatedReadLen * sizeof(char));

		c = (char)getc(readsFile);
		while (c != '\n' && !feof(readsFile))
		{
			if (readLen >= allocatedReadLen)
			{
				allocatedReadLen += READ_LENGTH_ALLOC;
				read->seq = (char *)realloc(read->seq, allocatedReadLen * sizeof(char));
				read->rc = (char *)realloc(read->rc, allocatedReadLen * sizeof(char));
				read->qual = (char *)realloc(read->qual, allocatedReadLen * sizeof(char));
			}
			read->seq[readLen] = nt4_table[(int)c];
			c = (char)getc(readsFile);
			readLen++;
		}
		read->len = readLen;
		if (feof(readsFile))
			fastq_error(readsFname);

		while (c != '+' && !feof(readsFile))
		{
			c = (char)getc(readsFile);
		}
		if (feof(readsFile))
			fastq_error(readsFname);

		// line 3 (+ ...)
		while (c != '\n' && !feof(readsFile))
		{
			c = (char)getc(readsFile);
		}
		if (feof(readsFile))
			fastq_error(readsFname);

		// line 4 (quality values)
		int qualLen = 0;
		c = (char)getc(readsFile);
		while (c != '\n' && !feof(readsFile))
		{
			if (qualLen <= readLen)
			{
				read->qual[qualLen] = c;
			}
			qualLen++;
			c = (char)getc(readsFile);
		}
		if (qualLen != readLen)
		{
			printf("Error: The number of quality score symbols does not match the length of the read sequence.\n");
			exit(1);
		}
		read->qual[qualLen] = '\0';

		// compute the reverse complement
		for (int i = 0; i < read->len; i++)
		{
			read->rc[read->len - 1 - i] = nt4_complement[(int)read->seq[i]];
		}

		if (read->len > reads->max_len)
		{
			reads->max_len = read->len;
		}
		reads->count++;

		if (reads->count == limit)
		{
			break;
		}
	}
	printf("Loaded %d reads from %s.\n", reads->count, readsFname);

	fclose(readsFile);
	return reads;
}

char *strdup(const char *str)
{
	int n = strlen(str) + 1;
	char *dup = malloc(n * sizeof(char));
	if (dup)
	{
		strcpy(dup, str);
	}
	return dup;
}

// Assumed format: @chr_lpos_rpos_strand_mpos1_mpos2_..._mposn
// Positions are 1-based
void parse_read_mapping(read_t *read)
{
	int _count = 0;
	int idx = 0;
	char c = read->name[idx];
	while (c != '\0')
	{
		if (c == '_')
		{
			_count++;
		}
		idx++;
		c = read->name[idx];
	}

	read->num_mref_pos = _count - 3;
	read->mref_pos = (bwtint_t *)malloc(read->num_mref_pos * sizeof(bwtint_t));

	int token_index = 0;
	const char delimiters[] = "_";
	char *read_name = strdup(read->name);
	char *token = strtok(read_name, delimiters);
	while (token != NULL)
	{
		if (token_index == 1)
		{
			sscanf(token, "%" SCNbwtint_t "", &read->ref_pos_l);
		}
		else if (token_index == 2)
		{
			sscanf(token, "%" SCNbwtint_t "", &read->ref_pos_r);
		}
		else if (token_index == 3)
		{
			read->strand = (strcmp(token, "nm") == 0) ? 0 : 1;
		}
		else if (token_index > 3)
		{
			sscanf(token, "%" SCNbwtint_t "", &read->mref_pos[token_index - 4]);
		}
		token = strtok(NULL, delimiters);
		token_index++;
	}
	free(read_name);
}

void free_reads(reads_t *reads)
{
	for (int i = 0; i < reads->count; i++)
	{
		if (reads->reads[i].seq)
			free(reads->reads[i].seq);
		if (reads->reads[i].rc)
			free(reads->reads[i].rc);
		if (reads->reads[i].qual)
			free(reads->reads[i].qual);
		if (reads->reads[i].aln_path)
			free(reads->reads[i].aln_path);
		if (reads->reads[i].mref_pos)
			free(reads->reads[i].mref_pos);
	}
	free(reads->reads);
	free(reads);
}

void print_read(read_t *read)
{
	printf("READ %s \n", read->name);
	printf("FWD: ");
	for (int i = 0; i < read->len; i++)
	{
		printf("%c", iupacChar[nt4_gray[(int)read->seq[i]]]);
	}
	printf("\n");
	printf("RC: ");
	for (int i = 0; i < read->len; i++)
	{
		printf("%c", iupacChar[nt4_gray[(int)read->rc[i]]]);
	}
	printf("\n");
}

/* Compression */

void pack_word(const unsigned char *input, uint32_t *output, const bwtint_t length)
{
	uint32_t c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_WORD; i++)
	{
		c = 0;
		j = i * CHARS_PER_WORD;
		for (k = 0; k < CHARS_PER_WORD; k++)
		{
			c = c | (input[j + k] << (BITS_IN_WORD - (k + 1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
	if (i * CHARS_PER_WORD < length)
	{
		c = 0;
		j = i * CHARS_PER_WORD;
		for (k = 0; j + k < length; k++)
		{
			c = c | (input[j + k] << (BITS_IN_WORD - (k + 1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
}

void unpack_word(const unsigned char *input, unsigned char *output, const bwtint_t length)
{
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_WORD; i++)
	{
		c = input[i];
		j = i * CHARS_PER_WORD;
		for (k = 0; k < CHARS_PER_WORD; k++)
		{
			output[j + k] = c >> (BITS_IN_WORD - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
	if (i * CHARS_PER_WORD < length)
	{
		c = input[i];
		j = i * CHARS_PER_WORD;
		for (k = 0; j + k < length; k++)
		{
			output[j + k] = c >> (BITS_IN_WORD - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
}

void pack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length)
{
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_BYTE; i++)
	{
		c = 0;
		j = i * CHARS_PER_BYTE;
		for (k = 0; k < CHARS_PER_BYTE; k++)
		{
			c = c | (unsigned char)(input[j + k] << (BITS_IN_BYTE - (k + 1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
	if (i * CHARS_PER_BYTE < length)
	{
		c = 0;
		j = i * CHARS_PER_BYTE;
		for (k = 0; j + k < length; k++)
		{
			c = c | (unsigned char)(input[j + k] << (BITS_IN_BYTE - (k + 1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
}

void unpack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length)
{
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_BYTE; i++)
	{
		c = input[i];
		j = i * CHARS_PER_BYTE;
		for (k = 0; k < CHARS_PER_BYTE; k++)
		{
			output[j + k] = c >> (BITS_IN_BYTE - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
	if (i * CHARS_PER_BYTE < length)
	{
		c = input[i];
		j = i * CHARS_PER_BYTE;
		for (k = 0; j + k < length; k++)
		{
			output[j + k] = c >> (BITS_IN_BYTE - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
}

#endif
