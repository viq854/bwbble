#ifndef IO_H_
#define IO_H_

#include "align.h"
#include "common.h"

#define SEQ_BATCH_ALLOC_LEN			262144
#define READ_LENGTH_ALLOC  150
#define ANN_ALLOC_LEN				256
#define MAX_SEQ_NAME_LEN 			256
#define NUM_READS_ALLOC 			1000000

// Compression
#define BITS_PER_CHAR               4 // number of bits per character in the BWT
#define CHARS_PER_BYTE              (BITS_IN_BYTE / BITS_PER_CHAR)
#define CHARS_PER_WORD              8//(BITS_IN_WORD / BITS_PER_CHAR)
#define CHARS_PER_128BITS           32
#define BYTES_IN_WORD 				4
#define BITS_IN_WORD 				32
#define BITS_IN_BYTE 				8
#define CHAR_MASK					15 //((1<<BITS_PER_CHAR)-1)
#define CHAR_COUNT_2BIT_MASK		3
#define LS_2BYTES_MASK				65535//((1<<(2*BITS_IN_BYTE))-1)
#define MS_2BYTES_MASK				(LS_2BYTES_MASK << (2*BITS_IN_BYTE)) //4294901760

// ALPHABET
#define ALPHABET_SIZE 16
static const unsigned char iupacChar[16] =  {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};
static const unsigned char grayVal[16] =    { 0,   1,   3,   2,   6,   7,   5,   4,   12,  13,  15,  14,  10,  11,  9,   8 };
static const unsigned char iupacGrayOrd[16]={ 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15};
static const unsigned char grayIupac[16] =  { 15,  0,   2,   1,   6,   5,   3,   4,   14,  13,  11,  12,  7,   8,   10,  9 };
static const unsigned char iupacCompl[16] = { 0,   15,  8,   7,   4,   11,  12,  3,   2,   13,  10,  5,   6,   9,   14,  1};
static const unsigned char is_snp[16] = 	{ 0,   0,   1,   0,   1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   0};
static const uint32_t char_masks[16] = { 0x0, 0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666,
											  0x77777777, 0x88888888, 0x99999999, 0xAAAAAAAA, 0xBBBBBBBB, 0xCCCCCCCC,
											  0xDDDDDDDD, 0xEEEEEEEE, 0xFFFFFFFF};
static const uint32_t char_masks128[16][4] __attribute__((aligned(16))) =
													  { {0x0, 		0x0, 		0x0, 		0x0},
														{0x11111111, 0x11111111, 0x11111111, 0x11111111},
														{0x22222222, 0x22222222, 0x22222222, 0x22222222},
														{0x33333333, 0x33333333, 0x33333333, 0x33333333},
														{0x44444444, 0x44444444, 0x44444444, 0x44444444},
														{0x55555555, 0x55555555, 0x55555555, 0x55555555},
														{0x66666666, 0x66666666, 0x66666666, 0x66666666},
														{0x77777777, 0x77777777, 0x77777777, 0x77777777},
														{0x88888888, 0x88888888, 0x88888888, 0x88888888},
														{0x99999999, 0x99999999, 0x99999999, 0x99999999},
														{0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA},
														{0xBBBBBBBB, 0xBBBBBBBB, 0xBBBBBBBB, 0xBBBBBBBB},
														{0xCCCCCCCC, 0xCCCCCCCC, 0xCCCCCCCC, 0xCCCCCCCC},
														{0xDDDDDDDD, 0xDDDDDDDD, 0xDDDDDDDD, 0xDDDDDDDD},
														{0xEEEEEEEE, 0xEEEEEEEE, 0xEEEEEEEE, 0xEEEEEEEE},
														{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}};
#define CHAR_MASK_A 0xFFFFFFFF
#define CHAR_MASK_G 0x33333333
#define CHAR_MASK_C 0x77777777
#define CHAR_MASK_T 0x11111111

static const uint32_t word_fractions_masks128[32][4] __attribute__((aligned(16))) =
													 { {0x0FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x000FFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x0000FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000FFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x000000FF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x0000000F, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x000FFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0000FFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000FFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x000000FF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0000000F, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0FFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00FFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x000FFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0000FFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000FFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x000000FF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0000000F, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0FFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x00FFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x000FFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0000FFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x00000FFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x000000FF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0000000F},
													   {0x00000000, 0x00000000, 0x00000000, 0x00000000},
													};

static const uint32_t word_fractions_masks[8] = { 0x0FFFFFFF, 0x00FFFFFF, 0x000FFFFF, 0x0000FFFF,
                                                  0x00000FFF, 0x000000FF, 0x0000000F, 0x0};


#define NUM_NUCLEOTIDES 4
#define BASES_PER_NUCLEOTIDE 7 // (N is ignored)

// bases are sorted by gray code order
static const unsigned char nucl_bases_table[NUM_NUCLEOTIDES][BASES_PER_NUCLEOTIDE] = {
			{8,  9, 11, 12, 13, 14, 15} /*A*/,
			{2,  3,  4,  5, 11, 12, 13} /*G*/,
			{4,  5,  6,  7, 8,  9,  11} /*C*/,
			{1,  2,  5,  6, 9,  13, 14} /*T*/};

static const unsigned char nt4_gray[5] = {15/*A*/, 3/*G*/, 7/*C*/, 1/*T*/, 10/*N*/};
static const unsigned char nt4_gray_val[5] = {8/*A*/, 2/*G*/, 4/*C*/, 1/*T*/, 15/*N*/};
static const unsigned char nt4_complement[5] = {3/*A*/, 2/*G*/, 1/*C*/, 0/*T*/, 4/*N*/};

// encoding: A=0, G=1, C=2, T=3, N=4
static const unsigned char nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*A*/,4,2/*C*/,4,4,4,1/*G*/,4,4,4,4,4,4,4/*N*/,4,
	4, 4, 4, 4,  3/*T*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*a*/,4,2/*c*/,4,4,4,1/*g*/,4,4,4,4,4,4,4/*n*/,4,
	4, 4, 4, 4,  3/*t*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static const unsigned char nt16_table[256] = {
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    0/*$*/,10,10,10, 10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 15/*A*/,5/*B*/,7/*C*/,  13/*D*/,10, 10,3/*G*/,      9/*H*/,10,10,2/*K*/,  10,8/*M*/,10/*N*/,10,
	10, 10, 12/*R*/,4/*S*/,     1/*T*/,10,11/*V*/,14/*W*/,  10,6/*Y*/,10,10,      10, 10, 10, 10,
	10, 15/*a*/,5/*b*/,7/*c*/,  13/*d*/,10, 10,3/*g*/,      9/*h*/,10,10,2/*k*/,  10,8/*m*/,10/*n*/,10,
	10, 10, 12/*r*/,4/*s*/,     1/*t*/,10,11/*v*/,14/*w*/,  10,6/*y*/,10,10,      10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
};

typedef struct {
	// read length
	int len;
	// read name
	char name[MAX_SEQ_NAME_LEN+1];
	// read sequence
	char* seq;
	// reverse complement sequence
	char* rc;
	// quality scores
	char* qual;

	// original mapping information
	int strand;
	bwtint_t ref_pos_l;
	bwtint_t ref_pos_r;
	bwtint_t* mref_pos;
	int num_mref_pos;

	// alignment information
	alns_t* alns;
	int aln_type;
	int aln_top1_count;
	int aln_top2_count;
	int mapQ; // single-end mapping quality
	int num_mm;
	int num_gapo;
	int num_gape;
	int aln_strand; // fwd or rev
	int aln_score; // aln score
	bwtint_t aln_pos; // position in the reference
	bwtint_t aln_sa;
	int aln_length;
	char* aln_path;
} read_t;

// collection of reads
typedef struct {
	// number of reads
	unsigned int count;
	unsigned int max_len;
	// compressed read sequence
	read_t* reads;
} reads_t;

typedef struct {
	char name[MAX_SEQ_NAME_LEN];
	// range in the combined genome
	bwtint_t start_index;
	bwtint_t end_index;
} seq_annotation_t;

typedef struct {
	// number of different sequences in the fasta
	int num_seq;
	seq_annotation_t* seq_anns;
} fasta_annotations_t;

void fasta2ref(const char *fastaFname, const char* refFname, const char* annFname, unsigned char** seq, bwtint_t *totalSeqLen);
void ref2seq(const char* refFname, unsigned char** seq, bwtint_t* seqLen);
void fasta2pac(const char *fastaFname, const char* pacFname, const char* annFname);
void pac2seq(const char *pacFname, unsigned char** seq, bwtint_t *seqLength);
reads_t* fastq2reads(const char *readsFname);
void seq2rev_compl(unsigned char* seq, bwtint_t seqLen, unsigned char** rcSeq);
void parse_read_mapping(read_t* read);

void print_read(read_t* read);
void free_reads(reads_t* reads);

fasta_annotations_t* annf2ann(const char *annFname);
void free_ann(fasta_annotations_t* annotations);

// Compression
void pack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length);
void unpack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length);
void pack_word(const unsigned char *input, unsigned int *output, const bwtint_t length);
void unpack_word(const unsigned char *input, unsigned char *output, const bwtint_t length);


#endif
