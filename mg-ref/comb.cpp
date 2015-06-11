///////////////////////////////////////////
// Combine SNPs and INDELs with a reference
// genome (fasta file)
// Read SNPs from SNP.extract.chrxx.data
// Read INDELs from INDEL.extract.chrxx.data
// Lin Huang <linhuang@cs.stanford.edu>, 24 May 2012
///////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <getopt.h>
using namespace std;

#define MAX_CHR_LENGTH 1000000000 
#define ALPHABET_SIZE 16
#define BASE_SIZE 4
#define BASE_SET_SIZE 8

static const int gray_code[] =      { 0,   1,   3,   2,   6,   7,   5,   4,  12,  13,  15,  14,  10,  11,   9,   8};
static const char abbr[] = {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};
static const int base_order[] = {8 /*A*/, 4 /*C*/, 2 /*G*/, 1 /*T*/};
static const char base[] = {'A', 'C', 'G', 'T'};
static const char base_set[BASE_SIZE][BASE_SET_SIZE] = {{'A', 'N', 'M', 'H', 'V', 'R', 'D', 'W'}, /*A*/
                                {'C', 'N', 'S', 'B', 'Y', 'M', 'H', 'V'}, /*C*/
                                {'G', 'N', 'K', 'S', 'B', 'V', 'R', 'D'}, /*G*/
                                {'T', 'N', 'K', 'B', 'Y', 'H', 'D', 'W'}}; /*T*/

typedef struct
{
	int window_size;
	int min_occ;
	int max_occ;
	int min_occ_specified;
        int max_occ_specified;
} pars_t;

long long get_min(long long a, long long b)
{
  	if(a > b) return b; else return a;
}

long long get_max(long long a, long long b)
{
  	if(a > b) return a; else return b;
}

int inSet(char test, char test_base)
{
  	int i, j;
  	for(i = 0; i < BASE_SIZE; i++)
  	{
    		if(test_base == base[i])
    		{
      			for(j = 0; j < BASE_SET_SIZE; j++)
      			{
        			if(test == base_set[i][j] || test == base_set[i][j] - 'A' + 'a')
        			{
          				return 1;
        			}
      			}
      			return 0;
    		}
  	}
}

void print_multigenome(string multifasta_filename, string bubble_filename, char *chromosome, long long start, string header, long long flag, long long & total_snp_number, long long & low_end_snp_number, long long & high_end_snp_number, pars_t* pars)
{
  	string chr;
  	stringstream line_stream(header);
  	line_stream >> chr;
  	chr.erase(chr.begin());

  	string extract_filename;
  	ifstream ext;
  	extract_filename = "mg-ref-output/SNP.extract.chr";
  	extract_filename += chr;
	extract_filename += ".data";
  	ext.open(extract_filename.c_str());

  	ofstream multifasta, bubble;
  	if(flag == 1)
  	{
    		multifasta.open(multifasta_filename.c_str());
    		bubble.open(bubble_filename.c_str());
  	}
  	else
  	{
    		multifasta.open(multifasta_filename.c_str(), ios::out | ios::app);
    		bubble.open(bubble_filename.c_str(), ios::out | ios::app);
  	}

  	multifasta << header << endl;
  	bubble << header << endl;

	long long i;
	if(ext.good())
	{
	  	long long pos, occ;
  		char ref, alt;
	  	int bit[BASE_SIZE], total;
  		while(ext >> pos >> ref >> alt >> occ)
	  	{
			if(pars->min_occ_specified && occ < pars->min_occ) 
			{
				low_end_snp_number++;
				continue;
			}

			if(pars->max_occ_specified && occ > pars->max_occ)
			{
				high_end_snp_number++;
				chromosome[pos] = alt;
				continue;
			}

    			total_snp_number++;
    			for(i = 0; i < BASE_SIZE; i++)
	    		{
      				if(inSet(chromosome[pos], base[i]) || inSet(ref, base[i]) || inSet(alt, base[i]))
      				{
        				bit[i] = 1;
	      			}
      				else
      				{
        				bit[i] = 0;
	      			}
    			}

    			total = 0;
	    		for(i = 0; i < BASE_SIZE; i++)
    			{
      				total += (bit[i] * base_order[i]);
	    		}

    			for(i = 0; i < ALPHABET_SIZE; i++)
    			{
      				if(gray_code[i] == total)
	      			{
        				chromosome[pos] = abbr[i];
      				}
    			}
	  	}
		ext.close();
	}

  	for(i = 1; i < start; i++)
 	{
    		multifasta << chromosome[i];
    		bubble << chromosome[i];
    		if(!(i % 60))
    		{
      			multifasta << endl;
      			bubble << endl;
    		}
  	}
  	if((start - 1) % 60)
  	{
    		multifasta << endl;
    		bubble << endl;
  	}

  	multifasta.close();
  	bubble.close();
}

void insert_SNP(string fasta_filename, string multifasta_filename, string bubble_filename, pars_t* pars)
{
  	ifstream fasta;
  	fasta.open(fasta_filename.c_str());

  	string header, line;
  	long long i, start, g = 0;
	long long total_snp_number = 0, low_end_snp_number = 0, high_end_snp_number = 0;
  	char *chromosome;

  	chromosome = (char*)malloc(MAX_CHR_LENGTH * sizeof(char));

  	while(getline(fasta, line))
  	{
    		if(line.length() >= 1 && line[0] == '>')
    		{
      			if(g)
      			{
        			print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, total_snp_number, low_end_snp_number, high_end_snp_number, pars);
      			}
      			g++;

      			start = 1;
      			header = line;
    		}
    		else
    		{
      			for(i = 0; i < line.length(); i++)
      			{
        			chromosome[start] = line[i], start++;
      			}
    		}
  	}
  	print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, total_snp_number, low_end_snp_number, high_end_snp_number, pars);
	printf("total snp number is %lld\n", total_snp_number);
	printf("low end snp number is %lld\n", low_end_snp_number);
	printf("high end snp number is %lld\n", high_end_snp_number);

  	fasta.close();
}

void print_bubble(string chr, string schr, string bubble_filename, string data_filename, char* chromosome, long long& indel_count, long long start, int &flag, long long & total_indel_number, long long & low_end_indel_number, pars_t* pars)
{
  	string extract_filename;
  	long long pos, occ;
  	int i;
  	string indel_ref, indel_alt;
	long long A, B_minus_A, C, D_minus_C, ref_len, alt_len;

  	ifstream ext;
  	ofstream bubble;
  	bubble.open(bubble_filename.c_str(), ios::out | ios::app);
	ofstream data;
	if(!flag)
	{
		data.open(data_filename.c_str());
		flag = 1;
	}
	else
	{
		data.open(data_filename.c_str(), ios::out | ios::app);
	}

  	extract_filename = "mg-ref-output/INDEL.extract.chr";
  	extract_filename += schr;
	extract_filename += ".data";
  	ext.open(extract_filename.c_str());

	if(ext.good())
	{
	  	while(ext >> pos >> indel_ref >> indel_alt >> occ)
  		{
			total_indel_number++;

	    		bubble << ">" << "bubble" << indel_count << " " << chr << " " << get_max(pos - pars->window_size, 1) << endl;
			A = get_max(pos - pars->window_size, 1);
			B_minus_A = get_min(pars->window_size, pos - 1);
			C = pos + indel_ref.length();
			D_minus_C = get_min(pars->window_size, start - pos - indel_ref.length()) - 1;
			if(indel_ref[0] != '.') ref_len = indel_ref.length(); else ref_len = 0;
			if(indel_alt[0] != '.') alt_len = indel_alt.length(); else alt_len = 0;

			data << chr << endl;
			data << A << "\t" << B_minus_A << "\t" << C << "\t" << D_minus_C << "\t" << ref_len << "\t" << alt_len << endl;

    			for(i = get_min(pars->window_size, pos - 1); i > 0; i--)
	    		{
      				bubble << chromosome[pos - i];
    			}
	    		if(indel_alt[0] != '.')
    			{
      				bubble << indel_alt;
	    		}
    			for(i = 0; i < get_min(pars->window_size, start - pos - indel_ref.length()); i++)
    			{
      				bubble << chromosome[pos + indel_ref.length() + i];
	    		}
    			bubble << endl;

    			indel_count++;
	  	}

  		ext.close();
	}
  	bubble.close();
	data.close();
}

void comp_bubble(string multifasta_filename, string bubble_filename, string data_filename, pars_t* pars)
{
  	ifstream multifasta;
  	multifasta.open(multifasta_filename.c_str());
  	string line, chr, schr;
  	int total, g = 0, flag = 0;
  	int i, j;
  	long long indel_count = 0, start;
	long long total_indel_number = 0, low_end_indel_number = 0;
  	char *chromosome;

  	chromosome = (char*)malloc(MAX_CHR_LENGTH * sizeof(char));

  	while(getline(multifasta, line))
  	{
    		if(line.length() >= 1 && line[0] == '>')
    		{
      			if(!g)
      			{
        			g++;
      			}
      			else
      			{
        			print_bubble(chr, schr, bubble_filename, data_filename, chromosome, indel_count, start, flag, total_indel_number, low_end_indel_number, pars);
      			}

      			start = 1;
			chr = line;
			chr.erase(chr.begin());

      			stringstream line_stream(line);
      			line_stream >> schr;
      			schr.erase(schr.begin());
    		}
    		else
    		{
      			for(i = 0; i < line.length(); i++)
      			{
        			chromosome[start] = line[i], start++;
      			}
    		}
  	}
  	print_bubble(chr, schr, bubble_filename, data_filename, chromosome, indel_count, start, flag, total_indel_number, low_end_indel_number, pars);
	printf("total indel number is %lld\n", total_indel_number);

  	multifasta.close();
}

static int usage()
{
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage: comb <input.fasta> <output.fasta> <output_bubble.fasta> <bubble.data>\n");
	fprintf(stderr, "Option:  -w INT  window size [default: 124]\n");
	fprintf(stderr, "         -i INT  minimum occurrence\n");
	fprintf(stderr, "         -a INT  maximum occurrence\n");
  	fprintf(stderr, "\n");
  	return 1;
}

void set_default_pars(pars_t* pars)
{
	pars->window_size = 124;
	pars->min_occ_specified = 0;
	pars->max_occ_specified = 0;
}

int main(int argc, char* argv[])
{
  	if(argc < 5) return usage();

	pars_t* pars = (pars_t*) calloc(1, sizeof(pars_t));
	set_default_pars(pars);
	int c;
	while ((c = getopt(argc, argv, "w:i:a:")) >= 0) 
	{
			switch (c) {
				case 'w': 	pars->window_size = atoi(optarg); 
					  	if(pars->window_size < 0)
					  	{
							fprintf(stderr, "window size shouldn't be negative.\n");
							return 1;
					  	}
					  	break;
				case 'i':       pars->min_occ_specified = 1;
						pars->min_occ = atoi(optarg);
                                                break;
				case 'a':       pars->max_occ_specified = 1;
						pars->max_occ = atoi(optarg);
                                                break;
				case '?': 	usage(); return 1;
				default: 	return 1;
			}
	}

  	string fasta_filename = argv[optind];
  	string multifasta_filename = argv[optind+1];
  	string bubble_filename = argv[optind+2];
	string data_filename = argv[optind+3];

  	insert_SNP(fasta_filename, multifasta_filename, bubble_filename, pars);
  	comp_bubble(multifasta_filename, bubble_filename, data_filename, pars);

	return 0;
}
