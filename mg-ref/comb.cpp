///////////////////////////////////////////
// Combine SNPs and INDELs with a reference
// genome (fasta file)
// Read SNPs from SNP.extract.chrxx.data
// Read INDELs from INDEL.extract.chrxx.data
// Lin Huang, 24 May 2012
///////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>
using namespace std;

#define MAX_CHR_LENGTH 1000000000 
#define WINDOW_SIZE 124
#define ALPHABET_SIZE 16
#define BASE_SIZE 4
#define BASE_SET_SIZE 8
#define MAX_OCC 1090
#define OCC_THRESHOLD 3

static const int gray_code[] =      { 0,   1,   3,   2,   6,   7,   5,   4,  12,  13,  15,  14,  10,  11,   9,   8};
static const char abbr[] = {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};
static const int base_order[] = {8 /*A*/, 4 /*C*/, 2 /*G*/, 1 /*T*/};
static const char base[] = {'A', 'C', 'G', 'T'};
static const char base_set[BASE_SIZE][BASE_SET_SIZE] = {{'A', 'N', 'M', 'H', 'V', 'R', 'D', 'W'}, /*A*/
                                {'C', 'N', 'S', 'B', 'Y', 'M', 'H', 'V'}, /*C*/
                                {'G', 'N', 'K', 'S', 'B', 'V', 'R', 'D'}, /*G*/
                                {'T', 'N', 'K', 'B', 'Y', 'H', 'D', 'W'}}; /*T*/

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

void print_multigenome(string multifasta_filename, string bubble_filename, char *chromosome, long long start, string header, long long flag, long long & total_snp_number, long long & low_end_snp_number, long long & high_end_snp_number)
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
			if(occ < OCC_THRESHOLD) 
			{
				low_end_snp_number++;
				continue;
			}

			if(occ > MAX_OCC - OCC_THRESHOLD)
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

void insert_SNP(string fasta_filename, string multifasta_filename, string bubble_filename)
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
        			print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, total_snp_number, low_end_snp_number, high_end_snp_number);
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
  	print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, total_snp_number, low_end_snp_number, high_end_snp_number);
	printf("total snp number is %lld\n", total_snp_number);
	printf("low end snp number is %lld\n", low_end_snp_number);
	printf("high end snp number is %lld\n", high_end_snp_number);

  	fasta.close();
}

void print_bubble(string chr, string bubble_filename, string data_filename, char* chromosome, long long& indel_count, long long start, int &flag, long long & total_indel_number, long long & low_end_indel_number)
{
  	string extract_filename;
  	long long pos, occ;
  	int i;
  	string indel_ref, indel_alt;

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
  	extract_filename += chr;
	extract_filename += ".data";
  	ext.open(extract_filename.c_str());

	if(ext.good())
	{
	  	while(ext >> pos >> indel_ref >> indel_alt >> occ)
  		{
			if(occ < OCC_THRESHOLD) 
			{
				low_end_indel_number++;
				continue;
			}

			total_indel_number++;

	    		bubble << ">" << "bubble" << indel_count << " " << chr << " " << get_max(pos - WINDOW_SIZE, 1) << endl;
			data << chr << " " << get_max(pos - WINDOW_SIZE, 1) << " " << get_min(WINDOW_SIZE, pos - 1) << " " << indel_alt.size() << " " << get_min(WINDOW_SIZE, start - pos - indel_ref.length()) << " " << indel_ref.size() << endl;
    			for(i = get_min(WINDOW_SIZE, pos - 1); i > 0; i--)
	    		{
      				bubble << chromosome[pos - i];
    			}
	    		if(indel_alt[0] != '.')
    			{
      				bubble << indel_alt;
	    		}
    			for(i = 0; i < get_min(WINDOW_SIZE, start - pos - indel_ref.length()); i++)
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

void comp_bubble(string multifasta_filename, string bubble_filename, string data_filename)
{
  	ifstream multifasta;
  	multifasta.open(multifasta_filename.c_str());
  	string line, chr;
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
        			print_bubble(chr, bubble_filename, data_filename, chromosome, indel_count, start, flag, total_indel_number, low_end_indel_number);
      			}

      			start = 1;
      			stringstream line_stream(line);
      			line_stream >> chr;
      			chr.erase(chr.begin());
    		}
    		else
    		{
      			for(i = 0; i < line.length(); i++)
      			{
        			chromosome[start] = line[i], start++;
      			}
    		}
  	}
  	print_bubble(chr, bubble_filename, data_filename, chromosome, indel_count, start, flag, total_indel_number, low_end_indel_number);
	printf("total indel number is %lld\n", total_indel_number);
	printf("low end indel number is %lld\n", low_end_indel_number);

  	multifasta.close();
}

static int usage()
{
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage: comb <input.fasta> <output.fasta> <output_bubble.fasta> <bubble.data>\n");
  	fprintf(stderr, "\n");
  	return 1;
}

int main(int argc, char* argv[])
{
  	if(argc != 5) return usage();

  	string fasta_filename = argv[1];
  	string multifasta_filename = argv[2];
  	string bubble_filename = argv[3];
	string data_filename = argv[4];

  	insert_SNP(fasta_filename, multifasta_filename, bubble_filename);
  	comp_bubble(multifasta_filename, bubble_filename, data_filename);

	return 0;
}
