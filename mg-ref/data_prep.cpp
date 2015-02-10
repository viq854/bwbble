///////////////////////////////////////////
// Extract SNPs and INDELs from vcf file
// Write SNPs into SNP.extract.chrxx.data
// Write INDELs into INDEL.extract.chrxx.data
// Lin Huang <linhuang@cs.stanford.edu>, 24 May 2012
///////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;

#define CHROM_POS 0
#define POS_POS 1
#define REF_POS 3
#define ALT_POS 4
#define FILTER_POS 6
#define INFO_POS 7
#define FORMAT_POS 8

int pass_screen(string chr, string pos, string ref, string alt, string filter, string info, int allele_freq)
{
  	if(info.find("VT=SNP") != string::npos || info.find("VT=INDEL") != string::npos)
  	{
    		return 1;
  	}
  	return 0;
}

int has(vector<string> chr_list, string chr)
{
  	return find(chr_list.begin(), chr_list.end(), chr) != chr_list.end();
}

vector<string> vcf_extract(string input_filename, int par_clear, vector<string> chr_list)
{
  	int i;
  	string line, attr, chr, chr_ori = "", pos, ref, alt, filter, info, alt_case;
  	int attr_count = 0, allele_freq;
  	vector<string> new_chr_list;

  	string SNP_filename, INDEL_filename;
  	ifstream vcf;
  	vcf.open(input_filename.c_str());
  	ofstream SNP, INDEL;

  	while(getline(vcf, line))
  	{
    		if(line[0] == '#' && line[1] == '#')
    		{
      			continue;
    		}
    		else
    		{
      			break;
    		}
  	}

  	stringstream line_stream(line);
  	while(getline(line_stream, attr, '\t'))
  	{
    		attr_count++;
  	}

  	while(getline(vcf, line))
  	{
		allele_freq = 0;
    		stringstream line_stream(line);
    		for(i = 0; i < attr_count; i++)
    		{
      			getline(line_stream, attr, '\t');
      			switch(i)
      			{
        			case CHROM_POS:
          				chr = attr;
          				break;
        			case POS_POS:
          				pos = attr;
          				break;
        			case REF_POS:
          				ref = attr;
          				break;
        			case ALT_POS:
          				alt = attr;
          				break;
        			case FILTER_POS:
          				filter = attr;
          				break;
        			case INFO_POS:
          				info = attr;
          				break;
        			default:
          				break;
      			}
			if(i > FORMAT_POS)
			{
				if(attr[0] == '1' || attr[2] == '1') allele_freq++;
			}
    		}

    		if(pass_screen(chr, pos, ref, alt, filter, info, allele_freq))
    		{
      			if(chr.compare(chr_ori) != 0)
      			{
        			SNP_filename = "mg-ref-output/SNP.extract.chr";
        			SNP_filename += chr;
        			SNP_filename += ".data";

        			INDEL_filename = "mg-ref-output/INDEL.extract.chr";
        			INDEL_filename += chr;
        			INDEL_filename += ".data";

        			if(SNP.is_open())
        			{
          				SNP.close();
        			}
        			if(INDEL.is_open())
        			{
          				INDEL.close();
        			}

        			if(par_clear && !has(chr_list, chr) && !has(new_chr_list, chr))
        			{
          				SNP.open(SNP_filename.c_str());
          				INDEL.open(INDEL_filename.c_str());
          				new_chr_list.push_back(chr);
        			}
        			else
        			{
          				SNP.open(SNP_filename.c_str(), ios::out | ios::app);
          				INDEL.open(INDEL_filename.c_str(), ios::out | ios::app);
        			}

        			chr_ori = chr;
      			}

      			stringstream alt_stream(alt);
      			while(getline(alt_stream, alt_case, ','))
      			{
        			if(ref.length() == 1 && alt_case.length() == 1 && alt_case[0] != '.')
        			{
          				SNP << pos << "\t" << ref << "\t" << alt_case << "\t" << allele_freq << endl;
        			}
        			else if(ref.length() != alt_case.length() || (ref.length() == 1 && alt_case.length() == 1 && alt_case[0] == '.'))
        			{
          				INDEL << pos << "\t" << ref << "\t" << alt_case << "\t" << allele_freq << endl;
        			}
      			}
    		}	
  	}

  	vcf.close();
  	SNP.close();
  	INDEL.close();

  	return new_chr_list;
}

static int usage()
{
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage:   data_prep [option] <input1.vcf> <input2.vcf> ... \n");
  	fprintf(stderr, "Option:  -c  clear all SNP.extract.chrxx.data and INDEL.extract.chrxx.data files before usage\n");
  	fprintf(stderr, "Example: data_prep -c ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf\n");
  	fprintf(stderr, "\n");
  	return 1;
}

int main(int argc, char* argv[])
{
  	if(argc < 2) return usage();

  	int par_clear;
  	if(!strcmp(argv[1], "-c"))
  	{
    		par_clear = 1;
  	}
  	else
  	{
  	  	par_clear = 0;
  	}
  
  	int i;
  	vector<string> chr_list, new_chr_list;
  	string input_filename;

  	for(i = 1 + par_clear; i < argc; i++)
  	{
    		input_filename = argv[i];
    		cout << input_filename << endl;
    		new_chr_list = vcf_extract(input_filename, par_clear, chr_list);
    		chr_list.insert(chr_list.end(), new_chr_list.begin(), new_chr_list.end());
    		vector<string>::iterator it;
    		for(it = chr_list.begin(); it < chr_list.end(); it++)
		{
      			cout << *it << " ";
		}
    		cout << endl;
  	}

  	return 0;
}
