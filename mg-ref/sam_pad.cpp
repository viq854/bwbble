///////////////////////////////////////////
// for the reads aligned to bubbles
// compute their positions in the reference genome
// Lin Huang <linhuang@cs.stanford.edu>, 9 Feb 2015
///////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <sstream>
using namespace std;

typedef struct
{
	string ann;
	long long A;
	long long B_minus_A;
	long long C;
	long long D_minus_C;
	int ref_len;
	int alt_len;
} bubble_t;

void read_in_bubbles(string input_filename, vector<bubble_t> &bubbles)
{
	ifstream input;
	input.open(input_filename.c_str());

	string ann, line;
	long long A, B_minus_A, C, D_minus_C;
	int ref_len, alt_len;
	while(getline(input, ann))
	{
		getline(input, line);
		stringstream line_stream(line);
		line_stream >> A >> B_minus_A >> C >> D_minus_C >> ref_len >> alt_len;
		bubble_t bubble = {ann, A, B_minus_A, C, D_minus_C, ref_len, alt_len};
		bubbles.push_back(bubble);
	}

	input.close();
}

void process_sam(string input_filename, vector<bubble_t> bubbles, string output_filename)
{
	ifstream input;
        input.open(input_filename.c_str());
	ofstream output;
	output.open(output_filename.c_str());

	string line, qname, flag, ann, pos, bubble;
	size_t p;
	long long which_bubble, locus;
	while(getline(input, line))
	{
		if(line.length() >= 1 && line[0] == '@') {output << line << endl; continue;}
		stringstream line_stream(line);
		getline(line_stream, qname, '\t');
		getline(line_stream, flag, '\t');
		getline(line_stream, ann, '\t');
		getline(line_stream, pos, '\t');

		output << line;
		p = ann.find("bubble");
		if(!p)
		{
			stringstream ann_stream(ann);
			ann_stream >> bubble;
			which_bubble = atoll(bubble.substr(6).c_str());
			output << "\tbC:Z:" << bubbles[which_bubble].ann << "\tbP:Z:";

			locus = atoll(pos.c_str());
			if(locus >= 1 && locus <= bubbles[which_bubble].B_minus_A)
			{
				output << bubbles[which_bubble].A + locus - 1;
			}
			else if(locus >= bubbles[which_bubble].B_minus_A + bubbles[which_bubble].alt_len + 1 && locus <= bubbles[which_bubble].B_minus_A + bubbles[which_bubble].alt_len + bubbles[which_bubble].D_minus_C + 1)
			{
				output << locus + bubbles[which_bubble].C - (bubbles[which_bubble].B_minus_A + bubbles[which_bubble].alt_len + 1);
			}
			else
			{
				output << bubbles[which_bubble].B_minus_A + bubbles[which_bubble].A << "-" << bubbles[which_bubble].B_minus_A + bubbles[which_bubble].A + bubbles[which_bubble].ref_len - 1;
			}
		}
		output << endl;
	}

	input.close();
	output.close();
}

static int usage()
{
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: sam_pad <bubble.data> <sam.input> <sam.output>\n");
        fprintf(stderr, "\n");
        return 1;
}

int main(int argc, char* argv[])
{
	if(argc < 4) return usage();

	string meta_filename = argv[1];
	string sam_filename = argv[2];
	string output_filename = argv[3];

	vector<bubble_t> bubbles;
	read_in_bubbles(meta_filename, bubbles);

	process_sam(sam_filename, bubbles, output_filename);

	return 0;
}
