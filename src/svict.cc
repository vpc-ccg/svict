#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <string>
#include <map>
#include <zlib.h>
#include <getopt.h>
#include "partition.h"
#include "common.h"
#include "assembler.h"
#include "svict_caller.h"
#include "genome.h"
#include "simulator.h"
#include "common.h"
#include "record.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "extractor.h"
#include "sort.h"

using namespace std;

char versionNumber[10]="1.0.1";
// 10.46
// S:[step]:END

/********************************************************************/
inline string space (int i) 
{
	return string(i, ' ');
}

/********************************************************************/
inline string itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}

inline bool check_input(string filename, bool resume, vector<string> &chromosomes){

	Parser* parser;
	FILE *fi = fopen(filename.c_str(), "rb");

	bool ok = true;
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	int ftype = 1; // 0 for BAM and 1 for SAM
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		ftype = 0;

	if ( !ftype )
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();

	if(comment == "" || comment[0] != '@'){
		E("Error: Input file is not a BAM/SAM or has a missing/malformed header.\n");
		ok = false;
	}
	else{
		
		vector<string> comments = split(comment, '\n');
		int sorted = 0;
		vector<string> hd = split( comments[0], '\t');
		for ( int id = 0; id < hd.size(); id++ )
		{
			if (hd[id] == "SO:coordinate")
			{ sorted = 1; break;}
		}

		if(  !sorted && !resume){
			E("Error: Input BAM/SAM does not appear to be coordinate sorted as required.\n");
			ok = false;


		}
		for(int id = 1; id < comments.size(); id++)
		{
			if ( !(strncmp( "@SQ", comments[id].c_str(), 3) ) )
			{
				chromosomes.push_back( split( comments[id], '\t')[1].substr(3) );
			}
		}
	}

	delete parser;
	
	return ok;
}

inline bool check_ref(string filename){

	string line;
	ifstream ref_file(filename);

	if (ref_file.is_open()) {
		getline(ref_file, line);
		if(line[0] != '>'){
			cerr << "Error: Specified reference file (" + filename + ") does not appear to be in fasta format." << endl;
			return false;
		}
	}
	else{
		cerr << "Error: Could not open reference file: " + filename << endl;
		return false;
	}

	return true;
}

inline bool check_gtf(string filename){

	string line;
	ifstream gtf_file(filename);

	if (gtf_file.is_open()) {
		getline(gtf_file, line);
		if(line[0] != '#'){
			cerr << "Error: Specified GTF file (" + filename + ") does not appear to be correctly formatted, or is missing a header." << endl;
			return false;
		}
	}
	else{
		cerr << "Error: Could not open GTF file: " + filename << endl;
		return false;
	}

	return true;
}

/********************************************************************/
void predict (const string &input_file, const string &reference, const string &gtf, const string &out_vcf, const string &print_fastq, const bool print_reads, const bool print_stats, 
					int k, int assembler_overlap, int anchor_len, int min_support, int max_support, int uncertainty, int min_length, int max_length, int sub_optimal, const bool LOCAL_MODE,
					int window_size, int min_sc, int max_fragment_size, double clip_ratio, bool use_indel, bool heuristic, const vector<string> &chromosomes, bool resume)
{
	svict_caller predictor(k, assembler_overlap, anchor_len, input_file, reference, gtf, print_reads, print_stats, chromosomes);
	if(resume){
		predictor.resume(out_vcf, input_file, print_fastq, uncertainty, min_length, max_length, sub_optimal, LOCAL_MODE);
	}
	else{ 
		predictor.run(out_vcf, print_fastq, min_support, max_support, uncertainty, min_length, max_length, sub_optimal, LOCAL_MODE, window_size, min_sc, max_fragment_size, clip_ratio, use_indel, heuristic);
	}
}

/********************************************************************/
void simulate_SVs (string in_file, string ref_file, bool from_bed) {

	simulator sim(ref_file);

	if(from_bed){
		sim.simulate(in_file);
	}
	else{
		sim.simulate_from_file(in_file);
	}
	
}

/********************************************************************/
void annotate (string gtf, string in_file, string out_file, bool genomic) {

	string line, chr, best_gene_s, best_name_s, best_trans_s, best_gene_e, best_name_e, best_trans_e, context_s, context_e, key, gene_id, trans_id;
	int start, end, t_start, t_end;
	ifstream infile (in_file);
	ofstream outfile (out_file);
	vector<uint32_t> vec_best_s, vec_best_e;
	vector<string> tokens;
	map <string,vector<isoform>> iso_gene_map;
	map <string,vector<gene_data>> gene_sorted_map;
	map <string,vector<uint32_t>> vec_all_s, vec_all_e;
	vector<string> contexts ={"intergenic", "intronic", "non-coding", "UTR", "exonic-CDS"};
	ensembl_Reader(gtf.c_str(), iso_gene_map, gene_sorted_map );

	if (infile.is_open()){
		while ( getline (infile,line)){

			tokens = split(line, '\t'); 
			chr = tokens[0];
			start = atoi(tokens[1].c_str());
			end = atoi(tokens[2].c_str());

			if(genomic){
				locate_interval(chr, start, start, gene_sorted_map[chr], 0, iso_gene_map, best_gene_s, best_name_s, best_trans_s, vec_best_s, vec_all_s);
				locate_interval(chr, end, end, gene_sorted_map[chr], 0, iso_gene_map, best_gene_e,  best_name_e, best_trans_e, vec_best_e, vec_all_e);
			}
			else{
				gene_id = tokens[3];
				trans_id = tokens[4];

				locate_interval(chr, start, start, gene_id, trans_id, gene_sorted_map[chr], 0, iso_gene_map, best_gene_s, best_name_s, best_trans_s, vec_best_s, vec_all_s);
				locate_interval(chr, end, end, gene_id, trans_id, gene_sorted_map[chr], 0, iso_gene_map, best_gene_e, best_name_e, best_trans_e, vec_best_e, vec_all_e);
			}
			
			for (auto & match : vec_all_s){
				key = match.first;
				if(!vec_all_e[key].empty()){

					t_start = match.second[1];
					t_end = vec_all_e[key][1];

					context_s = contexts[vec_best_s[0]];
					context_e = contexts[vec_best_e[0]];

					if(t_start < t_end){
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\t" << t_start << "\t" << t_end << "\t" << context_s << "\t" << context_e << endl;
					}
					else if(t_start > t_end){
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\t" << t_end << "\t" << t_start << "\t" << context_e << "\t" << context_s << endl;
					}
					else{
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\tmissing\tmissing\tNA\tNA" << endl;
					}
				}
				else{
					outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\tmissing\tmissing\tNA\tNA" << endl;
				}
			}

			vec_best_s.clear();
			vec_best_e.clear();
			vec_all_s.clear();
			vec_all_e.clear();
		}
		infile.close();
		outfile.close();
	}

	
}
/**********************************************/
void printHELP()
{
	HELP( "Description");
	HELP( "\tsvcit  -- Structural Variant in ctDNA Sequencing Data\n");
	HELP( "Usage:");
	HELP( "\tsvict -i [FILE] -r [FILE]");

	HELP( "Required Parameters:");
	HELP( "\t-i, --input [FILE]\n\t\t Input alignment file. This file should be a sorted SAM or BAM file.\n");
	HELP( "\t-r, --reference [FILE]\n\t\tReference geneome that the reads are algined to." );

	HELP( "\nMain Optional Parameters:");
	HELP( "\t-o, --output [STRING]\n\t\tPrefix of output file (default out)\n");
	HELP( "\t-g, --annotation [FILE]\n\t\tGTF file. Enables annotation of SV calls and fusion identification.\n" );
	HELP( "\t-s, --min_support [INT]\n\t\tThe minimum number of supporting reads required to be considered a SV (default 2).\n");
	HELP( "\t-S, --max_support [INT]\n\t\tThe maximum number of supporting reads required to be considered a SV, could be used to filter out germline calls (default unlimited).\n");
	HELP( "\t-m, --min_length [INT]\n\t\tMin SV length (default 30).\n");
	HELP( "\t-M, --max_length [INT]\n\t\tMax SV length (default 20000).");

	HELP( "\t\nAdditional Parameters:");
	HELP( "\t-h, --help\n\t\tShows help message.\n");
	HELP( "\t-v, --version\n\t\tShows current version.\n");
	HELP( "\t-p, --print_reads\n\t\tPrint all contigs and associated reads as additional output out.vcf.reads, out is the prefix specified in -o|--output, when activated.\n");
	HELP( "\t-P, --print_stats:\n\t\tPrint statistics of detected SV calls and fusions to stderr.\n");
	HELP( "\t-w, --window_size [INT]\n\t\tThe size of the sliding window collecting all reads with soft-clip/split positions in it to form the breakpoint specific cluster for contig assembly. \n\t\tLarger window size could assign a read to more clusters for potentially higher sensitivity with the cost of increased running time (default 3).\n");
	HELP( "\t-d, --min_sc [INT]\n\t\tMinimum soft clip length for a read to be considered as unmapped or incorrectly mapped to be extracted for contig assembly (default 10).\n");
	HELP( "\t-n, --no_indel\n\t\tIgnore indels in the input BAM/SAM (I and D in cigar) when extracting potential breakpoints.\n");
	HELP( "\t-O, --assembler_overlap [INT]\n\t\tThe minimum lenth of overlaps between 2 reads in overlap-layout-consensus contig assembly (default 50).\n");
	HELP( "\t-a, --anchor [INT]\n\t\tAnchor length. Intervals shorter than this value would be discarded in interval chaining procedure for locating contigs (default 30).\n");
	HELP( "\t-k, --kmer [INT]\n\t\tk-mer length to index and remap assembled contigs to reference genome (default 14).\n");
	HELP( "\t-u, --uncertainty [INT]\n\t\tUncertainty around the breakpoint position.\n\t\tSee \"Interval Chaining for Optimal Mapping\" in publication (default 8).\n");
	HELP( "\t-c, --sub_optimal [INT]\n\t\tFor a contig, SViCT will report all paths which are not worse than the optimal by c on the DAGs to achieve potentially better sensitivity. \n\t\tSee \"Interval Chaining for Optimal Mapping\" in publication (default 0 - co-optimals only, negative value disables).\n");
	HELP( "\t-H, --heuristic \n\t\tUse clustering heuristic (good for data with PCR duplicates).\n");
	HELP( "\t-D, --dump_contigs\n\t\tDump contigs in fastq format for mapping.\n");
	HELP( "\t-R, --resume:\n\t\tResume at the interval chaining stage with mapped contigs.\n");

	HELP( "Example:");
	HELP( "\tsvict -i input.bam -r human_genome.fa -o final\n\tThis command will generate prediction result final.vcf directly from input.sam.\n\n");
}
/********************************************************************/
int main(int argc, char *argv[])
{

	int opt;
	int opt_index;

	string	input_sam  = "" ,
			reference  = "" ,
			out_prefix = "out" ,
			annotation = "" ,
			contig_file = "";
	vector<string> chromosomes; chromosomes.reserve(128);

	int k = 14, O = 50, a = 30, s = 2, S = 5000, u = 8, c = 0, w = 3, m = 30, M = 20000, min_sc = 10, max_fragment_size = 200; 
	bool barcodes = false, print_reads = false, print_stats = false, dump_contigs = false, resume = false, indel = true, heuristic = false;
	double clip_ratio = 0.99;
	int ref_flag  = 0;

	static struct option long_opt[] =
	{
		{ "help", no_argument, 0, 'h' },
		{ "version", no_argument, 0, 'v' },
		{ "input", required_argument, 0, 'i' },
		{ "reference", required_argument, 0, 'r' },
		{ "output", required_argument, 0, 'o' },
		{ "annotation", required_argument, 0, 'g' },
		{ "barcodes", no_argument, 0, 'b' },
		{ "print_reads", no_argument, 0, 'p' },
		{ "print_stats", no_argument, 0, 'P' },
		{ "kmer", required_argument, 0, 'k' },
		{ "assembler_overlap", required_argument, 0, 'O' },
		{ "anchor", required_argument, 0, 'a' },
		{ "min_support", required_argument, 0, 's' },
		{ "max_support", required_argument, 0, 'S' },
		{ "uncertainty", required_argument, 0, 'u' },
		{ "sub_optimal", required_argument, 0, 'c' },
		{ "min_length", required_argument, 0, 'm' },
		{ "max_length", required_argument, 0, 'M' },
		{ "window_size", required_argument, 0, 'w' },
		{ "min_sc", required_argument, 0, 'd' },
		{ "no_indel", no_argument, 0, 'n' },
		{ "heuristic", no_argument, 0, 'H'},
		{ "dump_contigs", no_argument, 0, 'D'},
		{ "resume", no_argument, 0, 'R'},
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvi:r:o:g:bpPk:O:a:s:S:u:c:m:M:w:d:nHDR", long_opt, &opt_index )  ) )
	{
		switch(opt)
		{
			case 'h':
				printHELP();
				return 0;
			case 'v':
				fprintf(stdout, "SViCT v%s\n", versionNumber );
				return 0;
			case 'b':
				barcodes = true; 
				break;
			case 'p':
				print_reads = true; 
				break;
			case 'P':
				print_stats = true; 
				break;
			case 'i':
				input_sam.assign( optarg );
				break;
			case 'r':
				reference.assign( optarg );
				ref_flag = 1;
				break;
			case 'o':
				out_prefix.assign( optarg );
				break;
			case 'g':
				annotation.assign( optarg );
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'O':
				O = atoi(optarg);
				break;
			case 'a':
				a = atoi(optarg);
				break;
			case 's':
				s = atoi(optarg);
				break;
			case 'S':
				S = atoi(optarg);
				break;
			case 'u':
				u = atoi(optarg);
				break;
			case 'c':
				c = atoi(optarg);
				break;
			case 'm':
				m = atoi(optarg);
				break;
			case 'M':
				M = atoi(optarg);
				break;
			case 'w':
				w = atoi(optarg);
				break;
			case 'd':
				min_sc = atoi(optarg);
				break;
			case 'n':
				indel = false; 
				break;
			case 'H':
				heuristic = true; 
				break;
			case 'D':
				dump_contigs = true; 
				break;
			case 'R':
				resume = true; 
				break;
			case '?':
				fprintf(stderr, "Unknown parameter: %s\n", long_opt[opt_index].name);
				return 1;
			default:
				printHELP();
				return 0;
		}
	}

	int pass = 1;
	string msg = "";
	// sanity checking
	if( k < 5 ) {
		msg += "\tError: K must be greater than 5\n";
		pass = 0;
		}
	if( a < k ) {
		msg +=  "\tError: K must be less than or equal to anchor length\n";
		pass = 0;
		}
	if( s < 1 ) {
		msg += "\tError: Read support threshold must be an integer >= 1\n";
		pass = 0;
		}
	if( s > S ) {
		msg += "\tError: Min read support should be less than or equal to max read support\n";
		pass = 0;
		}
	if ( O < 0 ){
		msg += "\tError: Assembler overlap must be a positive integer";
		pass = 0;
	}
	if ( u < 0 ){
		msg += "\tError: Uncertainty must be a positive integer";
		pass = 0;
	}
	if ( w < 0 ){
		msg += "\tError: Window size must be a positive integer";
		pass = 0;
	}
	else if(w > 10){
		msg += "\tError: Window size must be <= 10";
		pass = 0;
	}
	if ( min_sc < 0 ){
		msg += "\tError: Minimum soft clip must be a positive integer";
		pass = 0;
	}
	if ( min_sc > m ){
		E("d > m. Adjusting d to equal m.\n");
		min_sc = m;
	}
	if ( m <= u ){
		msg += "\tError: Min SV length should be > uncertainty";
		pass = 0;
	}
	if ( m > M ){
		msg += "\tError: Min SV length should be <= max SV length\n";
		pass =  0;
	}


	if ( !ref_flag ){
		msg += "\tError: Reference Genome (specify by -r ) is required to predict SV\n";
		pass =  0;
	}

	
	if ( !pass )
	{
		E("SViCT does not accept the following parameter values:\n\n%s\n\n", msg.c_str() );
		E("Run ./svict -h to check help message for more information\n\n");
		return 0;
	}

	if(!check_input(input_sam, resume, chromosomes))return 0;

	if(!check_ref(reference))return 0;

	if(annotation != ""){
		if(!check_gtf(annotation))return 0;
	}

	if(dump_contigs || resume){
		contig_file = (out_prefix + ".fastq");
	}

	predict(input_sam, reference, annotation, (out_prefix + ".vcf"), contig_file, print_reads, print_stats,  k, O, a, s, S, u, m, M, c, false, w, min_sc, max_fragment_size, clip_ratio, indel, heuristic, chromosomes, resume);

	return 0;	
}
