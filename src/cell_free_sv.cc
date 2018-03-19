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
#include "assembler_old.h"
#include "assembler_ext.h"
#include "kmistrvar.h"
#include "genome.h"
#include "simulator.h"
#include "common.h"
#include "record.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "extractor.h"
#include "sort.h"

using namespace std;

char versionNumber[10]="1.0.0";
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

/********************************************************************/
void partify (const string &read_file, const string &mate_file, const string &out, int threshold) 
{
	gzFile gzf = gzopen(mate_file.c_str(), "rb");
	unordered_map<string, string> mymap;

	const int MAXB = 8096;
	char buffer[MAXB];
	char name[MAXB], read[MAXB];
	while (gzgets(gzf, buffer, MAXB)) {
		sscanf(buffer, "%s %s", name, read);
		string read_name = string(name+1); 
		if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/')
			read_name = read_name.substr(0, read_name.size() - 2);
		mymap[read_name] = string(read);
	}

	genome_partition pt(read_file, threshold, mymap); 

	FILE *fo = fopen(out.c_str(), "wb");
	FILE *fidx = fopen((out + ".idx").c_str(), "wb");
	while (pt.has_next()) {
		auto p = pt.get_next();
		size_t i = pt.dump(p, fo);
		fwrite(&i, 1, sizeof(size_t), fidx);
	}
	fclose(fo);
	fclose(fidx);
	printf("%d\n", pt.get_cluster_id());
}

/********************************************************************/
void predict (const string &partition_file, const string &reference, const string &gtf, const bool barcodes, const bool print_reads, const bool print_stats, const string &out_vcf,
					int k, int anchor_len, int min_support, int max_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE,
					int min_dist, int max_dist, int max_num_read, double clip_ratio)
{
	kmistrvar predictor(k, anchor_len, partition_file, reference, gtf, barcodes, print_reads, print_stats);
	predictor.run_kmistrvar(out_vcf, min_support, max_support, uncertainty, min_length, max_length, LOCAL_MODE, min_dist, max_dist, max_num_read, clip_ratio);
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
	LOG( "SVICT: Structural Variant In CTDNA Sequencing Data\n");
	LOG( "\t-h|--help:\tShows help message.");
	LOG( "\t-v|--version:\tShows current version.");
	LOG( "\t\nMandatory Parameters:");
	LOG( "\t-i|--input:\tInput file. (SAM/BAM)");
	LOG( "\t-o|--output:\tPrefix or output file");
	LOG( "\t\nParameters for Supplementary Information:");
	LOG( "\t-r|--reference:\tReference Genome. Required for SV detection." );
	LOG( "\t-g|--annotation:\tGTF file for gene annotation." );
	LOG( "\t\nOptional Parameters:");
	LOG( "\t-b|--barcode:\tInput reads contain barcodes.");
	LOG( "\t-p|--print_reads:\tPrint all contigs and associated reads as additional output.");
	LOG( "\t-P|--print_stats:\tPrint statistics as additional output.");
	LOG( "\t-c|--cluster:\tClustering threshold (default 1000).");
	LOG( "\t-k|--kmer:\tKmer length (default 14).");
	LOG( "\t-a|--anchor:\tAnchor length (default 40).");
	LOG( "\t-s|--min_support:\tMin Read Support (default 2).");
	LOG( "\t-S|--max_support:\tMax Read Support (default unlimited).");
	LOG( "\t-u|--uncertainty:\tUncertainty (default 8).");
	LOG( "\t-m|--min_length:\tMin SV length (default 60).");
	LOG( "\t-M|--max_length:\tMax SV length (default 20000).");
	LOG( "\t-d|--min_sc:\tMin soft clip to consider (default 5).");
	LOG( "\t-D|--max_dist:\tMax cluster distance (default 1000).");
	LOG( "\t-n|--max_reads:\tMax number of reads allowed in a cluster.\n\t\tA region with more than the number of OEA/clipped reads will not be considered in prediction. (default 200).\n");


	LOG( "\t\nExample Command:");
	LOG( "\tRunning SVICT from a SAM/BAM file can be done in single command:");
	LOG( "\t./SVICT -i input.sam -r human_genome.fa -o final\n\t\tThis command will generate prediction result final.vcf directly from input.sam.\n\n");
}
/********************************************************************/
int main(int argc, char *argv[])
{

	int opt;
	int opt_index;

	string	input_sam  = "" ,
			reference  = "" ,
			out_prefix = "out" ,
			annotation = "" ;

	int k = 14, a = 40, s = 2, S = 999999, u = 8, m = 60, M = 20000, min_sc = 10, max_dist = 1000, max_reads = 200; 
	bool barcodes = false, print_reads = false, print_stats = false;
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
		{ "anchor", required_argument, 0, 'a' },
		{ "min_support", required_argument, 0, 's' },
		{ "max_support", required_argument, 0, 'S' },
		{ "uncertainty", required_argument, 0, 'u' },
		{ "min_length", required_argument, 0, 'm' },
		{ "max_length", required_argument, 0, 'M' },
		{ "min_sc", required_argument, 0, 'd' },
		{ "max_dist", required_argument, 0, 'D' },
		{ "max_reads", required_argument, 0, 'n' },
	//	{ "both-mate", no_argument, 0, 'B' },
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvi:r:o:g:bpPk:a:s:S:u:m:M:d:D:n:", long_opt, &opt_index )  ) )
	{
		switch(opt)
		{
			case 'h':
				printHELP();
				return 0;
			case 'v':
				fprintf(stdout, "%s\n", versionNumber );
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
			case 'm':
				m = atoi(optarg);
				break;
			case 'M':
				M = atoi(optarg);
				break;
			case 'd':
				min_sc = atoi(optarg);
				break;
			case 'D':
				max_dist = atoi(optarg);
				break;
			case 'n':
				max_reads = atoi(optarg);
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
	if( max_reads <= 1 ){
		msg += "\tError: Max reads in a cluster must be larger than 1\n";
		pass = 0;
		}
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
	if ( u < 0 ){
		msg += "\tError:Uncertainty must be a positive integer";
		pass = 0;
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
		E("SVICT does not accept the following parameter values:\n\n%s\n\n", msg.c_str() );
		E("Check help message for more information\n\n");
		printHELP();
		return 0;
	}

	predict(input_sam, reference, annotation, barcodes, print_reads, print_stats, (out_prefix + ".vcf"), k, a, s, S, u, m, M, 0, min_sc, max_dist, max_reads, 0.99);

	return 0;
}
