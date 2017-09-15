#ifndef __KMISTRVAR__
#define __KMISTRVAR__

#include <cstring>
#include <bitset>
#include <limits>
#include "variant_caller.h"
#include "partition.h"
#include "common.h"
#include "assembler.h"
#include "assembler_old.h"
#include "assembler_ext.h"
#include "genome.h"
#include "annotation.h"

using namespace std;

class kmistrvar : public variant_caller
{
private:

	const int MAX_INTERVALS = 100000;
	const int MAX_INTERVAL_LEN = 10000;
	const short MASK = 6;
	const short MASK_RC = 4;
	const int CON_NUM_DEBUG = -10235;
	const int REPEAT_LIMIT1 = 5;   //intra-chromosome
	const int REPEAT_LIMIT2 = 50;  //inter-chromosome
	const int REPEAT_LIMIT3 = 400; //total
	const int CON_REPEAT_LIMIT = 2;
	const int ANCHOR_SIZE = 40;
	const bool USE_ANNO = true;
	const bool PRINT_READS = true;
	const bool PRINT_STATS = false;
	const bool USE_BARCODES = false;

	struct mapping{
		string seq;
		string chr;
		bool rc;
		long loc;
		int len;
		int con_loc;
		int con_id;
		int id;
	};

	vector<contig> all_contigs;
	vector<string> chromos = {"1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","MT","X","Y"};
	vector<string> contexts ={"intergenic", "intronic", "non-coding", "UTR", "exonic-CDS"};
	vector<pair<string,pair<int,int>>> regions;
	vector<vector<mapping>>** contig_mappings;
	vector<bool>** repeat;
	vector<vector<int>> contig_kmers;
	map <string,vector<isoform>> iso_gene_map;
	map <string,vector<gene_data>> gene_sorted_map;
	int* contig_kmers_index; 
	vector<mapping> all_intervals;
	vector<unordered_map<long long, vector<int>>> kmer_locations;
	long long kmer_mask;
	int num_kmer;
	int k;
	int num_intervals;

public:

	

private:

	void init();
	vector<string> split(string str, char sep = ' ');
	string itoa (int i);
	pair<int,pair<int,int>> compute_support(contig contig, int start, int end);
	vector<pair<pair<string, string>, int>> correct_reads(vector<pair<pair<string, string>, int>> reads);
	void print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id, mapping m1, mapping m2, string type);
	void print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id1, int id2, mapping m1, mapping m2, mapping m3, mapping m4, string type);
	void print_interval(string label, mapping& interval);
	void assemble(const string &range, int min_support, const bool LEGACY_ASSEMBLER, const bool LOCAL_MODE, int ref_flank);
	void generate_intervals(const bool LOCAL_MODE);
	void predict_variants(const string &out_vcf, const string &out_full, int uncertainty, int min_length, int max_length);
	bool bfs(const int DEPTH, int** rGraph, int s, int t, int parent[]);
	vector<vector<int>> fordFulkerson(const int DEPTH, int** rGraph, int s, int t);
	
public:
	
	kmistrvar(int kmer_len, int anchor_len, const string &partition_file, const string &reference, const string &gtf, const bool barcodes);
	~kmistrvar();
	void run_kmistrvar(const string &range, const string &out_vcf, const string &out_full, int min_support, int uncertainty, int min_length, int max_length, const bool LEGACY_ASSEMBLER, const bool LOCAL_MODE, int ref_flank);
};
#endif
