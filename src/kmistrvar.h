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

	const int MAX_READS_PER_PART = 10000;
	const int MAX_INTERVALS = 100000;
	//const int MAX_INTERVAL_LEN = 10000;
	const short MASK = 6;
	const short MASK_RC = 4;
	const int CON_NUM_DEBUG = -10235;
	const int REPEAT_LIMIT1 = 5;   //distant
	const int REPEAT_LIMIT2 = 50;  //close
	const int REPEAT_LIMIT3 = 400; //total 400
	const int CON_REPEAT_LIMIT = 2;
	const int ANCHOR_SIZE = 40;
	const bool USE_ANNO = true;
	const bool PRINT_READS = true;
	const bool PRINT_STATS = false;
	const bool USE_BARCODES = false;


	struct mapping{
		long loc;// : 29; 
		int len;// : 14;  // up to 16,383
		int con_loc;// : 14;
	};

	struct mapping_ext{
		char chr;// : 5;
		bool rc;
		long loc;// : 29;
		int len;// : 14;
		int con_loc;
		int con_id;
		int id;
	};

	struct sortable_mapping{
		int id;
		long loc;// : 29;

		bool operator<( const sortable_mapping& rhs){
			return this->loc < rhs.loc;
		}
	};

	struct interval_pair{
		bool visited;
		int id1;
		int id2;
	};

	vector<contig> all_contigs;
	vector<string> chromos = {"1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","MT","X","Y"};
	vector<string> contexts ={"intergenic", "intronic", "non-coding", "UTR", "exonic-CDS"};
	vector<pair<char,pair<int,int>>> regions;
	vector<vector<mapping>>** contig_mappings;
	vector<short>** repeat;
	vector<vector<int>> contig_kmers;
	vector<int> contig_kmer_counts;
	map <string,vector<isoform>> iso_gene_map;
	map <string,vector<gene_data>> gene_sorted_map;
	int* contig_kmers_index; 
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
	void print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id, mapping_ext m1, mapping_ext m2, string type);
	void print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id1, int id2, mapping_ext m1, mapping_ext m2, mapping_ext m3, mapping_ext m4, string type);
	mapping_ext copy_interval(char chr, bool rc, int con_id, mapping& interval);
	void print_interval(string label, mapping_ext& interval);
	void assemble(const string &range, int min_support, const bool LOCAL_MODE, int ref_flank);
	void generate_intervals(const bool LOCAL_MODE);
	void predict_variants(const string &out_vcf, const string &out_full, int uncertainty, int min_length, int max_length);
	bool bfs(const int DEPTH, int** rGraph, int s, int t, int parent[]);
	vector<vector<int>> fordFulkerson(const int DEPTH, int** rGraph, int s, int t);
	
public:
	
	kmistrvar(int kmer_len, int anchor_len, const string &partition_file, const string &reference, const string &gtf, const bool barcodes);
	~kmistrvar();
	void run_kmistrvar(const string &range, const string &out_vcf, const string &out_full, int min_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE, int ref_flank);
};
#endif
