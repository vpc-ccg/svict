#ifndef __KMISTRVAR__
#define __KMISTRVAR__
#define GCC_VERSION (__GNUC__ * 1000 \
                     + __GNUC_MINOR__ * 10 )

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
#include "extractor.h"
#include "../include/tsl/sparse_map.h"

using namespace std;

class kmistrvar : public variant_caller
{
private:

	const int MAX_INTERVALS = 100000;
	const short MASK = 6;
	const short MASK_RC = 4;
	const int CON_NUM_DEBUG = -1435;//411//2034;//1357;//8562;
	const int REPEAT_LIMIT = 1;
	const int CON_REPEAT_LIMIT = 100;//20;//2;
	const int ANCHOR_SIZE = 40;
	const bool USE_ANNO = true;
	const bool PRINT_READS = false;
	const bool PRINT_STATS = false;
	const bool USE_BARCODES = false;
	const vector<string> chromos = {"1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","MT","X","Y"};
	const vector<string> contexts = {"intergenic", "intronic", "non-coding", "UTR", "exonic-CDS"};
	const vector<string> sv_types = {"INV", "INS", "DEL", "DUP", "TRANS", "BND"};
	const short INV = 0; const short INS = 1; const short DEL = 2; const short DUP = 3; const short TRANS = 4; const short BND = 5; const short INSL = 6; const short INSR = 7; 

	struct contig_metrics{
		double num_reads;
		double len;
		double max_dist;
	};

	struct mapping{
		long loc;// : 29; 
		int len;// : 14;  // up to 16,383
		int con_loc;// : 14;
#if GCC_VERSION == 4080
		bool operator<( const mapping& rhs) const{  //To address known GCC 4.8 bug
			return this->len > rhs.len;
		}
#else 
		bool operator<( const mapping& rhs){
			return this->len > rhs.len;
		}
#endif
	};

	struct mapping_ext{
		char chr;// : 5;
		bool rc;
		long loc;// : 29;
		int len;// : 14;
		int con_loc;
		int con_id;
		int id;
#if GCC_VERSION == 4080
		bool operator<( const mapping_ext& rhs) const{  //To address known GCC 4.8 bug
			return this->len > rhs.len;
		}
#else 
		bool operator<( const mapping_ext& rhs){
			return this->len > rhs.len;
		}
#endif
	};

	struct sortable_mapping{
		int id;
		long loc;// : 29;

#if GCC_VERSION == 4080
		bool operator<( const sortable_mapping& rhs) const{  //To address known GCC 4.8 bug
			return this->loc < rhs.loc;
		}
#else 
		bool operator<( const sortable_mapping& rhs){
			return this->loc < rhs.loc;
		}
#endif
	};

	struct interval_pair{
		bool visited;
		int id1;
		int id2;
	};

	struct last_interval{
		long loc;				//64
		unsigned int len;		//16
		unsigned short repeat;	//8
	};

	struct result{

		string ref_seq;
		string alt_seq;
		string info;
		bool one_bp;
		long end;
		long clust;
		int con;
		int u_id;
		unsigned short sup;
		long pair_loc;
		char pair_chr;
	};

	vector<contig> all_contigs; // 10% of memory
	vector<compressed_contig> all_compressed_contigs; 
	vector<contig_metrics> all_contig_metrics;
	vector<pair<char,pair<int,int>>> regions;
	vector<vector<mapping>>** contig_mappings;
	vector<last_interval>** last_intervals;
	tsl::sparse_map<long, vector<result>>** results;
	vector<pair<int,char>> cluster_info;
	vector<vector<int>> contig_kmers;
	vector<int> contig_kmer_counts;
	vector<tsl::sparse_map<int, vector<int>>> kmer_locations;
	map <string,vector<isoform>> iso_gene_map;
	map <string,vector<gene_data>> gene_sorted_map;
	unsigned long long kmer_mask;
	unsigned long long num_kmer;
	double rate_param;
	bool** valid_mappings;
	int* contig_kmers_index;  // 35% of memeory at k=14
	int k;
	int num_intervals;
	int u_ids;

bool TP[14000];

public:

	

private:

	void init();
	vector<string> split(string str, char sep = ' ');
	string itoa (int i);
	string con_string(int& id, int start, int len);
	compressed_contig compress(contig& con);
	pair<unsigned short,pair<unsigned short,unsigned short>> compute_support(int& id, int start, int end);
	vector<pair<string, string>> correct_reads(vector<pair<string, string>> reads);
	long add_result(int id, mapping_ext& m1, mapping_ext& m2, short type, char pair_chr, long pair_loc);
	void print_results(FILE* fo_vcf, FILE* fr_vcf, int uncertainty);
	mapping_ext copy_interval(char chr, bool rc, int con_id, mapping& interval);
	void print_interval(string label, mapping_ext& interval);
	void assemble(int min_support, int max_support, int uncertainty, const bool LOCAL_MODE, int min_dist, int max_dist, int max_num_read, double clip_ratio);
	void probalistic_filter();
	void index();
	void generate_intervals(const string &out_vcf, const bool LOCAL_MODE);
	void generate_intervals_bwa(const string &sam_file, const bool LOCAL_MODE);
	void predict_variants(const string &out_vcf, int uncertainty, int min_length, int max_length);
	bool bfs(const int DEPTH, int** rGraph, int s, int t, int parent[]);
	vector<vector<int>> fordFulkerson(const int DEPTH, int** rGraph, int s, int t);
	
public:
	
	kmistrvar(int kmer_len, int anchor_len, const string &input_file, const string &reference, const string &gtf, const bool barcodes, const bool print_reads, const bool print_stats);
	~kmistrvar();
	void run_kmistrvar(const string &out_vcf, int min_support, int max_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE, int min_dist, int max_dist, int max_num_read, double clip_ratio = 0.99);
};
#endif
