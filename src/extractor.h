#ifndef __EXTRACTOR__
#define __EXTRACTOR__

#include <deque>
#include "common.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "record.h"

using namespace std;

class extractor 
{

private:

	struct read{
		string name;
		string seq;
	};

	struct sa_read{
		string readname;
		uint32_t flag;
	};

	struct sortable_read{
		string name;
		string seq;
		uint32_t flag;
		int sc_loc;

#if GCC_VERSION == 4080
		bool operator<( const sortable_read& rhs) const{  //To address known GCC 4.8 bug
			return this->sc_loc < rhs.sc_loc;
		}
#else 
		bool operator<( const sortable_read& rhs){
			return this->sc_loc < rhs.sc_loc;
		}
#endif
	};

public:

	struct cluster{
		vector<pair<string, string>> reads;
		vector<sa_read> sa_reads;
		int start;
		int end;
		string ref;
	};

private:

	Parser* parser;
	unordered_map<string,pair<string, string>> supply_dict;
	unordered_map<string, Record> map_oea;
	unordered_map<string, Record> map_read;
	unordered_map<string, Record> map_orphan;
	vector<sortable_read> sorted_soft_clips;
	vector<cluster> supple_clust;
	deque<sortable_read> local_reads; 
	cluster orphan_clust;
	string cur_ref;
	int min_sc;
	int max_dist;
	int max_num_read;
	double clip_ratio;
	int index = 0;
	const bool use_indel = true;
	const int INSERT_SIZE = 600;

private:
	int parse_sc( const char *cigar, int &match_l, int &read_l );
	bool has_supply_mapping( const char *attr );
	vector<pair<int, int>> extract_bp(string& cigar, int& mapped, int sc_loc, bool use_indel);
	int dump_oea( const Record &rc, read &tmp, vector<pair<int, int>> &bps, double clip_ratio );
	int dump_mapping( const Record &rc, read &tmp, vector<pair<int, int>> &bps, double clip_ratio );
	bool dump_supply( const string& readname, const int flag, read &tmp);
	void extract_reads();

public:
	extractor(string filename, int min_sc, int max_dist, int max_num_read, double clip_ratio = 0.99);
	~extractor();
	extractor::cluster& get_next_cluster_heuristic(int uncertainty, int min_support);
	extractor::cluster& get_next_cluster(int uncertainty, int min_support);
	bool has_next_cluster();
	void clear_maps();
};



#endif
