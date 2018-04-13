#ifndef __EXTRACTOR__
#define __EXTRACTOR__

#include <deque>
#include <set>
#include "common.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "record.h"

using namespace std;

class extractor 
{

private:

	struct breakpoint{
		int sc_loc;
		short len;
		short pos;
	};

	struct read{
		string name;
		string seq;
	};

	struct sa_read{
		string name;
		uint32_t flag;
	};

	struct sortable_read{
		string name;
		string seq;
		breakpoint bp;
		uint32_t flag;
		bool supple;

		bool operator<( const sortable_read& rhs) const{
			return (this->bp.sc_loc == rhs.bp.sc_loc) ? ((this->bp.pos == rhs.bp.pos) ?  this->seq < rhs.seq : this->bp.pos < rhs.bp.pos) : this->bp.sc_loc < rhs.bp.sc_loc;
		}
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
	set<sortable_read> sorted_soft_clips;
	vector<sortable_read> indexed_soft_clips;
	vector<cluster> supple_clust;
	deque<sortable_read> local_reads; 
	cluster orphan_clust;
	string cur_ref;
	int min_sc;
	int max_dist;
	int max_num_read;
	double clip_ratio;
	int cur_pos = -1;
	int cur_type = -1;
	int skip_pos = -1;
	int skip_count = 0;
	int index = 0;
	const bool use_indel = true;
	const int INSERT_SIZE = 600;
	const short DLEFT = 0; const short LEFT = 1; const short BOTH = 2; const short RIGHT = 3; const short DRIGHT = 4;

private:
	int parse_sc( const char *cigar, int &match_l, int &read_l );
	vector<breakpoint> extract_bp(string& cigar, int& mapped, int sc_loc, bool use_indel);
	int dump_oea( const Record &rc, read &tmp, vector<breakpoint> &bps, double clip_ratio );
	int dump_mapping( const Record &rc, read &tmp, vector<breakpoint> &bps, double clip_ratio );
	bool dump_supply( const string& readname, const int flag, read &tmp);
	void extract_reads();

public:
	extractor(string filename, int min_sc, int max_dist, int max_num_read, double clip_ratio = 0.99);
	~extractor();
	extractor::cluster& get_next_cluster(int uncertainty, int min_support, bool heuristic);
	bool has_next_cluster();
	void clear_maps();
};



#endif
