#ifndef __EXTRACTOR__
#define __EXTRACTOR__

#include "common.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "record.h"

using namespace std;

class extractor 
{


private:

	Parser* parser;
	//unordered_map<string,pair<string, string>> supply_dict;
	unordered_map<string, Record> map_oea;
	unordered_map<string, Record> map_read;
	unordered_map<string, Record> map_orphan;
	int min_dist;
	int max_dist;
	int max_num_read;
	double clip_ratio;
	bool both_mates;
	bool two_pass;



public:
	unordered_map<string,pair<string, string>> supply_dict;

	struct read{
		string name;
		string seq;
	};
	
	struct sa_read{
		string 	name;
		//char *name;
		int 	flag;
		size_t 	pos;
	};
	struct cluster{
		vector<pair<string, string>> reads;
		vector< sa_read > sa_reads;
		int start;
		int end;
		string ref;
		int sup;
	};


private:
	int md_length( char *md);
	int parse_supply( const char *attr, char *rname, int &loc, int &nm );
	int parse_sc( const char *cigar, int &match_l, int &read_l );
	int get_endpoint( const uint32_t pos, const uint32_t pair_pos, const int match_l, const int tlen, int &t_s, int &t_e );
	int process_orphan( const Record &rc, FILE *forphan, FILE *f_int, int ftype);
	int process_oea( const Record &rc, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, int &min_length);
	int examine_mapping( const Record &rc, FILE *f_map, FILE *f_unmap, FILE *f_int, int ftype, double clip_ratio, int &min_length  );
	int dump_oea( const Record &rc, read &tmp, int &anchor_pos, bool both_mates );
	int dump_mapping( const Record &rc, read &tmp, int &anchor_pos, double clip_ratio, bool both_mates );
	int parse_sa( const char *attr );
	bool has_supply_mapping( const char *attr );
	int scan_supply_mappings( const string filename, int ftype );
	int check_supply_mappings( const Record &rc );
	//int dump_supply( const char *readname, const int flag, const size_t pos, bool both_mates, read &tmp);
	//int dump_supply( const char *readname, const int flag, const size_t pos, bool both_mates, read &tmp, const char *ref_name, const int sa_pos);

public:
	extractor(string filename, int min_dist, int max_dist, int max_num_read, double clip_ratio = 0.99, bool both_mates = false, bool two_pass = true );
	~extractor();
	extractor::cluster get_next_cluster();
	int dump_supply( const char *readname, const int flag, const size_t pos, bool both_mates, read &tmp);
	bool has_next_cluster();
	void clear_maps();
};



#endif
