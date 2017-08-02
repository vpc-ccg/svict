#ifndef __ANNOTATION__
#define __ANNOTATION__
#include <map>
#include <vector>
#include <string>
#include "common.h"
using namespace std;

#define L(c,...) fprintf(stdout,c,##__VA_ARGS__)
#define E(c,...) fprintf(stderr,c,##__VA_ARGS__)

const int MAX_LINE_ANNO = 200000;
const int TOKEN_LENGTH_ANNO = 20000;
const bool NEW_GTF = true;

typedef struct
{
	int start;
	int end;
}Region;


typedef struct
{
	uint32_t s;
	uint32_t e;
	int f;
}CDS;

typedef struct
{
	uint32_t e_s;
	uint32_t e_e;
	uint32_t c_s;
	uint32_t c_e;
	int f;
	int cds; // 0 for UTR/non-protein coding gene, 2 for complete CDS, 1 for state-changing exon
}Exon;

typedef struct
{
	string id;
	string gid;
	string ref;
	string tname;
	string gname;
	string src;
	string fea;
	char strand;
	vector<Exon> exon;
	int cds_iso;
	int len;
}isoform;

typedef struct
{
	uint32_t start;
	uint32_t end;
	char strand;
	string chr;
	string gene_id;
	string gene_name;
	int cds_gene;
}gene_data;

bool comp_isoform_start(const isoform &is_1, const isoform &is_2);
bool comp_exon( const Exon &exon1, const Exon &exon2); 
bool comp_region( const Region &r1, const Region &r2); 
bool comp_gene_in_contig( const gene_data &g1, const gene_data &g2); 
uint32_t overlap_l( const uint32_t &s1, const uint32_t &e1, const uint32_t &s2, const uint32_t &e2);

void adjust_gene( map<string, gene_data> &map_gene, const string g_id, char *gene_name, const isoform &new_iso);

void ensembl_Reader(const char *gtf_file, map<string, vector<isoform> > &iso_gene_map, map<string, vector<gene_data> > &gene_sorted_map);
void bed_reader( char *bed_file, map< string, vector<Region> > &interval_map);
void locate_in_isoform( uint32_t s, uint32_t e, const isoform &iso, vector<uint32_t> &t_vec, bool genomic);
int  locate_interval( const string &ref, uint32_t s, uint32_t e, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map, string &best_gene, string &best_trans , vector<uint32_t> &vec_best);
int  locate_interval( const string &ref, uint32_t s, uint32_t e, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map, string &best_gene, string &best_trans , vector<uint32_t> &vec_best, map <string, vector<uint32_t>> &vec_all);
int  locate_interval( const string &ref, uint32_t s, uint32_t e, const string &gene_id, const string &trans_id, const vector<gene_data> &gene_vector, int pos, const map<string, vector< isoform> > &iso_gene_map, string &best_gene, string &best_trans , vector<uint32_t> &vec_best, map <string, vector<uint32_t>> &vec_all);
#endif
