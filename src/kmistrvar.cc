#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cfenv>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include "kmistrvar.h"
#pragma STDC FENV_ACCESS ON

using namespace std;
#define strequ(A,B) (strcmp(A.c_str(), B.c_str()) == 0)
kmistrvar::kmistrvar(int kmer_len, const int anchor_len, const string &input_file, const string &reference, const string &gtf, const bool barcodes, const bool print_reads, const bool print_stats) : 
	variant_caller(input_file, reference), ANCHOR_SIZE(anchor_len), USE_ANNO((gtf != "")), USE_BARCODES(barcodes), PRINT_READS(print_reads), PRINT_STATS(print_stats) {

	k = kmer_len;
	num_kmer = (unsigned long long)pow(4,k);
	kmer_mask = (unsigned long long)(num_kmer-1);
	num_intervals = 1;
	u_ids = 0;
	
	if(USE_ANNO) ensembl_Reader(gtf.c_str(), iso_gene_map, gene_sorted_map );

	contig_mappings = new vector<vector<mapping>>*[2];
	for(int rc=0; rc < 2; rc++){
		contig_mappings[rc] = new vector<vector<mapping>>[25];
	}

	last_intervals = new vector<last_interval>*[2];
	for(int rc=0; rc < 2; rc++){
		last_intervals[rc] = new vector<last_interval>[25];
	}

	results = new unordered_map<long, vector<result>>*[6];
	for(int sv=0; sv < sv_types.size(); sv++){
		results[sv] = new unordered_map<long, vector<result>>[25];
	}

	init();
}

kmistrvar::~kmistrvar(){

	delete[] contig_kmers_index;

	for(int rc=0; rc < 2; rc++){
		delete[] contig_mappings[rc];
	}
	delete[] contig_mappings;

	for(int rc=0; rc < 2; rc++){
		delete[] last_intervals[rc];
	}
	delete[] last_intervals;

	for(int sv=0; sv < sv_types.size(); sv++){
		delete[] results[sv];
	}
	delete[] results;

}

void kmistrvar::init(){

	for(int sv=0; sv < sv_types.size(); sv++){
		for(int chr=0; chr < chromos.size(); chr++){
			results[sv][chr] = unordered_map<long, vector<result>>();
		}
	}
}


//================================================================
// Helper functions						   
//================================================================

vector<string> kmistrvar::split(string str, char sep)
{
	vector<string> ret;

	istringstream stm(str);
	string token;
	while(getline(stm, token, sep)) ret.push_back(token);

	return ret;
}

inline string kmistrvar::itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}

inline double factorial(double n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline string kmistrvar::con_string(int& id, int start, int len)
{

	string out = "";

	if(PRINT_READS){

		if(all_contigs[id].data.empty())return out;
		out = all_contigs[id].data.substr(start, len);
	}
	else{

		vector<bool>& data = all_compressed_contigs[id].data;
		if(data.empty())return out;
		int adj_start = start*2;
		int adj_len = len*2;

		for(int i = adj_start; i < adj_start+adj_len; i+=2){
			bool bit1 = data[i];
			bool bit2 = data[i+1];

			if(bit1){
				if(bit2){
					out += 'G';
				}
				else{
					out += 'T';
				}
			}
			else{
				if(bit2){
					out += 'C';
				}
				else{
					out += 'A';
				}
			}
		}
	}

	return out;
}

compressed_contig kmistrvar::compress(contig& con){

	compressed_contig comp_con;
	vector<bool> data;
	vector<unsigned short> coverage(con.data.length());

	data.reserve(con.data.length()*2);

	for(char& n : con.data){
		switch (n) {
			case 'A': data.push_back(0);
					  data.push_back(0);
					  break;
			case 'T': data.push_back(1);
					  data.push_back(0);
					  break;
			case 'C': data.push_back(0);
					  data.push_back(1);
					  break;
			default:  data.push_back(1);  //Ns turn to G because I don't care! 
					  data.push_back(1);
		}
	}

	for (int j = 0; j < con.support(); j++){
		for (int k=0; k < con.read_information[j].seq.length(); k++){
			coverage[k+con.read_information[j].location_in_contig]++;
		}
	}

	comp_con.data = data;
	comp_con.coverage = coverage;

	return comp_con;
}

pair<unsigned short,pair<unsigned short,unsigned short>> kmistrvar::compute_support(int& id, int start, int end){

	unsigned short len = PRINT_READS ? (unsigned short)all_contigs[id].data.length() : (unsigned short)all_compressed_contigs[id].data.size()/2;
	end = min(end, len-1);
	unsigned short max = 0;
	unsigned short min = 0;
	unsigned short sum = 0;
	pair<unsigned short,pair<unsigned short,unsigned short>> result;
	vector<unsigned short> coverage;

	if(PRINT_READS){

		contig& con = all_contigs[id];

		for (int k=0; k < con.data.length(); k++){
			coverage.push_back((unsigned short)0);
		}

		for (int j = 0; j < con.support(); j++){
			for (int k=0; k < con.read_information[j].seq.length(); k++){
				coverage[k+con.read_information[j].location_in_contig]++;
			}
		}
	}
	else{
		coverage = all_compressed_contigs[id].coverage;
	}

	for (int k = start; k <= end; k++)
	{
		if (coverage[k]>max){
			max = coverage[k];
		}
		if (coverage[k]<min){
			min = coverage[k];
		}
		sum += coverage[k];
	}

	result = {sum/len,{min, max}};

	return result;
}

vector<pair<string, string>> kmistrvar::correct_reads(vector<pair<string, string>> reads){

	unordered_map<string, unordered_map<string, int>> consensus;
	vector<pair<string, string>> corrected_reads;
	int max_count = 0;
	int len = 0;
	string c = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";  //quick and dirty fix TODO  CCAGCGCCT-TATCGTGCA_
	string g = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
	string max_seq = "";

	for(auto & read : reads){
		len = read.second.length() < 50 ? read.second.length() : 50;
		if(read.second.substr(0,len) == c.substr(0,len) ||
			read.second.substr(read.second.length()-len,len) == c.substr(0,len) ||
			read.second.substr(0,len) == g.substr(0,len) ||
			read.second.substr(read.second.length()-len,len) == g.substr(0,len))continue; 
		consensus[read.first.substr((read.first.length()-20), 20)][read.second]++;
	}

	for(auto & barcode : consensus){
		
		max_count = 0;
		max_seq = "";
		
		for(auto & seq : barcode.second){
			if(seq.second > max_count){
				max_count = seq.second;
				max_seq = seq.first;
			}
		}

		corrected_reads.push_back({barcode.first, max_seq});
	}

	return corrected_reads;
}

kmistrvar::mapping_ext kmistrvar::copy_interval(char chr, bool rc, int con_id, mapping& interval){

	mapping_ext mapping_copy;
	mapping_copy.chr = chr;
	mapping_copy.rc = rc;
	mapping_copy.loc = interval.loc;
	mapping_copy.len = interval.len;
	mapping_copy.con_loc = interval.con_loc;
	mapping_copy.con_id = con_id;
	mapping_copy.id = -1;

	return mapping_copy;
}

void kmistrvar::print_interval(string label, mapping_ext& interval){
	cerr << label << "\tRef_Loc: " << interval.loc << "-" << (interval.loc+interval.len) << "\tCon_Loc: " << ((interval.con_loc+k)-interval.len) << "-" << (interval.con_loc+k) 
		 << "\tLen:" << interval.len << "\tChr: " << chromos[interval.chr] << "\tRC: " << interval.rc << "\tCon_ID: " << interval.con_id << endl; 
}

long kmistrvar::add_result(int id, mapping_ext& m1, mapping_ext& m2, short type, char pair_chr, long pair_loc){

	long loc = m1.rc ? m1.loc : (m1.loc + m1.len);
	long end = m2.rc ? (m2.loc + m2.len) : m2.loc;
	int contig_dist = (m2.con_loc-m2.len)-m1.con_loc; 
	int contig_len = PRINT_READS ? all_contigs[id].data.length() : all_compressed_contigs[id].data.size()/2;
	int contig_sup;
	int contig_loc = cluster_info[id].first;
	char contig_chr = cluster_info[id].second;
	char strand = '+';
	result new_result;
	new_result.info = "";
	new_result.one_bp = false;
	new_result.u_id = u_ids++;

	//TO think about:
	//end location for 2 contig
	//alt seq for 2 contig

	if(m1.rc && m2.rc)strand = '-';

	if(type == INS){
		new_result.ref_seq = con_string(id, m1.con_loc+k-1, 1);
		new_result.alt_seq = con_string(id, m1.con_loc+k-1, contig_dist+1);
		if(new_result.alt_seq == "")new_result.alt_seq = "<INS>"; //two contig cases where contigdist <= uncertainty
	}
	else if(type == INSL){
		type = INS;
		new_result.ref_seq = "N";
		new_result.alt_seq = con_string(id, 0, (m1.con_loc+k-m1.len));
		new_result.info = ":L";
		new_result.one_bp = true;
		loc = end;
	}
	else if(type == INSR){
		type = INS;
		new_result.ref_seq = con_string(id, m1.con_loc+k-1, 1);
		new_result.alt_seq = con_string(id, m1.con_loc+k-1, (contig_len-(m1.con_loc+k-1)+1));
		new_result.info = ":R";
		new_result.one_bp = true;
		end = loc;
	}
	else if(type == DEL){
		new_result.ref_seq = "<DEL>";
		new_result.alt_seq = con_string(id, m1.con_loc+k-1, 1);
	}
	else if(type == DUP){

		new_result.info = ":INTERSPERSED";
		new_result.ref_seq = "N";
		new_result.alt_seq = "<DUP>";

		if(strand == '+'){
			if(m1.loc < m2.loc && (m1.loc + m1.len)-m2.loc > MIN_SV_LEN_DEFAULT){
				new_result.info = ":TANDEM";
				loc = m2.loc;
				end = (m1.loc + m1.len);
				int len = ((m1.loc + m1.len)-m2.loc)+1;

				if(m1.con_loc+k-1-len >= 0){
					new_result.ref_seq = con_string(id, m1.con_loc+k-1-len, 1);
					new_result.alt_seq = con_string(id, m1.con_loc+k-1-len, len+1);
				}
				else{
					new_result.alt_seq = con_string(id,m1.con_loc+k-len, len);
				}
			}
			else if(m1.loc > m2.loc && (m1.loc + m1.len) > (m2.loc + m2.len)){
				new_result.info = ":TANDEM";
				new_result.ref_seq = con_string(id,m1.con_loc+k-1, 1);
				new_result.one_bp = true;
			}
			else{
				new_result.info = ":INTERSPERSED"; //TODO we have the alt sequence
				new_result.ref_seq = con_string(id,m1.con_loc+k-1, 1);
				new_result.one_bp = true;
			}
		}
		else{
			new_result.info = ""; //negative
		}
	}
	else if(type == INV){

		new_result.ref_seq = con_string(id,m1.con_loc+k-1, 1);
		new_result.alt_seq = "<INV>";
		new_result.info = ":LONG";
		new_result.one_bp = true;

		if(m1.rc){
			new_result.ref_seq = "N";
		}
		else if(!m2.rc){
			new_result.alt_seq = con_string(id,m1.con_loc+k-1, contig_dist+1);
			new_result.info = ":MICRO";
		}
	}
	else if(type == TRANS){

		bool is_anchor1 = ((m1.chr == contig_chr) && (abs(m1.loc - contig_loc) < MAX_ASSEMBLY_RANGE));
		bool is_anchor2 = ((m2.chr == contig_chr) && (abs(m2.loc - contig_loc) < MAX_ASSEMBLY_RANGE));

		if(!is_anchor1){
			new_result.ref_seq = con_string(id,m2.con_loc+k-1-m2.len, 1);;
			new_result.alt_seq = "<TRANS>";
			new_result.info = "right_bp";
			new_result.one_bp = true;
		}
		else if(!is_anchor2){
			new_result.ref_seq = con_string(id,m1.con_loc+k-1, 1);
			new_result.alt_seq = "<TRANS>";
			new_result.info = "left_bp";
			new_result.one_bp = true;
		}
		else{
			new_result.info = "small";
			new_result.ref_seq = con_string(id,m1.con_loc+k-1, 1);
			new_result.alt_seq = con_string(id,m1.con_loc+k-1, contig_dist+1);
		}
	}

	if(loc > end){
		if(type == INV || type == DEL || type == DUP){
			long tmp = loc;
			loc = end;
			end = tmp;
		}
	}

	//quick and dirty fix
	if(type == INS){

		int g_count = 0;

		for(int i = 0; i < new_result.alt_seq.length(); i++){
			if(new_result.alt_seq[i] == 'G'){
				g_count++;
			}
			else{
				g_count = 0; 
			}

			if(g_count > ANCHOR_SIZE){
				return 0;
			}
		}
	}

	pair<int,pair<int,int>> support = compute_support(id, min(m1.con_loc, m2.con_loc-m2.len+k), max(m1.con_loc, m2.con_loc-m2.len+k));
	contig_sup = support.second.second;

	new_result.end = end;
	new_result.clust = contig_loc;
	new_result.con = id;
	new_result.sup = contig_sup;
	new_result.pair_loc = pair_loc;
	new_result.pair_chr = pair_chr;

	results[type][m1.chr][loc].push_back(new_result);

	return loc;
}

void kmistrvar::print_results(FILE* fo_vcf, FILE* fr_vcf, int uncertainty){

	long loc, ploc, start, end, pend;
	char pchr;
	string best_gene1 = "", best_name1 = "", best_trans1 = "", context1 = "";
	string best_gene2 = "", best_name2 = "", best_trans2 = "", context2 = "";
	string best_gene3 = "", best_name3 = "", best_trans3 = "", context3 = "";
	string type, info;
	string cluster_ids, cluster_ids_pair;
	string contig_ids, contig_ids_pair;
	string filter = "PASS";
	string anno = "";
	vector<uint32_t> vec_best1, vec_best2, vec_best3;
	result c, cp; //consensus

	// DEL signal for interspersed duplications
	for(char chr = 0; chr < chromos.size(); chr++){
		if(!results[DEL][chr].empty()){
			for(auto& rp1 : results[DEL][chr]){
				for(auto& rp2 : results[DUP][chr]){
					if(abs(rp1.first - rp2.first) < uncertainty || abs(rp1.first - rp2.second[0].end) < uncertainty || abs(rp1.second[0].end - rp2.first) < uncertainty || abs(rp1.second[0].end - rp2.second[0].end) < uncertainty){ // 2bp would be safer
						for (auto& r : rp1.second){
							r.end = rp2.second[0].end; // Not great
							r.ref_seq = "N";
 							r.alt_seq = "N";
							rp2.second.push_back(r);
						}
						rp1.second.clear();
						break;
					}
				}
			}
		}
	}

	// DEL signal for translocations
	// for(char chr = 0; chr < chromos.size(); chr++){
	// 	if(!results[DEL][chr].empty()){
	// 		for(auto& rp1 : results[DEL][chr]){
	// 			if(rp1.second.empty())continue;
	// 			for(char chr2 = 0; chr2 < chromos.size(); chr2++){
	// 				for(auto& rp2 : results[TRANS][chr2]){

	// 					if(abs(rp2.second[0].sup - rp1.second[0].sup) < 50){
	// 						if(chr == chr2 && abs(rp2.first - rp1.second[0].end) < uncertainty ){
	//  							for (auto& r : rp1.second){

	//  								r.ref_seq = "N";
	//  								r.alt_seq = "N";

	// 	 							if(rp2.second[0].pair_loc == -1){
	// 									r.end = rp2.second[0].end; // Not great
	// 									r.pair_loc = -1;
	// 									rp2.second.push_back(r);
	// 								}
	// 								else if (rp2.second[0].pair_loc == 0){
	// 									rp2.second[0].pair_loc = -1;
	// 									r.pair_loc = rp2.first;
	// 									r.end = rp1.first;
	// 									results[TRANS][rp2.second[0].pair_chr][rp2.second[0].end-1].push_back(r);
	// 								}

	// 								for (auto& r2 : rp2.second){
	// 									r2.pair_loc = -1;
	// 								}
	// 							}
								
	// 							rp1.second.clear();
	// 							break;
	// 						}
	// 						else if(rp2.second[0].pair_chr == chr && ( abs(rp2.second[0].end - rp1.first) < uncertainty || abs(rp2.second[0].pair_loc - rp1.second[0].end) < uncertainty ) ){
	// 							for (auto& r : rp1.second){

	//  								r.ref_seq = "N";
	//  								r.alt_seq = "N";

	// 	 							if(rp2.second[0].pair_loc){
	// 									r.end = rp2.second[0].end; // Not great
	// 									r.pair_loc = rp2.second[0].pair_loc;
	// 									rp2.second.push_back(r);
	// 								}
	// 								else {
	// 									rp2.second[0].pair_loc = rp1.second[0].end;
	// 									r.pair_loc = -1;
	// 									r.pair_chr = chr2;
	// 									r.end = rp2.first+1;
	// 									results[TRANS][rp2.second[0].pair_chr][rp1.second[0].end].push_back(r);
	// 								}

	// 								for (auto& r2 : rp2.second){
	// 									r2.pair_loc = rp1.second[0].end;
	// 								}
	// 							}


	// 							rp1.second.clear();
	// 							break;
	// 						}
	// 					}
	// 				}
	// 				if(rp1.second.empty())break;
	// 			}
	// 		}
	// 	}
	// }

//STATS
int bnd_count = 0;
int not_bnd_count = 0;

	for(int sv = 0; sv < sv_types.size(); sv++){
		for(char chr = 0; chr < chromos.size(); chr++){
			if(results[sv][chr].empty())continue;
			for(auto& result_pair : results[sv][chr]){

				loc = result_pair.first;
				unordered_set<long> clusters;
				unordered_set<paired_result, paired_result_hash> pairs; 
				unordered_set<int> contigs;
				filter = "PASS";
				c.sup = 0;
				c.end = 0;

				if(result_pair.second.empty())continue;

				for(auto& r : result_pair.second){
					if(r.end != c.end && r.sup > c.sup){
						c = r;
						pairs.clear();
						pairs.insert({r.pair_loc, r.end, r.pair_chr});

						clusters.clear();
						clusters.insert(r.clust);

						contigs.clear();
						contigs.insert(r.con);
					}
					else if(r.end == c.end || r.pair_loc > 0 || r.pair_chr != c.pair_chr){
						clusters.insert(r.clust);
						contigs.insert(r.con);
						c.sup += r.sup;

					//	if(r.pair_loc != 0){
							pairs.insert({r.pair_loc, r.end, r.pair_chr});
							// c.pair_loc = r.pair_loc;
							// c.pair_chr = r.pair_chr;
					//	}					
					}

					if(PRINT_READS){
						contig con = all_contigs[r.con];
						fprintf(fr_vcf, ">Cluster: %d Contig: %d MaxSupport: %d BP: %d Reads: \n", r.clust, r.con, r.sup, 0); //TODO find way to get start loc in contig
						fprintf(fr_vcf, "ContigSeq: %s\n", con.data.c_str());
						for (auto &read: con.read_information)
							fprintf(fr_vcf, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());
					}
				}

				cluster_ids = "";
				contig_ids = "";
				type = sv_types[sv];

				if(clusters.empty() || contigs.empty()){
					cerr << "ERROR: Invalid call of type " << type << ": " << chromos[chr] << ":" << loc << " ploc " << chromos[c.pair_chr] << ":" << c.pair_loc  << endl;	
					continue;
				}

				if(PRINT_READS){
					for(auto& clust : clusters) {
						cluster_ids = cluster_ids + (to_string(clust) + ",");
					}  

					for(auto& con : contigs) {
						contig_ids = contig_ids + (to_string(con) + ",");
					}   

					cluster_ids.resize(cluster_ids.length()-1);
					contig_ids.resize(contig_ids.length()-1);
				}
				else{
					cluster_ids = to_string(clusters.size());
					contig_ids = to_string(contigs.size());
				}

				for(auto& p : pairs){

					ploc = p.loc;
					pend = p.end;
					pchr = p.chr;

					if(ploc == 0 || ploc == loc || (ploc > 0 && (results[sv][pchr].find(ploc) == results[sv][pchr].end()))){ //TODO deal with this last case

						if(c.one_bp)filter = "TWO_BP"; //Counter-intuitive I think, but rules are rules

						if(USE_ANNO){
							locate_interval(chromos[chr], loc, loc, gene_sorted_map[chromos[chr]], 0, iso_gene_map, best_gene1, best_name1, best_trans1, vec_best1);
							context1 = contexts[vec_best1[0]];
							locate_interval(chromos[pchr], pend, pend, gene_sorted_map[chromos[pchr]], 0, iso_gene_map, best_gene2, best_name2, best_trans2, vec_best2);
							context2 = contexts[vec_best2[0]];
							locate_interval(chromos[chr], loc, pend, gene_sorted_map[chromos[chr]], 0, iso_gene_map, best_gene3, best_name3, best_trans3, vec_best3);
							context3 = contexts[vec_best3[0]];


							anno = "ANNOL=" + context1 + "," + best_gene1 + "," + best_trans1 + "," + best_name1  + ";ANNOC=" + context3 + "," + best_gene3 + "," + best_trans3 + "," + best_name3 + ";ANNOR=" + context2 + "," + best_gene2 + "," + best_trans2 + "," + best_name2;  

							if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + ":FUSION";
							
						}

						if(sv == TRANS){

							if(c.info == "left_bp"){
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_1\t%s\t%s[%s:%d[\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_2;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), loc, c.u_id, c.ref_seq.c_str(), c.ref_seq.c_str(), chromos[pchr].c_str(), pend, filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_2\tN\t.N\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_1;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), (loc+1), c.u_id, filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
							}
							else if(c.info == "right_bp"){
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_1\tN\tN.\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_2;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), loc, c.u_id, filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_2\t%s\t]%s:%d]%s\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_1;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), (loc+1), c.u_id, c.ref_seq.c_str(), chromos[pchr].c_str(), pend, c.ref_seq.c_str(), filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
							}
							else{
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_1\t%s\t%s[%s:%d[\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_2;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), loc, c.u_id, c.ref_seq.c_str(), c.ref_seq.c_str(), chromos[pchr].c_str(), pend, filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
								fprintf(fo_vcf, "%s\t%d\tbnd_%d_2\t%s\t]%s:%d]%s\t.\t%s\tSVTYPE=BND;MATEID=bnd_%d_1;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
									chromos[chr].c_str(), (loc+1), c.u_id, c.ref_seq.c_str(), chromos[pchr].c_str(), (pend + c.alt_seq.length()), c.ref_seq.c_str(), filter.c_str(), c.u_id, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 
							}
						}
						else{

							start = loc;
							end = pend;

							// If it's more than uncertainty it's better to print it so we know something is wrong
							if(abs(pend-loc) <= uncertainty){
								start = min(loc, pend);
								end = max(loc, pend);
							}

							fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t.\t%s\tSVTYPE=%s%s;END=%d;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
								chromos[chr].c_str(), start, c.ref_seq.c_str(), c.alt_seq.c_str(), filter.c_str(), type.c_str(), c.info.c_str(), end, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str());
						}
					}
					else if(ploc > 0){

						unordered_set<long> clusters;
						unordered_set<int> contigs;
						cp.sup = 0;
						cp.end = 0;

						for(result& r : results[sv][pchr][ploc]){

							if(r.pair_loc < 0){

								if(r.end != cp.end && r.sup > cp.sup){

									cp = r;

									clusters.clear();
									clusters.insert(r.clust);

									contigs.clear();
									contigs.insert(r.con);
								}
								else if(r.end == cp.end){
									clusters.insert(r.clust);
									contigs.insert(r.con);
									cp.sup += r.sup;
									
								}
							}
						}

						cluster_ids_pair = "";
						contig_ids_pair = "";

						if(clusters.empty() || contigs.empty()){
							cerr << "ERROR: Invalid two contig call of type " << type << ": " << chromos[chr] << ":" << loc << " ploc " << chromos[pchr] << ":" << ploc  << endl;	
							continue;
						}

						if(PRINT_READS){
							for(auto& clust : clusters) {
								cluster_ids_pair = cluster_ids_pair + to_string(clust) + ",";
							}  

							for(auto& con : contigs) {
								contig_ids_pair = contig_ids_pair + to_string(con) + ",";
							}   

							cluster_ids_pair.resize(cluster_ids_pair.length()-1);
							contig_ids_pair.resize(contig_ids_pair.length()-1);
						}
						else{
							cluster_ids_pair = to_string(clusters.size());
							contig_ids_pair = to_string(contigs.size());
						}

						// int loc1 == m1.loc + m1.len  == loc
						// int loc2 == m2.loc			== c.end
						// int end1 == m4.loc			== cp.end
						// int end2 == m3.loc+m3.len    == c.pair_loc 

						if(sv == TRANS){
	bnd_count++;

							if(USE_ANNO){
								locate_interval(chromos[chr], loc, cp.end, gene_sorted_map[chromos[chr]], 0, iso_gene_map, best_gene1, best_name1, best_trans1, vec_best1);
								context1 = contexts[vec_best1[0]];
								locate_interval(chromos[pchr], c.end, ploc, gene_sorted_map[chromos[pchr]], 0, iso_gene_map, best_gene2, best_name2, best_trans2, vec_best2);
								context2 = contexts[vec_best2[0]];

								anno = "ANNOL=" + context1 + "," + best_gene1 + "," + best_trans1 + "," + best_name1 + ";ANNOR=" + context2 + "," + best_gene2 + "," + best_trans2 + "," + best_name2;  

								if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + ":FUSION";

							}

							if(abs(cp.end-loc) <= uncertainty){
								info = ":INS";
								if(abs(pend-ploc) <= uncertainty)continue; //1 FP
							}
							else{
								if(abs(pend-ploc) <= uncertainty){
									info = ":DEL";
								}
								else{
									info = ":SWAP";
								}
							}

							fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t%s[%s:%d[\t.\tPASS\tSVTYPE=BND%s;MATEID=bnd_%d;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
								chromos[chr].c_str(), loc, c.con, c.ref_seq.c_str(), c.ref_seq.c_str(), chromos[pchr].c_str(), pend, info.c_str(), cp.con, cluster_ids.c_str(), contig_ids.c_str(), c.sup, anno.c_str()); 

							if(USE_ANNO)anno = "ANNOL=" + context2 + "," + best_gene2 + "," + best_trans2 + "," + best_name2 + ";ANNOR=" + context1 + "," + best_gene1 + "," + best_trans1 + "," + best_name1; 

							fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t]%s:%d]%s\t.\tPASS\tSVTYPE=BND%s;MATEID=bnd_%d;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
								chromos[chr].c_str(), cp.end, cp.con, cp.ref_seq.c_str(), chromos[pchr].c_str(), ploc, cp.ref_seq.c_str(), info.c_str(), c.con, cluster_ids_pair.c_str(), contig_ids_pair.c_str(), cp.sup, anno.c_str()); 
						}
						else{
	not_bnd_count++;

							if(USE_ANNO){
								locate_interval(chromos[chr], loc, loc, gene_sorted_map[chromos[chr]], 0, iso_gene_map, best_gene1, best_name1, best_trans1, vec_best1);
								context1 = contexts[vec_best1[0]];
								locate_interval(chromos[pchr], ploc, ploc, gene_sorted_map[chromos[pchr]], 0, iso_gene_map, best_gene2, best_name2, best_trans2, vec_best2);
								context2 = contexts[vec_best2[0]];
								locate_interval(chromos[chr], loc, ploc, gene_sorted_map[chromos[chr]], 0, iso_gene_map, best_gene3, best_name3, best_trans3, vec_best3);
								context3 = contexts[vec_best3[0]];

								anno = "ANNOL=" + context1 + "," + best_gene1 + "," + best_trans1 + "," + best_name1  + ";ANNOC=" + context3 + "," + best_gene3 + "," + best_trans3 + "," + best_name3 + ";ANNOR=" + context2 + "," + best_gene2 + "," + best_trans2 + "," + best_name2;  
								
								if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + ":FUSION";
								
							}

							if(PRINT_READS){
								cluster_ids += "," + cluster_ids_pair;
								contig_ids += "," + contig_ids_pair;
							}
							else{
								cluster_ids = to_string(stoi(cluster_ids) + stoi(cluster_ids_pair));
								contig_ids = to_string(stoi(contig_ids) + stoi(contig_ids_pair));
							}

							if(sv == INV || sv == INS){
								start = loc;
								end = cp.end;

								if(abs(cp.end-loc) <= uncertainty){
									start = min(loc, ploc);
									end = max(loc, cp.end);
								}
							}
							else{
								start = loc;
								end = ploc;

								// If it's more than uncertainty it's better to print it so we know something is wrong
								if(abs(ploc-loc) <= uncertainty){
									start = min(loc, ploc);
									end = max(loc, ploc);
								}
							}
							
							fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t.\tPASS\tSVTYPE=%s%s;END=%d;CLUSTER=%s;CONTIG=%s;SUPPORT=%d;%s\n", 
								chromos[chr].c_str(), start, c.ref_seq.c_str(), c.alt_seq.c_str(), type.c_str(), c.info.c_str(),  end, cluster_ids.c_str(), contig_ids.c_str(), (c.sup+cp.sup), anno.c_str());
							
						}
					}// pair_loc > 0
				}//pairs
			}
		}
	}

if(PRINT_STATS){
cerr << "bnd_count " << bnd_count << endl;
cerr << "not_bnd_count " << not_bnd_count << endl;
}
}

//======================================================
// Main Functions
//======================================================


void kmistrvar::run_kmistrvar(const string &out_vcf, int min_support, int max_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE, int min_dist, int max_dist, int max_num_read, double clip_ratio)
{

clock_t begin = clock();
clock_t end;
double elapsed_secs;

	assemble(min_support, max_support, uncertainty, LOCAL_MODE, min_dist, max_dist, max_num_read, clip_ratio);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "ASSEMBLY TIME: " << elapsed_secs << endl;
cerr << endl;
begin = clock();
}

	//probalistic_filter();

// if(PRINT_STATS){
// end = clock();
// elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
// cerr << "PROBABILISTIC FILTER TIME: " << elapsed_secs << endl;
// cerr << endl;
// begin = clock();
// }

	index();

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "INDEX TIME: " << elapsed_secs << endl;
cerr << endl;
begin = clock();
}

	generate_intervals(out_vcf, LOCAL_MODE);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "MAPPING TIME: " << elapsed_secs << endl;
cerr << endl;
begin = clock();
}
	predict_variants(out_vcf, uncertainty, min_length, max_length);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "CALLING TIME: " << elapsed_secs << endl;
cerr << endl;
}

}

//======================================================
// Local Assembly Stage
//======================================================

void kmistrvar::assemble( int min_support, int max_support, int uncertainty, const bool LOCAL_MODE, int min_dist, int max_dist, int max_num_read, double clip_ratio)
{
	
	vector<contig> contigs;
	vector<contig> last_contigs;
	compressed_contig con;
	extractor ext(in_file, min_dist, max_dist, max_num_read, clip_ratio);
	double sum = 0;
	double count = 0; 

clock_t begin;
clock_t end;
double elapsed_secs;

//########### Read TP ##############
//FILE *writer = fopen("all.contigs.merged.reads", "wb");

// for(int i = 0; i < 14000; i++){
// 	TP[i] = false;
// }

// string line;
// ifstream myfile ("TP.contigs");

// if (myfile.is_open()) {
// 	while ( getline (myfile, line) ) { 
// 	  TP[stoi(line)] = true; 
// 	}
// 	myfile.close();
// }

//STATS
double part_count = 0;
double contig_count1 = 0;
double contig_count2 = 0;
double support_count = 0;
int read_count = 0;

	// If local mode enabled, use regions around contigs. Otherwise, add all chromosomes as regions
	if(!LOCAL_MODE){  
		for(char chr=0; chr < chromos.size(); chr++){
			regions.push_back({chr,{1, 300000000}});  //WARNING: This will not work with all reference files. Same with old version. 
		}
	}

	//Assemble contigs and build kmer index
	//=====================================
	cerr << "Assembling contigs..." << endl;

	while (1) {
		
		if(!ext.has_next_cluster()){
			break;
		}

		extractor::cluster p = ext.get_next_cluster(2, 2, false); //0,3 and 4,3 works for both!    5,2 for sim, 5,50 for 2     10,30  or 5, 200

		if (!p.reads.size() || p.reads.size() > 5000) //Temporary threshold until optimizations are finished. 
			continue;

		if(USE_BARCODES)p.reads = correct_reads(p.reads);

		if (!p.reads.size()) 
			continue;

		//Assemble contigs
		contigs = as.assemble(p.reads);

//STATS
if(PRINT_STATS){
	read_count += p.reads.size();
	part_count++;
	contig_count1 += contigs.size(); 	
}
	
		//Remove duplicates
		for(auto &contig: contigs){
			if (contig.support() >= min_support && contig.support() <= max_support && contig.data.size() <= MAX_ASSEMBLY_RANGE*2 ){
				for(auto &last_contig: last_contigs){
					if(contig.support() == last_contig.support() && contig.data == last_contig.data){
						contig.read_information.clear();   //investigate  chr9 5090790 and 5094130, tandem dup
					}
				}
			}
		}
		

		for (auto &contig: contigs){ 
			//TODO why are these long contigs happening? TP is in 6000bp contig :(
			if (contig.support() >= min_support && contig.support() <= max_support && contig.data.size() <= MAX_ASSEMBLY_RANGE*2 ){// && contig.data.size() <= (4*301)) { //Conservative threshold based on 4 times the longest read length of current short read sequencing technologies

//########### Dump Fastq ##############
// string qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
// cout << "@contig_" << contig_count2 << endl;
// cout << contig.data << endl;
// cout << "+" << endl;
// cout << qual.substr(0,contig.data.size()) << endl;

//STATS
if(PRINT_STATS){
	contig_count2++; 
	support_count += contig.support();	

	if(p.start >= 55174750 && p.start <= -55174800  ){
		cerr << "support: " << contig.support() << " " << contig.data.length() << " " << p.start << " " << cluster_info.size() << " " << p.ref << endl;
		cerr << contig.data << endl;
	}
}
				//if(LOCAL_MODE)regions.push_back({contig.cluster_chr, {(pt.get_start()-ref_flank), (pt.get_end()+ref_flank)}});

//########### Write Contigs with TP ##############
// cout << all_contig_metrics.size() << "\t" << contig.support() << "\t" << TP[all_contig_metrics.size()] << "\tchr" << p.ref << "\t" << p.start << endl;

//########### Write Clusters ##############
// fprintf(writer, ">Cluster: %d Contig: %d MaxSupport: %d TP: %d Reads: \n", p.start, all_contig_metrics.size(), contig.support(), TP[all_contig_metrics.size()]); //TODO find way to get start loc in contig
// fprintf(writer, "ContigSeq: %s\n", contig.data.c_str());
// for (auto &read: contig.read_information)
// 	fprintf(writer, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());

				
				// =================== Probablistic Filter =====================
				// set<double> positions;
				// double cur_pos = 0;
				// double last_pos = 0;
				// double max_dist = 0;
				// double num_reads = 1;


				// for (int j = 0; j < contig.support(); j++){
				// 	positions.insert((double)contig.read_information[j].location_in_contig);
				// }

				// for(double pos : positions){
				// 	if((pos-last_pos) > max_dist)max_dist = (pos-last_pos);
				// 	if(pos != 0){
				// 		sum += (pos-last_pos);
				// 		num_reads++; 
				// 	}
				// 	last_pos = pos;
				// }
				// count += num_reads;

				
				// double cov = (double)con.coverage[0];

				// all_contig_metrics.push_back({num_reads, (double)contig.data.length(), max_dist});
				// ==================================================================


				cluster_info.push_back({ p.start, (find(chromos.begin(), chromos.end(), p.ref) - chromos.begin()) });

				vector<bool> Ns = vector<bool>(contig.data.size(),false);

				//store N positions
				for (int i = 0; i < contig.data.length(); i++){
					if(contig.data[i] == 'N')Ns[i] = true;
				}

				all_contig_Ns.push_back(Ns);

				if(PRINT_READS){
					all_contigs.push_back(contig);
				}
				else{
					all_compressed_contigs.push_back(compress(contig));
				}

			}
		}
		//vector<contig>().swap(contigs);
		vector<contig>().swap(last_contigs);
		last_contigs = contigs;
	}

	if(all_contigs.empty() && all_compressed_contigs.empty()){
		cerr << "No contigs could be assembled. Exiting..." << endl;
		exit(1);
	}

	rate_param = 1/(sum/count);

//fclose(writer);

if(PRINT_STATS){
	cerr << "Read Count: " << read_count << endl;
	cerr << "Partition Count: " << part_count << endl;
	cerr << "Contigs: " << all_contigs.size() << " " << all_compressed_contigs.size() << endl;
	cerr << "Average Num Contigs Pre-Filter: " << (contig_count1/part_count) << endl;
	cerr << "Average Num Contigs Post-Filter: " << (contig_count2/part_count) << endl;
	cerr << "Average Contig Support: " << (support_count/contig_count2) << endl;
	cerr << "Rate Parameter Estimate: " << rate_param << endl;
}
}

void kmistrvar::probalistic_filter(){

	double pkall, pk, pk2;
	double pgk, pgk2;
	double contig_rate;
	int num_filtered = 0;

	for(int j = 0; j < all_contig_metrics.size(); j++){

		contig_metrics& m = all_contig_metrics[j];
		
		//Simple
		pk = (m.len-150) > 0 ? m.num_reads/(m.len-150) : 0;
		pgk = 1-pow((1-pow((1-pk),m.max_dist)),(m.num_reads-1));

		//Poisson
		pkall = (pow(rate_param*m.len, m.num_reads)*exp(-rate_param*m.len))/factorial(m.num_reads);

		//pk2 = (rate_param*m.max_dist)*exp(-rate_param*m.max_dist);
		//pgk2 = 1-pow((1-(1-pk2)),(m.num_reads-1));
		feclearexcept(FE_ALL_EXCEPT);
    	errno = 0;
		//pk2 = exp(-rate_param*m.max_dist);
		contig_rate = (m.len-150) > 0 ? m.num_reads/(m.len-150) : 0;
		pk2 = exp(-contig_rate*m.max_dist);
		if(errno == ERANGE)pk2 = 0;

		pgk2 = 1-pow((1-pk2),(m.num_reads-1));

		if(pgk2 < 0){  //0.001 for sim -2DUP -2FP, 0.00001 -2FP  0.0000001
		//if(pgk < 0.3){  //0.2 -4FP
			if(PRINT_READS){
				all_contigs[j].data.clear();
			}
			else{
				all_compressed_contigs[j].data.clear();
				num_filtered++;
			}
		}
	}

if(PRINT_STATS){
	cerr << "Filtered " << num_filtered << " contigs." << endl;
}

}


//====================================================== 
// Contig Indexing Stage
//======================================================  
													
void kmistrvar::index()
{
	
	int contig_id = 0;
	int max_hits = 0;
	int i, j, x;
	unsigned long long index = 0;
	vector<int> contig_list;
	contig_kmers.push_back(contig_list);

	//Assemble contigs and build kmer index
	//=====================================
	cerr << "Indexing contigs..." << endl;

	contig_kmers_index = new int[num_kmer]; //Dangerous thing to do, think of different solution

	for(int i = 0; i < num_kmer; i++){
		contig_kmers_index[i] = 0;
	} 
	
	//TODO: Contig merging
	for (contig_id = 0; contig_id < max(all_contigs.size(), all_compressed_contigs.size()); contig_id++){
		
		string con = PRINT_READS ? all_contigs[contig_id].data :  con_string(contig_id, 0, all_compressed_contigs[contig_id].data.size()/2); 
		tsl::sparse_map<int, vector<int>> kmer_location;
		kmer_location.reserve(con.size());
		kmer_locations.push_back(kmer_location);
		
		i = 0, j = 0, x = 0;

		// Per chromosome repeat flag for each contig
		for(int rc=0; rc < 2; rc++){
			for(int chr=0; chr < 25; chr++){
				last_intervals[rc][chr].push_back((last_interval){ 0, k, 0});
			}
		}

		if(con.empty())continue;

		// k-1 kmer for first round
		for(j = 0; j < k-1; j++){ 
			if(j > 0 && con[i+j] == 'N'){
				i += (j+1);
				j = -1;
				index = 0;
				continue;
			}
			index <<= 2;
			index |= ((con[i+j] & MASK) >> 1);
		}

		// Skip N nucleotides
		for (i; i < con.length()-k+1; i++){

			if(con[i+k-1] == 'N'){
				i+=k;
				while(con[i] == 'N'){
					i++;
				}

				index = 0;

				for(x = 0; x < k-1; x++){
					if(con[i+x] == 'N'){
						i += (x+1);
						x = -1;
						index = 0;
						continue;
					}
					index <<= 2;
					index |= ((con[i+x] & MASK) >> 1);
				}
				i--;
				continue;
			}



			// Compute next kmer via bit shift
			index <<= 2;
			index |= ((con[i+k-1] & MASK) >> 1); 
			index &= kmer_mask;

			//if(i % (k+2) != 0)continue;

			int& cur_index = contig_kmers_index[index];

			// Add to index
			if(!cur_index){
				contig_kmers_index[index] = contig_kmers.size();

				vector<int> contig_list;
				contig_list.reserve(1);

				contig_list.push_back(contig_id);
				contig_kmers.push_back(contig_list);
			}
			else{
				if(contig_kmers[cur_index].back() != contig_id){
					contig_kmers[cur_index].push_back(contig_id);
				}
			}
			kmer_locations[contig_id][cur_index].push_back(i); 
			
		}
	}

	// Initialize mappings vector based on number of contigs
	for(int rc=0; rc <= 1; rc++){
		for(int chr=0; chr < 25; chr++){
			contig_mappings[rc][chr] = vector<vector<mapping>>(contig_id, vector<mapping>(1, {0, k, -1}));
		}
	}

//STATS
if(PRINT_STATS){
	cerr << "Number of unique kmers: " << contig_kmers.size() << endl;
	cerr << "Max hits: " << max_hits << endl;
}	
}


//======================================================
// Interval Mapping Stage
//======================================================
void kmistrvar::generate_intervals_bwa(const string &sam_file, const bool LOCAL_MODE){

	for(int rc=0; rc <= 1; rc++){
		for(int chr=0; chr < 25; chr++){
			contig_mappings[rc][chr] = vector<vector<mapping>>(30826, vector<mapping>());
		}
	}

	SAMParser* parser = new SAMParser(sam_file);
	string comment = parser->readComment();
	

	while(parser->hasNext()){

		const Record &rec = parser->next();
		mapping interval;
		interval.loc = rec.getLocation();
		string cigar = string(rec.getCigar());
		string optional = string(rec.getOptional());
		string supple = "";	
		char chromo = (find(chromos.begin(), chromos.end(), string(rec.getChromosome())) - chromos.begin());
		char* con_name = strdup(rec.getReadName()); 
		strtok(con_name, "_");
		uint32_t flag = rec.getMappingFlag();
		int con_id = stoi(strtok(NULL, "_"));
		int len = 0;
		int con_loc = 0;
		int con_len = all_compressed_contigs[con_id].data.size()/2;
		bool neg = ( 0x10 == (flag&0x10));

		for(char c : cigar){

			if(isdigit(c)){
				len *= 10;
				len +=  int(c - '0');
			}
			else{
				if(c == 'S' || c == 'H' || c == 'I'){
					con_loc += len;
				}
				else if(c == 'M'){
					interval.len = len;
					interval.con_loc = neg ? con_len-con_loc-k : con_loc+len-k;
					if(len >= k)contig_mappings[neg][chromo][con_id].push_back(interval);
					con_loc += len;
					interval.loc += len;
				}
				else{ 
					interval.loc += len;
				}
				len = 0;
			}
		}

		free(con_name);

		parser->readNext();
	}

	delete parser;
}

void kmistrvar::generate_intervals(const string &out_vcf, const bool LOCAL_MODE)
{

	//Read reference and map contig kmers
	//=====================================
	cerr << "Generating intervals..." << endl;

	char chromo = -1;
	char last_chromo = chromo;
	int use_rc = 2; //1 no, 2 yes
	int jj = 0;
	int id = 1;  //zero reserved for sink
	unsigned long long index = 0;
	unsigned long long rc_index = 0;
	unsigned long long cur_index = 0;
	unsigned long long cur_index2 = 0;
	unsigned long long index_mask = 0;
	long iend;
	int shift = 2*(k-1)-1;
	int i, ii, ii1, ik, iik, iik1, j, js, je, x, y, z, rc, len, ilen, k1 = k+1, km1 = k-1;
	int MAX_CONTIG_LEN = 0; 
	int MAX_INTERVAL_LEN = 20000;
	int max_interval_len_k1;
	int num_contigs = PRINT_READS ? all_contigs.size() : all_compressed_contigs.size();
	int gen_start, gen_end, interval_len, contig_len, contig_lenk1;
	bool con_repeat;
	vector<pair<int,int>> starts;
	string gen;
	FILE *fo_vcf = fopen(out_vcf.c_str(), "wb");

	fprintf(fo_vcf, "##fileformat=VCFv4.2\n##source=SVICT-v1.0\n");

//STATS
int anchor_intervals = 0;
int non_anchor_intervals = 0;
int pop_back1 = 0;
int pop_back2 = 0;
int pop_back3 = 0;
int pop_back4 = 0;
int max_ints = 0;
int max_int_len = 0;
int max_con_len = 0;
int max_starts = 0;
int max_sub = 0;
int index_misses = 0;
long kmer_hits = 0;
long kmer_hits_con = 0;
double interval_count = 0;
int corrected_bp = 0;

	// Initialization 
	if(PRINT_READS){
		for(auto& contig : all_contigs){
			if(contig.data.length() > MAX_CONTIG_LEN) MAX_CONTIG_LEN = contig.data.length();
		}
	}
	else{
		for(auto& contig : all_compressed_contigs){
			if(contig.data.size()/2 > MAX_CONTIG_LEN) MAX_CONTIG_LEN = contig.data.size()/2;
		}
	}
	
	MAX_CONTIG_LEN += k1; //Rare case of SNP near the end of the longest contig. 
	MAX_INTERVAL_LEN = MAX_CONTIG_LEN+100;
	max_interval_len_k1 = MAX_INTERVAL_LEN-k-1;

	bool** valid_mappings = new bool*[MAX_INTERVAL_LEN]; //think of a better solution

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		valid_mappings[i] = new bool[MAX_CONTIG_LEN];
	}

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		for(int j=0; j < MAX_CONTIG_LEN; j++){
			valid_mappings[i][j] = false;
		}
	}


	// Traverse each chromosome or region (Local mode)
	for (auto &region: regions){

		chromo = region.first;
		gen = ref.extract(chromos[chromo], region.second.first, region.second.second); //get region sequence
		gen_start = max(0,region.second.first-1);
		gen_end = gen.length()-k+1;

		fprintf(fo_vcf, "##contig=<ID=%s,length=%d>\n", chromos[chromo].c_str(), gen.length());

		js = LOCAL_MODE ? jj : 0;
		je = LOCAL_MODE ? (jj+1) : num_contigs;

		// Skip leading Ns 
		i = 0;
		while(gen[i] == 'N'){
			i++;
		}

		index = 0;
		rc_index = 0;

		// k-1 kmer for first round
		for(x = 0; x < km1; x++){    
			if(x > 0 && gen[i+x] == 'N'){
				i += (x+1);
				x = -1;
				index = 0;
				rc_index = 0;
				continue;
			}
			index <<= 2;
			index |= ((gen[i+x] & MASK) >> 1);
			rc_index <<= 2;
			rc_index |= (((gen[(i+km1)-x] & MASK) ^ MASK_RC) >> 1);
		}

		// Skip N nucleotides
		for (i; i < gen_end; i++){ //-ANCHOR_SIZE
			if(gen[i+km1] == 'N'){
				i+=k;
				while(gen[i] == 'N'){
					i++;
				}

				index = 0;
				rc_index = 0;

				for(x = 0; x < km1; x++){
					if(gen[i+x] == 'N'){
						i += (x+1);
						x = -1;
						index = 0;
						rc_index = 0;
						continue;
					}
					index <<= 2;
					index |= ((gen[i+x] & MASK) >> 1);
					rc_index <<= 2;
					rc_index |= (((gen[(i+km1)-x] & MASK) ^ MASK_RC) >> 1);
				}
				i--;
				continue;
			}

			// Precompute some arithmetic 
			// index_mask = 2nd and 3rd bit of char, which is a unique value for each nucleotide 
			ii = i + gen_start;// + 1;
			ii1 = ii-1;
			iik1 = ii1+k;
			index_mask = (gen[i+km1] & MASK);

			// Forward and reverse strand
			for(rc=0; rc < use_rc; rc++){

				// Convert kmer string to index
				// Bit shift (reuse last kmer) -> shift mask to 1st and 2nd bit -> mask bits for char of last kmer
				// Slightly modified for reverse compliment to flip bits
				if(rc){
					rc_index = ((rc_index >> 2) | ((index_mask ^ MASK_RC) << shift));
					if(!contig_kmers_index[rc_index])continue;
					cur_index = rc_index;
				}
				else{
					index = (((index << 2) | (index_mask >> 1)) & kmer_mask);
					if(!contig_kmers_index[index])continue;
					cur_index = index;
				}

kmer_hits++;
//contig_kmer_counts[contig_kmers_index[cur_index]]++; //178956970

				for(const int& j : contig_kmers[contig_kmers_index[cur_index]]){

					// Local mode temporarily disabled
					//if( (!LOCAL_MODE || j == jj)){ 
kmer_hits_con++;
						last_interval& l_interval = last_intervals[rc][chromo][j];
						iend = l_interval.loc + l_interval.len;

						// Else end last interval and start new one
						if(iend < ii1){

							// Get last interval for this contig
 							vector<mapping>& intervals = contig_mappings[rc][chromo][j];

							// Check if sufficiently long
							if(l_interval.len >= ANCHOR_SIZE && l_interval.len < max_interval_len_k1){

								mapping& cur_interval = intervals.back();
								cur_interval.loc = l_interval.loc;
								cur_interval.len = l_interval.len;
								intervals.push_back((mapping){ii, k, -1}); 
							}

							l_interval.loc = ii;
							l_interval.len = k;

						}
						// Extend interval by one
						else if(iend == iik1){
							l_interval.len++;
						}
						// Extend interval by k+1 if a SNP/SNV is detected
						else if(iend == ii1){
							l_interval.len += k1;
						}

					//	if(LOCAL_MODE)break;
					//}
				}
			}
		}// gen, slowest loop 



		//Remove intervals without consectutive kmers in contigs
		//======================================================
		for(rc=0; rc < use_rc; rc++){
			for(j = js; j < je; j++){

				contig_len = PRINT_READS ? all_contigs[j].data.length() : all_compressed_contigs[j].data.size()/2;
				last_interval& l_interval = last_intervals[rc][chromo][j];
				mapping& back_interval = contig_mappings[rc][chromo][j].back();

				if(l_interval.loc == back_interval.loc){
					contig_mappings[rc][chromo][j].pop_back();
 					if(contig_mappings[rc][chromo][j].empty())continue;
				}
				else{
					if(l_interval.len >= ANCHOR_SIZE && l_interval.len < max_interval_len_k1){
						back_interval.loc = l_interval.loc;
						back_interval.len = l_interval.len;
					}
				}

				if(contig_len == 0)continue;
		 		
		 		pair<int,char>& bp = cluster_info[j];
		 		vector<bool>& Ns = all_contig_Ns[j];
				vector<mapping> valid_intervals;
				contig_lenk1 = contig_len+k1;

if(contig_mappings[rc][chromo][j].size() > max_ints)max_ints = contig_mappings[rc][chromo][j].size();

if(contig_len > max_con_len)max_con_len = contig_len;

				for(auto &interval: contig_mappings[rc][chromo][j]){

					if(interval.len > contig_lenk1)continue; 

//STATS
if(PRINT_STATS){
int& contig_loc = cluster_info[j].first;
char& contig_chr = cluster_info[j].second;
if(chromo == contig_chr && abs(interval.loc - contig_loc) < MAX_ASSEMBLY_RANGE){
anchor_intervals++;
}
else{
non_anchor_intervals++;
}
}					

if(interval.len > max_int_len)max_int_len = interval.len;

int sub_count = 0;

					// intialization
					starts.clear();
					
					ii = interval.loc - gen_start;// + 1;
					len = interval.len-km1;
					interval_len = interval.len+k1; //maybe off by one	
					con_repeat = false;					

					cur_index = 0;

					if(rc){
						// k-1 kmer for first round
						for(x = 0; x < km1; x++){    
							cur_index <<= 2;
							cur_index |= (((gen[(ii+interval.len-1)-x] & MASK) ^ MASK_RC) >> 1);
						}
					}
					else{
						// k-1 kmer for first round
						for(x = 0; x < km1; x++){    
							cur_index <<= 2;
							cur_index |= ((gen[ii+x] & MASK) >> 1);
						}
					}

					// Traverse the interval
					for (i = 0; i < len; i++){

						if(rc){
							cur_index = (((cur_index << 2) | (((gen[ii+(interval.len-1-i)-km1] & MASK) ^ MASK_RC) >> 1)) & kmer_mask);
						}
						else{
							cur_index = (((cur_index << 2) | ((gen[ii+i+km1] & MASK) >> 1)) & kmer_mask);
						}

						int& kmer_index = contig_kmers_index[cur_index];

						if(kmer_index && binary_search(contig_kmers[kmer_index].begin(), contig_kmers[kmer_index].end(), j)){

							for(auto &loc: kmer_locations[j][kmer_index]){

								if(!Ns[loc]){
									valid_mappings[i][loc] = true;

									// Record starts when interval no longer consecutive
									if(i == 0 || loc == 0){
										starts.push_back({i,loc});
									}
									else if((i < k1 || loc < k1) && !valid_mappings[i-1][loc-1]){
										starts.push_back({i,loc});
									}
									else if(!valid_mappings[i-1][loc-1] && !valid_mappings[i-k1][loc-k-1]){
										starts.push_back({i,loc});
									}
								}
							}	
						}

						if(starts.size() > CON_REPEAT_LIMIT){
							con_repeat = true;
							break;
						}
					}

if(starts.size() > max_starts)max_starts = starts.size();


					// Check each start point
					for(auto& start : starts){

						x = start.first;
						y = start.second;
						
						// Traverse diagonal
						while(x < interval.len && y < contig_len && valid_mappings[x][y]){
							valid_mappings[x++][y++] = false;
							if(y+k < contig_len){
								if(!valid_mappings[x][y]){
									if(valid_mappings[x+k][y+k]){
										x +=k;
										y +=k;
									}
								}
							}
						}

						x--;
						y--;

						ilen = (y - start.second)+k;

						//TODO: Use this table to identify duplications, although maybe no point since all can be found without this
						if(ilen > ANCHOR_SIZE && !con_repeat){

							//Add to final list of intervals
							mapping sub_interval;
							sub_interval.loc = interval.loc + start.first;
							sub_interval.len = ilen;
							sub_interval.con_loc = y;
							valid_intervals.push_back(sub_interval); 
sub_count++;							
if(sub_count > max_sub)max_sub = sub_count;
						}
					}
				}//intervals


				//SNP near BP error correction.
				if(!valid_intervals.empty() && valid_intervals.size() < 800){ // Arbitrary for now
corrected_bp++;
					sort(valid_intervals.begin(), valid_intervals.end(), con_interval_compare());
					for(i = 0; i < valid_intervals.size(); i++){
						mapping& bp_map = valid_intervals[i];
						if(i > 0 && abs(bp_map.loc - bp.first) < 3){
	 						for(x = i-1; x >= 0; x--){
	 							mapping& prior_map = valid_intervals[x];
	 							int e_con_loc = (prior_map.con_loc+k);
	 							int bp_con_loc = ((bp_map.con_loc+k)-bp_map.len);
	 							int missed_len = (bp_con_loc-e_con_loc)-1;
	 							
	 							if(bp_con_loc > e_con_loc+1 && (bp_con_loc-e_con_loc) <= k){
	 								string contig_seq = PRINT_READS ? all_contigs[j].data : con_string(j, 0, all_compressed_contigs[j].data.size()/2);
	 								string missed_con_seq_left = contig_seq.substr(e_con_loc+1, missed_len);
	 								string missed_con_seq_right = contig_seq.substr(e_con_loc, missed_len);
	 								string missed_ref_seq_left = gen.substr(prior_map.loc+prior_map.len+1,missed_len);
	 								string missed_ref_seq_right = gen.substr(bp_map.loc-missed_len,missed_len);

	 								if(missed_ref_seq_left == missed_con_seq_left){
	 									prior_map.con_loc += missed_len+1;
	 									prior_map.len += missed_len+1;
	 								}
	 								else if(missed_ref_seq_right == missed_con_seq_right){
	 									bp_map.loc -= missed_len+1;
	 									bp_map.len += missed_len+1;
	 								}
	 								break;
	 							}
	 						}
	 						break;	
	 					}
	 					else if(i < valid_intervals.size()-1 && abs((bp_map.loc+bp_map.len) - bp.first) < 3){
							for(x = i+1; x < valid_intervals.size(); x++){
	 							mapping& post_map = valid_intervals[x];
	 							int bp_con_loc = (bp_map.con_loc+k);
	 							int s_con_loc = ((post_map.con_loc+k)-post_map.len);
	 							int missed_len = (s_con_loc-bp_con_loc)-1;

	 							if(s_con_loc > bp_con_loc+1 && (s_con_loc-bp_con_loc) <= k){
	 								string contig_seq = PRINT_READS ? all_contigs[j].data : con_string(j, 0, all_compressed_contigs[j].data.size()/2);
	 								string missed_con_seq_left = contig_seq.substr(bp_con_loc+1, missed_len);
	 								string missed_con_seq_right = contig_seq.substr(bp_con_loc, missed_len);
	 								string missed_ref_seq_left = gen.substr(bp_map.loc+bp_map.len+1, missed_len);
	 								string missed_ref_seq_right = gen.substr(post_map.loc-missed_len, missed_len);

	 								if(missed_ref_seq_left == missed_con_seq_left){
	 									bp_map.con_loc += missed_len+1;
	 									bp_map.len += missed_len+1;
	 								}
	 								else if(missed_ref_seq_right == missed_con_seq_right){
	 									post_map.loc -= missed_len+1;
	 									post_map.len += missed_len+1;
	 								}
	 								break;
	 							}
	 						}
	 						break;
	 					}						
					}
				}

//STATS
interval_count += valid_intervals.size(); 

				num_intervals += valid_intervals.size();
				contig_mappings[rc][chromo][j] = move(valid_intervals);
			}//j
		}//rc


		jj++;

	}//regions

//STATS
if(PRINT_STATS){
	cerr << "Intervals: " << num_intervals << endl;
	cerr << "Anchor Intervals: " << anchor_intervals << endl;
	cerr << "Non-Anchor Intervals: " << non_anchor_intervals << endl;
	cerr << "Average Intervals per Contig: " << (interval_count/(double)num_contigs) << endl;
	cerr << "Pop backs Short: " << pop_back1 << endl;
	cerr << "Pop backs Repeat1: " << pop_back2 << endl;
	cerr << "Pop backs Repeat2: " << pop_back3 << endl;
	cerr << "Pop backs Last: " << pop_back4 << endl;

	cerr << "Max Intervals: " << max_ints << endl;
	cerr << "Max Interval Len: " << max_int_len << endl;
	cerr << "Max Contig Len: " << max_con_len<< endl;
	cerr << "Max Starts: " << max_starts << endl;
	cerr << "Max Subintervals: " << max_sub << endl;
	cerr << "Index Misses: " << index_misses << endl;
	cerr << "Kmer Hits: " << kmer_hits << endl;
	cerr << "Kmer Hits Total: " << kmer_hits_con << endl;
	cerr << "Corrected BP: " << corrected_bp << endl;
}

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		delete[] valid_mappings[i];
	}
	delete[] valid_mappings;

	fclose(fo_vcf);
}


//======================================================
// Predict Structural Variants Stage
//======================================================

void kmistrvar::predict_variants(const string &out_vcf, int uncertainty, int min_length, int max_length)
{
	//Predict structural variants
	//=====================================
	cerr << "Predicting short structural variants..." << endl;

//STATS
int normal_edges = 0, ins_edges = 0, source_edges = 0, sink_edges = 0, singletons = 0;
int path_count = 0;
int intc1 = 0, intc2 = 0, intc3 = 0, intc4 = 0, intc5 = 0, intc6 = 0;
int anchor_fails = 0, ugly_fails = 0;
int short_dup1 = 0, short_dup2 = 0, short_dup3 = 0, short_del = 0, short_inv1 = 0, short_inv2 = 0, short_ins1 = 0, short_ins2 = 0, short_trans = 0, mystery1 = 0, mystery2 = 0, mystery3 = 0;
int long_del = 0, long_inv1 = 0, long_inv2 = 0, long_inv3 = 0, long_ins = 0, long_trans = 0, long_dup = 0;
int one_bp_del = 0, one_bp_dup = 0, one_bp_inv = 0, one_bp_ins = 0, one_bp_trans = 0;
if(PRINT_STATS)num_intervals = 0;

	int interval_count = 0;
	int one_bp_id = 0;
	int pair_id = 0;
	int num_contigs = PRINT_READS ? all_contigs.size() : all_compressed_contigs.size();
	int rc = 0;
	int use_rc = 1;
	int contig_loc, contig_len, loc, id, ref_dist, con_dist;
	unsigned long long index = 0;
	bool found = false;
	char contig_chr;
	string contig_seq;
	mapping_ext cur_interval;
	mapping_ext source = {-1, 0, 0, -1, 0, 0, 0};
	mapping_ext sink = {-1, 0, 0, -1, 0, 0, 0};
	min_length = max(min_length, uncertainty+1);

	vector<sortable_mapping>* sorted_intervals = new vector<sortable_mapping>[25];	// Sorted list for traversal 
	unordered_map<int, vector<int>> interval_pair_ids;								// Mapping from interval -> interval_pairs
	vector<mapping_ext> all_intervals;												// Bulk data of unused intervals
	vector<interval_pair> interval_pairs; 											// Pairs of intervals

	FILE *fo_vcf = fopen(out_vcf.c_str(), "ab");
	FILE *fr_vcf = fopen((out_vcf + ".reads").c_str(), "wb");

	fprintf(fo_vcf, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(fo_vcf, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
	if(PRINT_READS){
		fprintf(fo_vcf, "##INFO=<ID=CLUSTER,Number=.,Type=Integer,Description=\"Cluster IDs supporting this SV\">\n");
		fprintf(fo_vcf, "##INFO=<ID=CONTIG,Number=.,Type=Integer,Description=\"Contigs IDs supporting this SV\">\n");
	}
	else{
		fprintf(fo_vcf, "##INFO=<ID=CLUSTER,Number=1,Type=Integer,Description=\"Number of clusters supporting this SV\">\n");
		fprintf(fo_vcf, "##INFO=<ID=CONTIG,Number=1,Type=Integer,Description=\"Number of contigs supporting this SV\">\n");
	}
	fprintf(fo_vcf, "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Maximum basepair read support at or nearby the breakpoint\">\n");
	if(USE_ANNO){
		fprintf(fo_vcf, "##INFO=<ID=ANNOL,Number=4,Type=String,Description=\"Annotation information for region left of the BP: genomic context, gene ID, transcript ID, gene name\">\n");
		fprintf(fo_vcf, "##INFO=<ID=ANNOC,Number=4,Type=String,Description=\"Annotation information for inserted region (if applicable): genomic context, gene ID, transcript ID, gene name\">\n");
		fprintf(fo_vcf, "##INFO=<ID=ANNOR,Number=4,Type=String,Description=\"Annotation information for region right of the BP: genomic context, gene ID, transcript ID, gene name\">\n");
	}
	fprintf(fo_vcf, "##FILTER=<ID=TWO_BP,Description=\"Support for both breakpoints\">\n");
	fprintf(fo_vcf, "#CHROM  POS ID  REF ALT QUAL  FILTER  INFO\n");

	// Traverse each contig
	for(int j = 0; j < num_contigs; j++){

		contig_loc = cluster_info[j].first;
		contig_chr = cluster_info[j].second;
		contig_seq = PRINT_READS ? all_contigs[j].data : con_string(j, 0, all_compressed_contigs[j].data.size()/2);
		contig_len = PRINT_READS ? all_contigs[j].data.length() : all_compressed_contigs[j].data.size()/2;

		if(contig_len == 0)continue;

		vector<vector<mapping_ext>> intervals(contig_len+1, vector<mapping_ext>());
		unordered_map<int, pair<int,int>> lookup;
		vector<mapping_ext> grouped_intervals; 
		vector<int> mapped = vector<int>(contig_len, 0); 
		vector<bool>& Ns = all_contig_Ns[j];
		bool visited[contig_len];
		memset(visited, 0, sizeof(visited));
		id = 1;
		interval_count = 1;

		// Handle repeats, sort by length and select any interval included an unmapped contig base. 
		for(char chr=0; chr < 25; chr++){  
			for(rc=0; rc <= use_rc; rc++){ 
				if(!contig_mappings[rc][chr][j].empty()){
					for(auto &interval: contig_mappings[rc][chr][j]){  

						grouped_intervals.push_back(copy_interval(chr, rc, j, interval));
					}
					contig_mappings[rc][chr][j].clear();
				}
			}
		}

		sort(grouped_intervals.begin(), grouped_intervals.end());

		for(auto& interval : grouped_intervals){
		
			for(int x = ((interval.con_loc+k)-interval.len); x < (interval.con_loc+k); x++){
				if(mapped[x] < REPEAT_LIMIT){
					for(int y = ((interval.con_loc+k)-interval.len); y < (interval.con_loc+k); y++){
						mapped[y]++;
					}

					loc = (interval.con_loc+k-interval.len);

					//Ns at breakpoint correction
					while(loc-1 >= 0 && Ns[loc-1]){
						loc--;
					}

					if(loc < 0){
						print_interval("ERROR Invalid Interval:", interval);
						continue;
					}

					intervals[loc].push_back(interval);
					id++;
if(j == CON_NUM_DEBUG)print_interval("Included", interval);
					break;
				}
			}

		}

		intervals[contig_len].push_back(source);
		intervals[contig_len].push_back(sink);
		lookup[0] = {contig_len,0};
		lookup[id] = {contig_len,id};

		for(int i = 0; i < contig_len; i++){
			for(int x = 0; x < intervals[i].size(); x++){
				intervals[i][x].id = interval_count++;
				lookup[intervals[i][x].id] = {i,x};
			}
		}

		interval_count++;
num_intervals += interval_count;

		// Initialization
		int** contig_graph = new int*[interval_count]; //think of a better solution

		for(int i=0; i < interval_count; i++){
			contig_graph[i] = new int[interval_count];
		}

		for(int i=0; i < interval_count; i++){
			for(int j=0; j < interval_count; j++){
				contig_graph[i][j] = 0;
			}
		}

		BUILD_DIST = ANCHOR_SIZE;

		//Build Edges
		//==================================
		for(int i = 0; i < contig_len; i++){

			if(!intervals[i].empty()){

				for(auto &interval1: intervals[i]){
					found = false; 

					//Build edge for balanced SVs and deletions
					for(int u = max(0, i+interval1.len-BUILD_DIST); u <= min(contig_len-1, i+interval1.len+BUILD_DIST); u++){ //was k
						if(!intervals[u].empty()){
							for(auto &interval2: intervals[u]){	

								if(!interval1.rc && !interval2.rc){
									ref_dist = interval2.loc-(interval1.loc+interval1.len);
								}
								else if(interval1.rc && !interval2.rc){
									ref_dist = (interval2.loc-interval1.loc);
								}
								else if(!interval1.rc && interval2.rc){
									ref_dist = (interval2.loc+interval2.len-(interval1.loc-interval1.len));
								}
								else{
									ref_dist = ((interval2.loc-interval2.len)-interval1.loc);
								}

								con_dist = (interval2.con_loc+k-interval2.len)-(interval1.con_loc+k);

								if(abs(abs(ref_dist)-abs(con_dist)) <= uncertainty || ref_dist <= uncertainty || con_dist <= BUILD_DIST){//uncertainty){ //Temporary until performance issues are resolved

									visited[u] = true;
									found = true;
									contig_graph[interval1.id][interval2.id] = -interval2.len;
normal_edges++;
								}
							}
						}
					}

					//Build edges for insertions
					if(!found){

						//to sink
						if(!found && ( ((interval1.con_loc+k-interval1.len) <= BUILD_DIST) || ((contig_len - (interval1.con_loc+k)) <= BUILD_DIST) || interval1.id == (id-1)) ){ 
							contig_graph[interval1.id][id] = -interval1.len;
sink_edges++;
						}
					}

					//from source
					if(!visited[i] && ( ((interval1.con_loc+k-interval1.len) <= BUILD_DIST) || ((contig_len - (interval1.con_loc+k)) <= BUILD_DIST)  || interval1.id == 1) ){ 
						contig_graph[0][interval1.id] = -interval1.len;
source_edges++;
					}
				}
			}
		}

cur_debug_id = j;

		// Run Ford Fulkerson for max flow
		//====================================

		//vector<vector<int>> paths2 = fordFulkerson(id+1, contig_graph, 0, id); 
		vector<pair<vector<int>, int>> paths = dijkstra(id+1, contig_graph, 0, id, 5);



		// Predict short variants
		//====================================

		if(!paths.empty()){

			int max_paths = PATH_LIMIT;
			int short_dist = paths[0].second;

			// Check each path and classify variant
			for(int p = 0; p < min(max_paths, (int)paths.size()); p++){  //top x paths

if(j == CON_NUM_DEBUG)cerr << paths[p].second << " - " << short_dist << " = " << (paths[p].second - short_dist) << endl;

				if(paths[p].second - short_dist > SUB_OPTIMAL)break; //only co-optimal

				vector<int>& path = paths[p].first;

if(j == CON_NUM_DEBUG){

	cerr << "All examined paths:" << endl;

	for(int i = path.size()-2; i > 0; i--){
		cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].loc << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].loc + intervals[lookup[path[i]].first][lookup[path[i]].second].len << "\t";
	}
	cerr << endl;
	for(int i = path.size()-2; i > 0; i--){
		cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].con_loc+k - intervals[lookup[path[i]].first][lookup[path[i]].second].len << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].con_loc+k << "\t";
	}
	cerr << "\n========================================================" << endl;
}
				mapping_ext cur_i, i1, i2;
				int chromo_dist, contig_dist;
				int a1 = -1;
				int a2 = -1;
				bool is_anchor;
				bool anchor_found = false;

				for(int i = path.size()-2; i > 0; i--){

					// Check if anchor (located near partition range)
					cur_i = intervals[lookup[path[i]].first][lookup[path[i]].second];

					is_anchor = ((cur_i.chr == contig_chr) && (abs(cur_i.loc+(cur_i.len/2) - contig_loc) < (MAX_ASSEMBLY_RANGE*2)));
					found = false;
					if(is_anchor)anchor_found = true;

					//Note: For a small variant to be called it requires two anchors flanking the variant
					// Otherwise "interval pairs" are created and passed on to long structural variant detection. 

					// If first anchor not set 
					if(a1 == -1){
						if(is_anchor){

							a1 = i;

							//leading non-anchor region
							if(a1 != path.size()-2){
								
				 				i1 = intervals[lookup[path[a1+1]].first][lookup[path[a1+1]].second];
				 				i2 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
				 				chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
				 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

								if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){ 
intc1++;
									i1.id = all_intervals.size();
									all_intervals.push_back(i1);
									i2.id = all_intervals.size();
									all_intervals.push_back(i2);
									sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
									sorted_intervals[i2.chr].push_back({i2.id,i2.loc});
									interval_pairs.push_back({false,i1.id,i2.id});
									interval_pair_ids[i1.id].push_back(pair_id);
						 		 	interval_pair_ids[i2.id].push_back(pair_id++);
								}
							}
							else{

								i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];

								if(i1.con_loc+k-i1.len > ((int)contig_len - (i1.con_loc+k)) || a1 != 1){
									//Insertion, Left Side
									if(i1.con_loc+k-i1.len > ANCHOR_SIZE && !i1.rc){
				intc2++;	
										i1.id = all_intervals.size();
										all_intervals.push_back(i1);
										sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
										interval_pairs.push_back({false,-1,i1.id});
							 		 	interval_pair_ids[i1.id].push_back(pair_id++);
									}
								}
								else{
									//Insertion, Right Side, singleton case
									if(((int)contig_len - (i1.con_loc+k)) > ANCHOR_SIZE && !i1.rc){
				intc3++;
										i1.id = all_intervals.size();
										all_intervals.push_back(i1);
										sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
										interval_pairs.push_back({false,i1.id,-1});
							 		 	interval_pair_ids[i1.id].push_back(pair_id++);
									}
								}	
							}
						}
						else{
anchor_fails++;
						}
					}
					// If second anchor not set
					else if(a2 == -1){

						// Interval is anchor with the same orientation as the other anchor
						if((is_anchor) && cur_i.rc == intervals[lookup[path[a1]].first][lookup[path[a1]].second].rc){// && !cur_i.rc){

							a2 = i;

							// No interval in between anchors
							if(a1 - a2 == 1){

					 			i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
				 				i2 = intervals[lookup[path[a2]].first][lookup[path[a2]].second];

					 			chromo_dist = i2.loc-(i1.loc+i1.len);
					 			contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

					 			if(contig_dist > 0){

				 					bool allNs = true;

				 					for(int c = i1.con_loc+k+1; c < (i2.con_loc+k-i2.len); c++){
				 						if(!Ns[c]){
				 							allNs = false;
				 							break;
				 						}
				 					}
				 					if(allNs)contig_dist = 0;
				 				}

					 			//Negative Strand (Rare)
					 			if(i1.rc && i2.rc){ 

					 				chromo_dist = i1.loc-(i2.loc+i2.len);

					 				//DUP
					 				if(i2.loc+i2.len-i1.loc >= min_length && i2.loc+i2.len-i1.loc <= max_length && abs(contig_dist) <= uncertainty){
						 				 // only 3 FP
short_dup1++;
						 			}
						 			else{
						 				//INS
										if(abs(chromo_dist) <= uncertainty && contig_dist > uncertainty){
											add_result(j, i1, i2, INS, i2.chr, 0);
short_ins1++;
										}//DEL
										else if(chromo_dist > uncertainty && abs(contig_dist) <= uncertainty){
											add_result(j, i1, i2, DEL, i2.chr, 0);  //+1 TP +2FP, probably lots of redundancy (careful, may inflate support)
intc4++;
						 					i1.id = all_intervals.size();
											all_intervals.push_back(i1);
											i2.id = all_intervals.size();
											all_intervals.push_back(i2);
											sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
											sorted_intervals[i2.chr].push_back({i2.id,i2.loc});
											interval_pairs.push_back({false,i1.id,i2.id});
											interval_pair_ids[i1.id].push_back(pair_id);
								 		 	interval_pair_ids[i2.id].push_back(pair_id++);
										}
										else if(abs(chromo_dist - contig_dist) <= uncertainty && min(abs(chromo_dist), abs(contig_dist)) > uncertainty){// max(min_length, uncertainty+1)){
											add_result(j, i1, i2, INV, i2.chr, 0);
short_inv1++;
										}//???
										else if(chromo_dist > uncertainty && contig_dist > uncertainty){
											//zero
//mystery1++;	
										}//SNP
										else{
//mystery2++;								//FP
										}
						 			}
					 			}
					 			//Positive Strand
					 			else{ 
					 				//DUP
					 				if(i1.loc+i1.len-i2.loc >= min_length && i1.loc+i1.len-i2.loc <= max_length && abs(contig_dist) <= uncertainty){
						 				add_result(j, i1, i2, DUP, i2.chr, 0);
short_dup2++;
						 			}
						 			else{
						 				//INS
										if(abs(chromo_dist) <= uncertainty && contig_dist > uncertainty){
											add_result(j, i1, i2, INS, i2.chr, 0);
short_ins2++;
										}//DEL
										else if(chromo_dist > uncertainty && abs(contig_dist) <= uncertainty){
intc4++;
if(j == CON_NUM_DEBUG)cerr << "Should be here" << endl;
											i1.id = all_intervals.size();
											all_intervals.push_back(i1);
											i2.id = all_intervals.size();
											all_intervals.push_back(i2);
											sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
											sorted_intervals[i2.chr].push_back({i2.id,i2.loc});
											interval_pairs.push_back({false,i1.id,i2.id});
											interval_pair_ids[i1.id].push_back(pair_id);
								 		 	interval_pair_ids[i2.id].push_back(pair_id++);
										}//INV
										else if(abs(chromo_dist - contig_dist) <= uncertainty && min(abs(chromo_dist), abs(contig_dist)) >= max(min_length, uncertainty+1)){ // min_length){
											add_result(j, i1, i2, INV, i2.chr, 0);
short_inv1++;
										}//??? Maybe SNPs/errors near breakpoint
										else if(chromo_dist > uncertainty && contig_dist > uncertainty){
											if(chromo_dist - abs(contig_dist) >= MIN_SV_LEN_DEFAULT){
												add_result(j, i1, i2, DEL, i2.chr, 0);
mystery1++;
											}
										}
										else{
mystery2++;
										}
						 			}
					 			}
							}
							// Interval(s) between anchors, specifically TRANS and INV
							else{

								// Inner non-anchor region
								for(int u = a1; u > a2; u--){

					 				i1 = intervals[lookup[path[u]].first][lookup[path[u]].second];
				 					i2 = intervals[lookup[path[u-1]].first][lookup[path[u-1]].second];
					 				chromo_dist = i2.loc-(i1.loc+i1.len);
					 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

					 				if(i1.rc && i2.rc)chromo_dist = i1.loc-(i2.loc+i2.len);

					 				if(abs(contig_dist) > uncertainty)break;

					 				if((u != a1 && u-1 != a2) && (i1.chr != i2.chr || i1.rc != i2.rc || abs(chromo_dist) > uncertainty))break;

					 				if(u-1 == a2)found = true;
								}

								if(found){

					 				i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
				 					i2 = intervals[lookup[path[a2]].first][lookup[path[a2]].second];

				 					if(intervals[lookup[path[a1-1]].first][lookup[path[a1-1]].second].rc != i1.rc  && intervals[lookup[path[a2+1]].first][lookup[path[a2+1]].second].rc != i2.rc){
				 					//	add_result(j, i1, i2, INV, i2.chr, 0);
short_inv2++;
									}
									else if(i1.chr != i2.chr){
										add_result(j, i1, i2, TRANS, i2.chr, 0);
short_trans++;
									}
								}
								else{
									if(contig_dist > uncertainty && i1.rc == i2.rc && i1.chr == i2.chr){
ugly_fails++;
										add_result(j, i1, i2, DEL, i2.chr, 0);
									}
								}
							}
							a1 = a2;
					 		a2 = -1;
						}
						// Interval is anchor with different orientation from the other anchor (one BP INV)
						else if(is_anchor && cur_i.rc != intervals[lookup[path[a1]].first][lookup[path[a1]].second].rc){ 

							for(int x = a1-1; x >= 0; x--){
								//The fuck is this?
								if(intervals[lookup[path[a1]].first][lookup[path[a1]].second].rc == cur_i.rc)continue;
							}

							a2 = i;

							if(a1 - a2 == 1){

					 			i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
				 				i2 = intervals[lookup[path[a2]].first][lookup[path[a2]].second];

					 			chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
				 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

								if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){
	intc5++;
									i1.id = all_intervals.size();
									all_intervals.push_back(i1);
									i2.id = all_intervals.size();
									all_intervals.push_back(i2);
									sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
									sorted_intervals[i2.chr].push_back({i2.id,i2.loc});
									interval_pairs.push_back({false,i1.id,i2.id});
									interval_pair_ids[i1.id].push_back(pair_id);
						 		 	interval_pair_ids[i2.id].push_back(pair_id++);
								}
					 		}

					 		a1 = a2;
					 		a2 = -1;

						}
						// trailing non-anchor region
						else if(i == 1){

			 				i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
				 			i2 = intervals[lookup[path[a1-1]].first][lookup[path[a1-1]].second];
			 				chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
			 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

							if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){
intc5++;
								i1.id = all_intervals.size();
								all_intervals.push_back(i1);
								i2.id = all_intervals.size();
								all_intervals.push_back(i2);
								sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
								sorted_intervals[i2.chr].push_back({i2.id,i2.loc});
								interval_pairs.push_back({false,i1.id,i2.id});
								interval_pair_ids[i1.id].push_back(pair_id);
					 		 	interval_pair_ids[i2.id].push_back(pair_id++);
							}
						}
						else{
if(!is_anchor)anchor_fails++;
						}

						// Insertion Right Side, after last anchor interval
						if(a1 == 1 || a2 == 1){

							i1 = (a2 == 1) ? intervals[lookup[path[a2]].first][lookup[path[a2]].second] : intervals[lookup[path[a1]].first][lookup[path[a1]].second];
							if(((int)contig_len - (i1.con_loc+k)) > ANCHOR_SIZE && !i1.rc){
intc6++;
								i1.id = all_intervals.size();
								all_intervals.push_back(i1);
								sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
								interval_pairs.push_back({false,i1.id,-1});
					 		 	interval_pair_ids[i1.id].push_back(pair_id++);
							}
						}
					}//anchor set
				}//path

				if(!anchor_found)max_paths++;

			}//paths	
		}//empty


//Debug graph
if(j == CON_NUM_DEBUG){

	for(int i=0; i <= id; i++){
		for(int j=0; j <= id; j++){
			cerr << contig_graph[i][j] << ", ";
		}
		cerr << endl;
	}

	cerr << "All generated paths:" << endl;

	for(auto &path_pair: paths){
		vector<int>& path = path_pair.first;
		for(int i = path.size()-2; i > 0; i--){
			cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].loc << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].loc + intervals[lookup[path[i]].first][lookup[path[i]].second].len << "\t";
		}
		cerr << endl;
		for(int i = path.size()-2; i > 0; i--){
			cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].con_loc+k - intervals[lookup[path[i]].first][lookup[path[i]].second].len << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].con_loc+k << "\t";
		}
		cerr << endl;
	}
}

		for(int i=0; i < interval_count; i++){
			delete[] contig_graph[i];
		}
		delete[] contig_graph;

	}//j


	
	// Long Structural Variants
	//============================================

	cerr << "Predicting long structural variants..." << endl;

	int w, w_id, z, z_id;
	int last_z_loc = 0;
	int last_w_loc = 0;
	int next_w_loc = 0;
	int chromo_dist, contig_dist, contig_dist1, contig_dist2;

	for(char chr=0; chr < 25; chr++){
		if(sorted_intervals[chr].empty())continue;
		sort(sorted_intervals[chr].begin(), sorted_intervals[chr].end());

		last_z_loc = 0;
		last_w_loc = 0;
		next_w_loc = 0;

		for(z=0; z < sorted_intervals[chr].size()-1; z++){  

			z_id = sorted_intervals[chr][z].id;
			found = false;
 
			if(!interval_pair_ids[z_id].empty()){

				w = (z-1 > 0) ? z-1 : 0;
				mapping_ext iz = all_intervals[z_id];
				mapping_ext ix, iy, iw;
				next_w_loc = all_intervals[sorted_intervals[chr][w].id].rc ? sorted_intervals[chr][w].loc : sorted_intervals[chr][w].loc + all_intervals[sorted_intervals[chr][w].id].len;

				while(w >= 0 && next_w_loc >= last_w_loc  && abs(iz.loc - sorted_intervals[chr][w].loc) < max_length){

					w_id = sorted_intervals[chr][w].id;

					if(!interval_pair_ids[w_id].empty()){

						iw = all_intervals[w_id];

						chromo_dist = (iz.rc && iw.rc) ? iw.loc-(iz.loc+iz.len) : iz.loc-(iw.loc+iw.len);

						// Try matching all interval pairs w,x with downstream interval pairs y,z within a user-specified max distance
						for(int x: interval_pair_ids[w_id]){
							for(int y: interval_pair_ids[z_id]){   //Shit, all combinations, TODO: Check if always small. 
								interval_pair& ip1 = interval_pairs[x];
								interval_pair& ip2 = interval_pairs[y];

								// If not visited and compatible 
								if(!ip1.visited && !ip2.visited && ip1.id1 == w_id && ip2.id2 == z_id && ip1.id2 != z_id && ip2.id1 != w_id && ip1.id1 != -1 && ip2.id2 != -1){

									//Both have empty interals on opposite sides, INS
									if(ip1.id2 == -1 && ip2.id1 == -1){
										if(abs(chromo_dist) <= uncertainty){
											ip1.visited = true;
											ip2.visited = true;
											add_result(iw.con_id, iw, iw, INS, iz.chr, add_result(iz.con_id, iz, iz, INS, iz.chr, -1));
											found = true;
long_ins++;
										}
									}
									else if(ip1.id2 != -1 && ip2.id1 != -1){

										ix = all_intervals[ip1.id2];
										iy = all_intervals[ip2.id1];
										contig_dist1 = abs(iw.con_loc - (ix.con_loc-ix.len));
										contig_dist2 = abs(iy.con_loc - (iz.con_loc-iz.len));		

										// Opposite sides are on the same chromosome 
										if(ix.chr == iy.chr && contig_dist1 <= uncertainty && contig_dist2 <= uncertainty){

											// If internal intervals are mapped to the opposite strand of external intervals INV
											if(iw.rc != ix.rc && iy.rc != iz.rc){ //if(ix.rc || iy.rc){

												if(ix.rc && iy.rc){

													if(ix.loc > iy.loc && (ix.loc+ix.len)-iy.loc < max_length  && abs(((ix.loc+ix.len)-iy.loc) - chromo_dist) <= uncertainty){
														ip1.visited = true;
														ip2.visited = true;
														add_result(iw.con_id, iw, ix, INV, iy.chr, add_result(iz.con_id, iy, iz, INV, iz.chr, -1));		
														found = true;												
long_inv1++;
													}
													else if(ix.loc < iy.loc && (iw.loc+iw.len)-iz.loc < max_length && abs(((iw.loc+iw.len)-iz.loc) - (iy.loc-(ix.loc+ix.len))) <= uncertainty){
														ip1.visited = true;
														ip2.visited = true;
														add_result(iw.con_id, iw, ix, INV, iy.chr, add_result(iz.con_id, iy, iz, INV, iz.chr, -1));
														found = true; //This causes some to be missed, 1 inv, 1 trans
long_inv2++;
													}
												}
												else if(iw.rc && iz.rc){
													if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < max_length  && abs(((iy.loc+iy.len)-ix.loc) - chromo_dist) <= uncertainty){
														ip1.visited = true;
														ip2.visited = true;
														add_result(iw.con_id, iw, ix, INV, iy.chr, add_result(iz.con_id, iy, iz, INV, iz.chr, -1));
														found = true;
long_inv3++;
													}
												}
												
											}
											else{
						
												if(iw.rc && ix.rc && iy.rc && iz.rc){
													// Rare, if ever. Ignore for now.
												}
												else{
													if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < max_length && (iw.chr != ix.chr || (ix.loc > iz.loc+iz.len || iy.loc+iy.len < iw.loc))){ //last check = 1 Trans
														ip1.visited = true;
														ip2.visited = true;
														add_result(iw.con_id, iw, ix, TRANS, iy.chr, add_result(iz.con_id, iy, iz, TRANS, iz.chr, -1));
														found = true;
long_trans++;
													}
													//else if(ix.loc > iy.loc && (ix.loc+ix.len)-iy.loc < max_length && abs(iw.loc-iy.loc) <= uncertainty && abs(iz.loc-ix.loc) <= uncertainty){
													else if(ix.loc > iw.loc && iz.loc > iy.loc && abs((iw.loc+iw.len)-ix.loc) < max_length  && abs((iy.loc+iy.len)-iz.loc) < max_length && abs(iw.loc-iy.loc) <= uncertainty && abs(iz.loc-ix.loc) <= uncertainty){
														// ip1.first = true;
														// ip2.first = true;
														// print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "DEL");
														// Does not work well
long_del++;
													}
													else if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < max_length && abs(iz.loc-(iw.loc+iw.len)) <= uncertainty){
														// Mostly false positives?
long_dup++;
													}
												}
											}
										}
									}//ins or not
								}//not used yet and correct side
								if(found)break;
							}//ip2
							if(found)break;
						}//ip1
					}
					if(found){

						int cur_z_loc = iz.rc ? iz.loc + iz.len : iz.loc;

						if(abs(cur_z_loc-last_z_loc) > uncertainty){
							last_w_loc = last_z_loc; //So not really "w" but ok
							last_z_loc = cur_z_loc;
						}
						else{
							last_w_loc = iw.rc ? iw.loc : iw.loc + iw.len;
						}
						break;
					}
					else{
						w--;
					}
				} //w
				if(found)continue;
			}
		}//z
	}


	//Single interval pair case (WARNING: may cause many false positives)
	for(char chr=0; chr < 25; chr++){
		if(sorted_intervals[chr].empty())continue;
		sort(sorted_intervals[chr].begin(), sorted_intervals[chr].end());

		for(w=0; w < sorted_intervals[chr].size()-1; w++){

			w_id = sorted_intervals[chr][w].id;

			if(!interval_pair_ids[w_id].empty()){

				mapping_ext iw, ix;

				//Same as above but with only w,x
				for(int x: interval_pair_ids[w_id]){

					interval_pair& ip1 = interval_pairs[x];

					if(!ip1.visited && (ip1.id1 == -1 || ip1.id2 == -1)){
						int contig_id = (ip1.id1 == -1) ? all_intervals[ip1.id2].con_id : all_intervals[ip1.id1].con_id;
						string contig_seq = PRINT_READS ? all_contigs[contig_id].data : con_string(contig_id, 0, all_compressed_contigs[contig_id].data.size()/2);
						vector<bool>& Ns = all_contig_Ns[contig_id];
						bool Nfix = true;
						for(int i = 0; i < 10; i++){
							if(!Ns[i] && !Ns[contig_seq.size()-i-1])Nfix = false;
						}
						if(Nfix)continue;
					}

					if(!ip1.visited && ip1.id1 != -1 && ip1.id2 != -1){

						iw = all_intervals[ip1.id1];
						ix = all_intervals[ip1.id2];
						contig_dist = abs(iw.con_loc - (ix.con_loc-ix.len));
						chromo_dist = (ix.rc && iw.rc) ? iw.loc-(ix.loc+ix.len) : ix.loc-(iw.loc+iw.len);
						contig_len = PRINT_READS ? all_contigs[ix.con_id].data.length() : all_compressed_contigs[ix.con_id].data.size()/2;

						if(contig_dist <= uncertainty){
							if(( (iw.chr != ix.chr || abs(ix.loc - (iw.loc+iw.len)) > max_length/2) ) ){// && (iw.len + ix.len) > contig_len/4){	 //was: 2,3,3  //iw.rc == ix.rc &&  real events are transversions! 
								add_result(iw.con_id, iw, ix, TRANS, ix.chr, 0);
								ip1.visited = true;
		one_bp_trans++;
							}
							else if(abs(ix.loc - iw.loc) <= max_length && abs(ix.loc - iw.loc) >= min_length){
								if(iw.rc != ix.rc){
									add_result(iw.con_id, iw, ix, INV, ix.chr, 0); //ALL FP come from here
		one_bp_inv++;
								}
								else if(chromo_dist >= min_length){
									add_result(iw.con_id, iw, ix, DEL, ix.chr, 0);
		one_bp_del++;
								}
								else if(iw.loc+iw.len-ix.loc >= min_length){
									//DUP
		one_bp_dup++;
								}
								else{
		mystery3++;
								}
								ip1.visited = true;
							}
						}
					}
					else if(!ip1.visited && ip1.id1 == -1){    //Two more possible with extra FP
						ix = all_intervals[ip1.id2];
						contig_len = PRINT_READS ? all_contigs[ix.con_id].data.length() : all_compressed_contigs[ix.con_id].data.size()/2;

						if(ix.len >= contig_len/3 && (ix.con_loc+k-ix.len) > contig_len/4){
							add_result(ix.con_id, ix, ix, INSL, ix.chr, 0);
		one_bp_ins++;
						}
						ip1.visited = true;
					}
					else if(!ip1.visited && ip1.id2 == -1){
						iw = all_intervals[ip1.id1];
						contig_len = PRINT_READS ? all_contigs[iw.con_id].data.length() : all_compressed_contigs[iw.con_id].data.size()/2;

						if(iw.len >= contig_len/3 && (contig_len-iw.con_loc+k) > contig_len/4){
							add_result(iw.con_id, iw, iw, INSR, iw.chr, 0);
		one_bp_ins++;
						}
						ip1.visited = true;
					}
				}//ip1
			}
		}//i
	}//chr

if(PRINT_STATS){

	cerr << "Used intervals: " << num_intervals << endl;
	cerr << "Normal Edges: " << normal_edges << endl;
	cerr << "INS Edges: " << ins_edges << endl;
	cerr << "Source Edges: " << source_edges << endl;
	cerr << "Sink Edges: " << sink_edges << endl;
	cerr << "Singletons: " << singletons << endl;

	cerr << "Paths: " << path_count << endl;

	cerr << "Short_INV1: " << short_inv1 << endl;
	cerr << "Short_INV2: " << short_inv2 << endl;
	cerr << "Short_INS1: " << short_ins1 << endl;
	cerr << "Short_INS2: " << short_ins2 << endl;
	cerr << "Short_DEL: " << short_del << endl;
	cerr << "Short_DUP1: " << short_dup1 << endl;
	cerr << "Short_DUP2: " << short_dup2 << endl;
	cerr << "Short_DUP3: " << short_dup3 << endl;
	cerr << "Short_TRANS: " << short_trans << endl;

	cerr << "Long_INV1: " <<  long_inv1 << endl;
	cerr << "Long_INV2: " <<  long_inv2 << endl;
	cerr << "Long_INV3: " <<  long_inv3 << endl;
	cerr << "Long_INS: "  <<  long_ins << endl;
	cerr << "Long_DEL: "  <<  long_del << endl;
	cerr << "Long_TRANS: "<<  long_trans << endl;
	cerr << "Long_DUP: "<<  long_dup << endl;

	cerr << "1bp_INV: "  <<  one_bp_inv << endl;
	cerr << "1bp_DEL: "  <<  one_bp_del << endl;
	cerr << "1bp_DUP: "  <<  one_bp_dup << endl;
	cerr << "1bp_TRANS: "<<  one_bp_trans << endl;
	cerr << "1bp_INS: " << one_bp_ins << endl;

	cerr << "Interval_pairs: " << interval_pairs.size() << endl;

	cerr << "Leading: " << intc1 << endl;
	cerr << "INS IP left: " << intc2 << endl;
	cerr << "INS IP right: " << intc3 << endl;
	cerr << "DEL IP: " << intc4 << endl;
	cerr << "Mystery 1: " << mystery1 << endl;
	cerr << "Mystery 2: " << mystery2 << endl;
	cerr << "Mystery 3: " << mystery3 << endl;
	cerr << "Trailing: " << intc5 << endl;
	cerr << "INS IP right: " << intc6 << endl;

	cerr << "Anchor fails: " << anchor_fails << endl;
	cerr << "Ugly Fails: " << ugly_fails << endl;
}

	print_results(fo_vcf, fr_vcf, uncertainty);

	cerr << "Done." << endl;

	fclose(fo_vcf);
	fclose(fr_vcf);

	delete[] sorted_intervals;

}

//Introduction to Algorithms 3rd Edition by Clifford Stein, Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest

/* Returns true if there is a path from source 's' to sink 't' in
  residual graph. Also fills parent[] to store the path */
bool kmistrvar::bfs(const int DEPTH, int** rGraph, int s, int t, int parent[])
{
	// Create a visited array and mark all vertices as not visited
	bool visited[DEPTH];
	memset(visited, 0, sizeof(visited));

	// Create a queue, enqueue source vertex and mark source vertex
	// as visited
	queue <int> q;
	q.push(s);
	visited[s] = true;
	parent[s] = -1;

	// Standard BFS Loop
	while (!q.empty())
	{
		int u = q.front();
		q.pop();

		for (int v=0; v<DEPTH; v++)
		{
			if (visited[v]==false && rGraph[u][v] > 0)
			{
				q.push(v);
				parent[v] = u;
				visited[v] = true;
			}
		}
	}

	// If we reached sink in BFS starting from source, then return
	// true, else false
	return (visited[t] == true);
}
 
// Returns the maximum flow from s to t in the given graph
vector<vector<int>> kmistrvar::fordFulkerson(const int DEPTH, int** rGraph, int s, int t)
{
	int u, v;
	vector<vector<int>> paths;

	// Create a residual graph and fill the residual graph with
	// given capacities in the original graph as residual capacities
	// in residual graph

	int parent[DEPTH];	// This array is filled by BFS and to store path

	int max_flow = 0;	// There is no flow initially

	// Augment the flow while there is path from source to sink
	while (bfs(DEPTH, rGraph, s, t, parent))
	{	
		vector<int> path;
		// Find minimum residual capacity of the edges along the
		// path filled by BFS. Or we can say find the maximum flow
		// through the path found.
		int path_flow = 999999999;
		for (v=t; v!=s; v=parent[v])
		{
			u = parent[v];
			//cerr << v << "\t";
			if(v != t)path.push_back(v);
			path_flow = min(path_flow, rGraph[u][v]);
		}
		//cerr << endl;

		// update residual capacities of the edges and reverse edges
		// along the path
		for (v=t; v != s; v=parent[v])
		{
			u = parent[v];
			rGraph[u][v] -= path_flow;
			rGraph[v][u] += path_flow;
		}

		//if(path.size() > 1)
		paths.push_back(path);

		// Add path flow to overall flow
		max_flow += path_flow;
	}

	// Return the overall flow
	return paths;
}
 
int kmistrvar::minDistance(const int DEPTH, int dist[], bool sptSet[])
{
	// Initialize min value
	int min = 9999999, min_index;

	for (int v = 0; v < DEPTH; v++){
		if (sptSet[v] == false && dist[v] <= min){
		 	min = dist[v];
		 	min_index = v;
		}
	}

	return min_index;
}
 
// Function to print shortest path from source to j
// using parent array
void kmistrvar::getPath(vector<int>& path, int dist[], int parent[], int j)
{

	path.push_back(j);

	// Base Case : If j is source
	if (parent[j]==-1)return;

	getPath(path, dist, parent, parent[j]);
}

pair<vector<pair<int,int>>,int> kmistrvar::findNextPath(const int DEPTH, int** graph, treeNode* node, vector<pair<int,int>> ancestors, int dist[], int parent[]){

	pair<vector<pair<int,int>>,int> min_sub_path;
	pair<vector<pair<int,int>>,int> cur_sub_path;
	pair<int,int> min_edge;
	int min_dist = 99999;
	int cur_dist;
	node->ancestors = ancestors;

	for (int v=node->v; v != -1; v=parent[v]){
		for(int w = 0; w < DEPTH; w++){
			if(graph[w][v] && w != parent[v]){

				if(!node->children.empty() && node->children.find(w) != node->children.end()){
					cur_sub_path = findNextPath(DEPTH, graph, node->children[w], ancestors, dist, parent);
					if(!cur_sub_path.first.empty()){
						cur_dist = ((cur_sub_path.second-dist[parent[v]])+dist[node->v]);
						if(cur_dist < min_dist){
							min_sub_path = cur_sub_path;
							min_sub_path.second = cur_dist;
							min_edge = {parent[v],w};
							min_dist = cur_dist;
						}
					}
				}
				else{
					cur_dist = ((dist[w]-dist[parent[v]])+dist[node->v]);
					if(cur_dist < min_dist){
						min_edge = {parent[v],w};
						min_dist = cur_dist;
					}
				}
			}
		}
	}

	if((node->children.empty() || node->children.find(min_edge.second) == node->children.end()) && min_dist != 99999){
		ancestors.push_back(min_edge);
		min_sub_path.first = ancestors;
		min_sub_path.second = min_dist;
	}

	return min_sub_path;
}
 
// Funtion that implements Dijkstra's single source shortest path
// algorithm for a graph represented using adjacency matrix
// representation
vector<pair<vector<int>, int>> kmistrvar::dijkstra(const int DEPTH, int** graph, int s, int t, int max_paths)
{
	int dist[DEPTH];  // The output array. dist[i] will hold
		   // the shortest distance from src to i

	// sptSet[i] will true if vertex i is included / in shortest
	// path tree or shortest distance from src to i is finalized
	bool sptSet[DEPTH];

	// Parent array to store shortest path tree
	int parent[DEPTH];
	
	vector<pair<vector<int>, int>> paths;

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < DEPTH; i++)
	{
		parent[i] = -1;
		dist[i] = 9999999;
		sptSet[i] = false;
	}
 
    // Distance of source vertex from itself is always 0
	dist[s] = 0;
 
    // Find shortest path for all vertices
	for (int count = 0; count < DEPTH-1; count++)
	{
		// Pick the minimum distance vertex from the set of
		// vertices not yet processed. u is always equal to src
		// in first iteration.
		int u = minDistance(DEPTH, dist, sptSet);

		// Mark the picked vertex as processed
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the
		// picked vertex.
		for (int v = 0; v < DEPTH; v++){

			// Update dist[v] only if is not in sptSet, there is
			// an edge from u to v, and total weight of path from
			// src to v through u is smaller than current value of
			// dist[v]
			if (!sptSet[v] && graph[u][v] && (dist[u] + graph[u][v]) <= dist[v]){
				parent[v]  = u;
				dist[v] = dist[u] + graph[u][v];
			}  
		}
	}

	if(parent[t] != -1){

		vector<int> path;
		vector<treeNode> nodes;
		unordered_map<string, bool> visted;
		int path_tree[max_paths];

		//shortest path
		getPath(path, dist, parent, t);
		paths.push_back({path, dist[t]});//parent[t]]});
		treeNode sink;
		sink.v = t;
		max_paths--;

		while(max_paths > 0){

			vector<int> next_path;
			pair<vector<pair<int,int>>,int> min_sub_path = findNextPath(DEPTH, graph, &sink, sink.ancestors, dist, parent);

			if(min_sub_path.first.empty())break;
		
			int v = t;
			int w;
			treeNode* cur_node = &sink;

//if(cur_debug_id == CON_NUM_DEBUG)cerr << "t: " << t  << " size: " << min_sub_path.first.size() << " replace: " << min_sub_path.first.back().first << " with: " << min_sub_path.first.back().second << " dist: " << min_sub_path.second << endl;

			for(int i = 0; i < min_sub_path.first.size()+1; i++){
				w = (i == min_sub_path.first.size()) ? s : min_sub_path.first[i].first;
				for (v; v != s && v != w; v=parent[v]){
					next_path.push_back(v);
				}
				if(v == s)break;

				if(i == min_sub_path.first.size()-1){
					treeNode new_node;
					new_node.v = min_sub_path.first.back().second;
					new_node.ancestors = min_sub_path.first;
					nodes.push_back(new_node);
					cur_node->children[new_node.v] = &nodes.back();
				}

				cur_node = cur_node->children[min_sub_path.first[i].second];
				v = cur_node->v;
			}
			next_path.push_back(v);

			paths.push_back({next_path, min_sub_path.second});
			max_paths--;
		}
	}
 
	return paths;
}