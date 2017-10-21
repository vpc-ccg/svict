#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include "kmistrvar.h"

using namespace std;
#define strequ(A,B) (strcmp(A.c_str(), B.c_str()) == 0)
kmistrvar::kmistrvar(int kmer_len, const int anchor_len, const string &partition_file, const string &reference, const string &gtf, const bool barcodes) : 
	variant_caller(partition_file, reference), ANCHOR_SIZE(anchor_len), USE_ANNO((gtf != "")), USE_BARCODES(barcodes){

	k = kmer_len;
	num_kmer = (int)pow(4,k);
	kmer_mask = (long long)(num_kmer-1);
	num_intervals = 1;
	contig_kmers_index = new int[num_kmer]; 
	if(USE_ANNO) ensembl_Reader(gtf.c_str(), iso_gene_map, gene_sorted_map );

	contig_mappings = new vector<vector<mapping>>*[2];
	for(int rc=0; rc < 2; rc++){
		contig_mappings[rc] = new vector<vector<mapping>>[25];
	}

	repeat = new vector<short>*[2];
	for(int rc=0; rc < 2; rc++){
		repeat[rc] = new vector<short>[25];
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
		delete[] repeat[rc];
	}
	delete[] repeat;

}

void kmistrvar::init(){

	for(int i = 0; i < num_kmer; i++){
		contig_kmers_index[i] = 0;
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

pair<int,pair<int,int>> kmistrvar::compute_support(contig contig, int start, int end){

	string con = contig.data;
	pair<int,pair<int,int>> result;
	int max = 0;
	int min = 0;
	int sum = 0;
	int coverage[con.length()];

	for (int k=0; k < con.length(); k++){
		coverage[k]=0; 
	}

	for (int j = 0; j < contig.support(); j++){
		for (int k=0; k < contig.read_information[j].seq.length(); k++){
			coverage[k+contig.read_information[j].location_in_contig]++;
		}
	}

	for (int k = start; k < end; k++)
	{
		if (coverage[k]>max){
			max = coverage[k];
		}
		if (coverage[k]<min){
			min = coverage[k];
		}
		sum += coverage[k];
	}

	result = {sum/con.length(),{min, max}};

	return result;
}

vector<pair<pair<string, string>, int>> kmistrvar::correct_reads(vector<pair<pair<string, string>, int>> reads){

	unordered_map<string, unordered_map<string, int>> consensus;
	vector<pair<pair<string, string>, int>> corrected_reads;
	int max_count = 0;
	int len = 0;
	string c = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";  //quick and dirty fix TODO
	string g = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
	string max_seq = "";

	for(auto & read : reads){
		len = read.first.second.length() < 50 ? read.first.second.length() : 50;
		if(read.first.second.substr(0,len) == c.substr(0,len) ||
			read.first.second.substr(read.first.second.length()-len,len) == c.substr(0,len) ||
			read.first.second.substr(0,len) == g.substr(0,len) ||
			read.first.second.substr(read.first.second.length()-len,len) == g.substr(0,len))continue; 
		consensus[read.first.first.substr((read.first.first.length()-20), 20)][read.first.second]++;
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

		corrected_reads.push_back({{barcode.first, max_seq}, 0});
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
	cout << label << "\tRef_Loc: " << interval.loc << "-" << (interval.loc+interval.len) << "\tCon_Loc: " << ((interval.con_loc+k)-interval.len) << "-" << (interval.con_loc+k) 
		 << "\tLen:" << interval.len << "\tChr: " << chromos[interval.chr] << "\tRC: " << interval.rc << "\tCon_ID: " << interval.con_id << endl; 
}

void kmistrvar::print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id, mapping_ext m1, mapping_ext m2, string type){

	contig con = all_contigs[id];
	int loc = m1.rc ? m1.loc : (m1.loc + m1.len);
	int end = m2.rc ? (m2.loc + m2.len) : m2.loc;
	int chromo_dist = m2.loc - (m1.loc + m1.len);	
	int contig_dist = (m2.con_loc-m2.len)-m1.con_loc; 
	int contig_loc = con.cluster_loc;
	int contig_sup = con.support();
	int pos = 0;
	float score;
	char strand = '+';
	string ref, alt;
	string info = "";
	string contig_seq = con.data;
	char contig_chr = con.cluster_chr;
	string best_gene1 = "", best_trans1 = "", context1 = "";
	string best_gene2 = "", best_trans2 = "", context2 = "";
	string best_gene3 = "", best_trans3 = "", context3 = "";
	vector<uint32_t> vec_best1, vec_best2, vec_best3;

	if(m1.rc && m2.rc)strand = '-';

	try{
		if(type == "INS"){
			score = abs(chromo_dist); //wrong with negative strand
			ref = contig_seq.substr(m1.con_loc+k-1, 1);
			alt = contig_seq.substr(m1.con_loc+k-1, contig_dist+1);
		}
		else if(type == "INSL"){
			type = "INS";
			score = abs(chromo_dist);
			ref = "N";
			alt = contig_seq.substr(0, (m1.con_loc+k-m1.len));
			loc = end;
		}
		else if(type == "INSR"){
			type = "INS";
			score = abs(chromo_dist);
			ref = contig_seq.substr(m1.con_loc+k-1, 1);
			alt = contig_seq.substr(m1.con_loc+k-1, (contig_seq.length()-(m1.con_loc+k-1)+1));
			end = loc;
		}
		else if(type == "DEL"){
			score = abs(contig_dist);
			ref = "<DEL>";
			alt = contig_seq.substr(m1.con_loc+k-1, 1);
		}
		else if(type == "DUP"){

			score = abs(chromo_dist);
			info = "Interspersed";

			if(strand == '+'){
				if(m1.loc < m2.loc && (m1.loc + m1.len)-m2.loc > MIN_SV_LEN_DEFAULT){
					info = "Tandem";
					loc = m2.loc;
					end = (m1.loc + m1.len);
					int len = ((m1.loc + m1.len)-m2.loc)+1;
					if(m1.con_loc+k-1-len >= 0){
						ref = contig_seq.substr(m1.con_loc+k-1-len, 1);
						alt = contig_seq.substr(m1.con_loc+k-1-len, len+1);
					}
					else{
						ref = "N";
						alt = contig_seq.substr(m1.con_loc+k-len, len);
					}
				}
				else if(m1.loc > m2.loc && (m1.loc + m1.len) > (m2.loc + m2.len)){
					info = "Tandem";
					ref = contig_seq.substr(m1.con_loc+k-1, 1);
					alt = "<DUP>";
				}
				else{
					info = "Interspersed"; //TODO we have the alt sequence
					ref = contig_seq.substr(m1.con_loc+k-1, 1);
					alt = "<DUP>";
				}
			}
			else{
				info = "Negative";
				ref = "N";
				alt = "<DUP>"; //TODO
			}
		}
		else if(type == "INV"){
			score = abs(chromo_dist);
			ref = contig_seq.substr(m1.con_loc+k-1, 1);
			alt = "<INV>";
			info = "OneBP";

			if(m1.rc){
				ref = "N";
			}
			else if(!m2.rc){
				alt = contig_seq.substr(m1.con_loc+k-1, contig_dist+1);
				info = "Micro";
			}
		}
		else if(type == "TRANS"){

			score = abs(chromo_dist);
			bool is_anchor1 = ((m1.chr == contig_chr) && (abs(m1.loc - contig_loc) < MAX_ASSEMBLY_RANGE));
			bool is_anchor2 = ((m2.chr == contig_chr) && (abs(m2.loc - contig_loc) < MAX_ASSEMBLY_RANGE));

			if(!is_anchor1){
				ref = "N";
				alt = "<TRANS>";
			}
			else if(!is_anchor2){
				ref = contig_seq.substr(m1.con_loc+k-1, 1);
				alt = "<TRANS>";
			}
			else{
				ref = contig_seq.substr(m1.con_loc+k-1, 1);
				alt = contig_seq.substr(m1.con_loc+k-1, contig_dist+1);
			}
		}
		else if(type == "ALL"){
			//FOR TESTING ONLY
			score = 1;
			ref = "N";
			alt = "N";
		}
	}
	catch(const exception& e){
		print_interval("m1",m1);
		print_interval("m2",m2);
		cerr << "Invalid call of type: " << type << ", " << info << ", skipping..." << endl;
		cerr << e.what() << endl;
		return;
	}

	if(loc > end){
		if(type == "INV" || type == "DEL"){
			int tmp = loc;
			loc = end;
			end = tmp;
		}
	}

	score = -10*log(score);
	pair<int,pair<int,int>> support = compute_support(con, min(m1.con_loc, m2.con_loc-m2.len+k), max(m1.con_loc, m2.con_loc-m2.len+k));
	contig_sup = support.second.second;
	//if(contig_sup > 50)return;

	if(type == "TRANS"){

		if(USE_ANNO){
			locate_interval(chromos[m1.chr], loc, loc, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(chromos[m2.chr], m2.loc, (m2.loc + m2.len), gene_sorted_map[chromos[m2.chr]], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];

			if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;Source=%s,%d,%d,+\n", 
			chromos[m1.chr].c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
			info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), chromos[m2.chr].c_str(), m2.loc, (m2.loc + m2.len));
	}
	else{

		if(USE_ANNO){
			locate_interval(chromos[m1.chr], loc, loc, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(chromos[m2.chr], end, end, gene_sorted_map[chromos[m2.chr]], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];
			locate_interval(chromos[m1.chr], loc, end, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene3, best_trans3, vec_best3);
			context3 = contexts[vec_best3[0]];

			if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		if(info == "Interspersed"){
			fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s;Source=%s,%d,%d,+\n", 
				chromos[m1.chr].c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
				info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str(), chromos[m2.chr].c_str(), m2.loc, (m2.loc + m2.len));
		}
		else{
			fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s\n", 
				chromos[m1.chr].c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
				info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str());
		}
	}
	 
	fprintf(fo_full, "%s\t%s\t%d\t%d\t%s\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\n", 
		type.c_str(), chromos[m1.chr].c_str(), m1.loc,  (m1.loc + m1.len), chromos[m2.chr].c_str(), m2.loc, (m2.loc + m2.len), id, 
		chromos[m1.chr].c_str(), (m1.con_loc-(m1.len-k)), (m1.con_loc+k), chromos[m2.chr].c_str(), (m2.con_loc-(m2.len-k)), (m2.con_loc+k), contig_sup, info.c_str());

	if(PRINT_READS){
		fprintf(fr_vcf, ">Cluster: %d Contig: %d MaxSupport: %d BP: %d Reads: \n", contig_loc, id, contig_sup, (m1.con_loc+k-1));
		fprintf(fr_vcf, "ContigSeq: %s\n", contig_seq.c_str());
		for (auto &read: con.read_information)
			fprintf(fr_vcf, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());
	}
} 

void kmistrvar::print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id1, int id2, mapping_ext m1, mapping_ext m2, mapping_ext m3, mapping_ext m4, string type){

	contig con1 = all_contigs[id1];
	contig con2 = all_contigs[id2];
	int loc1 = m1.loc + m1.len;
	int loc2 = m2.loc;
	int end1 = m4.loc;
	int end2 = m3.loc+m3.len;
	int chromo_dist = m4.loc - (m1.loc + m1.len);	
	int contig_loc1 = con1.cluster_loc;
	int contig_loc2 = con2.cluster_loc;
	int contig_sup1 = con1.support();
	int contig_sup2 = con2.support();
	float score;
	string ref, ref2, alt;
	string contig_seq1 = con1.data;
	string contig_seq2 = con2.data;
	string best_gene1 = "", best_trans1 = "", context1 = "";
	string best_gene2 = "", best_trans2 = "", context2 = "";
	string best_gene3 = "", best_trans3 = "", context3 = "";
	vector<uint32_t> vec_best1, vec_best2, vec_best3;

	if(type != "INV" && (m1.rc && m2.rc && m3.rc && m4.rc)){
		loc1 = m2.loc + m2.len;
		loc2 = m1.loc;
		end1 = m3.loc;
		end2 = m4.loc+m4.len;
	}

	try{
		if(type == "INS"){
			score = abs(chromo_dist);
			ref = contig_seq1.substr(m1.con_loc+k-1, 1);
			alt = "<INS>";
			end2 = end1;
		}
		else if(type == "DEL"){
			score = 1;
			ref = "<DEL>";
			alt = contig_seq1.substr(m1.con_loc+k-1, 1);
			end2 = loc2;
		}
		else if(type == "DUP"){
			score = abs(chromo_dist);
			ref = contig_seq1.substr(m1.con_loc+k-1, 1);
			alt = "<DUP>";
		}
		else if(type == "INV"){
			score = abs(chromo_dist-(end2-loc2));
			ref = contig_seq1.substr(m1.con_loc+k-1, 1);
			alt = "<INV>";
			if(m2.rc && m3.rc){
				end2 = end1;
			}
			else if(m1.rc && m4.rc){
				loc1 = m2.loc;
				end2 = m3.loc + m3.len;
			}
		}
		else if(type == "TRANS"){
			score = abs(chromo_dist);
			ref = contig_seq1.substr(m1.con_loc+k-1, 1);
			ref2 = contig_seq2.substr(m4.con_loc+k-1-m4.len, 1);
			alt = "<TRANS>";
		}
	}
	catch(const exception& e){
		cerr << "Invalid call of type: " << type << ", skipping..." << endl;
		cerr << e.what() << endl;
		return;
	}

	if(loc1 > end2 && type == "INV"){
		int tmp = loc1;
		loc1 = end2;
		end2 = tmp;
	}

	score = -10*log(score);
	pair<int,pair<int,int>> support = compute_support(con1, min(m1.con_loc, m2.con_loc-m2.len+k), max(m1.con_loc, m2.con_loc-m2.len+k));// m1.con_loc, m2.con_loc-m2.len+k);
	contig_sup1 = support.second.second;
	support = compute_support(con2, min(m3.con_loc, m4.con_loc-m4.len+k), max(m3.con_loc, m4.con_loc-m4.len+k));// m3.con_loc, m4.con_loc-m4.len+k);
	contig_sup2 = support.second.second;
	// if(contig_sup1 > 50 || contig_sup2 > 50)return;

	if(type == "TRANS"){

		if(USE_ANNO){
			locate_interval(chromos[m1.chr], loc1, end1, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(chromos[m4.chr], loc2, end2, gene_sorted_map[chromos[m2.chr]], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];

			if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t%s[%s:%d[\t%f\tPASS\tSVTYPE=BND;MATEID=bnd_%d;Cluster=%d;Contig=%d;Support=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s\n", 
			chromos[m1.chr].c_str(), loc1, m1.id, ref.c_str(), ref.c_str(), chromos[m2.chr].c_str(), loc2, score, m4.id, contig_loc1, id1, contig_sup1,
			context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str()); 
		fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t]%s:%d]%s\t%f\tPASS\tSVTYPE=BND;MATEID=bnd_%d;Cluster=%d;Contig=%d;Support=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s\n", 
			chromos[m4.chr].c_str(), end1, m4.id, ref2.c_str(), chromos[m2.chr].c_str(), end2, ref2.c_str(), score, m1.id, contig_loc2, id2, contig_sup2,
			context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str()); 
	}
	else{

		if(USE_ANNO){
			locate_interval(chromos[m1.chr], loc1, loc1, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(chromos[m4.chr], end2, end2, gene_sorted_map[chromos[m2.chr]], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];
			locate_interval(chromos[m1.chr], loc1, end2, gene_sorted_map[chromos[m1.chr]], 0, iso_gene_map, best_gene3, best_trans3, vec_best3);
			context3 = contexts[vec_best3[0]];

			if(vec_best1[0] > 0 && vec_best2[0] > 0 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster1=%d;Contig1=%d;Support1=%d;Cluster2=%d;Contig2=%d;Support2=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s\n", 
			chromos[m1.chr].c_str(), loc1, ref.c_str(), alt.c_str(), score, type.c_str(), end2, contig_loc1, id1, contig_sup1, contig_loc2, id2, contig_sup2,
			 context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str());
	}

	fprintf(fo_full, "%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n", 
		type.c_str(), chromos[m1.chr].c_str(), m1.loc,  (m1.loc + m1.len), chromos[m4.chr].c_str(), m4.loc, (m4.loc + m4.len),
		type.c_str(), chromos[m2.chr].c_str(), m2.loc,  (m2.loc + m2.len), chromos[m3.chr].c_str(), m3.loc, (m3.loc + m3.len), 
		id1, chromos[m1.chr].c_str(), (m1.con_loc-(m1.len-k)), (m1.con_loc+k), chromos[m4.chr].c_str(), (m4.con_loc-(m4.len-k)), (m4.con_loc+k), contig_sup1,
		id2, chromos[m2.chr].c_str(), (m2.con_loc-(m2.len-k)), (m2.con_loc+k), chromos[m3.chr].c_str(), (m3.con_loc-(m3.len-k)), (m3.con_loc+k), contig_sup2);

	if(PRINT_READS){
		fprintf(fr_vcf, ">Cluster: %d Contig: %d MaxSupport: %d Reads: \n", contig_loc1, id1, contig_sup1);
		fprintf(fr_vcf, "ContigSeq: %s\n", contig_seq1.c_str());
		for (auto &read: con1.read_information)
			fprintf(fr_vcf, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());
		fprintf(fr_vcf, ">Cluster: %d Contig: %d MaxSupport: %d Reads: \n", contig_loc2, id2, contig_sup2);
		fprintf(fr_vcf, "ContigSeq: %s\n", contig_seq2.c_str());
		for (auto &read: con2.read_information)
			fprintf(fr_vcf, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());
	}
} 







//======================================================
// Main Functions
//======================================================


void kmistrvar::run_kmistrvar(const string &range, const string &out_vcf, const string &out_full, int min_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE, int ref_flank)
{


clock_t begin = clock();
clock_t end;
double elapsed_secs;

	assemble(range, min_support, LOCAL_MODE, ref_flank);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "Assembly Time: " << elapsed_secs << endl;
begin = clock();
}

	generate_intervals(LOCAL_MODE);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "Mapping Time: " << elapsed_secs << endl;
begin = clock();
}
	predict_variants(out_vcf, out_full, uncertainty, min_length, max_length);

if(PRINT_STATS){
end = clock();
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cerr << "Calling Time: " << elapsed_secs << endl;
}

}


//======================================================
// Local Assembly and Contig Indexing Stage
//======================================================

void kmistrvar::assemble(const string &range, int min_support, const bool LOCAL_MODE, int ref_flank)
{
	
	int contig_id = 0;
	int cur_test = 0;
	int max_hits = 0;
	int loc, i, j, x;
	long long index = 0;
	string con;
	string cluster_chr;
	vector<int> contig_list;
	contig_kmers.push_back(contig_list);
	contig_kmer_counts.push_back(0);

//STATS
double part_count = 0;
double contig_count1 = 0;
double contig_count2 = 0;
double support_count = 0;

	// If local mode enabled, use regions around contigs. Otherwise, add all chromosomes as regions
	if(!LOCAL_MODE){  
		for(char chr=0; chr < chromos.size(); chr++){//string & chr : chromos){
			regions.push_back({chr,{1, 300000000}});  //WARNING: This will not work with all reference files. Same with old version. 
		}
	}

	//Assemble contigs and build kmer index
	//=====================================
	cerr << "Assembling contigs..." << endl;
	while (1) {
		
		vector<pair<pair<string, string>, int>> p;
		vector<contig> contigs;

		//Read in partition file
		p = pt.read_partition(part_file, range, min_support, MAX_READS_PER_PART);

		if(USE_BARCODES)p = correct_reads(p);

		if (!p.size()) 
			break;

		//Assemble contigs
		contigs = as.assemble(p); 

//STATS
if(PRINT_STATS){
	part_count++;
	contig_count1 += contigs.size(); 	
}
		
		//TODO: Contig merging
		for (auto &contig: contigs){ 
			if (contig.support() >= min_support) {

//STATS
if(PRINT_STATS){
	contig_count2++; 
	support_count += contig.support();	

	if(pt.get_start() == 13322962){
		cout << "support: " << contig.support() << " " << contig.data.length() << " " << pt.get_start() << " " << contig_id << endl;
		cout << contig.data << endl;
	}
}
				
				contig.cluster_chr = find(chromos.begin(), chromos.end(), pt.get_reference()) - chromos.begin(); //pt.get_reference();
				contig.cluster_loc = pt.get_start();
				if(LOCAL_MODE)regions.push_back({contig.cluster_chr, {(pt.get_start()-ref_flank), (pt.get_end()+ref_flank)}});
				
				con = contig.data;
				unordered_map<long long, vector<int>> kmer_location;
				kmer_locations.push_back(kmer_location);
				all_contigs.push_back(contig);
				i = 0, j = 0, x = 0;

				// Per chromosome repeat flag for each contig
				for(int rc=0; rc < 2; rc++){
					for(int chr=0; chr < 25; chr++){
						repeat[rc][chr].push_back(-1);
					}
				}

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

					// Add to index
					if(!contig_kmers_index[index]){
						vector<int> contig_list = vector<int>(1);
						contig_kmers_index[index] = contig_kmers.size();
						contig_list.push_back(contig_id);
						contig_kmers.push_back(contig_list);
//contig_kmer_counts.push_back(0);
					}
					else{
						if(contig_kmers[contig_kmers_index[index]].back() != contig_id){
							contig_kmers[contig_kmers_index[index]].push_back(contig_id);
							if(contig_kmers[contig_kmers_index[index]].size() > max_hits)max_hits = contig_kmers[contig_kmers_index[index]].size();
						}
					}
					kmer_locations[contig_id][index].push_back(i); //TODO get rid of this. Use the lookup table dummy.

				}
				contig_id++;
			}
		}
	}

	//sketchy
	//contig_kmers_index[178956970] = 0;
	//contig_kmers_index[0] = 0;

	if(all_contigs.empty()){
		cerr << "No contigs could be assembled. Exiting..." << endl;
		exit(1);
	}

	// Initialize mappings vector based on number of contigs
	for(int rc=0; rc <= 1; rc++){
		for(int chr=0; chr < 25; chr++){
			contig_mappings[rc][chr] = vector<vector<mapping>>(contig_id, vector<mapping>(1, {0, k, -1}));
		}
	}

//STATS
if(PRINT_STATS){
	cerr << "Partition Count: " << part_count << endl;
	cerr << "Average Num Contigs Pre-Filter: " << (contig_count1/part_count) << endl;
	cerr << "Average Num Contigs Post-Filter: " << (contig_count2/part_count) << endl;
	cerr << "Average Contig Support: " << (support_count/contig_count2) << endl;
	cerr << "Number of unique kmers: " << contig_kmers.size() << endl;
	cerr << "Max hits: " << max_hits << endl;
}
	
}


//======================================================
// Interval Mapping Stage
//======================================================

void kmistrvar::generate_intervals(const bool LOCAL_MODE)
{

	//Read reference and map contig kmers
	//=====================================
	cerr << "Generating intervals..." << endl;

	int use_rc = 2; //1 no, 2 yes
	//int chromo = 0;
	int jj = 0;
	int id = 1;  //zero reserved for sink
	long long index = 0;
	long long rc_index = 0;
	long long cur_index = 0;
	long long cur_index2 = 0;
	long long index_mask = 0;
	int shift = 2*(k-1)-1;
	int i, ii, ii1, ik, iik1, j, js, je, x, y, z, rc, len, ilen;
	int MAX_CONTIG_LEN = 0; 
	int MAX_INTERVAL_LEN = 20000;
	int num_contigs = all_contigs.size();
	int gen_start, gen_end, interval_len, contig_len;
	bool con_repeat;

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

	vector<pair<int,int>> starts;
	string gen;
	char chromo = -1;
	char last_chromo = chromo;


	// Initialization 
	for(auto& contig : all_contigs){
		if(contig.data.length() > MAX_CONTIG_LEN) MAX_CONTIG_LEN = contig.data.length();
	}
	MAX_CONTIG_LEN += (k+1); //Rare case of SNP near the end of the longest contig. 
	MAX_INTERVAL_LEN = MAX_CONTIG_LEN+100;
	int ref_counts[MAX_INTERVAL_LEN];

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
		for(x = 0; x < k-1; x++){    
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
			rc_index |= (((gen[(i+k-1)-x] & MASK) ^ MASK_RC) >> 1);
		}

		// Skip N nucleotides
		for (i; i < gen_end; i++){ //-ANCHOR_SIZE
			if(gen[i+k-1] == 'N'){
				i+=k;
				while(gen[i] == 'N'){
					i++;
				}

				index = 0;
				rc_index = 0;

				for(x = 0; x < k-1; x++){
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
					rc_index |= (((gen[(i+k-1)-x] & MASK) ^ MASK_RC) >> 1);
				}
				i--;
				continue;
			}

			// Precompute some arithmetic 
			// index_mask = 2nd and 3rd bit of char, which is a unique value for each nucleotide 
			ii = i + gen_start;// + 1;
			ii1 = ii-1;
			iik1 = ii1+k;
			index_mask = (gen[i+k-1] & MASK);

			// Forward and reverse strand
			for(rc=0; rc < use_rc; rc++){

				// Convert kmer string to index
				// Bit shift (reuse last kmer) -> shift mask to 1st and 2nd bit -> mask bits for char of last kmer
				// Slightly modified for reverse compliment to flip bits
				if(rc){
					rc_index = ((rc_index >> 2) | ((index_mask ^ MASK_RC) << shift));
					if(!contig_kmers_index[rc_index])continue;
					cur_index = rc_index;

					// cur_index2 = 0;

					// for(x = 0; x < k; x++){    
					// 	cur_index2 <<= 2;
					// 	cur_index2 |= (((gen[(i+ANCHOR_SIZE-1)-x] & MASK) ^ MASK_RC) >> 1);
					// }
				}
				else{
					index = (((index << 2) | (index_mask >> 1)) & kmer_mask);
					if(!contig_kmers_index[index])continue;
					cur_index = index;

					// cur_index2 = 0;

					// for(x = 0; x < k; x++){    
					// 	cur_index2 <<= 2;
					// 	cur_index2 |= ((gen[i+ANCHOR_SIZE-k+x] & MASK) >> 1);
					// }
				}

kmer_hits++;
//contig_kmer_counts[contig_kmers_index[cur_index]]++; //178956970

				for(const int& j : contig_kmers[contig_kmers_index[cur_index]]){

					// Local mode temporarily disabled
					//if( (!LOCAL_MODE || j == jj)){ 
kmer_hits_con++;

						// Get last interval for this contig
 						vector<mapping>& intervals = contig_mappings[rc][chromo][j];
 						
 						// Start new interval if empty
						// if(intervals.empty()){
						// 	if(!contig_kmers_index[cur_index2])continue;
						// 	intervals.push_back((mapping){ii, k, -1});
						// 	continue;
						// }

						mapping& last_interval = intervals.back();

						// Else end last interval and start new one
						//if((ii-(last_interval.loc + last_interval.len) > 1)){// && contig_kmers_index[cur_index2]){ //-1 dup, -1 trans
						if((ii1 > (last_interval.loc + last_interval.len))){

							// Check if sufficiently long
							if(last_interval.len < ANCHOR_SIZE){
								intervals.pop_back();
pop_back1++;
							}
							else if(repeat[rc][chromo][j] > -1){

								//NOTE: This is not perfect since if the new minimum is large, later intervals larger than one of the initial 4 but less than the minimum, would be skipped.
								if((chromo != all_contigs[j].cluster_chr || (abs(last_interval.loc+(last_interval.len/2) - all_contigs[j].cluster_loc) > (MAX_ASSEMBLY_RANGE*2))) || intervals.size() > REPEAT_LIMIT2){

									if(last_interval.len > intervals[repeat[rc][chromo][j]].len){
										intervals[repeat[rc][chromo][j]].loc = last_interval.loc;
										intervals[repeat[rc][chromo][j]].len = last_interval.len;
									}
									intervals.pop_back();
pop_back2++;
								}
								else{
									if(last_interval.len < intervals[repeat[rc][chromo][j]].len){
										repeat[rc][chromo][j] = intervals.size()-1;
									}
								}
							}

							// Set repeat flag
							if(intervals.size() > REPEAT_LIMIT1 && repeat[rc][chromo][j] == -1){
								repeat[rc][chromo][j] = 0;
								for(x = 1; x < intervals.size(); x++){
									if(intervals[x].len < intervals[repeat[rc][chromo][j]].len){
										repeat[rc][chromo][j] = x;
									}
								}
							}

							// Create new interval
							intervals.push_back((mapping){ii, k, -1}); 
						}
						// Extend interval by one
						else if(last_interval.loc == iik1-last_interval.len){  
							last_interval.len++;
						}
						// Extend interval by k+1 if a SNP/SNV is detected
						else if((ii1 == (last_interval.loc + last_interval.len))){ 
							last_interval.len += (k+1);
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

				if(contig_mappings[rc][chromo][j].back().len < ANCHOR_SIZE){
pop_back4++;
					contig_mappings[rc][chromo][j].pop_back();
					if(contig_mappings[rc][chromo][j].empty())continue;
				}
		 
				vector<mapping> valid_intervals;
				string& contig_seq = all_contigs[j].data;
				contig_len = (contig_seq.length()+k+1);

if(contig_mappings[rc][chromo][j].size() > max_ints)max_ints = contig_mappings[rc][chromo][j].size();

if(contig_len > max_con_len)max_con_len = contig_len;

				for(auto &interval: contig_mappings[rc][chromo][j]){

					if(interval.len+k+1 >= MAX_INTERVAL_LEN)continue;

//STATS
if(PRINT_STATS){
int& contig_loc = all_contigs[j].cluster_loc;
char& contig_chr = all_contigs[j].cluster_chr;
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
					len = interval.len-k+1;
					interval_len = interval.len+k+1; //maybe off by one						

					for(x=0; x < interval_len; x++){
					 	ref_counts[x] = 0;
						for(y=0; y < contig_len; y++){
							valid_mappings[x][y] = false;
						}
					}

					cur_index = 0;

					if(rc){
						// k-1 kmer for first round
						for(x = 0; x < k-1; x++){    
							cur_index <<= 2;
							cur_index |= (((gen[(ii+interval.len-1)-x] & MASK) ^ MASK_RC) >> 1);
						}
					}
					else{
						// k-1 kmer for first round
						for(x = 0; x < k-1; x++){    
							cur_index <<= 2;
							cur_index |= ((gen[ii+x] & MASK) >> 1);
						}
					}

					// Traverse the interval
					for (i = 0; i < len; i++){

						if(rc){
							cur_index = (((cur_index << 2) | (((gen[ii+(interval.len-1-i)-k+1] & MASK) ^ MASK_RC) >> 1)) & kmer_mask);
						}
						else{
							cur_index = (((cur_index << 2) | ((gen[ii+i+k-1] & MASK) >> 1)) & kmer_mask);
						}


						if(!contig_kmers_index[cur_index])continue;
						if(binary_search(contig_kmers[contig_kmers_index[cur_index]].begin(), contig_kmers[contig_kmers_index[cur_index]].end(), j)){
							//ref_counts[x] = 0;

								// Build diagonals
							for(auto &loc: kmer_locations[j][cur_index]){

								valid_mappings[i][loc] = true;
								valid_mappings[i+1][loc+1] = false;
								valid_mappings[i+k+1][loc+k+1] = false;
								ref_counts[i]++;

								// Record starts when interval no longer consecutive
								if(i == 0 || loc == 0){
									starts.push_back({i,loc});
								}
								else if(i <= k+1 && loc <= k+1 && !valid_mappings[i-1][loc-1]){
									starts.push_back({i,loc});
								}
								else if(i > k+1 && loc > k+1 && (!valid_mappings[i-1][loc-1] && !valid_mappings[i-k-2][loc-k-2])){
									starts.push_back({i,loc});
								}
							}	
						}
					}

if(starts.size() > max_starts)max_starts = starts.size();

					//if(starts.size() > CON_REPEAT_LIMIT)continue;

					// Check each start point
					for(auto& start : starts){

						x = start.first;
						y = start.second;
						con_repeat = false;

						// Traverse diagonal
						while(x < interval.len && y < contig_seq.length() && valid_mappings[x++][y++]){
							if(y+k < contig_seq.length()){
								if(!valid_mappings[x][y] && valid_mappings[x+k][y+k]){
									x +=k;
									y +=k;
								}
								else if(!valid_mappings[x][y] && contig_seq[y+k] == 'N'){
									valid_mappings[x][y] = true;
								}
							}
						}
						x-=2;
						y-=2;

						ilen = (y - start.second)+k;


						//TODO: Use this table to identify duplications, although maybe no point since all can be found without this


						if(ilen > ANCHOR_SIZE){

							// Ignore repetitive sequences within contig
							for(z = start.first; z <= x; z++){
								if(ref_counts[z] > CON_REPEAT_LIMIT){
									con_repeat = true;
								}
								else{
									con_repeat = false;
									break;
								}
							}

							//Add to final list of intervals
							if(!con_repeat){
sub_count++;
								mapping sub_interval;
								sub_interval.loc = interval.loc + start.first;
								sub_interval.len = ilen;
								sub_interval.con_loc = y;
								valid_intervals.push_back(sub_interval); 
							}
if(sub_count > max_sub)max_sub = sub_count;
						}
					}
				}//intervals

//STATS
interval_count += valid_intervals.size(); 


				num_intervals += valid_intervals.size();
				contig_mappings[rc][chromo][j] = move(valid_intervals);
			}
		}


		jj++;

	}//regions

//STATS
if(PRINT_STATS){
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

	// for(int& i : contig_kmer_counts){
	// 	cout << i << endl;
	// }

}

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		delete[] valid_mappings[i];
	}
	delete[] valid_mappings;
}


//======================================================
// Predict Structural Variants Stage
//======================================================

void kmistrvar::predict_variants(const string &out_vcf, const string &out_full, int uncertainty, int min_length, int max_length)
{
	//Predict structural variants
	//=====================================
	cerr << "Predicting structural variants..." << endl;
	cerr << "Building Graph..." << endl;

//STATS
int anchor_intervals = 0;
int non_anchor_intervals = 0;
int normal_edges = 0, ins_edges = 0, source_edges = 0, sink_edges = 0, singletons = 0;
int path_count = 0;
int intc1 = 0, intc2 = 0, intc3 = 0, intc4 = 0, intc5 = 0, intc6 = 0;
int anchor_fails = 0, ugly_fails = 0;
int short_dup1 = 0, short_dup2 = 0, short_dup3 = 0, short_del = 0, short_inv1 = 0, short_inv2 = 0, short_ins1 = 0, short_ins2 = 0, short_trans = 0, mystery1 = 0, mystery2 = 0, mystery3 = 0;
int long_del = 0, long_inv1 = 0, long_inv2 = 0, long_ins = 0, long_trans = 0, long_dup = 0;
int one_bp_del = 0, one_bp_dup = 0, one_bp_inv = 0, one_bp_ins = 0, one_bp_trans = 0;


if(PRINT_STATS){
	cerr << "Intervals: " << num_intervals << endl;
	cerr << "Contigs: " << all_contigs.size() << endl;
	num_intervals = 0;
	// for(auto & interval : all_intervals){
	// 	int contig_loc = all_contigs[interval.con_id].cluster_loc;
	// 	char contig_chr = all_contigs[interval.con_id].cluster_chr;
	// 	bool is_anchor = ((interval.chr == contig_chr) && (abs(interval.loc+(interval.len/2) - contig_loc) < (MAX_ASSEMBLY_RANGE*2)));
	// 	if(interval.chr == contig_chr && abs(interval.loc - contig_loc) < MAX_ASSEMBLY_RANGE){
	// 		anchor_intervals++;
	// 	}
	// 	else{
	// 		non_anchor_intervals++;
	// 	}
	// }
	// cerr << "Anchor Valid Intervals: " << anchor_intervals << endl;
	// cerr << "Non-Anchor Valid Intervals: " << non_anchor_intervals << endl;
}
// END STATS

	int interval_count = 0;
	int one_bp_id = 0;
	int pair_id = 0;
	int num_contigs = all_contigs.size();
	int rc = 0;
	int use_rc = 1;
	int contig_loc, contig_len, loc, id;
	long long index = 0;
	bool found = false;
	char contig_chr;
	string contig_seq;
	contig cur_contig;
	mapping_ext cur_interval;
	mapping_ext source = {-1, 0, 0, -1, 0, 0, 0};
	mapping_ext sink = {-1, 0, 0, -1, 0, 0, 0};

	vector<sortable_mapping>* sorted_intervals = new vector<sortable_mapping>[25];	// Sorted list for traversal 
	unordered_map<int, vector<int>> interval_pair_ids;								// Mapping from interval -> interval_pairs
	vector<mapping_ext> all_intervals;												// Bulk data of unused intervals
	vector<interval_pair> interval_pairs; 											// Pairs of intervals

	FILE *fo_vcf = fopen(out_vcf.c_str(), "wb");
	FILE *fr_vcf = fopen((out_vcf + ".reads").c_str(), "wb");
	FILE *fo_full = fopen(out_full.c_str(), "wb");

	// Traverse each contig
	for(int j = 0; j < num_contigs; j++){
		
		cur_contig = all_contigs[j];
		contig_seq = cur_contig.data;
		contig_chr = cur_contig.cluster_chr;
		contig_loc = cur_contig.cluster_loc;
		contig_len = contig_seq.length();
		vector<vector<mapping_ext>> intervals(contig_len+1, vector<mapping_ext>());
		unordered_map<int, pair<int,int>> lookup; 
		bool visited[contig_len];
		memset(visited, 0, sizeof(visited));
		id = 1;
		interval_count = 1;

		// Sort intervals by contig location
		for(char chr=0; (id <= REPEAT_LIMIT3) && chr < 25; chr++){  //(interval_count <= REPEAT_LIMIT3) && 
			for(rc=0;  (id <= REPEAT_LIMIT3) && rc <= use_rc; rc++){ //(interval_count <= REPEAT_LIMIT3) &&
				if(!contig_mappings[rc][chr][j].empty()){
					for(auto &interval: contig_mappings[rc][chr][j]){          
						cur_interval = copy_interval(chr, rc, j, interval);
						loc = (cur_interval.con_loc+k-cur_interval.len);
						if(loc < 0){
							print_interval("ERROR Invalid Interval:", cur_interval);
							continue;
						}
						// Sort ties by length
						if(!intervals[loc].empty()){
							found = false;
							for(int i = 0; i < intervals[loc].size(); i++){
								if(cur_interval.len > intervals[loc][i].len){
									intervals[loc].insert(intervals[loc].begin()+i,cur_interval);
									found = true;
									break;
								}
							}
							if(!found)intervals[loc].push_back(cur_interval);
						}
						else{
							intervals[loc].push_back(cur_interval);
						}
						id++;
						if(id > REPEAT_LIMIT3)break;
if(j == CON_NUM_DEBUG)print_interval("Included", cur_interval);
					}
					contig_mappings[rc][chr][j].clear();
				}
			}
		}

		if(id > REPEAT_LIMIT3)continue;

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

		//Build Edges
		//==================================
		for(int i = 0; i < contig_len; i++){

			if(!intervals[i].empty()){

				for(auto &interval1: intervals[i]){
					found = false;

					//Build edge for balanced SVs and deletions
					for(int u = max(0, i+interval1.len-ANCHOR_SIZE); u <= min(contig_len-1, i+interval1.len+ANCHOR_SIZE); u++){ //was k
						if(!intervals[u].empty()){
							for(auto &interval2: intervals[u]){	
								//if(abs(abs(interval1.loc+interval1.len-interval2.loc) - abs(interval1.con_loc+k-(interval2.con_loc+k-interval2.len))) <= uncertainty ||   //TODO Solve this mystery
									//(interval1.con_loc+k-(interval2.con_loc+k-interval2.len)) <= uncertainty){
									//((interval2.con_loc+k-interval2.len)-interval1.con_loc+k) <= uncertainty){
									visited[u] = true;
									found = true;
									contig_graph[interval1.id][interval2.id] = 1;
normal_edges++;
								//}
							}
						}
					}

					//Build edges for insertions
					if(!found){
						for(int u = min((int)contig_seq.length()-ANCHOR_SIZE-1, i+interval1.len+MIN_SV_LEN_DEFAULT); u < contig_seq.length()-ANCHOR_SIZE; u++){
							if(!intervals[u].empty()){
								for(auto &interval2: intervals[u]){	
									if(abs(interval1.loc+interval1.len-interval2.loc) <= uncertainty){
										visited[u] = true;
										found = true;
										contig_graph[interval1.id][interval2.id] = 1;
ins_edges++;
									}
								}
								if(found)break;
							}
						}
						//to sink
						if(!found){
							if(contig_graph[0][interval1.id] > 0){
								contig_graph[0][interval1.id] = 0;   //remove singletons (doesn't happen apparently)
singletons++;
							}														// It happens, but apparently this is not stopping it. Consider all contig coverage.
							else{
								contig_graph[interval1.id][id] = 2;
							}
sink_edges++;
						}
					}

					//from source
					if(!visited[i]){
						contig_graph[0][interval1.id] = 2;
source_edges++;
					}
				}
			}
		}


		// Run Ford Fulkerson for max flow
		//====================================

		vector<vector<int>> paths = fordFulkerson(id+1, contig_graph, 0, id); //TODO!!!!! do not generate paths that do not contain anchors


		// Predict short variants
		//====================================

path_count += paths.size();

		if(!paths.empty()){

			int max_paths = 1;

			// Check each path and classify variant
			for(int p = 0; p < min(max_paths, (int)paths.size()); p++){  //top x paths

				vector<int>& path = paths[p];

				mapping_ext cur_i, i1, i2;
				int chromo_dist, contig_dist;
				int a1 = -1;
				int a2 = -1;
				bool is_anchor;
				bool first = true;
				bool anchor_found = false;

				for(int i = path.size()-1; i >= 0; i--){

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
							if(first && a1 != path.size()-1){
								
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
								first = false;
							}
							//Insertion, Left Side
							else if(a1 == path.size()-1){ 
								i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
								if(i1.con_loc+k-i1.len > ANCHOR_SIZE){
intc2++;	
									i1.id = all_intervals.size();
									all_intervals.push_back(i1);
									sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
									interval_pairs.push_back({false,-1,i1.id});
						 		 	interval_pair_ids[i1.id].push_back(pair_id++);
								}
							}

							//Insertion, Right Side
							if(a1 == 0){
								i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
								if(((int)contig_seq.length() - (i1.con_loc+k)) > ANCHOR_SIZE){
intc3++;
									i1.id = all_intervals.size();
									all_intervals.push_back(i1);
									sorted_intervals[i1.chr].push_back({i1.id,i1.loc});
									interval_pairs.push_back({false,i1.id,-1});
						 		 	interval_pair_ids[i1.id].push_back(pair_id++);
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

					 			//Negative Strand (Rare)
					 			if(i1.rc && i2.rc){ 

					 				chromo_dist = i1.loc-(i2.loc+i2.len);

					 				//DUP
					 				if(i2.loc+i2.len-i1.loc >= max(ANCHOR_SIZE, uncertainty+1) && i2.loc+i2.len-i1.loc < max_length && abs(contig_dist) <= uncertainty){
						 				print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP");
short_dup1++;
						 			}
						 			else{
						 				//INS
										if(abs(chromo_dist) <= uncertainty && contig_dist > uncertainty){
											print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INS");
short_ins1++;
										}//DEL
										else if(chromo_dist > uncertainty && abs(contig_dist) <= uncertainty){
											print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DEL"); //+1 TP +2FP, probably lots of redundancy
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
										else if(abs(chromo_dist - contig_dist) <= uncertainty){
											print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INV");
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
					 				if(i1.loc+i1.len-i2.loc >= max(ANCHOR_SIZE, uncertainty+1) && i1.loc+i1.len-i2.loc < max_length && abs(contig_dist) <= uncertainty){
						 				print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP");
short_dup2++;
						 			}
						 			else{
						 				//INS
										if(abs(chromo_dist) <= uncertainty && contig_dist > uncertainty){
											print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INS");
short_ins2++;
										}//DEL
										else if(chromo_dist > uncertainty && abs(contig_dist) <= uncertainty){
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
										}//INV
										else if(abs(chromo_dist - contig_dist) <= uncertainty && min(abs(chromo_dist), abs(contig_dist)) >= MIN_SV_LEN_DEFAULT){
										//	print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INV"); //All FP!!!!
//short_inv1++;
										}//??? Maybe SNPs/errors near breakpoint
										else if(chromo_dist > uncertainty && contig_dist > uncertainty){
											if(chromo_dist - abs(contig_dist) >= MIN_SV_LEN_DEFAULT){
												print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DEL");
mystery1++;
											}
										}
										else if(i1.loc+i1.len-i2.loc >= max(MIN_SV_LEN_DEFAULT, uncertainty+1) && (i1.con_loc+k-i1.len) == 0 && (i2.con_loc+k) == contig_len){
											//print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP"); 
short_dup3++;
										}//SNP
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
								//		print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INV");
short_inv2++;
									}
									else if(i1.chr != i2.chr){
										print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "TRANS");
short_trans++;
									}
								}
								else{
									if(abs(contig_dist) > uncertainty && i1.rc == i2.rc && i1.chr == i2.chr){
ugly_fails++;
										print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DEL");
									}
								}
							}
							a1 = a2;
					 		a2 = -1;
						}
						// Interval is anchor with different orientation from the other anchor (one BP INV)
						else if(is_anchor && cur_i.rc != intervals[lookup[path[a1]].first][lookup[path[a1]].second].rc){ 

							for(int x = a1-1; x >= 0; x--){

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
						else if(i == 0){

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

						// Insertion Right Side
						if(a1 == 0){

							i1 = intervals[lookup[path[a1]].first][lookup[path[a1]].second];
							if(((int)contig_seq.length() - (i1.con_loc+k)) >= ANCHOR_SIZE){
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

	for(auto &path: paths){
		for(int i = path.size()-1; i >= 0; i--){
			cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].loc << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].loc + intervals[lookup[path[i]].first][lookup[path[i]].second].len << "\t";
		}
		cerr << endl;
		for(int i = path.size()-1; i >= 0; i--){
			cerr << intervals[lookup[path[i]].first][lookup[path[i]].second].loc << "-" << intervals[lookup[path[i]].first][lookup[path[i]].second].loc + intervals[lookup[path[i]].first][lookup[path[i]].second].len << "\t";
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

	cerr << "Long structural variants..." << endl;

	int w, w_id, z, z_id;
	int chromo_dist, contig_dist, contig_dist1, contig_dist2;

	for(char chr=0; chr < 25; chr++){
		if(sorted_intervals[chr].empty())continue;
		sort(sorted_intervals[chr].begin(), sorted_intervals[chr].end());

		for(w=0; w < sorted_intervals[chr].size()-1; w++){  

			w_id = sorted_intervals[chr][w].id;
 
			if(!interval_pair_ids[w_id].empty()){

				z = w+1;
				mapping_ext iw = all_intervals[w_id];
				mapping_ext ix, iy, iz;

				while(z < sorted_intervals[chr].size() && (all_intervals[sorted_intervals[chr][z].id].loc - iw.loc) < max_length){

					z_id = sorted_intervals[chr][z].id;

					if(!interval_pair_ids[z_id].empty()){

						iz = all_intervals[z_id];
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
											print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, iw, iz, iz, "INS");
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
	//cerr << "Got here " << ix.loc << " " << iy.loc << " " << ((ix.loc+ix.len)-iy.loc) << " " << abs(((ix.loc+ix.len)-iy.loc) - chromo_dist) << endl;
													if(ix.loc > iy.loc && (ix.loc+ix.len)-iy.loc < max_length  && abs(((ix.loc+ix.len)-iy.loc) - chromo_dist) <= uncertainty){
														ip1.visited = true;
														ip2.visited = true;
														print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "INV");
			long_inv1++;
													}
												}
												else if(iw.rc && iz.rc){
													if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < max_length  && abs(((iy.loc+iy.len)-ix.loc) - chromo_dist) <= uncertainty){
														ip1.visited = true;
														ip2.visited = true;
														print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "INV");
			long_inv2++;
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
														print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "TRANS");
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
							}//ip2
						}//ip1
					}
					z++;
				} //z
			}
		}//w
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

					if(!ip1.visited && ip1.id1 != -1 && ip1.id2 != -1){

						iw = all_intervals[ip1.id1];
						ix = all_intervals[ip1.id2];
						contig_dist = abs(iw.con_loc - (ix.con_loc-ix.len));
						chromo_dist = (ix.rc && iw.rc) ? iw.loc-(ix.loc+ix.len) : ix.loc-(iw.loc+iw.len);
						contig_len = all_contigs[ix.con_id].data.length();

						if(contig_dist <= uncertainty){
							if((iw.chr != ix.chr || abs(ix.loc - iw.loc) > max_length) && contig_len-(iw.len + ix.len) < contig_len/2){	
								print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "TRANS");
								ip1.visited = true;
		one_bp_trans++;
							}
							else if(abs(ix.loc - iw.loc) <= max_length ){
								if(iw.rc != ix.rc){
									print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "INV");
		one_bp_inv++;
								}
								else if(chromo_dist >= min_length){
									print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "DEL");
		one_bp_del++;
								}
								else if(iw.loc+iw.len-ix.loc >= max(MIN_SV_LEN_DEFAULT, uncertainty+1)){
									//print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "DUP");
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
						contig_len = all_contigs[ix.con_id].data.length();

						if(ix.len >= contig_len/2 && (ix.con_loc+k-ix.len) > contig_len/3){
							print_variant(fo_vcf, fr_vcf, fo_full, ix.con_id, ix, ix, "INSL");
		one_bp_ins++;
						}
						ip1.visited = true;
					}
					else if(!ip1.visited && ip1.id2 == -1){
						iw = all_intervals[ip1.id1];
						contig_len = all_contigs[iw.con_id].data.length();

						if(iw.len >= contig_len/2 && (contig_len-iw.con_loc+k) > contig_len/3){
							print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, iw, "INSR");
		one_bp_ins++;
						}
						ip1.visited = true;
					}
				}//ip1
			}
		}//i
	}



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

	cerr << "Done." << endl;

	fclose(fo_vcf);
	fclose(fr_vcf);
	fclose(fo_full);

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
 
