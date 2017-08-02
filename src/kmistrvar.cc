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

kmistrvar::kmistrvar(int kmer_len, const string &partition_file, const string &reference, const string &gtf) : 
	variant_caller(partition_file, reference){

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

	repeat = new vector<bool>*[2];
	for(int rc=0; rc < 2; rc++){
		repeat[rc] = new vector<bool>[25];
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
		contig_kmers_index[i] = -1;
	}

	all_intervals.push_back((mapping){"","", 0, 0, -1, 0, 0, 0});
}

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

void kmistrvar::print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id, mapping m1, mapping m2, string type){

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
	string contig_chr = con.cluster_chr;
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
				else if(m1.dup || m2.dup){
					info = "Interspersed";
					if(m1.dup){
						loc = m2.loc;
						if(m1.con_loc+k-1-m1.len >= 0){
							ref = contig_seq.substr(m1.con_loc+k-1-m1.len, 1);
							alt = contig_seq.substr(m1.con_loc+k-1-m1.len, m1.len+1);
						}
						else{
							ref = "N";
							alt = contig_seq.substr(m1.con_loc+k-m1.len, m1.len);
						}
						m2 = m1;
					}
					else if(m2.dup){
						if(m2.con_loc+k-1-m2.len >= 0){
							ref = contig_seq.substr(m2.con_loc+k-1-m2.len, 1);
							alt = contig_seq.substr(m2.con_loc+k-1-m2.len, m2.len+1);
						}
						else{
							ref = "N";
							alt = contig_seq.substr(m2.con_loc+k-m2.len, m2.len);
						}
						
					}
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
		
		//cerr << type << " What?? " << info << " " << m1.rc << " " << m2.rc << " " << (end-loc) << endl;
	}

	score = -10*log(score);
	//pair<int,pair<int,int>> support = compute_support(con, loc, loc+1); //TODO
	//contig_sup = support.first;
	//if(contig_sup > 50)return;

	if(type == "TRANS"){

		if(USE_ANNO){
			locate_interval(m1.chr, loc, loc, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(m2.chr, m2.loc, (m2.loc + m2.len), gene_sorted_map[m2.chr], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];

			if(vec_best1[0] > 1 && vec_best2[0] >1 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;Source=%s,%d,%d,+\n", 
			m1.chr.c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
			info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), m2.chr.c_str(), m2.loc, (m2.loc + m2.len));
	}
	else{

		if(USE_ANNO){
			locate_interval(m1.chr, loc, loc, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(m2.chr, end, end, gene_sorted_map[m2.chr], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];
			locate_interval(m1.chr, loc, end, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene3, best_trans3, vec_best3);
			context3 = contexts[vec_best3[0]];

			if(vec_best1[0] > 1 && vec_best2[0] >1 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		if(info == "Interspersed"){
			fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s;Source=%s,%d,%d,+\n", 
				m1.chr.c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
				info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str(), m2.chr.c_str(), m2.loc, (m2.loc + m2.len));
		}
		else{
			fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;Contig=%d;Support=%d;Info=%s;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s\n", 
				m1.chr.c_str(), loc, ref.c_str(), alt.c_str(), score, type.c_str(), end, contig_loc, id, contig_sup, 
				info.c_str(), context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str());
		}
	}
	 
	fprintf(fo_full, "%s\t%s\t%d\t%d\t%s\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\n", 
		type.c_str(), m1.chr.c_str(), m1.loc,  (m1.loc + m1.len), m2.chr.c_str(), m2.loc, (m2.loc + m2.len), id, 
		m1.chr.c_str(), (m1.con_loc-(m1.len-k)), (m1.con_loc+k), m2.chr.c_str(), (m2.con_loc-(m2.len-k)), (m2.con_loc+k), contig_sup, info.c_str());

	if(PRINT_READS){
		fprintf(fr_vcf, ">Cluster: %d Contig: %d MaxSupport: %d Reads: \n", contig_loc, id, contig_sup);
		fprintf(fr_vcf, "ContigSeq: %s\n", contig_seq.c_str());
		for (auto &read: con.read_information)
			fprintf(fr_vcf, "+ %d %d %s %s\n", read.location_in_contig, read.seq.size(), read.name.c_str(), read.seq.c_str());
	}
} 

void kmistrvar::print_variant(FILE* fo_vcf, FILE* fr_vcf, FILE* fo_full, int id1, int id2, mapping m1, mapping m2, mapping m3, mapping m4, string type){

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
		//cerr << "2contig " << type << " " << (end2-loc1) << endl;
	}

	score = -10*log(score);
	// pair<int,pair<int,int>> support = compute_support(con1, loc1, loc1+1);
	// contig_sup1 = support.first;
	// support = compute_support(con2, loc2, loc2+1);
	// contig_sup2 = support.first;
	// if(contig_sup1 > 50 || contig_sup2 > 50)return;

	if(type == "TRANS"){

		if(USE_ANNO){
			locate_interval(m1.chr, loc1, end1, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(m4.chr, loc2, end2, gene_sorted_map[m2.chr], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];

			if(vec_best1[0] > 1 && vec_best2[0] >1 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t%s[%s:%d[\t%f\tPASS\tSVTYPE=BND;MATEID=bnd_%d;Cluster=%d;Contig=%d;Support=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s\n", 
			m1.chr.c_str(), loc1, m1.id, ref.c_str(), ref.c_str(), m2.chr.c_str(), loc2, score, m4.id, contig_loc1, id1, contig_sup1,
			context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str()); 
		fprintf(fo_vcf, "%s\t%d\tbnd_%d\t%s\t]%s:%d]%s\t%f\tPASS\tSVTYPE=BND;MATEID=bnd_%d;Cluster=%d;Contig=%d;Support=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s\n", 
			m4.chr.c_str(), end1, m4.id, ref2.c_str(), m2.chr.c_str(), end2, ref2.c_str(), score, m1.id, contig_loc2, id2, contig_sup2,
			context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str()); 
	}
	else{

		if(USE_ANNO){
			locate_interval(m1.chr, loc1, loc1, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene1, best_trans1, vec_best1);
			context1 = contexts[vec_best1[0]];
			locate_interval(m4.chr, end2, end2, gene_sorted_map[m2.chr], 0, iso_gene_map, best_gene2, best_trans2, vec_best2);
			context2 = contexts[vec_best2[0]];
			locate_interval(m1.chr, loc1, end2, gene_sorted_map[m1.chr], 0, iso_gene_map, best_gene3, best_trans3, vec_best3);
			context3 = contexts[vec_best3[0]];

			if(vec_best1[0] > 1 && vec_best2[0] >1 && best_gene1 != best_gene2)type = type + "-FUSION";
		}

		fprintf(fo_vcf, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster1=%d;Contig1=%d;Support1=%d;Cluster2=%d;Contig2=%d;Support2=%d;ContextL=%s;GeneL=%s;TransL=%s;ContextR=%s;GeneR=%s;TransR=%s;ContextC=%s;GeneC=%s;TransC=%s\n", 
			m1.chr.c_str(), loc1, ref.c_str(), alt.c_str(), score, type.c_str(), end2, contig_loc1, id1, contig_sup1, contig_loc2, id2, contig_sup2,
			 context1.c_str(), best_gene1.c_str(), best_trans1.c_str(), context2.c_str(), best_gene2.c_str(), best_trans2.c_str(), context3.c_str(), best_gene3.c_str(), best_trans3.c_str());
	}

	fprintf(fo_full, "%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\tCONTIG:\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n", 
		type.c_str(), m1.chr.c_str(), m1.loc,  (m1.loc + m1.len), m4.chr.c_str(), m4.loc, (m4.loc + m4.len),
		type.c_str(), m2.chr.c_str(), m2.loc,  (m2.loc + m2.len), m3.chr.c_str(), m3.loc, (m3.loc + m3.len), 
		id1, m1.chr.c_str(), (m1.con_loc-(m1.len-k)), (m1.con_loc+k), m4.chr.c_str(), (m4.con_loc-(m4.len-k)), (m4.con_loc+k), contig_sup1,
		id2, m2.chr.c_str(), (m2.con_loc-(m2.len-k)), (m2.con_loc+k), m3.chr.c_str(), (m3.con_loc-(m3.len-k)), (m3.con_loc+k), contig_sup2);

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

void kmistrvar::print_interval(string label, mapping& interval){
	cerr << label << "\tRef_Loc: " << interval.loc << "-" << (interval.loc+interval.len) << "\tCon_Loc: " << ((interval.con_loc+k)-interval.len) << "-" << (interval.con_loc+k) 
		 << "\tLen:" << interval.len << "\tChr: " << interval.chr << "\tRC: " << interval.rc << "\tDup: " << interval.dup << "\tCon_ID: " << interval.con_id << endl; 
}									

void kmistrvar::run_kmistrvar(const string &range, const string &out_vcf, const string &out_full, int min_support, int uncertainty, const bool LEGACY_ASSEMBLER, const bool LOCAL_MODE, int ref_flank)
{
	assemble(range, min_support, LEGACY_ASSEMBLER, LOCAL_MODE, ref_flank);
	generate_intervals(LOCAL_MODE);
	predict_variants(out_vcf, out_full, uncertainty);
}

void kmistrvar::assemble(const string &range, int min_support, const bool LEGACY_ASSEMBLER, const bool LOCAL_MODE, int ref_flank)
{
	
	int contig_id = 0;
	int cur_test = 0;
	int loc;
	long long index = 0;
	string con;
	string cluster_chr;

//STATS
double part_count = 0;
double contig_count1 = 0;
double contig_count2 = 0;
double support_count = 0;


	if(!LOCAL_MODE){  //WARNING: This will not work with all reference files. Same with old version. 
		if(TESTING){
			regions.push_back({"1",{1, 310000}});
			regions.push_back({"2",{1, 310000}});
		}
		else{
			regions.push_back({"1",{1, 300000000}});
			regions.push_back({"10",{1, 300000000}});
			regions.push_back({"11",{1, 300000000}});
			regions.push_back({"12",{1, 300000000}});
			regions.push_back({"13",{1, 300000000}});
			regions.push_back({"14",{1, 300000000}});
			regions.push_back({"15",{1, 300000000}});
			regions.push_back({"16",{1, 300000000}});
			regions.push_back({"17",{1, 300000000}});
			regions.push_back({"18",{1, 300000000}});
			regions.push_back({"19",{1, 300000000}});
			regions.push_back({"2",{1, 300000000}});
			regions.push_back({"20",{1, 300000000}});
			regions.push_back({"21",{1, 300000000}});
			regions.push_back({"22",{1, 300000000}});
			regions.push_back({"3",{1, 300000000}});
			regions.push_back({"4",{1, 300000000}});
			regions.push_back({"5",{1, 300000000}});
			regions.push_back({"6",{1, 300000000}});
			regions.push_back({"7",{1, 300000000}});
			regions.push_back({"8",{1, 300000000}});
			regions.push_back({"9",{1, 300000000}});
			regions.push_back({"MT",{1, 300000000}});
			regions.push_back({"X",{1, 300000000}});
			regions.push_back({"Y",{1, 300000000}});
		}
	}

	//Build contigs and check for kmers
	//==================================
	cerr << "Assembling contigs..." << endl;
	while (1) {
		
		vector<string> reads;
		vector<pair<pair<string, string>, int>> p;
		vector<contig> contigs;

		if(TESTING){

			if(cur_test >= TEST_PARTS)break;

			string line;

			bool start = false;
			vector<string> tokens;
			
			ifstream myfile(part_file);
			if(myfile.is_open()){
					while(getline(myfile,line)){
						tokens = split(line);
						if(tokens[0] == to_string(cur_test)){
							cluster_chr = tokens[4];
							loc = atoi(tokens[2].c_str());
							start = true;
							continue;
						}
						if(start && tokens[0] == to_string(cur_test+1)){
								break;
						}
						if(start){
							reads.push_back(tokens[1]);
						}
					}
				myfile.close();
			}
			cur_test++;
		}
		else{

			p = pt.read_partition(part_file, range);

			if (!p.size()) 
				break;
			if(p.size() < min_support)
				continue;
			if (p.size() > 10000) {
				fprintf(stderr, "Skipping a gigantic cluster %d\n", p.size());
				continue;
			}

			for (int i = 0; i < p.size(); i++) {
				reads.push_back(p[i].first.second);
			}
		}

		if(LEGACY_ASSEMBLER){
			contigs = as_old.assemble(p); 
		}
		else{
			contigs = as.assemble(reads); 
		}

//STATS
if(PRINT_STATS){
	part_count++;
	contig_count1 += contigs.size(); 	
}

		reads.clear();
		
		//TODO: Contig merging
		for (auto &contig: contigs){ 
			if (contig.support() >= min_support) {

//STATS
if(PRINT_STATS){
	contig_count2++; 
	support_count += contig.support();	

	if(pt.get_start() == 67685799){
		cout << "support: " << contig.support() << " " << contig.data.length() << " " << pt.get_start() << " " << contig_id << endl;
		cout << contig.data << endl;
	}
}
			
				if(TESTING){
					contig.cluster_chr = cluster_chr;
					contig.cluster_loc = loc;
					if(LOCAL_MODE)regions.push_back({cluster_chr, {(loc-ref_flank), (loc+ref_flank+300)}});
				}
				else{
					contig.cluster_chr = pt.get_reference();
					contig.cluster_loc = pt.get_start();
					if(LOCAL_MODE)regions.push_back({pt.get_reference(), {(pt.get_start()-ref_flank), (pt.get_end()+ref_flank)}});
				}
				
				con = contig.data;
				unordered_map<long long, vector<int>> kmer_location;
				kmer_locations.push_back(kmer_location);
				all_contigs.push_back(contig);
				int i = 0, x = 0, j = 0;

				for(int rc=0; rc < 2; rc++){
					for(int chr=0; chr < 25; chr++){
						repeat[rc][chr].push_back(false);
					}
				}

				for(j = 0; j < k-1; j++){    // k-1 kmer for first round
					if(j > 0 && con[i+j] == 'N'){
						i += (j+1);
						j = -1;
						index = 0;
						continue;
					}
					index <<= 2;
					index |= ((con[i+j] & MASK) >> 1);
				}

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

				
					index <<= 2;
					index |= ((con[i+k-1] & MASK) >> 1); 
					index &= kmer_mask;

//cerr << index << " " << con.substr(i, k) << " 2" << endl;

					if(contig_kmers_index[index] == -1){
						vector<int> contig_list;
						contig_kmers_index[index] = contig_kmers.size();
						contig_list.push_back(contig_id);
						contig_kmers.push_back(contig_list);
					}
					else{
						if(contig_kmers[contig_kmers_index[index]].back() != contig_id)contig_kmers[contig_kmers_index[index]].push_back(contig_id);
					}
					kmer_locations[contig_id][index].push_back(i);
				}
				contig_id++;
			}
		}
	}

	if(all_contigs.empty()){
		cerr << "No contigs could be assembled. Exiting..." << endl;
		exit(1);
	}

	for(int rc=0; rc <= 1; rc++){
		for(int chr=0; chr < 25; chr++){
			contig_mappings[rc][chr] = vector<vector<mapping>>(contig_id);
		}
	}

//STATS
if(PRINT_STATS){
	cerr << "Partition Count: " << part_count << endl;
	cerr << "Average Num Contigs Pre-Filter: " << (contig_count1/part_count) << endl;
	cerr << "Average Num Contigs Post-Filter: " << (contig_count2/part_count) << endl;
	cerr << "Average Contig Support: " << (support_count/contig_count2) << endl;
}

	
}

void kmistrvar::generate_intervals(const bool LOCAL_MODE)
{

	//Read reference and map contig kmers
	//=====================================
	cerr << "Generating intervals..." << endl;

	int use_rc = 2; //1 no, 2 yes
	int chr_num = 0;
	int jj = 0;
	int id = 1;  //zero reserved for sink
	long long index = 0;
	long long rc_index = 0;
	long long cur_index = 0;
	int shift = 2*(k-1)-1;
	int i, ii, j, js, je, x, rc;
	int MAX_INTERVAL_LEN = 10000;
	int MAX_CONTIG_LEN = 0; 
	int num_contigs = all_contigs.size();
	int ref_counts[MAX_INTERVAL_LEN];
	int gen_start, gen_end;

//STATS
int anchor_intervals = 0;
int non_anchor_intervals = 0;
int pop_back = 0;
double interval_count = 0;

	vector<mapping>* intervals;
	vector<pair<int,int>> starts;
	mapping* last_interval; 
	string gen;
	string chromo = "first";
	string last_chromo = chromo;

	for(auto& contig : all_contigs){
		if(contig.data.length() > MAX_CONTIG_LEN)MAX_CONTIG_LEN = contig.data.length();
	}
	MAX_CONTIG_LEN += (k+1); //Rare case of SNP near the end of the longest contig.  

	bool** valid_mappings = new bool*[MAX_INTERVAL_LEN]; //think of a better solution

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		valid_mappings[i] = new bool[MAX_CONTIG_LEN];
	}

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		for(int j=0; j < MAX_CONTIG_LEN; j++){
			valid_mappings[i][j] = false;
		}
	}

	for (auto &region: regions){

		chr_num = find(chromos.begin(), chromos.end(), region.first) - chromos.begin();
		chromo = region.first;
		gen = ref.extract(chromo, region.second.first, region.second.second); 
		gen_start = max(0,region.second.first-1);
		gen_end = gen.length()-k+1;

		js = LOCAL_MODE ? jj : 0;
		je = LOCAL_MODE ? (jj+1) : num_contigs;

		i = 0;
		while(gen[i] == 'N'){
			i++;
		}

		index = 0;
		rc_index = 0;

		for(x = 0; x < k-1; x++){    // k-1 kmer for first round
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

		for (i; i < gen_end; i++){
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

			ii = i + gen_start;// + 1;

			for(rc=0; rc < use_rc; rc++){
				if(rc){
					rc_index >>= 2;
					rc_index |= (((gen[i+k-1] & MASK) ^ MASK_RC) << shift);
					cur_index = rc_index;
				}
				else{
					index <<= 2;
					index |= ((gen[i+k-1] & MASK) >> 1); 
					index &= kmer_mask;
					cur_index = index;
				}

				if(contig_kmers_index[cur_index] == -1)continue;

				for(const int& j : contig_kmers[contig_kmers_index[cur_index]]){

					if( (!LOCAL_MODE || j == jj)){ 

						intervals = &contig_mappings[rc][chr_num][j];
						last_interval = &intervals->back();

						if(intervals->empty()){
							intervals->push_back((mapping){"",chromo, rc, false, ii, k, -1, j, -1});
						}
						else if(last_interval->chr == chromo && last_interval->loc == ii-(last_interval->len-k)-1){
							last_interval->len++;
						}
						else if(last_interval->chr == chromo && (ii-(last_interval->loc + last_interval->len) == 1)){ //Allows for SNPs
							last_interval->len += (k+1);
						}
						else if(last_interval->chr != chromo || (ii-(last_interval->loc + last_interval->len) > 1)){

							if(last_interval->len < ANCHOR_SIZE){
								intervals->pop_back();
							}
							else if(repeat[rc][chr_num][j] && (chromo != all_contigs[j].cluster_chr || (abs(last_interval->loc+(last_interval->len/2) - all_contigs[j].cluster_loc) > (MAX_ASSEMBLY_RANGE*2))) && last_interval->len < intervals->at(intervals->size()-2).len){
								intervals->pop_back();
pop_back++;
							}
							else{
								for(int x = intervals->size()-1; x >= 0; x--){
									if(last_interval->len > intervals->at(x).len){
										intervals->insert(intervals->begin()+x,*last_interval);
										intervals->pop_back();
										if(repeat[rc][chr_num][j] && (chromo != all_contigs[j].cluster_chr || (abs(last_interval->loc+(last_interval->len/2) - all_contigs[j].cluster_loc) > (MAX_ASSEMBLY_RANGE*2))) || intervals->size() > REPEAT_LIMIT2){
pop_back++;
											intervals->pop_back();
										}
										break;
									}
								}
							}

							if(intervals->size() > REPEAT_LIMIT1){
								repeat[rc][chr_num][j] = true;
							}
							intervals->push_back((mapping){"",chromo, rc, false, ii, k, -1, j, -1});
						}

						if(LOCAL_MODE)break;
					}
				}
			}
		}// gen, slowest loop


		//Remove intervals without consectutive kmers in contigs
		//======================================================
		for(rc=0; rc < use_rc; rc++){
			for(j = js; j < je; j++){

				if(!contig_mappings[rc][chr_num][j].empty()){
					if(contig_mappings[rc][chr_num][j].back().len < ANCHOR_SIZE){
pop_back++;
						contig_mappings[rc][chr_num][j].pop_back();
						if(contig_mappings[rc][chr_num][j].empty())continue;
					}
				 
					vector<mapping> valid_intervals;
					string contig_seq = all_contigs[j].data;

					for(auto &interval: contig_mappings[rc][chr_num][j]){

//STATS
if(PRINT_STATS){
	int contig_loc = all_contigs[interval.con_id].cluster_loc;
	string contig_chr = all_contigs[interval.con_id].cluster_chr;
	if(interval.chr == contig_chr && abs(interval.loc - contig_loc) < MAX_ASSEMBLY_RANGE){
		anchor_intervals++;
	}
	else{
		non_anchor_intervals++;
	}
}
	
						string match = ref.extract(chromo, interval.loc+1, (interval.loc + interval.len));
						if(rc)match = reverse_complement(match);

						starts.clear();

						int ilen;
						int len = match.length()-k+1;
						int x, y;

						for(x=0; x < (match.length()+k+1); x++){
							ref_counts[x] = 0;
							for(y=0; y < (contig_seq.length()+k+1); y++){
								valid_mappings[x][y] = false;
							}
						}

//cerr << "Contig " << j << " of " << je << endl;
//cerr << all_contigs[j].data << endl;
//cerr << "Length: " << match.length() << " " << interval.loc <<	endl;
//cerr << " =================================== " << endl;
						for (i = 0; i < len; i++){

							cur_index = 0;

							for(x = 0; x < k; x++){
								cur_index <<= 2;
								cur_index |= ((match[i+x] & MASK) >> 1);
							}

							for(auto &loc: kmer_locations[j][cur_index]){

								valid_mappings[i][loc] = true;
								ref_counts[i]++;

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

						for(auto& start : starts){

							x = start.first;
							y = start.second;
							bool con_repeat = false;

							while(valid_mappings[x++][y++]){
								if(!valid_mappings[x][y] && valid_mappings[x+k][y+k]){
									x +=k;
									y +=k;
								}
								else if(!valid_mappings[x][y] && y+k < contig_seq.length() && contig_seq[y+k] == 'N'){
									valid_mappings[x][y] = true;
								}
							}
							x-=2;
							y-=2;

							ilen = (y - start.second)+k;

							//TODO: right side bias, incomplete intervals? Repeats?
							// for(int z=start.first; z <= x; z++){
							// 	if(ref_counts[z] > 1){
							// //		dup = true;
							// 	}
							// 	else{
							// 		dup = false;
							// 		break;
							// 	}
							// }


							if(ilen > ANCHOR_SIZE){
								for(int z=start.first; z <= x; z++){
	//cerr << ref_counts[z] << ",";
									if(ref_counts[z] > 2){
										con_repeat = true;
									}
									else{
										con_repeat = false;
										break;
									}
								}
//cerr << endl;	
								if(!con_repeat){
									mapping sub_interval = interval;
									sub_interval.seq = contig_seq.substr(start.second, ilen);
									sub_interval.loc += start.first;
									sub_interval.len = ilen;
									sub_interval.con_loc = y;
									//if(ref_counts[start.first] > 1 && ref_counts[x] > 1)sub_interval.dup = true; //&& ref_counts[(start.first + (x-start.first)/2)] > 1
									sub_interval.dup = false;//dup;
									sub_interval.id = id++;
									valid_intervals.push_back(sub_interval); 
								}
							}
						}
					}//intervals

//STATS
interval_count += valid_intervals.size(); 

					num_intervals += valid_intervals.size();
					contig_mappings[rc][chr_num][j] = valid_intervals;

					if(!valid_intervals.empty()){
						for(auto &interval: valid_intervals){
							all_intervals.push_back(interval);
						}
					}
				}
			}
		}
		jj++;
	}

//STATS
if(PRINT_STATS){
	cerr << "Anchor Intervals: " << anchor_intervals << endl;
	cerr << "Non-Anchor Intervals: " << non_anchor_intervals << endl;
	cerr << "Average Intervals per Contig: " << (interval_count/(double)num_contigs) << endl;
	cerr << "Pop backs: " << pop_back << endl;
}

	for(int i=0; i < MAX_INTERVAL_LEN; i++){
		delete[] valid_mappings[i];
	}
	delete[] valid_mappings;
}

void kmistrvar::predict_variants(const string &out_vcf, const string &out_full, int uncertainty)
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
	cerr << "Intervals: " << all_intervals.size() << endl;
	cerr << "Contigs: " << all_contigs.size() << endl;

	for(auto & interval : all_intervals){
		int contig_loc = all_contigs[interval.con_id].cluster_loc;
		string contig_chr = all_contigs[interval.con_id].cluster_chr;
		bool is_anchor = ((interval.chr == contig_chr) && (abs(interval.loc+(interval.len/2) - contig_loc) < (MAX_ASSEMBLY_RANGE*2)) && !interval.dup);
		if(interval.chr == contig_chr && abs(interval.loc - contig_loc) < MAX_ASSEMBLY_RANGE){
			anchor_intervals++;
		}
		else{
			non_anchor_intervals++;
		}
	}
	cerr << "Anchor Valid Intervals: " << anchor_intervals << endl;
	cerr << "Non-Anchor Valid Intervals: " << non_anchor_intervals << endl;
}
// END STATS

	int interval_count = 0;
	int pair_id = 0;
	int num_contigs = all_contigs.size();
	int rc = 0;
	int use_rc = 1;
	long long index = 0;
	bool found = false;

	vector<int>* interval_pair_ids = new vector<int>[all_intervals.size()];
	vector<pair<bool,pair<int,int>>> interval_pairs; 

	FILE *fo_vcf = fopen(out_vcf.c_str(), "wb");
	FILE *fr_vcf = fopen((out_vcf + ".reads").c_str(), "wb");
	FILE *fo_full = fopen(out_full.c_str(), "wb");

	for(int j = 0; j < num_contigs; j++){
		
		contig cur_contig = all_contigs[j];
		string contig_seq = cur_contig.data;
		string contig_chr = cur_contig.cluster_chr;
		int contig_loc = cur_contig.cluster_loc;
		int contig_len = contig_seq.length();
		vector<vector<mapping>> intervals(contig_len, vector<mapping>());
		unordered_map<int,int> interval_to_id;
		unordered_map<int,int> id_to_interval;
		bool visited[contig_len];
		memset(visited, 0, sizeof(visited));
		int loc;
		int id = 1;
		interval_count = 2;

		//sort by contig location
		for(int chr=0; chr < 25; chr++){
			for(rc=0; rc <= use_rc; rc++){
				if(!contig_mappings[rc][chr][j].empty()){
					for(auto &interval: contig_mappings[rc][chr][j]){
						loc = (interval.con_loc+k-interval.len);
						if(loc < 0){
							print_interval("ERROR Invalid Interval:", interval);
							continue;
						}
						//Sort by length
						if(!intervals[loc].empty()){
							found = false;
							for(int i = 0; i < intervals[loc].size(); i++){
								if(interval.len > intervals[loc][i].len){
									intervals[loc].insert(intervals[loc].begin()+i,interval);
									found = true;
									break;
								}
							}
							if(!found)intervals[loc].push_back(interval);
						}
						else{
							intervals[loc].push_back(interval);
						}
if(j == CON_NUM_DEBUG)print_interval("Included", interval);
					}
				}
			}
		}

		for(int i = 0; i < contig_seq.length(); i++){
			if(!intervals[i].empty()){
				for(auto &interval: intervals[i]){
					interval_to_id[interval.id] = id;
					id_to_interval[id++] = interval.id;
					interval_count++;
				}
			}	
		}

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
		for(int i = 0; i < contig_seq.length(); i++){
			if(!intervals[i].empty()){
				for(auto &interval1: intervals[i]){
					found = false;

					//Build edge for balanced SVs and deletions
					for(int u = max(0, i+interval1.len-ANCHOR_SIZE); u <= min((int)contig_seq.length()-1, i+interval1.len+ANCHOR_SIZE); u++){ //was k
						if(!intervals[u].empty()){
							for(auto &interval2: intervals[u]){	
								//if(abs(abs(interval1.loc+interval1.len-interval2.loc) - abs(interval1.con_loc+k-(interval2.con_loc+k-interval2.len))) <= uncertainty ||   //TODO Solve this mystery
									//(interval1.con_loc+k-(interval2.con_loc+k-interval2.len)) <= uncertainty){
									//((interval2.con_loc+k-interval2.len)-interval1.con_loc+k) <= uncertainty){
									visited[u] = true;
									found = true;
									contig_graph[interval_to_id[interval1.id]][interval_to_id[interval2.id]] = 1;
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
										contig_graph[interval_to_id[interval1.id]][interval_to_id[interval2.id]] = 1;
ins_edges++;
									}
								}
								if(found)break;
							}
						}
						//to sink
						if(!found){
							if(contig_graph[0][interval_to_id[interval1.id]] > 0){
								contig_graph[0][interval_to_id[interval1.id]] = 0;   //remove singletons (doesn't happen apparently)
singletons++;
							}														// It happens, but apparently this is not stopping it. Consider all contig coverage.
							else{
								contig_graph[interval_to_id[interval1.id]][id] = 2;
							}
sink_edges++;
						}
					}

					//from source
					if(!visited[i]){
						contig_graph[0][interval_to_id[interval1.id]] = 2;
source_edges++;
					}
				}
			}
		}

		vector<vector<int>> paths = fordFulkerson(id+1, contig_graph, 0, id); //TODO!!!!! do not generate paths that do not contain anchors

		//Predict short variants
		//====================================

path_count += paths.size();

		if(!paths.empty()){

			int max_paths = 1;

			for(int p = 0; p < min(max_paths, (int)paths.size()); p++){  //top x paths

				vector<int>& path = paths[p];

				mapping cur_i, i1, i2;
				int chromo_dist, contig_dist;
				int a1 = -1;
				int a2 = -1;
				bool is_anchor;
				bool first = true;
				bool anchor_found = false;

				for(int i = path.size()-1; i >= 0; i--){

					cur_i = all_intervals[id_to_interval[path[i]]];
					is_anchor = ((cur_i.chr == contig_chr) && (abs(cur_i.loc+(cur_i.len/2) - contig_loc) < (MAX_ASSEMBLY_RANGE*2)) && !cur_i.dup);
					found = false;
					if(is_anchor)anchor_found = true;

					if(a1 == -1){
						if(is_anchor){
							a1 = i;
							if(first && a1 != path.size()-1){
								//leading non-anchor region
								i1 = all_intervals[id_to_interval[path[a1+1]]];
				 				i2 = all_intervals[id_to_interval[path[a1]]];
				 				chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
				 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

								if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){ 
intc1++;
					 				interval_pairs.push_back({false,{i1.id,i2.id}}); 
					 				interval_pair_ids[i1.id].push_back(pair_id);
						 			interval_pair_ids[i2.id].push_back(pair_id++);
								}
								first = false;
							}
							else if(a1 == path.size()-1){ //Insertion, Left Side
								i1 = all_intervals[id_to_interval[path[a1]]];
								if(i1.con_loc+k-i1.len > ANCHOR_SIZE){
intc2++;
									interval_pairs.push_back({false,{-1,i1.id}});
									interval_pair_ids[i1.id].push_back(pair_id++);
								}
							}

							//Insertion, Right Side
							if(a1 == 0){
								i1 = all_intervals[id_to_interval[path[a1]]];
								if(((int)contig_seq.length() - (i1.con_loc+k)) > ANCHOR_SIZE){
intc3++;
									interval_pairs.push_back({false,{i1.id, -1}});
									interval_pair_ids[i1.id].push_back(pair_id++);
								}
							}
						}
						else{
anchor_fails++;
						}
					}
					else if(a2 == -1){

						if((is_anchor) && cur_i.rc == all_intervals[id_to_interval[path[a1]]].rc){// && !cur_i.rc){

							a2 = i;

							if(a1 - a2 == 1){

								i1 = all_intervals[id_to_interval[path[a1]]];
					 			i2 = all_intervals[id_to_interval[path[a2]]];

					 			chromo_dist = i2.loc-(i1.loc+i1.len);
					 			contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

					 			if(i1.rc && i2.rc){ //Negative Strand (Rare)

					 				chromo_dist = i1.loc-(i2.loc+i2.len);

					 				if(i2.loc+i2.len-i1.loc >= max(ANCHOR_SIZE, uncertainty+1) && i2.loc+i2.len-i1.loc < MAX_SV_LEN && abs(contig_dist) <= uncertainty){
						 				print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP");
						 				//interval_pairs.push_back({false,{i1.id, i2.id}});
						 				//interval_pair_ids[i1.id].push_back(pair_id);
						 				//interval_pair_ids[i2.id].push_back(pair_id++);
short_dup1++;
						 			}
						 			else{
						 				//INS
										if(abs(chromo_dist) <= uncertainty && contig_dist > uncertainty){
											print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INS");
short_ins1++;
										}//DEL
										else if(chromo_dist > uncertainty && abs(contig_dist) <= uncertainty && !i1.dup && !i2.dup){
intc4++;									print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DEL"); //+1 TP +2FP, probably lots of redundancy
											interval_pairs.push_back({false,{i1.id, i2.id}});
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
										}
										else{
//mystery2++;									//FP
											//SNP
										}
						 			}
					 			}
					 			else{ //Positive Strand
					 				if(i1.loc+i1.len-i2.loc >= max(ANCHOR_SIZE, uncertainty+1) && i1.loc+i1.len-i2.loc < MAX_SV_LEN && abs(contig_dist) <= uncertainty){
						 				print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP");
						 				//interval_pairs.push_back({false,{i1.id, i2.id}});
						 				//interval_pair_ids[i1.id].push_back(pair_id);
						 				//interval_pair_ids[i2.id].push_back(pair_id++);
						 			//	print_interval("interval1", i1);
						 			//	print_interval("interval2", i2);
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
											interval_pairs.push_back({false,{i1.id, i2.id}});
						 					interval_pair_ids[i1.id].push_back(pair_id);
						 					interval_pair_ids[i2.id].push_back(pair_id++);
										}//INV
										else if(abs(chromo_dist - contig_dist) <= uncertainty && min(abs(chromo_dist), abs(contig_dist)) >= MIN_SV_LEN_DEFAULT){
										//	print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INV"); //All FP!!!!
//short_inv1++;
										}//???
										else if(chromo_dist > uncertainty && contig_dist > uncertainty){
											if(chromo_dist - abs(contig_dist) >= MIN_SV_LEN_DEFAULT){
												print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DEL");
mystery1++;
											}
										}
										else if(i1.loc+i1.len-i2.loc >= max(MIN_SV_LEN_DEFAULT, uncertainty+1) && (i1.con_loc+k-i1.len) == 0 && (i2.con_loc+k) == contig_len){
											//print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP"); 
short_dup3++;
										}
										else{
mystery2++;
											//SNP
										}
						 			}
					 			}
							}
							else{

								//Inner non-anchor region
								for(int u = a1; u > a2; u--){

									i1 = all_intervals[id_to_interval[path[u]]];
					 				i2 = all_intervals[id_to_interval[path[u-1]]];
					 				chromo_dist = i2.loc-(i1.loc+i1.len);
					 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

					 				if(i1.rc && i2.rc)chromo_dist = i1.loc-(i2.loc+i2.len);

					 				if(abs(contig_dist) > uncertainty)break;

					 				if((u != a1 && u-1 != a2) && (i1.chr != i2.chr || i1.rc != i2.rc || i1.dup != i2.dup || abs(chromo_dist) > uncertainty))break;

					 				if(u-1 == a2)found = true;
								}

								if(found){
									i1 = all_intervals[id_to_interval[path[a1]]];
					 				i2 = all_intervals[id_to_interval[path[a2]]];

					 				//This is a little shady, only one anchor checked? Maybe FP
									if(all_intervals[id_to_interval[path[a1-1]]].rc != i1.rc  && all_intervals[id_to_interval[path[a2+1]]].rc != i2.rc){
								//		print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "INV");
short_inv2++;
									}
									else if(all_intervals[id_to_interval[path[a1-1]]].dup){
										print_variant(fo_vcf, fr_vcf, fo_full, j, i1, i2, "DUP"); //Why is this not working?

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
						else if(is_anchor && cur_i.rc != all_intervals[id_to_interval[path[a1]]].rc){ 

							for(int x = a1-1; x >= 0; x--){
								if(all_intervals[id_to_interval[path[x]]].rc == cur_i.rc)continue;
							}

							a2 = i;

							if(a1 - a2 == 1){

								i1 = all_intervals[id_to_interval[path[a1]]];
					 			i2 = all_intervals[id_to_interval[path[a2]]];

					 			chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
				 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;

								if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){
	intc5++;
					 				interval_pairs.push_back({false,{i1.id,i2.id}}); 
					 				interval_pair_ids[i1.id].push_back(pair_id);
						 			interval_pair_ids[i2.id].push_back(pair_id++);
								}
					 		}

					 		a1 = a2;
					 		a2 = -1;

						}
						else if(i == 0){

							//trailing non-anchor region
							i1 = all_intervals[id_to_interval[path[a1]]];
			 				i2 = all_intervals[id_to_interval[path[a1-1]]];
			 				chromo_dist = (i1.rc && i2.rc) ? i1.loc-(i2.loc+i2.len) : i2.loc-(i1.loc+i1.len);
			 				contig_dist = (i2.con_loc-i2.len)-i1.con_loc;


							if(abs(chromo_dist) <= uncertainty || abs(contig_dist) <= uncertainty){
intc5++;
				 				interval_pairs.push_back({false,{i1.id,i2.id}}); 
				 				interval_pair_ids[i1.id].push_back(pair_id);
					 			interval_pair_ids[i2.id].push_back(pair_id++);
							}
						}
						else{
if(!is_anchor)anchor_fails++;
						}

						// Insertion Right Side
						if(a1 == 0){
							i1 = all_intervals[id_to_interval[path[a1]]];
							if(((int)contig_seq.length() - (i1.con_loc+k)) >= ANCHOR_SIZE){
intc6++;
								interval_pairs.push_back({false,{i1.id, -1}});
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
			cerr << all_intervals[id_to_interval[path[i]]].loc << "-" << all_intervals[id_to_interval[path[i]]].loc + all_intervals[id_to_interval[path[i]]].len << "\t";
		}
		cerr << endl;
		for(int i = path.size()-1; i >= 0; i--){
			cerr << all_intervals[id_to_interval[path[i]]].con_loc+k - all_intervals[id_to_interval[path[i]]].len << "-" << all_intervals[id_to_interval[path[i]]].con_loc+k  << "\t";
		}
		cerr << endl;
	}
}

		for(int i=0; i < interval_count; i++){
			delete[] contig_graph[i];
		}
		delete[] contig_graph;
	}//j


	cerr << "Long structural variants..." << endl;
	// Long Structural Variants
	//============================================

	for(int w=0; w < num_intervals-1; w++){

		int z = w+1;
		int chromo_dist, contig_dist1, contig_dist2;
		mapping iw = all_intervals[w];
		mapping ix, iy, iz;

		//if(iw.rc)continue;

		while(z < num_intervals && (all_intervals[z].chr == iw.chr) && (all_intervals[z].loc - iw.loc) < MAX_SV_LEN){

			if(!interval_pair_ids[w].empty() && !interval_pair_ids[z].empty()){
				iz = all_intervals[z];
				chromo_dist = (iz.rc && iw.rc) ? iw.loc-(iz.loc+iz.len) : iz.loc-(iw.loc+iw.len);

				// if(iz.rc){
				// 	z++;
				// 	continue;
				// }

				for(int x: interval_pair_ids[w]){
					for(int y: interval_pair_ids[z]){   //Shit, all combinations, TODO: Check if always small. 

						pair<bool,pair<int,int>>& ip1 = interval_pairs[x];
						pair<bool,pair<int,int>>& ip2 = interval_pairs[y];

						if(!ip1.first && !ip2.first && ip1.second.first == w && ip2.second.second == z && ip1.second.second != z && ip2.second.first != w && ip1.second.first != -1 && ip2.second.second != -1){
							if(ip1.second.second == -1 && ip2.second.first == -1){
								if(abs(chromo_dist) <= uncertainty){
									ip1.first = true;
									ip2.first = true;
									print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, iw, iz, iz, "INS");
long_ins++;
								}
							}
							else if(ip1.second.second != -1 && ip2.second.first != -1){

								ix = all_intervals[ip1.second.second];
								iy = all_intervals[ip2.second.first];
								contig_dist1 = abs(iw.con_loc - (ix.con_loc-ix.len));
								contig_dist2 = abs(iy.con_loc - (iz.con_loc-iz.len));

								if(ix.chr == iy.chr && contig_dist1 <= uncertainty && contig_dist2 <= uncertainty){ 
									if(iw.rc != ix.rc && iy.rc != iz.rc){ //if(ix.rc || iy.rc){
										if(ix.rc && iy.rc){
											if(ix.loc > iy.loc && (ix.loc+ix.len)-iy.loc < MAX_SV_LEN  && abs(((ix.loc+ix.len)-iy.loc) - chromo_dist) <= uncertainty){
												ip1.first = true;
												ip2.first = true;
												print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "INV");
long_inv1++;
											}
										}
										else if(iw.rc && iz.rc){
											if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < MAX_SV_LEN  && abs(((iy.loc+iy.len)-ix.loc) - chromo_dist) <= uncertainty){
												ip1.first = true;
												ip2.first = true;
												print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "INV");
long_inv2++;
											}
										}
										
									}
									else{
										if(iw.rc && ix.rc && iy.rc && iz.rc){
//cerr << "this actually happens2 " << ix.loc << " " << iy.loc << " " << iw.loc << " " << iz.loc << endl;
										}
										else{
											if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < MAX_SV_LEN && (iw.chr != ix.chr || (ix.loc > iz.loc+iz.len || iy.loc+iy.len < iw.loc))){ //last check = 1 Trans
												ip1.first = true;
												ip2.first = true;
												print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "TRANS");
long_trans++;
											}
											//else if(ix.loc > iy.loc && (ix.loc+ix.len)-iy.loc < MAX_SV_LEN && abs(iw.loc-iy.loc) <= uncertainty && abs(iz.loc-ix.loc) <= uncertainty){
											else if(ix.loc > iw.loc && iz.loc > iy.loc && abs((iw.loc+iw.len)-ix.loc) < MAX_SV_LEN  && abs((iy.loc+iy.len)-iz.loc) < MAX_SV_LEN && abs(iw.loc-iy.loc) <= uncertainty && abs(iz.loc-ix.loc) <= uncertainty){
												// ip1.first = true;
												// ip2.first = true;
												// print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iz.con_id, iw, ix, iy, iz, "DEL");
long_del++;
											}
											else if(ix.loc < iy.loc && (iy.loc+iy.len)-ix.loc < MAX_SV_LEN && abs(iz.loc-(iw.loc+iw.len)) <= uncertainty){
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
		} //x
	}//i


	//Single interval pair case (WARNING: may cause many false positives)
	for(int w=0; w < num_intervals-1; w++){

		int x = w+1;
		int chromo_dist, contig_dist, contig_len;
		mapping iw, ix;

		if(!interval_pair_ids[w].empty()){

			for(int x: interval_pair_ids[w]){

				pair<bool,pair<int,int>>& ip1 = interval_pairs[x];

				if(!ip1.first && ip1.second.first != -1 && ip1.second.second != -1){

					iw = all_intervals[ip1.second.first];
					ix = all_intervals[ip1.second.second];
					contig_dist = abs(iw.con_loc - (ix.con_loc-ix.len));
					chromo_dist = (ix.rc && iw.rc) ? iw.loc-(ix.loc+ix.len) : ix.loc-(iw.loc+iw.len);
					contig_len = all_contigs[ix.con_id].data.length();

					if(contig_dist <= uncertainty){
						if((iw.chr != ix.chr || abs(ix.loc - iw.loc) > MAX_SV_LEN) && contig_len-(iw.len + ix.len) < contig_len/2){	
							print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "TRANS");
							ip1.first = true;
one_bp_trans++;
						}
						else if(abs(ix.loc - iw.loc) <= MAX_SV_LEN ){
							if(iw.rc != ix.rc){
								print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, ix, "INV");
one_bp_inv++;
							}
							else if(chromo_dist >= ANCHOR_SIZE){
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
							ip1.first = true;
						}
					}
				}
				else if(!ip1.first && ip1.second.first == -1){    //Two more possible with extra FP
					ix = all_intervals[ip1.second.second];
					contig_len = all_contigs[ix.con_id].data.length();

					if(ix.len >= contig_len/2 && (ix.con_loc+k-ix.len) > contig_len/3){
						print_variant(fo_vcf, fr_vcf, fo_full, ix.con_id, ix, ix, "INSL");
one_bp_ins++;
					}
					ip1.first = true;
				}
				else if(!ip1.first && ip1.second.second == -1){
					iw = all_intervals[ip1.second.first];
					contig_len = all_contigs[iw.con_id].data.length();

					if(iw.len >= contig_len/2 && (contig_len-iw.con_loc+k) > contig_len/3){
						print_variant(fo_vcf, fr_vcf, fo_full, iw.con_id, iw, iw, "INSR");
one_bp_ins++;
					}
					ip1.first = true;
				}
			}//ip1
		}
	}//i



if(PRINT_STATS){

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

	delete[] interval_pair_ids;

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
 
