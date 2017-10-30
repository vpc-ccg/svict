#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <string>
#include <map>
#include <zlib.h>
#include "../include/cxxopts.hpp"
#include "partition.h"
#include "common.h"
#include "assembler.h"
#include "assembler_old.h"
#include "assembler_ext.h"
#include "kmistrvar.h"
#include "genome.h"
#include "simulator.h"
//#include "IntervalTree.h"
//#include <DeeZ.h>
#include "common.h"
#include "record.h"
#include "sam_parser.h"
#include "bam_parser.h"
#include "sort.h"

using namespace std;


// 10.46
// S:[step]:END

/********************************************************************/
inline string space (int i) 
{
	return string(i, ' ');
}

/********************************************************************/
inline string itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}

/********************************************************************/
void partify (const string &read_file, const string &mate_file, const string &out, int threshold) 
{
	gzFile gzf = gzopen(mate_file.c_str(), "rb");
	unordered_map<string, string> mymap;

	const int MAXB = 8096;
	char buffer[MAXB];
	char name[MAXB], read[MAXB];
	while (gzgets(gzf, buffer, MAXB)) {
		sscanf(buffer, "%s %s", name, read);
		string read_name = string(name+1); 
		if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/')
			read_name = read_name.substr(0, read_name.size() - 2);
		mymap[read_name] = string(read);
	}

	genome_partition pt(read_file, threshold, mymap); 

	FILE *fo = fopen(out.c_str(), "wb");
	FILE *fidx = fopen((out + ".idx").c_str(), "wb");
	while (pt.has_next()) {
		auto p = pt.get_next();
		size_t i = pt.dump(p, fo);
		fwrite(&i, 1, sizeof(size_t), fidx);
	}
	fclose(fo);
	fclose(fidx);
	printf("%d\n", pt.get_cluster_id());
}

/********************************************************************/
void predict (const string &partition_file, const string &reference, const string &gtf, const bool barcodes, const string &range, const string &out_vcf, const string &out_full, 
					int k, int anchor_len, int min_support, int max_support, int uncertainty, int min_length, int max_length, const bool LOCAL_MODE, int ref_flank)
{
	kmistrvar predictor(k, anchor_len, partition_file, reference, gtf, barcodes);
	predictor.run_kmistrvar(range, out_vcf, out_full, min_support, max_support, uncertainty, min_length, max_length, LOCAL_MODE, ref_flank);
}

/********************************************************************/
bool isOEA (uint32_t flag)
{
	return (/*!flag ||*/ ((flag & 0x8) && (flag & 0x4)) 
				      || ((flag & 0x8) && !(flag & 0x4)) 
				      || (!(flag & 0x8) && (flag & 0x4)));
}

/********************************************************************/
bool isOEA_mrsFAST (uint32_t flag)
{
	return (/*!flag ||*/ ((flag & 0x8) && (flag & 0x4)) 
				      || ((flag & 0x8) && !(flag & 0x4)) 
				      || (!(flag & 0x8) && (flag & 0x4)));
}

/********************************************************************/
void extractOEA (const string &path, const string &result, bool fasta = false) 
{
	FILE *fi = fopen(path.c_str(), "rb");
	char mc[2];
	fread(mc, 1, 2, fi);
	fclose(fi);

	Parser *parser;
	if (mc[0] == char(0x1f) && mc[1] == char(0x8b)) 
		parser = new BAMParser(path);
	else
		parser = new SAMParser(path);

	string comment = parser->readComment();

	gzFile output = gzopen(result.c_str(), "wb6");
	//if (!fasta)
	//	gzwrite(output, comment.c_str(), comment.size());
	while (parser->hasNext()) {
		const Record &rc = parser->next();
		if (isOEA(rc.getMappingFlag())) {
			string record;
			if (fasta) 
				record = S("@%s\n%s\n+\n%s\n", rc.getReadName(), rc.getSequence(), rc.getQuality());
			else {
				char *readName = (char*)rc.getReadName();
				int l = strlen(readName);
				if (l > 1 && readName[l - 2] == '/')
					readName[l - 2] = '\0';
				record = S("%s %d %s %d %d %s %s %d %d %s %s %s\n",
            		readName,
            		rc.getMappingFlag(),
            		rc.getChromosome(),
            		rc.getLocation(),
            		rc.getMappingQuality(),
            		rc.getCigar(),
            		rc.getPairChromosome(),
            		rc.getPairLocation(),
            		rc.getTemplateLength(),
            		rc.getSequence(),
            		rc.getQuality(),
            		rc.getOptional()
        		);
			}
			gzwrite(output, record.c_str(), record.size());
		}
		parser->readNext();
	}

	gzclose(output);
	delete parser;
}

/********************************************************************/
void mask (const string &repeats, const string &path, const string &result, int pad = 0, bool invert = false)
{
	const int MAX_BUF_SIZE = 2000;
	size_t b, e, m = 0;
	char ref[MAX_BUF_SIZE], name[MAX_BUF_SIZE], l[MAX_BUF_SIZE];
	
	map<string, vector<pair<size_t, size_t>>> masks;
	FILE *fi = fopen(repeats.c_str(), "r");
	int prev = 0; string prev_ref = "";
	while (fscanf(fi, "%s %lu %lu %s", ref, &b, &e, name) != EOF) {
		if (e - b < 2 * pad) continue;
		masks[ref].push_back({b + pad, e - pad});
	}
	for (auto &r: masks) 
		sort(r.second.begin(), r.second.end());
	if (invert) for (auto &r: masks) 
	{
		vector<pair<size_t, size_t>> v;
		size_t prev = 0;
		for (auto &p: r.second) {
			if (prev < p.first)
				v.push_back({prev, p.first});
			prev = p.second;
		}
		v.push_back({prev, 500 * 1024 * 1024});
		r.second = v;
	}
	fclose(fi);
	
	const int MAX_CHR_SIZE = 300000000;
	char *x = new char[MAX_CHR_SIZE];

	fi = fopen(path.c_str(), "r");
	FILE *fo = fopen(result.c_str(), "w");
	fgets(l, MAX_BUF_SIZE, fi);
	
	while (true) {
		if (feof(fi))
			break;
		
		fprintf(fo, "%s", l);
		char *p = l; while (*p && !isspace(*p)) p++; *p = 0;
		string n = string(l + 1);
		printf("Parsed %s\n", n.c_str());
		
		size_t xlen = 0;
		while (fgets(l, MAX_BUF_SIZE, fi)) {
			if (l[0] == '>') 
				break;
			memcpy(x + xlen, l, strlen(l) - 1);
			xlen += strlen(l) - 1;
		}
		x[xlen] = 0;

		for (auto &r: masks[n]) {
			if (r.first >= xlen) continue;
			if (r.second >= xlen) r.second = xlen;
			memset(x + r.first, 'N', r.second - r.first);
		}

		for (size_t i = 0; i < xlen; i += 120) 
			fprintf(fo, "%.120s\n", x + i);
	}

	fclose(fi);
	fclose(fo);
	delete[] x;
}

/********************************************************************/
void removeUnmapped (const string &path, const string &result)
{
	ifstream fin(path.c_str());
	FILE *fo = fopen(result.c_str(), "w");
	
	string l, a, b, aa, bb; 
	while (getline(fin, a)) {
		getline(fin, b);
		getline(fin, l);
		getline(fin, l);

		int ll = a.size();
		if (ll > 1 && a[ll - 2] == '/')
			a = a.substr(0, ll - 2);
		if (a == aa)
			continue;
		else {
			fprintf(fo, "%s %s\n", aa.c_str(), bb.c_str());
			aa = a, bb = b;
		}
	}
	fprintf(fo, "%s %s\n", aa.c_str(), bb.c_str());

	fin.close();
	fclose(fo);
}

/********************************************************************/
void sortSAM (const string &path, const string &result) 
{
	sortFile(path, result, 2 * GB);
}

/********************************************************************/
void copy_string_reverse( char *src, char *dest)
{
	int limit=strlen(src);
	for( int i = 0; i < limit; i ++)
	{
		dest[i]=src[limit-1-i];
	}
	dest[limit]='\0';
}

/**********************************************/
void copy_string_rc( char *src, char *dest)
{
	int limit=strlen(src);
	for( int i = 0; i < limit; i ++)
	{
		switch( src[limit-1-i])
		{
			case 'A':
				dest[i]='T';
				break;
			case 'C':
				dest[i]='G';
				break;
			case 'G':
				dest[i]='C';
				break;
			case 'T':
				dest[i]='A';
				break;
			default:
				dest[i]='N';
				break;
		}
	}
	dest[limit]='\0';
}

/********************************************************************/
void extractMrsFASTOEA (const string &sam, const string &out) {
	string s;
	ifstream fin(sam.c_str());
	char name1[300], seq1[300], qual1[300];
	char name2[300], seq2[300], qual2[300];
	char seq1_new[300], qual1_new[300];
	char seq2_new[300], qual2_new[300];
	int flag1, flag2;
	FILE *fm = fopen(string(out + ".mapped.fq").c_str(), "wb");
	FILE *fu = fopen(string(out + ".unmapd.fq").c_str(), "wb");

	while (getline(fin, s)) {
		if (s[0] == '@') continue;
		sscanf(s.c_str(), "%s %d %*s %*d %*s %*s %*s %*s %*s %s %s",
				name1, &flag1, seq1, qual1);

		getline(fin, s);
		sscanf(s.c_str(), "%s %d %*s %*d %*s %*s %*s %*s %*s %s %s",
				name2, &flag2, seq2, qual2);

		if (flag1 == 0 && flag2 == 4)  {
			if (checkNs(seq2)) wo(fm, name1, seq1, qual1), wo(fu, name2, seq2, qual2);
		}
		if (flag1 == 16 && flag2 == 4) {
			copy_string_reverse(qual1, qual1_new);
			copy_string_rc(seq1, seq1_new);
			if (checkNs(seq2)) wo(fm, name1, seq1_new, qual1_new), wo(fu, name2, seq2, qual2);
		}
		if (flag1 == 4 && flag2 == 0)  {
			if (checkNs(seq1)) wo(fu, name1, seq1, qual1), wo(fm, name2, seq2, qual2);
		}
		if (flag1 == 4 && flag2 == 16) {
			copy_string_reverse(qual2, qual2_new);
			copy_string_rc(seq2, seq2_new);
			if (checkNs(seq1)) wo(fu, name1, seq1, qual1), wo(fm, name2, seq2_new, qual2_new);
		}
	}
	fclose(fm), fclose(fu);
	fin.close();
}

/********************************************************************/
void simulate_SVs (string in_file, string ref_file, bool from_bed) {

	simulator sim(ref_file);

	if(from_bed){
		sim.simulate(in_file);
	}
	else{
		sim.simulate_from_file(in_file);
	}
	
}

/********************************************************************/
void annotate (string gtf, string in_file, string out_file, bool genomic) {

	string line, chr, best_gene_s, best_trans_s, best_gene_e, best_trans_e, context_s, context_e, key, gene_id, trans_id;
	int start, end, t_start, t_end;
	ifstream infile (in_file);
	ofstream outfile (out_file);
	vector<uint32_t> vec_best_s, vec_best_e;
	vector<string> tokens;
	map <string,vector<isoform>> iso_gene_map;
	map <string,vector<gene_data>> gene_sorted_map;
	map <string,vector<uint32_t>> vec_all_s, vec_all_e;
	vector<string> contexts ={"intergenic", "intronic", "non-coding", "UTR", "exonic-CDS"};
	ensembl_Reader(gtf.c_str(), iso_gene_map, gene_sorted_map );

	if (infile.is_open()){
		while ( getline (infile,line)){

			tokens = split(line, '\t'); 
			chr = tokens[0];
			start = atoi(tokens[1].c_str());
			end = atoi(tokens[2].c_str());

			if(genomic){
				locate_interval(chr, start, start, gene_sorted_map[chr], 0, iso_gene_map, best_gene_s, best_trans_s, vec_best_s, vec_all_s);
				locate_interval(chr, end, end, gene_sorted_map[chr], 0, iso_gene_map, best_gene_e, best_trans_e, vec_best_e, vec_all_e);
			}
			else{
				gene_id = tokens[3];
				trans_id = tokens[4];

				locate_interval(chr, start, start, gene_id, trans_id, gene_sorted_map[chr], 0, iso_gene_map, best_gene_s, best_trans_s, vec_best_s, vec_all_s);
				locate_interval(chr, end, end, gene_id, trans_id, gene_sorted_map[chr], 0, iso_gene_map, best_gene_e, best_trans_e, vec_best_e, vec_all_e);
			}
			
			for (auto & match : vec_all_s){
				key = match.first;
				if(!vec_all_e[key].empty()){

					t_start = match.second[1];
					t_end = vec_all_e[key][1];

					context_s = contexts[vec_best_s[0]];
					context_e = contexts[vec_best_e[0]];

					if(t_start < t_end){
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\t" << t_start << "\t" << t_end << "\t" << context_s << "\t" << context_e << endl;
					}
					else if(t_start > t_end){
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\t" << t_end << "\t" << t_start << "\t" << context_e << "\t" << context_s << endl;
					}
					else{
						outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\tmissing\tmissing\tNA\tNA" << endl;
					}
				}
				else{
					outfile << chr << "\t" << start << "\t" << end << "\t" << key << "\tmissing\tmissing\tNA\tNA" << endl;
				}
			}

			vec_best_s.clear();
			vec_best_e.clear();
			vec_all_s.clear();
			vec_all_e.clear();
		}
		infile.close();
		outfile.close();
	}

	
}

/********************************************************************/

int main(int argc, char* argv[])
{
	try{

		cxxopts::Options options(argv[0], "CellFreeSV: Structural Variant Calling in cfDNA Sequencing Data");
		options.positional_help("[optional args]");
		string input_sam, reference, out_prefix, annotation;
		int threshold, k, a, s, S, u, m, M;
		bool barcodes = false;

		//[kmer-length] [min-support] [uncertainty] [local-assembly] [local-mode] [reference-flank]

		options.add_options()
			("i,input", "Input SAM file (required)", cxxopts::value<std::string>(), "FILE")
			("r,reference", "Reference file (required)", cxxopts::value<std::string>(), "FILE")
			("o,output", "Output prefix", cxxopts::value<std::string>()->default_value("out"), "PREFIX")
			("g,annotation", "GTF annotation file", cxxopts::value<std::string>()->default_value(""), "FILE")  //TODO handle old and new
			("b,barcodes", "Input reads contain barcodes", cxxopts::value<bool>(barcodes))
			("c", "Clustering threshold (default 1000)", cxxopts::value<int>()->default_value("1000"), "INT")
			("k", "Kmer length (default 14)", cxxopts::value<int>()->default_value("14"), "INT")
			("a", "Anchor length (default 40)", cxxopts::value<int>()->default_value("40"), "INT")
			("s", "Min Read Support (default 2)", cxxopts::value<int>()->default_value("2"), "INT")
			("S", "Max Read Support (default unlimited)", cxxopts::value<int>()->default_value("999999"), "INT")
			("u", "Uncertainty (default 8)", cxxopts::value<int>()->default_value("8"), "INT")
			("m", "Min SV length (default 10)", cxxopts::value<int>()->default_value("40"), "INT")
			("M", "Max SV length (default 20000)", cxxopts::value<int>()->default_value("20000"), "INT")
			("h,help", "Print help");


		options.parse(argc, argv);

		if (options.count("help")){
			std::cout << options.help({""}) << std::endl;
			exit(0);
		}

		if (options.count("input")){
			input_sam = options["input"].as<std::string>();

			ifstream f(input_sam.c_str());
			if(!f.good()){
				throw cxxopts::OptionException("Input file does not exist");
			}
		}
		else{
			throw cxxopts::OptionException("No input file specified");
		}

		if (options.count("reference")){
			reference = options["reference"].as<std::string>();

			ifstream f(reference.c_str());
			if(!f.good()){
				throw cxxopts::OptionException("Reference file does not exist");
			}
		}
		else{
			throw cxxopts::OptionException("No Reference file specified");
		}

		out_prefix = options["output"].as<std::string>();

		if (options.count("annotation")){
			annotation = options["annotation"].as<std::string>();

			ifstream f(annotation.c_str());
			if(!f.good()){
				throw cxxopts::OptionException("Annotation file does not exist");
			}
		}
		

		threshold = options["c"].as<int>();
		k = options["k"].as<int>();
		a = options["a"].as<int>();
		s = options["s"].as<int>();
		S = options["S"].as<int>();
		u = options["u"].as<int>();
		m = options["m"].as<int>();
		M = options["M"].as<int>();

		if(threshold < 0){
			throw cxxopts::OptionException("Cluster threshold must be a positive integer");
		}

		if(k < 5){
			throw cxxopts::OptionException("K must be greater than 5");
		}
		else if(k > a){
			throw cxxopts::OptionException("K must be <= anchor length");
		}
		
		if(s < 2){
			throw cxxopts::OptionException("Read support threshold must be >= 2");
		}

		if(s > S){
			throw cxxopts::OptionException("Max read support should be <= max read support");
		}

		if(u < 0){
			throw cxxopts::OptionException("Uncertainty must be a positive integer");
		}

		if(m <= u){
			throw cxxopts::OptionException("Min SV length should be > uncertainty");
		}

		if(m > M){
			throw cxxopts::OptionException("Min SV length should be <= max SV length");
		}

		std::cout << "Arguments remain = " << argc << std::endl;

		cout << "Running with parameters: k=" << k << " a=" << a << " s=" << s << " u=" << u << " m=" << m << " M=" << M << endl;

		predict(input_sam, reference, annotation, barcodes, "0-999999999", (out_prefix + ".vcf"), (out_prefix + ".out"), k, a, s, S, u, m, M, 0, 0);

	} catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}


// int main(int argc, char **argv)
// {
// 	try {
// 		if (argc < 2) throw "Usage:\tCellFreeSV [mode=(?)]";

// 		string mode = argv[1];
// 		if (mode == "fastq") {
// 			if (argc < 4) throw "Usage:\tCellFreeSV fastq [sam-file] [output]";
// 			extractOEA(argv[2], argv[3], argc == 4 ? true : false);
// 		}
// 		else if (mode == "oea") {
// 			if (argc != 4) throw "Usage:\tCellFreeSV oea [sam-file] [output]";
// 			extractMrsFASTOEA(argv[2], argv[3]);
// 		}
// 		else if (mode == "mask" || mode == "maski") {
// 			if (argc != 6) throw "Usage:\tCellFreeSV mask/maski [repeat-file] [reference] [output] [padding]";
// 			mask(argv[2], argv[3], argv[4], atoi(argv[5]), mode == "maski");
// 		}
// 		else if (mode == "sort") {
// 			if (argc != 4) throw "Usage:\tCellFreeSV sort [sam-file] [output]";
// 			sortSAM(argv[2], argv[3]);
// 		}
// 		else if (mode == "rm_unmap") {
// 			if (argc != 4) throw "Usage:\tCellFreeSV rm_unmap [fq-file] [output]";
// 			removeUnmapped(argv[2], argv[3]);
// 		}
// 		else if (mode == "partition") {
// 			if (argc != 6) throw "Usage:\tCellFreeSV partition [read-file] [mate-file] [output-file] [threshold]";
// 			partify(argv[2], argv[3], argv[4], atoi(argv[5]));
// 		}
// 		else if (mode == "predict") {
// 			if (argc != 14) throw "Usage:\tCellFreeSV predict [partition-file] [reference] [gtf] [range] [output-file-vcf] [output-file-full] [kmer-length] [min-support] [uncertainty] [local-assembly] [local-mode] [reference-flank]"; 
// 			predict(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]));
// 		}
// 		else if (mode == "get_cluster") {
// 			if (argc != 4) throw "Usage:\tCellFreeSV get_cluster [partition-file] [range]";
// 			genome_partition pt;
// 			pt.output_partition( argv[2], argv[3]);
// 		}
// 		else if (mode == "simulate_SVs") {
// 			if (argc != 5) throw "Usage:\tCellFreeSV simulate_SVs [bed-or-SV-file] [ref-file] [is-bed]";
// 			simulate_SVs(argv[2], argv[3], atoi(argv[4]));
// 		}
// 		else if (mode == "annotate") {
// 			if (argc != 6) throw "Usage:\tCellFreeSV annotate [gtf-file] [bed-file] [out-file] [is-genomic]";
// 			annotate(argv[2], argv[3], argv[4], atoi(argv[5]));
// 		}
// 		else {
// 			throw "Invalid mode selected";
// 		}
// 	}
// 	catch (const char *e) {
// 		ERROR("Error: %s\n", e);
// 		exit(1);
// 	}
		
// 	return 0;
// }
