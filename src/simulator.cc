#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include "simulator.h"

using namespace std;

simulator::simulator(const string& ref_file) :
	ref(ref_file.c_str())
{	
	reference = ref_file;
}

simulator::~simulator()
{
}

void simulator::shuffle(int *arr, size_t n)
{
	if (n > 1) 
	{
		size_t i;
		srand(time(NULL));
		for (i = 0; i < n - 1; i++) 
		{
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = arr[j];
			arr[j] = arr[i];
			arr[i] = t;
		}
	}
}

void simulator::simulate_from_file(const string& sv_file)
{
	create_reference(reference, sv_file, 0);
}

void simulator::simulate(const string& bed_file, int min_small, int max_small, int min_large, int max_large, int min_offset, int max_offset)
{
	generate_SVs(bed_file, min_small, max_small, min_large, max_large, min_offset, max_offset);
	create_reference(reference, "", 0);
	create_reference(reference, "", 1);
}

void simulator::generate_SVs(const string& bed_file, int min_small, int max_small, int min_large, int max_large, int min_offset, int max_offset)
{
	string line, type, cur_type, chr, chr2, seq, strand, strand2, cur_seq, new_gen, gen, homo_tag;
	string cur_chr = "NA";
	int start, start2, end, end2, len, len2, sv_len, bp1, bp2, offset;
	int cur_SV = 0;
	int sv_id = 0;
	int last_bp2 = 0;
	int small, tandem, homo;
	ifstream infile(bed_file);
	vector<SV> small_trans;
	vector<SV> large_trans;
	vector<string> tokens;

	while (getline(infile,line)){

		tokens = split(line, '\t'); 
		chr = tokens[0];
		start = stoi(tokens[1]);
		end = stoi(tokens[2]);
		strand = tokens[3];
		len = end-start;
		if(chr != cur_chr)last_bp2 = 0;
		cur_chr = chr;

		string type = sv_types[(rand() % 7)];
		string sub_type = "INTERSPERSED";
		small = (rand() % 2);
		tandem = (rand() % 2);
		homo = (rand() % 10);
		homo_tag = homo ? "HOMO" : "HETRO";

		if(small){
			sv_len = (rand() % (max_small-min_small) + min_small);
			offset = (rand() % (max_offset-min_offset) + min_offset);
			bp1 = (rand() % len + start);
			//bp1 = (bp1 < last_bp2) ? last_bp2+50 : bp1;
			if(bp1 < last_bp2){
				cerr << "small " << bp1 << " " << last_bp2 << endl;
				continue;
			}

			if(bp1-start > end-bp1 && type != "INS" && type != "DUP"){
				bp1 -= sv_len;
				bp1 = (bp1 < last_bp2) ? last_bp2+50 : bp1;
			}
			
			bp2 = bp1+sv_len;
			last_bp2 = bp2;

			if(type == "INS"){
				seq = random_seq(sv_len);
				bp2 = bp1;
				cout << type << "\tMICRO\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << seq << "\t" << sv_id++ << endl;
			}
			else if(type == "DUP"){
				bp2 = bp1;
				if(tandem){
					sub_type = "TANDEM";
					offset = 0;
				}

				if(bp1-start > end-bp1){
					seq = ref.extract(chr,(bp1-sv_len-offset),(bp1-offset));
					cout << type << "\t" << sub_type << "\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp1 << "\t" << seq << "\t" << chr << "\t" << (bp1-sv_len-offset) << "\t" << (bp1-offset) << "\t" << sv_id++ << endl;
				}
				else{
					seq = ref.extract(chr,(bp1+offset),(bp1+sv_len+offset));
					cout << type << "\t" << sub_type << "\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp1 << "\t" << seq << "\t" << chr << "\t" << (bp1+offset) << "\t" << (bp1+sv_len+offset) << "\t" << sv_id++ << endl;
				}
			}
			else if(type == "TRANS"){
				small_trans.push_back({chr, bp1, bp2, homo, "TRANS", ref.extract(chr,bp1,bp2)});
				continue;
			}
			else if(type == "DEL"){
				seq = "N";
				cout << type << "\tMICRO\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << seq << "\t" << sv_id++ << endl;
			}
			else if(type == "INV"){
				seq = reverse_complement(ref.extract(chr,bp1,bp2));
				cout << type << "\tMICRO\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << seq << "\t" << sv_id++ << endl;
			
			}
			else{
				continue;
			}

			last_bp2 = bp2;
			SVs[chr].push_back({chr, bp1, bp2, homo, type, seq}); 
		}
		else{

			if(!getline(infile,line))break;

			tokens = split(line, '\t'); 
			chr2 = tokens[0];
			start2 = stoi(tokens[1]);
			end2 = stoi(tokens[2]);
			strand2 = tokens[3];
			len2 = end-start;

			sv_len = (rand() % (max_large-min_large) + min_large);
			bp1 = (rand() % len + start);
			//bp1 = (bp1 < last_bp2) ? last_bp2+50 : bp1;
			if(bp1 < last_bp2){
			cerr << "large " << bp1 << " " << last_bp2 << endl;
				continue;
			}

			if(chr == chr2 && (start2-end) < max_large){
				bp2 = (rand() % len2 + start2);
			}
			else{
				bp2 = bp1 + sv_len;
			}
			last_bp2 = bp2;

			seq = ref.extract(chr,bp1,bp2);

			if(type == "INS"){
				large_trans.push_back({chr, bp1, bp2, homo, "INS", seq});
				continue;
			}
			else if(type == "DUP"){
				if(tandem){
					cout << type << "\tTANDEM\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp1 << "\t" << seq << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << sv_id++ << endl;
					bp2 = bp1;
					last_bp2 = bp2;
				}
				else{
					large_trans.push_back({chr, bp1, bp2, homo, "DUP", seq});
					continue;
				}
			}
			else if(type == "TRANS"){
				large_trans.push_back({chr, bp1, bp2, homo, "TRANS", seq});
				continue;
			}
			else if(type == "DEL"){
				//large_trans.push_back({chr, bp1, bp2, homo, "DEL", seq});
				cout << type << "\tLONG\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << seq << "\t" << sv_id++ << endl;
			}
			else if(type == "INV"){
				seq = reverse_complement(seq);
				cout << type << "\tLONG\t" << homo_tag << "\t" << chr << "\t" << bp1 << "\t" << bp2 << "\t" << seq << "\t" << sv_id++ << endl;
			}
			else{
				continue;
			}

			SVs[chr].push_back({chr, bp1, bp2, homo, type, seq});
		}
	}

	type = "TRANS";

	//Resolve translocations
	for(small = 0; small <= 1; small++){

		int num_elem = small ? small_trans.size() : large_trans.size();
		int ind[num_elem];
		SV trans1;
		SV trans2;
		string type1, type2, seq1;

		for (int i=0; i < num_elem; i++){
			ind[i] = i;
		}
		shuffle(ind, num_elem);

		for(int i=0; i < num_elem-1; i+=2){

			trans1 = small ? small_trans[ind[i]] : large_trans[ind[i]];
			trans2 = small ? small_trans[ind[i+1]] : large_trans[ind[i+1]];
			type1 = trans1.type;
			type2 = trans2.type;
			seq1 = trans1.seq;
			homo_tag = (trans1.homo + trans2.homo) ? "HOMO" : "HETRO";
			vector<SV> to_create;
			vector<SV>::iterator it;

			trans1.seq = trans2.seq;
			trans1.homo += trans2.homo;
			
			trans2.seq = seq1;
			trans2.homo += trans1.homo;

			to_create.push_back(trans1);
			to_create.push_back(trans2);
			
			if(type1 == "INS"){
				to_create[1].type = "DEL";
				cout << type << "\t" << "INS-DEL" << "\t" << homo_tag << "\t" << trans1.chr << "\t" << trans1.bp1 << "\t" << trans1.bp1 << "\t" << trans1.seq 
						<< "\t" << trans2.chr << "\t" << trans2.bp1 << "\t" << (trans2.bp1+trans1.seq.length()) << "\t" << sv_id++ << endl;
			}
			else if(type1 == "DUP"){
				to_create.pop_back();
				cout << type << "\t" << "DUP" << "\t" << homo_tag << "\t" << trans1.chr << "\t" << trans1.bp1 << "\t" << trans1.bp1 << "\t" << trans1.seq 
						<< "\t" << trans2.chr << "\t" << trans2.bp1 << "\t" << (trans2.bp1+trans1.seq.length()) << "\t" << sv_id++ << endl;
			}
			else if(type1 == "TRANS"){
				to_create[1].type = "TRANS";
				cout << type << "\t" << "SWAP" << "\t" << homo_tag << "\t" << trans1.chr << "\t" << trans1.bp1 << "\t" << trans1.bp2 << "\t" << trans1.seq 
						<< "\t" << trans2.chr << "\t" << trans2.bp1 << "\t" << (trans2.bp1+trans1.seq.length()) << "\t" << sv_id++ << endl;
				cout << type << "\t" << "SWAP" << "\t" << homo_tag << "\t" << trans2.chr << "\t" << trans2.bp1 << "\t" << trans2.bp2 << "\t" << trans2.seq 
						<< "\t" << trans1.chr << "\t" << trans1.bp1 << "\t" << (trans1.bp1+trans2.seq.length()) << "\t" << sv_id++ << endl;
			}
			
			for(auto& trans : to_create){
				for(it=SVs[trans.chr].begin() ; it < SVs[trans.chr].end(); it++) {
					if(trans.bp1 < (*it).bp1) {
						SVs[trans.chr].insert(it,trans);
						break;
					}
				}

				if(it == SVs[trans.chr].end()){
					SVs[trans.chr].push_back(trans);
				}
			}
		}
	}

	cur_chr = "NA"; 
	ref.reset();

	infile.close();

}


void simulator::create_reference(const string& ref_file, const string& sv_file, bool hetro)
{

	string line, type, chr, seq, new_gen, gen;
	int len;
	int cur_SV = 0;
	string suffix = hetro ? ".allele2.fa" : ".allele1.fa";
	ofstream outfile((ref_file + suffix));

	if(sv_file != ""){
		ifstream infile(sv_file);
		vector<string> tokens;

		while (getline(infile,line)){
			tokens = split(line, ' '); 
			SVs[tokens[3]].push_back({tokens[3], atoi(tokens[4].c_str()), atoi(tokens[5].c_str()), atoi(tokens[2].c_str()), tokens[0], tokens[6]});
		}

		infile.close();
	}

	for(auto& chr : chromos){

		vector<SV> cur_SVs = SVs[chr];
		new_gen = "";
		gen = ref.extract(chr, 1, 300000000);

		for(int i = 0; i < gen.length(); i++){

			if(cur_SV < cur_SVs.size() && i == cur_SVs[cur_SV].bp1){ //TODO check for off-by-one
//if(!hetro)cerr << i << "/" << gen.length() << " "	<< cur_SV << "/" << cur_SVs.size() << " cur: " << cur_SVs[cur_SV].bp1 << endl;	
				if((cur_SVs[cur_SV].homo || !hetro)){	
					type = cur_SVs[cur_SV].type;
					seq = cur_SVs[cur_SV].seq;
					len = (cur_SVs[cur_SV].bp2-cur_SVs[cur_SV].bp1)-1;

//if(!hetro)cout << "SANITY\t" << type << "\t" << chr << "\t" << cur_SVs[cur_SV].bp1 << "\t" << cur_SVs[cur_SV].bp2 << "\t" << len << "\t" << seq << endl;
				
					if(type == "INS" || type == "DUP"){
						new_gen += seq;
						i--;
					}
					else if(type == "DEL"){
						i += len;
					}
					else if(type == "INV" || type == "TRANS"){
						new_gen += seq;
						i += len;
					}
				}
				else{
					new_gen += gen[i];
				}
				cur_SV++;
//if(!hetro)cerr << i << "/" << gen.length() << " next: " << cur_SVs[cur_SV].bp1 << endl;	
			}
			else{
				new_gen += gen[i];
			}
		}

		int pos = 0;

		outfile << ">" << chr << endl;
		while(pos < new_gen.length()){
			len = (pos+60 < new_gen.length()) ? 60 : (new_gen.length()-pos);
			outfile << new_gen.substr(pos, len) << endl;
			pos += 60;
		}

		cur_SV = 0;		
	}

	outfile.close();
	ref.reset();
}