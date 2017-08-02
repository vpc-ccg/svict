#include<iostream>
#include<set>
#include<vector>
#include<utility>
#include "assembler_old.h"
#include<cassert>
using namespace std;

inline char DNA(char c) {
	return (c=='T'?3:(c=='G'?2:(c=='C'?1:0)));
}


assembler_old::assembler_old (int mcs = 2000, int mgs = 10):
	max_contig_size (mcs), min_glue_size(mgs)
{
	initialize(max_contig_size, prime_seed);
}


assembler_old::~assembler_old (void) {
	for (int i = 0; i < max_contig_size; i++) {
		delete[] hash_p[i];
		delete[] hash_s[i];
	}
	delete[] hash_p;
	delete[] hash_s;
}


void assembler_old::initialize (int mcs, int ps) {
	hash_p = new set<sread*>*[mcs];
	hash_s = new set<sread*>*[mcs];
	for (int i = 0; i < mcs; i++) {
		hash_p[i] = new set<sread*>[ps];
		hash_s[i] = new set<sread*>[ps];
	}
}


void assembler_old::update (sread *r, char mode) {
	int h = 0, p;
	for (int i = 0; i < r->data.size(); i++) {
		h = ((h * 4) % prime_seed + DNA(r->data[i])) % prime_seed;
		p = i;
		assert(i < max_contig_size) ;
		if (!mode)
			hash_p[p][h].insert(r);
		else
			hash_p[p][h].erase(r);
	}

	int exp = 1;
	h = 0;
	for (int i = r->data.size() - 1; i >= 0; i--) {
		h = ((DNA(r->data[i]) * exp) % prime_seed + h) % prime_seed;
		p = r->data.size() - 1 - i;

		assert(p >= 0 && p < max_contig_size);
		if (!mode)
			hash_s[p][h].insert(r);
		else
			hash_s[p][h].erase(r);
		exp = (exp * 4) % prime_seed;
	}
}


bool assembler_old::full_compare (const string &a, const string &b, int glue_sz) {
	for (int i = 0; i < glue_sz + 1; i++)
		if (a[i] != b[b.size() - 1 - glue_sz + i]) 
			return 0;
	return 1;
}

bool assembler_old::assemble_single (int sz, int ha) {
	bool ret = 0;
	vector<sread*> rm;
	for (auto ix = hash_p[sz][ha].begin(); ix != hash_p[sz][ha].end(); ix++) if (!(*ix)->used) {
		auto iy = hash_s[sz][ha].begin();
		for (; iy != hash_s[sz][ha].end(); iy++) if (!(*iy)->used) {
			if (*iy == *ix) 
				continue;
			if (full_compare((*ix)->data, (*iy)->data, sz) 
					&& (*iy)->data.size() + (*ix)->data.size() - sz - 2 < max_contig_size)
				break;
		}
		if (iy == hash_s[sz][ha].end())
			continue;

		sread *nr = new sread();
		nr->data = (*iy)->data + string((*ix)->data.c_str() + sz + 1);
		nr->reads.insert(nr->reads.end(), (*iy)->reads.begin(), (*iy)->reads.end());
		nr->reads.insert(nr->reads.end(), (*ix)->reads.begin(), (*ix)->reads.end());
		for (int i = (*iy)->reads.size(); i < (*iy)->reads.size() + (*ix)->reads.size(); i++)
			nr->reads[i].first.second += (*iy)->data.size() - sz - 1;
		reads.insert(nr);
		update(nr, 0);

		rm.push_back(*ix);
		rm.push_back(*iy);
		(*ix)->used = 1;
		(*iy)->used = 1;
		ret = 1;
	}

	for (int i = 0; i < rm.size(); i++) {
		update(rm[i], 1);
		reads.erase(rm[i]);
		delete rm[i];
	}

	return ret;
}


void assembler_old::assemble (void) {
	while (1) {
		bool no_br = 0;
		// all prefix/suffix lengths
		for (int i = max_contig_size - 1; i >= min_glue_size; i--) { 
			// check pref=sfx
			for (int h = 0; h < prime_seed; h++) 
				if (hash_p[i][h].size() && hash_s[i][h].size()) {
					while (assemble_single(i, h))
						no_br = 1;
				}
		}
		if (!no_br) 
			break;
	}
}


vector<contig> assembler_old::assemble (const vector<pair<pair<string, string>, int>> &input) {
	for (int i = 0; i < input.size(); i++) {
		sread *r = new sread();
		reads.insert(r);
		r->data = input[i].first.second;
		r->reads.push_back({{input[i].second, 0}, input[i].first});
		update(r, 0);
	}
	assemble();

	vector<contig> result(reads.size());
	int ri = 0;
	for (auto it = reads.begin(); it != reads.end(); it++) {
		result[ri].data = (*it)->data;
		result[ri].read_information.resize((*it)->reads.size());
		for (int i = 0; i < (*it)->reads.size(); i++) {
			result[ri].read_information[i].location = (*it)->reads[i].first.first;
			result[ri].read_information[i].location_in_contig = (*it)->reads[i].first.second;
			result[ri].read_information[i].name = (*it)->reads[i].second.first;
			result[ri].read_information[i].seq = (*it)->reads[i].second.second;
		}
		ri++;

	//	result[ri++] = make_pair((*it)->data, (*it)->reads);
		update(*it, 1);
		delete *it;
	}
	reads.clear();

	return result;
}

