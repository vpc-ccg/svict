#ifndef __ASSEMBLEROLD__
#define __ASSEMBLEROLD__

#include<iostream>
#include<string>
#include<vector>
#include<utility>
#include<set>
#include "assembler.h"

using namespace std;
/*
struct Read {
	string name;
	string seq;
	int location;
	int location_in_contig;
};

struct contig {
	string data;
	vector<Read> read_information;

	contig (void) {
	}

	int support (void) const { 
		return read_information.size();
	}

	int start = -1;
	int get_start () {
		if (start > -1) return start;
		start = read_information[0].location_in_contig;
		for (int i = 1; i < read_information.size(); i++)
			start = min(start,  read_information[i].location_in_contig);
		return start;
	}

	int end = -1;
	int get_end () {
		if (end > -1) return end;
		end = read_information[0].location_in_contig + read_information[0].seq.size();
		for (int i = 1; i < read_information.size(); i++)
			end = max(end, read_information[i].location_in_contig + (int)read_information[i].seq.size());
		return end;
	}
};*/

class assembler_old {
private:
	struct sread {
		// first one real loc, second one in-contig loc
		vector<pair<pair<int, int>, pair<string, string>>> reads;
		string data;
		bool used;

		sread (void): 
			used(0), data("") {}
	};

private:
	const int max_contig_size;
	const int min_glue_size;
	static const int prime_seed=157;

	set<sread*> **hash_p;
	set<sread*> **hash_s;
	set<sread*> reads;

private:
	void initialize (int mcs, int ps);
	void update (sread *r, char mode);
	bool full_compare (const string &a, const string &b, int glue_sz);
	bool assemble_single (int sz, int ha);
	void assemble (void);
	void print ();

public:
	assembler_old (int, int);
	~assembler_old (void);

public:
	vector<contig> assemble (const vector<pair<pair<string, string>, int>> &input);
};

#endif
