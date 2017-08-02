#ifndef __SIMULATOR__
#define __SIMULATOR__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "common.h"
#include "genome.h"

using namespace std;

struct SV {
	string chr;
	int bp1;
	int bp2;
	int homo;
	string type;
	string seq;
};


class simulator
{
private:

	const string sv_types[7] = { "DEL", "INS", "INV", "DUP", "TRANS", "NONE", "NONE"};
	const vector<string> chromos = {"1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"}; //Venter has no MT

	unordered_map<string,vector<SV>> SVs;
	genome ref;
	string reference;

	void shuffle(int *arr, size_t n);
	void generate_SVs(const string& bed_file, int min_small, int max_small, int min_large, int max_large, int min_offset, int max_offset);
	void create_reference(const string& ref_file, const string& sv_file, bool hetro);

public:
	simulator(const string& ref_file);
	~simulator();
	void simulate_from_file(const string& sv_file);
	void simulate(const string& bed_file, int min_small = 10, int max_small = 1000, int min_large = 1000, int max_large = 10000, int min_offset = 10, int max_offset = 200);
};

#endif