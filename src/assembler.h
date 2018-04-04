#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <set>
using namespace std;

struct Read {
	string name;
	string seq;
	int location;
	int location_in_contig;
};

struct compressed_contig {
	vector<bool> data;
	vector<unsigned short> coverage;
};

struct contig {
	string data;
	vector<Read> read_information;

	int support (void) const { 
		return read_information.size();
	}
};

struct Node {
	string name;
	string seq;
	int indegree;
	vector<pair<int, int>> neighbors; // Read, Weight
	Node(): indegree(0) {}
};

struct compare {
    bool operator()(const string& first, const string& second) {
        return (first.size() == second.size()) ? first > second : first.size() >= second.size();
    }
   // !(a < b) && !(b < a)
};

typedef vector<Node> Graph;

class assembler {
private:
	static const int SEED = 92821; // must be prime
	static const int PRIME = 92821;//787;
	static const int K = 23;//787;
	static const int MAX_READS = 5000;//787;
	static const int READ_LEN = 301;//787;
	const int min_glue_size;
	long long power;
	long long** hashes;

	Graph graph;

private:
	bool validate(const string &a, const string &b, int glue_sz);
	bool validate2(const string &a, const string &b, int pos1, int pos2, int sz);
	bool validate3(const string &a, const string &b, int pos1, int sz);
	long long pow(int x, int y);
	long long rolling_hash(const string& str, int start, int end);

public:
	assembler();
	assembler(int);
	~assembler();
	vector<contig> path();
	vector<int> topsort();

public:
	vector<contig> assemble(vector<pair<string, string>>& reads);
};
