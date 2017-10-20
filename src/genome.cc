#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "genome.h"

using namespace std;

genome::genome(string filename)
{	
	file = filename;
	fin.open(file.c_str());
	reference.reserve(300000000);
	reference_name = "";
	char ch;

	fin.get(ch);
	load_next();
}

genome::~genome()
{
}

void genome::load_next(void)
{	
	reference.clear();
	if (fin.eof())
		return;

	getline(fin, reference_name);

	iss.str(reference_name);
    getline(iss, reference_name, ' ');
    iss.clear();

	string tmp;
	char ch;

	fin.get(ch);
	
	while (ch != '>' && !fin.eof())
	{
		reference+=ch;
		getline(fin, tmp);
		reference+=tmp;
		fin.get(ch);
	}
	transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
}

void genome::reset(void)
{	
	fin.close();
	fin.open(file.c_str());
	reference.clear();
	reference.reserve(300000000);
	reference_name = "";
	char ch;

	fin.get(ch);
	load_next();
}

string genome::extract(const string &rname, int start, int end)
{
	while (rname != reference_name) {
		fprintf(stderr, "Chr %s done\n",  reference_name.c_str(), rname.c_str());
		string last = reference_name;
		load_next();
		if(last == reference_name){
			fprintf(stderr, "Reached end of reference file.\n");
			return "";
		}
	}
	
	start = max(start, 1);
	end = max(end, 1);
	end = min(end, (int)reference.size());
	return reference.substr(start-1, end-start+1);
}

int genome::get_ref_length(const string &rname){
	
	while (rname != reference_name) {
		fprintf(stderr, "%s != %s\n",  reference_name.c_str(), rname.c_str());
		string last = reference_name;
		load_next();
		if(last == reference_name){
			fprintf(stderr, "Reached end of reference file.\n");
			return -1;
		}
	}

	return (int)reference.size();
}

