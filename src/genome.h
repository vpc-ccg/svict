#ifndef __GENOME__
#define __GENOME__

#include <iostream>
#include <sstream>

using namespace std;

class genome
{
private:
	ifstream fin;
	stringstream iss;
	string file;
	string reference_name;
	string reference;

public:
	genome(string);
	~genome();
	void load_next(void);
	void reset(void);
	string extract(const string&, int, int);
	int get_ref_length(const string&);
	string get_ref_name();
};

#endif
