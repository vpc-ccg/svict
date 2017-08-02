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
	void load_next(void);

public:
	genome(string);
	~genome();
	void reset(void);
	string extract(const string&, int, int);
	int get_ref_length(const string&);
};

#endif
