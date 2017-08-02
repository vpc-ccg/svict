#ifndef __ASSEMBLER_EXT__
#define __ASSEMBLER_EXT__

#include<iostream>
#include<string>
#include<vector>
#include<utility>
#include<set>
#include "assembler.h"

using namespace std;

class assembler_ext {
public:
	vector<contig> assemble (const vector<pair<string, string>> &input);
};

#endif
