#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdarg>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "assembler_ext.h"

using namespace std;

inline std::string S (const char* fmt, ...) {
	char *ptr = 0;
    va_list args;
    va_start(args, fmt);
    vasprintf(&ptr, fmt, args);
    va_end(args);
    std::string s = ptr;
    free(ptr);
    return s;
}

vector<contig> assembler_ext::assemble (const vector<pair<string, string>> &input) {
	FILE *fo = fopen("_data.fq", "w");
	for (int i = 0; i < input.size(); i++) {
		fprintf(fo, "@%s\n%s\n+\n%s\n", input[i].first.c_str(), input[i].second.c_str(), string(input[i].second.size(), 'A').c_str());
	}
	fclose(fo);

	mkdir("tmp", 0777);
	string command = S("python spades.py --only-assembler -s %s -o spades_tmp", "_data.fq");

	system("velveth tmp 29 -short -fastq _data.fq >/dev/null");
	system("velvetg tmp >/dev/null");
	
	vector<contig> result;
	ifstream fin("tmp/contigs.fa");
	string l;
	if (!getline(fin, l))
		return result;
	result.push_back(contig());
	result.back().read_information.resize(1);
	while (getline(fin, l)) {
		if (l[0] == '>') {
			result.push_back(contig());
			result.back().read_information.resize(1);
		}
		else
			result.back().data += l;
	}
	fin.close();

	return result;
}

