#ifndef __COMMON__
#define __COMMON__
#include <string>
#include <vector>

#define KB  1024LL
#define MB  KB * 1024LL
#define GB  MB * 1024LL

#define VERSION	0x10
#define MAGIC 	(0x07445A00 | VERSION)

#define ERROR(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define LOG(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define MAX_CHAR 1000000

using namespace std;

void wo (FILE *, char *, char *, char *);
char checkNs (char *);
string reverse_complement ( const string & );
string reverse ( const string & );
string S (const char* fmt, ...);
string random_seq(int length);
int check_AT_GC(const string &, const double &);
vector<string> split(const string &s, char delim);

inline char getDNAValue (char ch) {
	#define _ 0
	static char c[128] = {
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 0 15
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 16 31
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 32 47
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 48 63
		_,0,_,1,_,_,_,2,_,_,_,_,_,_,_,_,// 64 79
		_,_,_,_,3,_,_,_,_,_,_,_,_,_,_,_,// 80 95
		_,0,_,1,_,_,_,2,_,_,_,_,_,_,_,_,// 96 111
		_,_,_,_,3,_,_,_,_,_,_,_,_,_,_,_,// 112 127
	};	
	#undef _
	return c[ch];
}

#endif
