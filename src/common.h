#ifndef __COMMON__
#define __COMMON__
#include <assert.h>
#include <inttypes.h>
#include <limits>
#include <string>
#include <sstream>
#include <cstdio>
#include <vector>
#include <sys/time.h>
#include "logger.h"

#define KB  1024LL
#define MB  KB * 1024LL
#define GB  MB * 1024LL

#define VERSION	0x10
#define MAGIC 	(0x07445A00 | VERSION)

#define ERROR(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define LOG(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define HELP(c,...)\
	fprintf(stdout, c"\n", ##__VA_ARGS__)
#define MAX_CHAR 1000000

using namespace std;

void wo (FILE *, char *, char *, char *);
char checkNs (char *);
string reverse_complement ( const string & );
string reverse ( const string & );
string S (const char* fmt, ...);
string random_seq(int length);
int check_AT_GC(const string &, const double &);
//vector<string> split(const string &s, char delim);

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

/*************************************************/
inline std::vector<std::string> split (std::string s, char delim) 
{
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> elems;
	while (std::getline(ss, item, delim)) 
		elems.push_back(item);
	return elems;
}

/*************************************************/
inline uint64_t zaman() 
{
	struct timeval t;
	gettimeofday(&t, 0);
	return (t.tv_sec * 1000000ll + t.tv_usec);
}

/*************************************************/
// #define ZAMAN
#ifdef ZAMAN
	class __zaman__ { // crazy hack for gcc 5.1
	public:
		std::map<std::string, uint64_t> times;
		std::string prefix;

	public:
		static std::map<std::string, uint64_t> times_global;
		static std::string prefix_global;
		static std::mutex mtx;
	};
	extern thread_local __zaman__ __zaman_thread__;

	#define ZAMAN_VAR(s) \
		__zaman_time_##s
	#define ZAMAN_START(s) \
		int64_t ZAMAN_VAR(s) = zaman(); \
		__zaman_thread__.prefix += std::string(#s) + "_"; 
	#define ZAMAN_END(s) \
		__zaman_thread__.times[__zaman_thread__.prefix] += (zaman() - ZAMAN_VAR(s)); \
		__zaman_thread__.prefix = __zaman_thread__.prefix.substr(0, __zaman_thread__.prefix.size() - 1 - strlen(#s)); 
	#define ZAMAN_THREAD_JOIN() \
		{ 	std::lock_guard<std::mutex> __l__(__zaman__::mtx); \
			for (auto &s: __zaman_thread__.times) \
				__zaman__::times_global[__zaman__::prefix_global + s.first] += s.second; \
			__zaman_thread__.times.clear(); \
		}
	#define ZAMAN_START_P(s) \
		int64_t ZAMAN_VAR(s) = zaman(); \
		__zaman__::prefix_global += std::string(#s) + "_"; 
	#define ZAMAN_END_P(s) \
		__zaman__::times_global[__zaman__::prefix_global] += (zaman() - ZAMAN_VAR(s)); \
		__zaman__::prefix_global = __zaman__::prefix_global.substr(0, __zaman__::prefix_global.size() - 1 - strlen(#s)); \

	inline void ZAMAN_REPORT() 
	{ 
		using namespace std;
		ZAMAN_THREAD_JOIN();
		vector<string> p;
		for (auto &tt: __zaman__::times_global) { 
			string s = tt.first; 
			auto f = split(s, '_');
			s = "";
			for (int i = 0; i < f.size(); i++) {
				if (i < p.size() && p[i] == f[i])
					s += "  ";
				else
					s += "/" + f[i];
			}
			p = f;
			LOG("  %-40s: %'8.1lfs", s.c_str(), tt.second/1000000.0);
		}
	}
#else
	#define ZAMAN_VAR(s)	
	#define ZAMAN_START(s) 
	#define ZAMAN_END(s)
	#define ZAMAN_THREAD_JOIN()
	#define ZAMAN_START_P(s)
	#define ZAMAN_END_P(s)
	#define ZAMAN_REPORT()
#endif


#endif
