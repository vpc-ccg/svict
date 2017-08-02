#ifndef __LOGGER__
#define __LOGGER__
#include<string>

using namespace std;

class logger
{
private:
	FILE *output;

public:
	logger(string);
	~logger();
	void log(const char*, ...);
};


#endif
