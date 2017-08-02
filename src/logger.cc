#include <iostream>
#include <string>
#include <cstdio>
#include <cctype>
#include <cstdarg>
#include "logger.h"

using namespace std;

logger::logger(string name = "")
{
	if (name != "")
		output = fopen("name", "w");
	else
		output = stdout;
}

logger::~logger()
{
	if (output != stdout)
		fclose(output);
}

void logger::log(const char* format, ...)
{
	va_list args;
	va_start (args, format);
	vfprintf(output, format, args);
	va_end(args);
	fflush(output);
}


