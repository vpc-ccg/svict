#include "common.h"
#include "record.h"
#include "parser.h"
#include "sam_parser.h"


#include <assert.h>
using namespace std;

SAMParser::SAMParser (const string &filename) {
	input = fopen(filename.c_str(), "r");
	if (input == NULL)	
		exit(1);

	fseek(input, 0L, SEEK_END);
	file_size = ftell(input);
	fseek(input, 0L, SEEK_SET);
}

SAMParser::~SAMParser () {
	fclose(input);
}

string SAMParser::readComment (void)  {
	string s;
	while (fgets(currentRecord.line, MAXLEN, input)) 
		if (currentRecord.line[0] != '@') {
			parse();
			break;
		}
		else s += currentRecord.line;
	return s;
}

bool SAMParser::readNext (void)  {
	if (fgets(currentRecord.line, MAXLEN, input)) {
		assert(currentRecord.line[0] != '@');
		parse();
		return true;
	}
	return false;
}


bool SAMParser::hasNext (void) {
	return !feof(input);
}

size_t SAMParser::fpos (void) {
	return ftell(input);
}

size_t SAMParser::fsize (void) {
	return file_size;
}

void SAMParser::parse (void) {
	int l = strlen(currentRecord.line) - 1;
	while (l && (currentRecord.line[l] == '\r' || currentRecord.line[l] == '\n'))
		currentRecord.line[l--] = 0;
	char *c = currentRecord.strFields[0] = currentRecord.line;
	int f = 1, sfc = 1, ifc = 0;
	while (*c) {
		if (*c == '\t') {
			if (f == 1 || f == 3 || f == 4 || f == 7 || f == 8)
				currentRecord.intFields[ifc++] = atoi(c + 1);
			else
				currentRecord.strFields[sfc++] = c + 1;
		f++;
			*c = 0;
			if (f == 12) break;
		}
		c++;
	}	
	if (f == 11)
		currentRecord.strFields[sfc++] = c;
}

const Record &SAMParser::next (void) {
	return currentRecord;
}

string SAMParser::head (void) {
	return currentRecord.getChromosome();
}


