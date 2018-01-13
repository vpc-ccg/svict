#include "common.h"
#include "record.h"
#include "sam_parser.h"

#include <assert.h>
using namespace std;

SAMParser::SAMParser (const string &filename)
{
	Parser::fname = filename;

	input = fopen(filename.c_str(), "r");
	if (input == NULL)	
		{	exit(1);	}

	fseek(input, 0L, SEEK_END);
	file_size = ftell(input);
	fseek(input, 0L, SEEK_SET);
}

SAMParser::~SAMParser (void) 
{
	fclose(input);
}

string SAMParser::readComment (void)  
{
	string s;
	while ((currentRecord.lineLength = getline(&currentRecord.line, &currentRecord.lineSize, input)) != -1) {
		if (currentRecord.line[0] != '@') {
			parse(currentRecord);
			break;
		} else {
			s += currentRecord.line;
		}
	}
	return s;
}

bool SAMParser::readNext ()  
{
	if ((currentRecord.lineLength = getline(&currentRecord.line, &currentRecord.lineSize, input)) != -1) {
		assert(currentRecord.line[0] != '@');
		parse(currentRecord);
		return true;
	}
	return false;
}

bool SAMParser::hasNext (void) 
{
	return !feof(input);
}

size_t SAMParser::fpos (void) 
{
	return ftell(input);
}

size_t SAMParser::fsize (void) 
{
	return file_size;
}

void SAMParser::parse (Record &record) 
{
	char *line = &record.line[0];
	int l = strlen(line) - 1;
	while (l && (line[l] == '\r' || line[l] == '\n'))
		line[l--] = 0;
	record.strFields[0] = 0;
	char *c = line;
	int f = 1, sfc = 1, ifc = 0;
	while (*c) {
		if (*c == '\t') {
			if (f == 1 || f == 3 || f == 4 || f == 7 || f == 8)
				record.intFields[ifc++] = atoi(c + 1);
			else
				record.strFields[sfc++] = (c + 1) - line;
			f++;
			*c = 0;
			if (f == 12) {
				c++;
				while (*c) {
					//if (*c == '\t') *c = 0;
					c++;
				}
				break;
			}
		}
		c++;
	}	
	if (f == 11)
		record.strFields[sfc++] = c - line + 1;
	
	record.intFields[Record::IntField::LOC]--;
	record.intFields[Record::IntField::P_LOC]--;
}

Record SAMParser::next (void) 
{
	return std::move(currentRecord);
}

string SAMParser::head (void) 
{
	return currentRecord.getChromosome();
}
