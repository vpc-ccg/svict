#ifndef SAMParser_H
#define SAMParser_H

#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>

#include "common.h"
#include "parser.h"
#include "record.h"

class SAMParser: public Parser {
	FILE *input;
	
    size_t file_size;

    Record currentRecord;

public:
	SAMParser (const std::string &filename);
	~SAMParser (void);

public:
	std::string readComment (void);
	bool readNext ();
	bool hasNext (void);
	size_t fpos (void);
	size_t fsize (void);

public:
	void parse (Record &line);
	std::string head (void);
	Record next (void) ;
};

#endif
