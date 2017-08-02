#ifndef Parser_H
#define Parser_H

#include <stdio.h>
#include <string>

#include "record.h"

class Parser {
public:
	virtual ~Parser() {};
	virtual std::string readComment (void) = 0;
	virtual bool readNext (void) = 0;
	virtual bool hasNext (void) = 0;
	virtual size_t fpos (void) = 0;
	virtual size_t fsize (void) = 0;
	virtual const Record &next (void) = 0;
	virtual std::string head (void) = 0;
};

#endif
