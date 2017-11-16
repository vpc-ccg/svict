#ifndef Parser_H
#define Parser_H

#include <stdio.h>
#include <string>

#include "common.h"
#include "record.h"

class Parser {
protected:
	std::string fname;
public:
	virtual ~Parser() {};

	std::string fileName() const { return this->fname; }

	virtual std::string readComment (void) = 0;
	virtual bool hasNext (void) = 0;
	virtual size_t fpos (void) = 0;
	virtual size_t fsize (void) = 0;
	virtual std::string head (void) = 0;
	virtual Record next (void) = 0;
	virtual bool readNext () = 0;

};

#endif
