#ifndef __VARIANT_CALLER__
#define __VARIANT_CALLER__

#include <string>
#include "partition.h"
#include "common.h"
#include "assembler.h"
#include "genome.h"

using namespace std;

class variant_caller
{
protected:

	const int MAX_INDEL = 0;
	const int MIN_SV_LEN_DEFAULT = 10;
	const int MAX_ASSEMBLY_RANGE = 5000;
	const int OVERLAP_NEW = 50; //50
	const int DEBUG = 1410;//16034742;

	assembler as; 
	genome ref;
	genome_partition pt;
	string in_file;

public:

	

protected:

	void  init();

public:
	
	variant_caller(const int assembler_overlap, const string &input_file, const string &reference);
	~variant_caller();

};
#endif
