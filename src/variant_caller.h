#ifndef __VARIANT_CALLER__
#define __VARIANT_CALLER__

#include <string>
#include "partition.h"
#include "common.h"
#include "assembler.h"
#include "assembler_old.h"
#include "assembler_ext.h"
#include "genome.h"

using namespace std;

class variant_caller
{
protected:

	const int MAX_INDEL = 0;
	const int MIN_SV_LEN_DEFAULT = 10;
	const int MAX_ASSEMBLY_RANGE = 5000;
	const int OVERLAP_NEW = 50; //50
	const int OVERLAP_OLD = 20; //20
	const int DEBUG = 1410;//16034742;

	assembler as; 
	assembler_old as_old; 
	genome ref;
	genome_partition pt;
	string part_file;

public:

	

protected:

	void  init();

public:
	
	variant_caller(const string &partition_file, const string &reference);
	~variant_caller();

};
#endif
