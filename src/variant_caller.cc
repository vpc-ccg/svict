#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <cctype>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include "variant_caller.h"


using namespace std;

variant_caller::variant_caller(const string &partition_file, const string &reference) : 
		as(OVERLAP_NEW), ref(reference.c_str()){

	part_file = partition_file;
}

variant_caller::~variant_caller(){


}

void init(){


}

