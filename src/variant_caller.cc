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

variant_caller::variant_caller(const string &input_file, const string &reference) : 
		as(OVERLAP_NEW), ref(reference.c_str()){

	in_file = input_file;
}

variant_caller::~variant_caller(){


}

void init(){


}

