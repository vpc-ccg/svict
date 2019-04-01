#ifndef __SMOOTHER__
#define __SMOOTHER__

#include "fiboheap.h"
#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#define L(c,...) fprintf(stdout,c,##__VA_ARGS__)
#define E(c,...) fprintf(stderr,c,##__VA_ARGS__)

class smoother 
{

public:

	struct Xread {
		vector<int> clusters;
	};

	struct cluster {
		int id;
		int orig_id;
		int support, unresolved;
		vector<int> reads;
		bool operator< (const cluster& other) const {
			return make_pair(support, make_pair(unresolved, id)) > make_pair(other.support, make_pair(other.unresolved, other.id));
		}
	};

public:

	smoother();
	~smoother(void);
	void set_cover(vector<smoother::cluster>& clusters, vector<smoother::Xread>& reads, unordered_map<string, int>& read_to_id);

};

#endif
