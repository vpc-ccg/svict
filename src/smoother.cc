/// 786
#include "smoother.h"
#include <bits/stdc++.h>

using namespace std;

smoother::smoother(){
}

smoother::~smoother(void){	
}


void smoother::set_cover(vector<smoother::cluster>& clusters, vector<smoother::Xread>& reads, unordered_map<string, int>& read_to_id) 
{
	FiboHeap<cluster*> heap;
	unordered_map<int, FiboNode<cluster*>* > heap_handles;

	// inserting them all
	unordered_map<int,bool> x;
	for (auto &c: clusters) {
		if (c.unresolved == 0) continue;
		heap_handles[c.id] = heap.insert(&c);
		x[c.id] = false;
	}

	// assert that all shared clusters are in set cover
	int re = 0;
	for (auto &r: reads) {
		assert(r.clusters.size() >= 1);
		if (r.clusters.size() > 1) {
			re++;
			for (auto &c: r.clusters) {
				assert(heap_handles.find(c) != heap_handles.end());
			}
		}
	}

	unordered_set<int> W;

	while (!heap.isEmpty()) {
		cluster *c = heap.getMinimum(); heap.removeMinimum();
		//cout << c->id << '\t' << c->support <<'\t' << c->unresolved << '\t' << c->reads.size()<<'\n';
		
		if (!c->unresolved && !c->support) break;
		//E("{}-- {} {}", c->id, c->support, c->unresolved);

		c->unresolved = 0; // we will resolve all reads!
		heap_handles.erase(c->id);
		x[c->id] = true;

		//if (!c->unresolved) continue;

		for (auto &r: c->reads) {
			if (reads[r].clusters.size() <= 1) continue;

			for (auto &cr: reads[r].clusters) {
				unordered_map<int, FiboNode<cluster*>* >::const_iterator it = heap_handles.find(cr);
				if (it == heap_handles.end()) continue; // not in heap anymore

				assert(((it->second)->getValue())->id == cr);
				assert(((it->second)->getValue())->unresolved > 0);

				cluster* newc = it->second->getValue();
				newc->unresolved--;
				heap.update(it->second);
			}
			reads[r].clusters.clear();
			reads[r].clusters.push_back(c->id);
		}
	}

	for (auto &c: clusters) {
		c.reads.clear();
	}

	unordered_map<int, string> read_names;
	for (auto &r: read_to_id)
		read_names[r.second] = r.first;
	for (auto &r: reads) {
		assert(r.clusters.size() == 1);
		clusters[r.clusters[0]].reads.push_back(&r - &reads[0]);
	}
	int i = 0, rx = 0;
	for (auto &c: clusters) {
		if (c.reads.size() == 0) {
			//L("Removed: %d\n", c.orig_id);
			continue;
		}

		//L("%d %u\n", c.orig_id, c.reads.size());
		for (auto &r: c.reads) {
			//L("{}", read_names[r]);
		}
		
		rx += c.reads.size(), i++;
	}
	assert(rx == reads.size());
	//E("%d / %u sets resolved, %d / %d reads resolved\n", i, clusters.size(), re, rx);
}
