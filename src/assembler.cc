#include <iostream>
#include <set>
#include <vector>
#include <utility>
#include <cassert>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#include "common.h"
#include "assembler.h"

using namespace std;

assembler::assembler(): 
	assembler(20) 
{
}

assembler::assembler(int min_glue_size):
	 min_glue_size(min_glue_size)
{
}

assembler::~assembler(void) 
{	
}

// prefix of A == suffix of B
bool assembler::validate(const string &a, const string &b, int sz) 
{
	for (int i = 0; i < sz + 1; i++) {
		if (a[i] != b[b.size() - 1 - sz + i]) 
			return false;
	}
	return true;
}

vector<contig> assembler::assemble(vector<pair<string, string>>& reads) 
{

	auto lens = set<int>();
	for (auto &r: reads) {
		lens.insert(r.second.size());
	}

	auto s = set<string>(), R = set<string>();
	unordered_map<string, string> RN = unordered_map<string, string>();
/*	for (auto &r: reads) {
		for (int i = 0; i < r.size(); i++) {
			for (auto rl: lens) {
				if (i + rl >= r.size()) break;
				s.insert(r.substr(i, rl));
			}
		}
	}*/

	for (auto &r: reads) {
		//if (s.find(r) == s.end()) {
			R.insert(r.second);
			RN[r.second] = r.first;	
		//	s.insert(r);
		//}
	}
	vector<string> read_seqs = vector<string>(R.begin(), R.end());

	int min_len = *lens.begin();
	int max_len = *lens.rbegin();

	vector<unordered_map<int, vector<int>>> phash(max_len), shash(max_len);
	graph = Graph(read_seqs.size());

	for (int ri = 0; ri < read_seqs.size(); ri++) {
		auto &r = read_seqs[ri];
		for (int i = 0, hp = 0, hs = 0, exp = 1; i < r.size(); i++, exp = (exp << 2) % SEED) {
			hp = ((hp * 4) % SEED + getDNAValue(r[i])) % SEED;
			hs = ((getDNAValue(r[r.size() - 1 - i]) * exp) % SEED + hs) % SEED;
	
			if (i < min_glue_size) 
				continue;
			//if (i >= read_len) // assume only one contig is crazy
			//	break;
	
			for (auto ni: phash[i][hs]) {
				if (validate(read_seqs[ni], r, i)) 
					graph[ri].neighbors.push_back({ni, i});
			}
			for (auto ni: shash[i][hp]) {
				if (validate(r, read_seqs[ni], i)) 
					graph[ni].neighbors.push_back({ri, i});
			}
			phash[i][hp].push_back(ri);
			shash[i][hs].push_back(ri);
		}
		graph[ri].name = RN[read_seqs[ri]];
		graph[ri].seq = r;
	}

	auto result = path();
	//E("After assembly");

	return result;
}

vector<int> assembler::topsort() 
{
	vector<int> top;

	int no_edges = 0;
	for (auto &v: graph) {
		v.indegree = 0;
	}
	for (auto &v: graph) {
		int outdegree = 0;
		for (auto &n: v.neighbors) if (n.first != -1) {
			graph[n.first].indegree++, outdegree++;
		}
		no_edges += outdegree;
	}

	vector<bool> processed(graph.size(), false);
	int deleted = 0;
	while (no_edges > 0) {
		//E("in while");
		vector<int> zeros;
		for (auto &v: graph) {
			int vi = &v - &graph[0];
			if (!processed[vi] && v.indegree == 0) {
				processed[vi] = true;
				zeros.push_back(vi);
			}
		}
		//E("after first for");
		for (int zi = 0; zi < zeros.size(); zi++) {
			int z = zeros[zi];
			top.push_back(z);

			for (auto &n: graph[z].neighbors) {
				int ni = n.first;
				if (ni == -1 || processed[ni]) continue;
				graph[ni].indegree--;
				if (graph[ni].indegree == 0) {
					processed[ni] = true;
					zeros.push_back(ni);
				}
				no_edges--;
			}
		}
		//E("after second for");
		if (no_edges != 0) {
			//E("Not DAG, remaining {} edges", no_edges);
			tuple<int, int, int, int> minw(99999999, 99999999, 99999999, 99999999);
			//E("before for graph");
			for (auto &v: graph) {
				int vi = &v - &graph[0];
				//printf("vi: %d\n",vi);
				if (processed[vi]) continue;
				// delete for 
				for (auto &n: v.neighbors) {
					int ni = n.first;
					if(ni == -1)continue;
					int nix = &n - &v.neighbors[0];
					auto yy = make_tuple(  graph[ni].indegree, n.second, vi, nix);
					if (ni != -1 && !processed[ni] && yy < minw) {
						minw = yy;
					}
				}
			}
			//E("after for graph");
			int fr = get<2>(minw);
			int to = graph[fr].neighbors[get<3>(minw)].first;
			//E("  [{:2}] deleted: {}->{}, weight {}, indegree {}",
			//	++deleted, fr, to, get<1>(minw), get<0>(minw));

			graph[to].indegree--;
			no_edges--;
			graph[fr].neighbors[get<3>(minw)] = {-1, -1};
			//E("Before topsort in if");
		}
		//E("Before topsort iafter if");
		end:;
	}
	//E("After topsort");
	return top;
}

vector<vector<pair<int, int>>> max_path(const Graph &graph, vector<int> top) 
{
	vector<bool> visited(graph.size(), false);
	vector<vector<pair<int, int>>> paths;
	while (true) {
		vector<int> len(graph.size(), 0);
		vector<int> parent(graph.size(), -1);

		// maximize number of reads in the thing?!
		for (auto &vi: top) if (!visited[vi]) {
			for (auto &n: graph[vi].neighbors) {
				int ni = n.first;
				if (ni == -1) continue;
				int val = n.second ; //+ graph[ni].seq.size();
				if (!visited[ni] && len[ni] <= len[vi] + val) {
					len[ni] = len[vi] + val;
					parent[ni] = vi;
				}
			}
		}

		int max_len = -1, max_ni = -1;
		for (int ni = 0; ni < len.size(); ni++) {
			if (!visited[ni] && len[ni] > max_len) {
				max_len = len[ni], max_ni = ni;
			}
		}
		if (max_ni == -1) break;

		vector<pair<int, int>> path;
		while (parent[max_ni] != -1) {
			path.push_back({max_ni, len[max_ni] - len[parent[max_ni]]});// - graph[max_ni].seq.size()});
			visited[max_ni] = true;
			max_ni = parent[max_ni];
		}
		visited[max_ni] = true;
		path.push_back({max_ni, 0});
		reverse(path.begin(), path.end());
		paths.push_back(path);
	}

	return paths;
}

vector<contig> assembler::path() 
{
	// E("{}", graph.size());
	auto top = topsort();
	auto paths = max_path(graph, top);

	vector<contig> result;

	int maxpath = 1000;
	for (auto &path: paths) {
		contig c;

		string pad = "", seq = graph[path[0].first].seq.substr(0, 1);
		int prev_len = 1;
		for (auto &t: path) {
			pad += string(prev_len - t.second - 1, ' ');
			seq += graph[t.first].seq.substr(t.second + 1);
			//E("{:3} {:5}: {}{}", t.first, t.second, pad, graph[t.first].seq);
			prev_len = graph[t.first].seq.size();
			c.read_information.push_back({graph[t.first].name, graph[t.first].seq, 0, (int)pad.size()});
		}

		c.data = seq;
		//E("{:3} reads, {:5}: {}", path.size(), seq.size(), seq.substr(0, 100) + (seq.size() > 100 ? "..." : ""));
		//fmt::print("{} {}\n", seq.size(), seq);
		//E("");

		result.push_back(c);
	}

	return result;
}
