#include <string>
#include <cstdio>
#include <cstdlib>
#include <locale.h>
#include <getopt.h>
#include <libgen.h>
#include <sys/time.h>
#include <queue>
#include <utility>
#include <functional>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <zlib.h>
#include <climits>
#include <unordered_map>
#include <unistd.h>
#include <cassert>
#include "common.h"
#include "file.h"
#include "array.h"
#include "sort.h"
using namespace std;


unordered_map<string, int> chromosomes;
struct SAMNode {
	int chr;
	size_t pos;
	char *data;
	size_t data_sz;

	SAMNode(): chr(0), pos(0), data(0), data_sz(0) {
	}

	static bool sortComp (const SAMNode &x, const SAMNode &y) {
		return x.chr < y.chr || (x.chr == y.chr && x.pos < y.pos);
	}
	bool operator< (const SAMNode &s) const {
		return chr > s.chr || (chr == s.chr && pos > s.pos);
	}

	ssize_t readSAM(char *data, size_t sz) {
		this->data = data;
		data_sz = 0;
		
		while (data_sz < sz && data[data_sz] != '\t') data_sz++; data_sz++;
		while (data_sz < sz && data[data_sz] != '\t') data_sz++; data_sz++;
		size_t e = data_sz;
		while (data_sz < sz && data[data_sz] != '\t') data_sz++; 
		if (data_sz >= sz) return -1;
		string chrs = string(data + e, data_sz - e);
		e = ++data_sz;
		while (data_sz < sz && data[data_sz] != '\t') data_sz++;
		if (data_sz >= sz) return -1;
		data[data_sz] = 0; pos = atoi(data + e); data[data_sz] = '\t'; data_sz++;

		while (data_sz < sz && data[data_sz] != '\n') data_sz++;
		if (data_sz >= sz || data[data_sz] != '\n') return -1;
		
		if (string(chrs) == "*")
			chr = INT_MAX;
		else {
			unordered_map<string, int>::iterator it = chromosomes.find(chrs);
			if (it != chromosomes.end())
				chr = it->second;
			if (it == chromosomes.end()) {
				chr = chromosomes.size();
				chromosomes[chrs] = chr;
			}
		}	
		return ++data_sz;
	}

	ssize_t readBAM(char *data, size_t sz) {
		this->data = data;
		data_sz = *(uint32_t*)data + 4;
		if (data_sz > sz) return -1;
		int32_t *di = (int32_t*)(data + 4);
	 	chr = di[0] == -1 ? INT_MAX : di[0];
	 	pos = di[1] + 1;
		return data_sz;
	}
};

static bool isBAM;
static char *buffer;
static size_t bufsz;

size_t mergeSort (file **f, size_t fsz, file *fo, char *buffer, size_t bufsz) {
	LOG("Merging %lu files ...", fsz);
	size_t cnt = 0;
	size_t bsz = bufsz / fsz;
	vector<size_t> counts(fsz, 1);
	vector<size_t> offsets(fsz, 0);
	vector<size_t> lenghts(fsz, 0);
	
	priority_queue<pair<SAMNode, int> > pq;
	for (int fi = 0; fi < fsz; fi++) 
		pq.push(make_pair(SAMNode(), fi));

	while (!pq.empty()) {
		pair<SAMNode, int> p = pq.top(); pq.pop();
		fo->write(p.first.data, p.first.data_sz);
		
		int fi = p.second;
		if (--counts[fi] == 0) {
			memmove(buffer + fi * bsz, 
				buffer + fi * bsz + lenghts[fi] - offsets[fi], 
				offsets[fi]);
			lenghts[fi] = offsets[fi] + f[fi]->read(
				buffer + fi * bsz + offsets[fi], 
				bsz - offsets[fi]);
			ssize_t i;
			for (i = 0; i < lenghts[fi]; ) {
				SAMNode n;
				ssize_t p;
				if (isBAM)
					p = n.readBAM(buffer + fi * bsz + i, lenghts[fi] - i);
				else
					p = n.readSAM(buffer + fi * bsz + i, lenghts[fi] - i);
				if (p == -1) 
					break;
				pq.push(make_pair(n, fi));
				i += p, counts[fi]++, cnt++;
			}
			offsets[fi] = lenghts[fi] - i;
		}
	}
	return cnt;
} 

int detectFileType (const string &path) {
	FILE *fx = fopen(path.c_str(), "rb");
	char mc[4];
	fread(mc, 1, 4, fx);
	fclose(fx);
	if (mc[0] == char(0x1f) && mc[1] == char(0x8b)) 
		return 1;
	else if (((*(uint32_t*)mc) & 0xffffff00) == (MAGIC & 0xffffff00))
		return 2;
	else 
		return 0;
}

void sortFile (const string &path, const string &pathNew, size_t memLimit) {
	int ft = detectFileType(path);
	if (ft == 1) isBAM = true;
	if (ft == 2)
		exit (1);
		//throw DZException("File %s is DZ file, and it is already sorted", path.c_str());

	file *finput;
	vector<char> header;
	if (isBAM) {
		finput = new gzfile(path.c_str(), "rb");

		size_t ho = 0;

		header.resize(8);
		finput->read(header.data(), 4), ho += 4;
		finput->read(header.data() + ho, 4), ho += 4;
		int32_t len = *(int32_t*)(header.data() + ho - 4);
		header.resize(12 + len);
		finput->read(header.data() + ho, len), ho += len;
		
		finput->read(header.data() + ho, 4), ho += 4;
		int32_t chromosomesCount = *(int32_t*)(header.data() + ho - 4);
		for (int i = 0; i < chromosomesCount; i++) {
			int32_t clen;
			finput->read(&clen, 4);
			header.resize(ho + 8 + clen);
			copy((char*)&clen, (char*)&clen + 4, header.data() + ho), ho += 4;
			finput->read(header.data() + ho, clen), ho += clen;
			finput->read(header.data() + ho, 4), ho += 4;
		}
	}
	else {
		rawfile *f = new rawfile(path.c_str(), "rb");
		char *s = 0; size_t slen = 0;
		ssize_t len = 0, fpos = 0;
		while ((len = getline(&s, &slen, f->fh())) != -1) {
			if (s[0] != '@') {
				fseek(f->fh(), fpos, SEEK_SET);
				break;
			}
			header.insert(header.end(), s, s + len);
			fpos += len;
		}	
		finput = f;
	}
	bufsz = memLimit;
	buffer = (char*)malloc(memLimit + 1);

	vector<file*> files;
	vector<string> fileNames;

	int fi = 0, fp = 0;
	size_t offset = 0;
	Array<SAMNode> nodes(0, MB);
	while (!finput->eof()) {
		size_t sz = finput->read(buffer + offset, bufsz - offset) + offset;
		//DEBUG(">>>> %'llu %'llu", sz, bufsz - offset);

		ssize_t i;
		for (i = 0; i < sz; ) {
			SAMNode n;
			ssize_t p;
			if (isBAM)
				p = n.readBAM(buffer + i, sz - i);
			else
				p = n.readSAM(buffer + i, sz - i);
			if (p == -1) break;
			nodes.add(n);
			i += p;
		}

		//ZAMAN_START();
		sort(nodes.data(), nodes.data() + nodes.size(), SAMNode::sortComp);
		//radixSort(nodes, 0, 0, nodes.size());
		//ZAMAN_END("SORT");

		char fn[100];
		//ZAMAN_START();
		snprintf(fn, 100, "%s_%d.%02d", path.c_str(), fp, fi++);

		file *f;
		if (isBAM)
			f = new gzfile(fn, "wb1+");
		else
			f = new rawfile(fn, "wb+");
		for (int i = 0; i < nodes.size(); i++) 
			f->write(nodes.data()[i].data, nodes.data()[i].data_sz);
		f->close();
		if (isBAM)
			f = new gzfile(fn, "rb");
		else
			f = new rawfile(fn, "rb");

		files.push_back(f);
		fileNames.push_back(fn);
		//ZAMAN_END("FLUSH");

		memmove(buffer, buffer + i, offset = sz - i);
		
		LOG("Created %s with %'lu records", fn, nodes.size());
		nodes.resize(0);
	}

	finput->close();
	delete finput;

	int fsz = files.size();
	while (files.size() > 1) {
		fp++;
		fi = 0;
		vector<file*> nf;
		vector<string> nfn;
		for (int i = 0; i < files.size(); i += fsz) {
			char fn[100];
			snprintf(fn, 100, "%s_%d.%02d", path.c_str(), fp, fi++);
			file *f;
			if (isBAM)
				f = new gzfile(fn, "wb1+");
			else
				f = new rawfile(fn, "wb+");
			//if (files.size() <= fsz)
			//	f->write(header.data(), header.size());
			size_t sz;
			//ZAMAN_START();
			sz = mergeSort(files.data() + i, 
				min(fsz, (int)files.size() - i), 
				f, buffer, bufsz);
			f->close();
			if (isBAM)
				f = new gzfile(fn, "rb");
			else
				f = new rawfile(fn, "rb");
			//ZAMAN_END("SORT");
			nf.push_back(f);
			nfn.push_back(fn);

			//LOG("> %s ... %d", fn, sz);
		}
		for (int i = 0; i < files.size(); i++) {
			files[i]->close();
			delete files[i];
			unlink(fileNames[i].c_str());
		}
		files = nf;
		fileNames = nfn;
	}

	free(buffer);

	assert(files.size() == 1);
	files[0]->close();
	rename(fileNames[0].c_str(), pathNew.c_str());
}
