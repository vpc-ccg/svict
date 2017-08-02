#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include "common.h"
#include "record.h"
#include "bam_parser.h"


using namespace std;

BAMParser::BAMParser (const string &filename):
	dataSize(0), chromosomesCount(0), data(0), chromosomes(0)
{
	fd = fopen(filename.c_str(), "rb");
	fseek(fd, 0L, SEEK_END);
	file_size = ftell(fd);
	fseek(fd, 0L, SEEK_SET);

	input = gzdopen(fileno(fd), "rb");
	if (input == Z_NULL)	
		exit(1);

	char magic[5] = {0};
	gzread(input, magic, 4);
	assert(!strcmp(magic, "BAM\x1"));
}

BAMParser::~BAMParser () {
	for (int i = 0; i < chromosomesCount; i++)
		free(chromosomes[i]);
	chromosomes--;
	free(chromosomes);
	gzclose(input);
}

string BAMParser::readComment (void)  {
	int32_t len;
	gzread(input, &len, sizeof(int32_t));
	char *c = (char*)calloc(len + 1, 1);
	gzread(input, c, len);
	strncpy(currentRecord.line, c, len + 1);
	string s = c;
	free(c);

	readChromosomeInformation();
	readNext();

	return s;
}

void BAMParser::readChromosomeInformation (void) {
	gzread(input, &chromosomesCount, sizeof(int32_t));
	chromosomes = (char**)malloc((chromosomesCount + 1) * sizeof(char*));
	chromosomes[0] = (char*)"*";
	chromosomes++;
	for (int i = 0; i < chromosomesCount; i++) {
		int32_t len, chrlen;
		gzread(input, &len, sizeof(int32_t));
		chromosomes[i] = (char*)calloc(len + 1, 1);
		gzread(input, chromosomes[i], len);
		gzread(input, &chrlen, sizeof(int32_t));
	}
}

bool BAMParser::readNext (void) {
	int32_t bsize, cc;
	char *buf = currentRecord.line;
	if (gzread(input, &bsize, 4) != 4) 
		return false;

	assert(bsize < MAXLEN);
	if (bsize > dataSize)
		data = (char*)realloc(data, dataSize = bsize + KB);
	gzread(input, data, bsize);

	int32_t *di = (int32_t*)data;

	int bin = di[2] >> 16;
	int l_read_name = di[2] & 0xff;
	int n_cigar_op = di[3] & 0xffff;
	int l_seq = di[4];

	// rn
	strncpy(buf, data + 8 * 4, l_read_name);
	currentRecord.strFields[RN] = buf;
	buf += l_read_name;

	// flag
 	currentRecord.intFields[MF] = di[3] >> 16;

 	// chr
 	string chr = chromosomes[di[0]];
	strncpy(buf, chr.c_str(), chr.size() + 1);
	currentRecord.strFields[CHR] = buf;
	buf += chr.size() + 1;	 	

 	// loc
 	currentRecord.intFields[LOC] = di[1] + 1;

 	// mq
 	currentRecord.intFields[MQ] = (di[2] >> 8) & 0xff;
 	
 	// cigar
 	uint32_t *op = (uint32_t*)(data + 8 * 4 + l_read_name);
 	cc = 0;
	for (int i = 0; i < n_cigar_op; i++)
		cc += sprintf(buf + cc, "%d%c", op[i] >> 4, "MIDNSHP=X"[op[i] & 0xf]);
	if (n_cigar_op == 0)
		buf[cc++] = '*';
	n_cigar_op *= 4;
	buf[cc++] = 0;
 	currentRecord.strFields[CIGAR] = buf;
	buf += cc;

 	// p_chr
 	string pe_chr = chromosomes[di[5]];
 	if (pe_chr != "*" && pe_chr == chr)
 		pe_chr = "=";
 	strncpy(buf, pe_chr.c_str(), pe_chr.size() + 1);
	currentRecord.strFields[P_CHR] = buf;
	buf += pe_chr.size() + 1;

	// p_loc
 	currentRecord.intFields[P_LOC] = di[6] + 1;
 	
 	// tlen
 	currentRecord.intFields[TLEN] = di[7];

 	// seq
 	char *sq = data + 8 * 4 + l_read_name + n_cigar_op;
 	cc = 0;
 	for (int i = 0; i < l_seq; i++)
 		buf[cc++] = "=ACMGRSVTWYHKDBN"[(sq[i / 2] >> ((1 - i % 2) * 4)) & 0xf];
 	if (l_seq == 0)
 		buf[cc++] = '*';
 	buf[cc++] = 0;
 	currentRecord.strFields[SEQ] = buf;
	buf += cc;
 	
 	// qual
 	char *q = data + 8 * 4 + l_read_name + n_cigar_op + (l_seq + 1) / 2;
 	cc = 0;
 	for (int i = 0; i < l_seq; i++)
 		buf[cc++] = q[i] + 33;
 	if (l_seq == 0)
 		buf[cc++] = '*';
 	buf[cc++] = 0;
 	currentRecord.strFields[QUAL] = buf;
	buf += cc;

	currentRecord.strFields[OPT] = buf;

 	// optional data ...
 	int pos = 8 * 4 + l_read_name + n_cigar_op + (l_seq + 1) / 2 + l_seq;
 	if (pos >= bsize)
 		currentRecord.strFields[OPT][0] = 0;
 	else 
 		currentRecord.strFields[OPT]++; // avoid \t
	while (pos < bsize) {
		*buf = '\t', buf++;
		char t = data[pos + 2]; 
		buf += sprintf(buf, "%c%c:", data[pos], data[pos+1]); 
		pos += 3;
		switch (t) {
			case 'A': buf += sprintf(buf, "A:%c", data[pos]); pos++; break;
			case 'c': buf += sprintf(buf, "i:%d", *(int8_t*)(data + pos)); pos += 1; break;
			case 'C': buf += sprintf(buf, "i:%u", *(int8_t*)(data + pos)); pos += 1; break;
			case 's': buf += sprintf(buf, "i:%u", *(int16_t*)(data + pos)); pos += 2; break;
			case 'S': buf += sprintf(buf, "i:%d", *(int16_t*)(data + pos)); pos += 2; break;
			case 'i': buf += sprintf(buf, "i:%u", *(int32_t*)(data + pos)); pos += 4; break;
			case 'I': buf += sprintf(buf, "i:%d", *(int32_t*)(data + pos)); pos += 4; break;
			case 'f': buf += sprintf(buf, "f:%g", *(float*)(data + pos)); pos += 4; break;
			case 'd': buf += sprintf(buf, "d:%lg", *(double*)(data + pos)); pos += 8; break;
			case 'Z': 
			case 'H': {
				buf += sprintf(buf, "%c:", t); 
				while (data[pos]) 
					buf += sprintf(buf, "%c", data[pos++]); 
				pos++; 
				break;
			}
			case 'B': {
				uint8_t type = *(data + pos); pos++;
				buf += sprintf(buf, "B:%c", type);
				int i = *(int32_t*)(data + pos); pos += 4;
				while (i-- > 0) {
					buf += sprintf(buf, ",");
					switch(type) {
						case 'c': buf += sprintf(buf, "i:%d", *(int8_t*)(data + pos)); pos += 1; break;
						case 'C': buf += sprintf(buf, "i:%u", *(int8_t*)(data + pos)); pos += 1; break;
						case 's': buf += sprintf(buf, "i:%u", *(int16_t*)(data + pos)); pos += 2; break;
						case 'S': buf += sprintf(buf, "i:%d", *(int16_t*)(data + pos)); pos += 2; break;
						case 'i': buf += sprintf(buf, "i:%u", *(int32_t*)(data + pos)); pos += 4; break;
						case 'I': buf += sprintf(buf, "i:%d", *(int32_t*)(data + pos)); pos += 4; break;
						case 'f': buf += sprintf(buf, "f:%g", *(float*)(data + pos)); pos += 4; break;
					}
				}
				break;
			}
			default: assert(0);
		}
	}
	*buf = 0;

	return true;
}

bool BAMParser::hasNext (void) {
	return !gzeof(input);
}

size_t BAMParser::fpos (void) {
	return lseek(fileno(fd), 0L, SEEK_CUR);
}

size_t BAMParser::fsize (void) {
	return file_size;
}

const Record &BAMParser::next (void) {
	return currentRecord;
}

string BAMParser::head (void) {
	return currentRecord.getChromosome();
}

