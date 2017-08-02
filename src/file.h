#ifndef __MFILE__
#define __MFILE__

class file {
	public:
		file() {}
		virtual ~file(){}
		virtual void open(const char* fn, const char* m) = 0;
		virtual void close() = 0;
		virtual ssize_t read(void* d, size_t s) = 0;
		virtual ssize_t write(void* d, size_t s) = 0;
		virtual bool eof() = 0;
};

class rawfile: public file {
	FILE *f;
	public:
		rawfile(const char* fn, const char* m) { open(fn, m); }
		~rawfile(){}
		virtual void open(const char* fn, const char* m) { f = fopen(fn, m); }
		virtual void close() { fclose(f); }
		virtual ssize_t read(void* d, size_t s) { return fread(d, 1, s, f); }
		virtual ssize_t write(void* d, size_t s) { return fwrite(d, 1, s, f); }
		virtual bool eof() { return feof(f); }
		FILE *fh () { return f; }
};

class gzfile: public file {
	gzFile f;
	public:
		gzfile(const char* fn, const char* m) { open(fn, m); }
		~gzfile(){}
		virtual void open(const char* fn, const char* m) { f = gzopen(fn, m); /*gzbuffer(f, 128 * 1024);*/ }
		virtual void close() { gzclose(f); }
		virtual ssize_t read(void* d, size_t s) { 
			const size_t offset = 1 * (size_t)GB;
			if (s > offset) {
				return gzread(f, d, offset) + this->read((char*)d + offset, s - offset); 
			}
			else {
				return gzread(f, d, s);
			}
		}
		virtual ssize_t write(void* d, size_t s) { return gzwrite(f, d, s); }
		virtual bool eof() { return gzeof(f); }
};



#endif
