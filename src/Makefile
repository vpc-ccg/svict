CC=g++
CFLAGS= -w -c -O3 -std=gnu++0x  -D__USE_XOPEN2K8
SOURCES=partition.cc cell_free_sv.cc assembler.cc assembler_old.cc genome.cc variant_caller.cc kmistrvar.cc simulator.cc annotation.cc logger.cc assembler_ext.cc common.cc bam_parser.cc sam_parser.cc record.cc sort.cc extractor.cc
LDFLAGS=-lm -lz -lpthread -pthread
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=SVICT

all:	$(SOURCES) $(EXECUTABLE)
	rm -rf *.o

profile: PROFILE_FLAGS all

PROFILE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LDFLAGS = $(LDFLAGS) -pg -g)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	find . -name '*.o' -delete
