CXX?=g++
CXXFLAGS= -w -c -O3 -std=c++0x
SOURCES=partition.cc svict.cc assembler.cc genome.cc variant_caller.cc svict_caller.cc simulator.cc annotation.cc logger.cc common.cc bam_parser.cc sam_parser.cc record.cc sort.cc extractor.cc smoother.cc
LDFLAGS=-lm -lz -lpthread -pthread
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=svict

all:	$(SOURCES) $(EXECUTABLE)
	rm -rf *.o

profile: PROFILE_FLAGS all

PROFILE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -pg -g)
	$(eval LDFLAGS = $(LDFLAGS) -pg -g)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	find . -name '*.o' -delete
