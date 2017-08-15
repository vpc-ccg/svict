EXECUTABLE=CellFreeSV

all:
	$(MAKE) -C src all
	mv src/CellFreeSV .

profile:
	$(MAKE) -C src profile
	mv src/CellFreeSV .

clean:
	$(MAKE) -C src clean
	rm -rf ${EXECUTABLE}
