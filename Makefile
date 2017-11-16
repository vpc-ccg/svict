EXECUTABLE=SVICT

all:
	$(MAKE) -C src all
	mv src/${EXECUTABLE} .

profile:
	$(MAKE) -C src profile
	mv src/${EXECUTABLE} .

clean:
	$(MAKE) -C src clean
	rm -rf ${EXECUTABLE}
