LIBS = parsa/libparsa.a
DIRS = parsa tests
.PHONY: all tests $(LIBS) clean

all: $(LIBS) tests

tests: $(LIBS)
	cd $@ && $(MAKE)

parsa/libparsa.a:
	cd parsa && $(MAKE)

clean:
	for dir in $(DIRS) ; do\
		(cd $$dir && $(MAKE) clean); \
	done
