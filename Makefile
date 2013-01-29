LIBS = parsa/libparsa.a
DIRS = parsa tests
export CXXFLAGS += `xml2-config --cflags` -pg
export LDLIBS += `xml2-config --libs` -pg
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
