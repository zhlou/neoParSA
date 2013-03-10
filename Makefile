LIBS = parsa/libparsa.a
DIRS = parsa tests
export CXXFLAGS += `xml2-config --cflags` -pg
export LDLIBS += `xml2-config --libs` -pg
.PHONY: all tests clean libparsa rastrigin fly

all: libparsa tests

libparsa:
	cd parsa && $(MAKE)
	
rastrigin: libparsa
	cd tests && $(MAKE) rastrigin
	
fly: libparsa
	cd tests && $(MAKE) fly

tests: libparsa
	cd $@ && $(MAKE)

parsa/libparsa.a:
	cd parsa && $(MAKE)

clean:
	for dir in $(DIRS) ; do\
		(cd $$dir && $(MAKE) clean); \
	done
