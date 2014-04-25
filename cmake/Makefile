LIBS = parsa/libparsa.a
DIRS = parsa tests
export CXXFLAGS += `xml2-config --cflags`
export LDLIBS += `xml2-config --libs`
export MPICXX = mpicxx
export OMP_CXX = $(CXX)
export MPICH_CXX = $(CXX)
.PHONY: all tests clean libparsa rastrigin fly debug

all: libparsa tests

debug: CXXFLAGS += -g
debug: all

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
