LIBS = parsa/libparsa.a
.PHONY: all tests $(LIBS)

all: $(LIBS) tests

tests: $(LIBS)
	cd $@ && $(MAKE)

parsa/libparsa.a:
	cd parsa && $(MAKE)
