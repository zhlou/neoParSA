LIBS = parsa/libparsa.a
.PHONY: tests $(LIBS) 

tests: $(LIBS)
	cd $@ && $(MAKE)

parsa/libparsa.a:
	cd parsa && $(MAKE)
