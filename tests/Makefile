DIRS = tsp rastrigin fly
.PHONY: all clean $(DIRS)

all: $(DIRS)
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE)); \
	done
	
rastrigin:
	cd $@ && $(MAKE)
	
fly:
	cd $@ && $(MAKE)

tsp:
	cd $@ && $(MAKE)

clean:
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE) clean); \
	done
