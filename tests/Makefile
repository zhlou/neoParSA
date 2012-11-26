DIRS = tsp rastrigin
.PHONY: all clean

all: $(DIRS)
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE)); \
	done

clean:
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE) clean); \
	done
