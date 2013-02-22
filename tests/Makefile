DIRS = tsp rastrigin fly
.PHONY: all clean

all: $(DIRS)
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE)); \
	done

clean:
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE) clean); \
	done
