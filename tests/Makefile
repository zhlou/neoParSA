DIRS = tsp
.PHONY: $(DIRS) clean

$(DIRS):
	cd $@ && $(MAKE)

clean:
	for dir in $(DIRS); do\
		(cd $$dir && $(MAKE) clean); \
	done
