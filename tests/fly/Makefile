CXXFLAGS += `xml2-config --cflags` -I$(HOME)/include -g
LDLIBS += `xml2-config --libs` -L$(HOME)/lib -lgsl -lgslcblas
FLYOBJECTS = DataLists.o flyData.o maternal.o scoring.o solvers.o \
			zygotic.o fly.o
TARGET = printscore
.PHONY: all clean

all: $(TARGET)

printscore: printscore.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ printscore.o $(FLYOBJECTS) $(LDLIBS)
	
clean:
	rm -f $(FLYOBJECTS) $(TARGET) printscore.o