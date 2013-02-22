CXXFLAGS += `xml2-config --cflags`
LDLIBS += `xml2-config --libs`
FLYOBJECTS = DataLists.o flyData.o maternal.o scoring.o solvers.o \
			zygotic.o fly.o
TARGET = printscore
.PHONY: all clean

all: $(TARGET)

printscore: printscore.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ printscore.o $(FLYOBJECTS) $(LDLIBS)