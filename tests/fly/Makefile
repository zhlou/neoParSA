CXXFLAGS += -I$(HOME)/include -I../../parsa -g
LDLIBS += -L$(HOME)/lib -lgsl -lgslcblas -L../../parsa -lparsa
FLYOBJECTS = DataLists.o flyData.o maternal.o scoring.o solvers.o \
			zygotic.o fly.o
TARGET = printscore fly_sa
LIBPARSA = ../../parsa/libparsa.a
.PHONY: all clean

all: $(TARGET)

printscore: printscore.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_sa: fly_sa.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(FLYOBJECTS) $(TARGET) printscore.o fly_sa.o