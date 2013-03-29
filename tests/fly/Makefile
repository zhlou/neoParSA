CXXFLAGS += -I$(HOME)/include -I../../parsa
LDLIBS += -L$(HOME)/lib -lgsl -lgslcblas -L../../parsa -lparsa
FLYOBJECTS = DataLists.o flyData.o maternal.o scoring.o solvers.o \
			zygotic.o fly.o
TARGET = printscore fly_sa fly_sa-mpi scramble
LIBPARSA = ../../parsa/libparsa.a
.PHONY: all clean

all: $(TARGET)

printscore: printscore.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
scramble: scramble.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_sa: fly_sa.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_sa-mpi.o: fly_sa-mpi.cpp
	$(MPICXX) $(CXXFLAGS) -c $^
	
fly_sa-mpi: fly_sa-mpi.o $(FLYOBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(FLYOBJECTS) $(TARGET) printscore.o fly_sa.o fly_sa-mpi.o