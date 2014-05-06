all: tsp

CXX=mpic++
CXXFLAGS=-O3 -std=c++0x
LDFLAGS=-lpthread -lmpi -lmpi_cxx

tsp: tsp.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

clean:
	rm tsp 

run280:
	mpirun -np 4 ./tsp ../benchmark/a280.tsp

run100:
	mpirun -np 4 ./tsp ../benchmark/kroA100.tsp

run493:
	mpirun -np 4 ./tsp ../benchmark/d493.tsp
