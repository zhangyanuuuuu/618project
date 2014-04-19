
EXECUTABLE := tsp 

CU_FILES   := tsp.cu 

CC_FILES   :=  

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')
OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O3 -Wall

LIBS       :=
FRAMEWORKS := 

# Building on Linux
NVCCFLAGS=-O3 -m64 -arch compute_20 
LIBS += GL glut cudart
LDFLAGS=-L/usr/local/cuda/lib64/ -lcudart -lcurand

LDFRAMEWORKS := $(addprefix -framework , $(FRAMEWORKS))

NVCC=nvcc

OBJS=$(OBJDIR)/tsp.o 

.PHONY: dirs clean

default: $(EXECUTABLE)

dirs:
		mkdir -p $(OBJDIR)/

clean:
		rm -rf $(OBJDIR) *~ $(EXECUTABLE) 

$(EXECUTABLE): dirs $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDFRAMEWORKS)

$(OBJDIR)/%.o: %.cu
		$(NVCC) $< $(NVCCFLAGS) -c -o $@
