
EXECUTABLE := tsp 

CU_FILES   := tsp.cu 

CC_FILES   :=  

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')
OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O0 -Wall -g

LIBS       :=
FRAMEWORKS := 

ifeq ($(ARCH), Darwin)
# Building on mac
NVCCFLAGS=-O0 -m64 -arch compute_10
FRAMEWORKS += OpenGL GLUT
LDFLAGS=-L/usr/local/cuda/lib/ -lcudart
else
# Building on Linux
NVCCFLAGS=-O0 -m64 -g -G -arch compute_20
LIBS += GL glut cudart
LDFLAGS=-L/usr/local/cuda/lib64/ -lcudart
endif

LDLIBS  := $(addprefix -l, $(LIBS))
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
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

$(OBJDIR)/%.o: %.cpp
		$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/%.o: %.cu
		$(NVCC) $< $(NVCCFLAGS) -c -o $@
