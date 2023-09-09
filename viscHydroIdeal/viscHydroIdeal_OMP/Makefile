
#ROOTCFLAGS   := $(shell root-config --cflags)
#ROOTLIBS     := $(shell root-config --libs)

GSLLIBS      := $(shell gsl-config --libs)
OMPFLAGS := -fopenmp
OMPLIBS := -lgomp

#EXTRA_FLAGS = -D SIMPLE # EoS p=e/3
EXTRA_FLAGS  = -D TABLE # Laine EoS, tabulated

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native $(OMPFLAGS)
LD            = g++
LDFLAGS       = -O3 -march=native $(OMPFLAGS)

CXXFLAGS     += $(EXTRA_FLAGS)
LIBS          = $(SYSLIBS) $(GSLLIBS) $(OMPLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp master.cpp grid.cpp cell.cpp cnvrt.cpp init.cpp hydro.cpp eos0.cpp 
             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = vischydro
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

cfiles : 
	rm -rf hydro_output/*
	rm -rf output_after/* 

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
