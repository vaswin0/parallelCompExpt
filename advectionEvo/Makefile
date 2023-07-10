GSLLIBS      := $(shell gsl-config --libs)
OMPFLAGS	 := -fopenmp
OMPLIBS		 := -lgomp

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native $(OMPFLAGS)
LD            = g++
LDFLAGS       = -O3 -march=native $(OMPFLAGS)

LIBS          = $(SYSLIBS) $(GSLLIBS) $(OMPLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp cell.cpp grid.cpp init.cpp eos.cpp evolve.cpp
             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = Evo.exe
#--------------------------------------------------

$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@




