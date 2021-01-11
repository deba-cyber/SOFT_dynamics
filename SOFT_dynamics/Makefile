
CXX = g++
CXXFLAGS = -std=c++17 -O3 -malign-double   

#############################

INCDIR := include
SRCDIR := src
TESTDIR := test
OBJDIR := obj
BINDIR := bin
OUTPUTDIR := output_data
LDFLAGS := -lfftw3

#############################

OBJDIRINIT := $(OBJDIR)/objinit
OBJDIRDYN := $(OBJDIR)/objdyn
OBJDIRTEST := $(OBJDIR)/objtest

############################

# Header files for different executions

INITHEADER = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp $(INCDIR)/SOFT/tunneling_init_state.hpp  
DYNHEADER  = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp $(INCDIR)/SOFT/TROTTER_prop.hpp
TESTHEADER = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp

#############################

# Initial state preparation

init.exe: $(OBJDIRINIT)/constants.o $(OBJDIRINIT)/fileops.o $(OBJDIRINIT)/tunneling_init_state.o $(OBJDIRINIT)/get_init_state.o
	@echo "Generating Initial State ..."
	$(CXX) -o $@	$^

$(OBJDIRINIT)/constants.o: $(SRCDIR)/constants.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/fileops.o: $(SRCDIR)/fileops.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/tunneling_init_state.o: $(SRCDIR)/tunneling_init_state.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/get_init_state.o: $(SRCDIR)/get_init_state.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 


#############################

# Running SOFT dynamics

dyn.exe: $(OBJDIRDYN)/constants.o $(OBJDIRDYN)/fileops.o $(OBJDIRDYN)/TROTTER_prop.o $(OBJDIRDYN)/dynamics.o
	@echo "Running Split-Operator Fourier Transform dynamics ..."
	$(CXX) -o $@	$^	$(LDFLAGS)

$(OBJDIRDYN)/constants.o: $(SRCDIR)/constants.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/fileops.o: $(SRCDIR)/fileops.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/TROTTER_prop.o: $(SRCDIR)/TROTTER_prop.cpp $(DYNHEADER) 
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/dynamics.o: $(SRCDIR)/dynamics.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 


#############################

# Running analytical dynamics

testanal.exe: $(OBJDIRTEST)/constants.o $(OBJDIRTEST)/fileops.o $(OBJDIRTEST)/dynamics_analytical.o	
	@echo "Testing .."
	@echo "Running analytical time-evolution for initial state from linear combination of eigenstates"
	$(CXX)	-o $@	$^


$(OBJDIRTEST)/constants.o: $(TESTDIR)/constants.cpp $(TESTHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRTEST)/fileops.o: $(TESTDIR)/fileops.cpp $(TESTHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRTEST)/dynamics_analytical.o: $(TESTDIR)/dynamics_analytical.cpp 
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

.PHONY: clean


move:
	@echo "Moving execs after a run to bin directory ..."
	mv -f *.exe ./bin


clean:
	@echo "Cleaning all exec, objects, binaries in output directory ..."
	rm -f $(OBJDIRINIT)/* 
	rm -f $(OBJDIRDYN)/* 
	rm -f $(OBJDIRTEST)/* 
	rm -f $(BINDIR)/*
	rm -f $(OUTPUTDIR)/*.bin
	
