# WHAT
PROGRAM = fminer

# OPTIONS
VPATH         = ../libfminer
INCLUDE       = -I$(VPATH)/include/openbabel-2.0/
INCLUDE      +=  -I$(VPATH)

CC            = g++
CXXFLAGS      = -g -Wall -O3 $(INCLUDE)
LIBS	      = -lm -ldl -lopenbabel -lgsl -lgslcblas
LDFLAGS       = -L$(VPATH)/lib
RPATH         = -Wl,-rpath,$(VPATH)

# TARGETS
.PHONY:
all: $(PROGRAM) 
$(PROGRAM): -lfminer $(PROGRAM).o
	$(CC) $(CXXFLAGS) $(LIBS) $(LDFLAGS) $(RPATH) -o $@ $(OBJ) $^ 

.PHONY:
clean:
	-rm -rf *.o $(PROGRAM) 
