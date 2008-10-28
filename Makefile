# WHAT
PROGRAM = fminer

# OPTIONS
LIBDIR        = ../libfminer
INCLUDE       = -I$(LIBDIR)/include/openbabel-2.0/
INCLUDE      += -I$(LIBDIR)
LDFLAGS       = -L$(LIBDIR)/lib
LDFLAGS      += -L$(LIBDIR)

CC            = g++
CXXFLAGS      = -g -Wall -O3
LIBS	      = -lopenbabel -lfminer -lgsl -lgslcblas
RPATH         = -Wl,-rpath,$(LIBDIR):$(LIBDIR)/lib

# TARGETS
all: lib $(PROGRAM) 
lib: 
	$(MAKE) -C $(LIBDIR)
$(PROGRAM): main2.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) \
	      $(LIBS) \
	      $(LDFLAGS) \
	      $(RPATH) \
	      -o $@ $<

.PHONY:
clean:
	-rm -rf *.o $(PROGRAM) 
