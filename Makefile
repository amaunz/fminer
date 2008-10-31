# WHAT
PROGRAM = fminer

# OPTIONS
LIBDIR        = ../libfminer
INCLUDE       = -I$(LIBDIR)
INCLUDE       += -I/usr/local/include/openbabel-2.0/

LDFLAGS       = -L$(LIBDIR)

CC            = g++
CXXFLAGS      = -g -Wall -O3
LIBS	      = -lopenbabel -lfminer -lgsl -lgslcblas
RPATH         = -Wl,-rpath,$(LIBDIR)

# TARGETS
all: lib $(PROGRAM) 
lib: 
	$(MAKE) -C $(LIBDIR)
$(PROGRAM): main.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) \
	      $(LIBS) \
	      $(LDFLAGS) \
	      $(RPATH) \
	      -o $@ $<
doc: README
	rdoc

.PHONY:
clean:
	-rm -rf $(PROGRAM) 
