# WHAT
PROGRAM = fminer

# OPTIONS
LIBDIR        = libfminer
INCLUDE       = -I$(LIBDIR)
INCLUDE       += -I/usr/include/openbabel-2.0/

LDFLAGS       = -L$(LIBDIR)

CC            = g++
CXXFLAGS      = -g -Wall -O3
ifeq ($(OS), Windows_NT) # assume MinGW/Windows
LIBS	      = -lm -llibopenbabel-3 -llibgsl -llibgslcblas -llibfminer
else
LIBS	      = -lopenbabel -lfminer -lgsl -lgslcblas
endif
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

.PHONY:
clean:
	-rm -rf $(PROGRAM) $(PROGRAM).exe
