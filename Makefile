PROGRAM = fminer
FMINER_DIR = ../libfminer
INCLUDES = -I$(FMINER_DIR)/include/openbabel-2.0/
INCLUDES +=  -I$(FMINER_DIR)
LIBDIR = ../libfminer/lib

OBJ = $(FMINER_DIR)/database.o $(FMINER_DIR)/patterntree.o $(FMINER_DIR)/legoccurrence.o $(FMINER_DIR)/closeleg.o $(FMINER_DIR)/graphstate.o $(FMINER_DIR)/patterngraph.o $(FMINER_DIR)/path.o $(FMINER_DIR)/constraints.o

CC            = g++
CXXFLAGS      = -g -Wall -O3 $(INCLUDES)
LIBS	      = -lm -ldl -lopenbabel -lgsl -lgslcblas
LDFLAGS       = -L$(LIBDIR)
RPATH         = -Wl,-rpath,$(LIBDIR)

.PHONY:
all: $(PROGRAM) 

$(PROGRAM): $(OBJ) $(PROGRAM).o
	$(CC) $(CXXFLAGS) $(LIBS) $(LDFLAGS) -o $(PROGRAM) $(OBJ) $(PROGRAM).o 

$(FMINER_DIR)/database.o: $(FMINER_DIR)/database.h
$(FMINER_DIR)/patterntree.o: $(FMINER_DIR)/patterntree.cpp $(FMINER_DIR)/patterntree.h $(FMINER_DIR)/patterngraph.h $(FMINER_DIR)/graphstate.h
$(FMINER_DIR)/legoccurrence.o: $(FMINER_DIR)/legoccurrence.h $(FMINER_DIR)/legoccurrence.cpp $(FMINER_DIR)/closeleg.h $(FMINER_DIR)/database.h $(FMINER_DIR)/graphstate.h
$(FMINER_DIR)/closeleg.o: $(FMINER_DIR)/closeleg.cpp $(FMINER_DIR)/closeleg.h $(FMINER_DIR)/misc.h
$(FMINER_DIR)/graphstate.o: $(FMINER_DIR)/graphstate.cpp $(FMINER_DIR)/graphstate.h $(FMINER_DIR)/database.h $(FMINER_DIR)/misc.h
$(FMINER_DIR)/patterngraph.o: $(FMINER_DIR)/patterngraph.cpp $(FMINER_DIR)/patterngraph.h $(FMINER_DIR)/graphstate.h
$(FMINER_DIR)/path.o: $(FMINER_DIR)/path.cpp $(FMINER_DIR)/path.h $(FMINER_DIR)/patterntree.h $(FMINER_DIR)/patterngraph.h $(FMINER_DIR)/graphstate.h
$(FMINER_DIR)/constraints.o: $(FMINER_DIR)/constraints.cpp $(FMINER_DIR)/constraints.h

.PHONY:
clean:
	-rm -rf *.o $(PROGRAM) 
