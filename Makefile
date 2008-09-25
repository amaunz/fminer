# Edit according to your needs
PROGRAM = fminer

OBJ = database.o patterntree.o legoccurrence.o closeleg.o graphstate.o patterngraph.o path.o constraints.o

CC            = g++
CXXFLAGS      = -g -Wall -O3 -I./include/openbabel-2.0/
LIBS	      = -lm -ldl -lopenbabel -lgsl -lgslcblas
LDFLAGS       = -L./lib
RPATH         = -Wl,-rpath,./lib

.PHONY:
#all: graphstate.o
all: $(PROGRAM) 

$(PROGRAM): $(OBJ) $(PROGRAM).o
	$(CC) $(CXXFLAGS) $(LIBS) $(LDFLAGS) -o $(PROGRAM) $(OBJ) $(PROGRAM).o 

database.o: database.h
patterntree.o: patterntree.cpp patterntree.h patterngraph.h graphstate.h
legoccurrence.o: legoccurrence.h legoccurrence.cpp closeleg.h database.h graphstate.h
closeleg.o: closeleg.cpp closeleg.h misc.h
graphstate.o: graphstate.cpp graphstate.h database.h misc.h
patterngraph.o: patterngraph.cpp patterngraph.h graphstate.h
path.o: path.cpp path.h patterntree.h patterngraph.h graphstate.h
constraints.o: constraints.cpp constraints.h

.PHONY:
clean:
	-rm -rf *.o $(PROGRAM) 
