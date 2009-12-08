# Makefile
# © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

# This file is part of Fminer (fminer).
#
# Fminer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Fminer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Fminer.  If not, see <http://www.gnu.org/licenses/>.

# PATH TO LIBFMINER SOURCE DIRECTORY:
LIBDIR        = ../libfminer

INCLUDE_OB    =  
INCLUDE_GSL   =   
LDFLAGS_OB    = 
LDFLAGS_GSL   =   
INCLUDE_RB    =

# EXAMPLES:
# ADJUST COMPILER PATH TO OPENBABEL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_OB    = -I/usr/local/include/openbabel-2.0
INCLUDE_OB   += -I/home/openbabel-2.2.1/include
# ADJUST COMPILER PATH TO GSL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_GSL   = -I/usr/include
INCLUDE_GSL  += -I/c/Program\ Files/GnuWin32/include/

# ADJUST LINKER PATH TO OPENBABEL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_OB    = -L/usr/local/lib
LDFLAGS_OB   += -L/c/Users/am/fminer
# ADJUST LINKER PATH TO GSL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_GSL   =
LDFLAGS_GSL  += -L/c/Program\ Files/GnuWin32/bin/



# FOR LINUX: INSTALL TARGET DIRECTORY
DESTDIR       = /usr/local/bin


# NORMALLY NO ADJUSTMENT NECESSARY BELOW THIS LINE. Exit and try 'make' now.
# WHAT
PROGRAM = fminer
# OPTIONS
CC            = g++
INCLUDE       = -I$(LIBDIR) $(INCLUDE_OB) $(INCLUDE_GSL)
LDFLAGS       = -L$(LIBDIR) $(LDFLAGS_OB) $(LDFLAGS_GSL)
CXXFLAGS      = -Wall -O3 -g
ifeq ($(OS), Windows_NT) # assume MinGW/Windows
LIBS	      = -lm -llibopenbabel-3 -llibgsl -llibgslcblas -llibfminer
else
LIBS	      = -lopenbabel -lfminer -lgsl -lgslcblas
endif

# TARGETS
all: $(PROGRAM) 
$(PROGRAM): main.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) \
	      $(LIBS) \
	      $(LDFLAGS) \
	      -o $@ $<

install: $(PROGRAM)
	cp -P fminer $(DESTDIR)

.PHONY:
clean:
	-rm -rf $(PROGRAM) $(PROGRAM).exe
