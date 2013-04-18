#    This file is part of MLML
#
#    Copyright (C) 2013 University of Southern California and
#             		Meng Zhou, Jianghan Qu, Andrew D. Smith
#
#    Authors: Meng Zhou, Jianghan Qu and Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

PROGS = mlml
SRC_ROOT = $(MLML_ROOT)

CXX = g++
CFLAGS = -Wall -fmessage-length=50
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
SMITHLAB_CPP = smithlab_cpp

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

INCLUDEDIRS = $(SMITHLAB_CPP)

LIBS = -lgsl -lgslcblas 

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ -I $(INCLUDEDIRS) $(LIBS)

clean:
	@-make -C smithlab_cpp clean
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
