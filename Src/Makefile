# Makefile
#

#objects := $(patsubst %.cpp, %.o, $(wildcard *.cpp))
objects := calign.o smallapp.o hashapp.o alignapp.o
objects2 := cleanseq.o

CXX := g++

CXXFLAGs :=

BUILD_TYPE :=

ifeq ($(BUILD_TYPE), debug)
	CXXFLAGS += -g
else
	CXXFLAGS += -O2 -w
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
else
	CXXFLAGS +=
endif

INC := -I ./

clus: AbundantOTU.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o AbundantOTU+ AbundantOTU.cpp $(objects) $(INC) 
$(objects) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)

CleanMap: cleanmap.cpp $(objects2)
	$(CXX) $(CXXFLAGS) -o CleanMap cleanmap.cpp $(objects2) $(INC) 
$(objects2) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)

clean:
	$(RM) *.o AbundantOTU+ 
