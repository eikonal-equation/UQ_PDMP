## Compiler flags
CXXFLAGS = -std=c++14 -O2

## Directories
BOOSTPATH = /usr/local/include/boost
INCLUDE = -I $(BOOSTPATH) -I ./include/ 
LIB = -L $(BOOSTPATH)
SRC = ./src/
CPPFILES := $(wildcard $(SRC)*.cpp)
HPPFILES := $(wildcard ./include/*.hpp)

## Command line Arguments
TEST =
NX =
NS =

## Commands
main: $(CPPFILES) $(HPPFILES)
	 $(CXX) $(CXXFLAGS) $(INCLUDE) $(CPPFILES) $(LIB) -o $@

run: main
	./main $(TEST) $(NX) $(NS)

clean:
	$(RM) main

realclean:
	rm output/*
	clean

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@