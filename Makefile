SRCDIR = src
OBJDIR = obj
BINDIR = bin

CXX = g++
CXXFLAGS = -g -I/opt/local/bin/seqan/include -Wall -pedantic -ansi --std=c++14 

EXENAME = main

OBJS = $(OBJDIR)/TCC_Matrix.o $(OBJDIR)/file_io.o $(OBJDIR)/kallisto_util.o $(OBJDIR)/util.o 

all: main

main: $(OBJDIR)/main.o $(OBJS) 
	$(CXX) -o $(BINDIR)/$(EXENAME) $^ -pthread

debug: $(OBJDIR)/debug_util.o $(OBJDIR)/util.o
	$(CXX) -o $(BINDIR)/debug $^

# File not on github. Used to decide whether multithreading would be worth it.
# (It is.)
#timing: $(OBJDIR)/timing_tests.o $(OBJS)
#	$(CXX) -o $(BINDIR)/timing $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/*
