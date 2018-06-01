SRCDIR = src
OBJDIR = obj
BINDIR = bin

CXX = g++
CXXFLAGS = -g -Wall -pedantic -ansi --std=c++11 

EXENAME = main

OBJS = $(OBJDIR)/TCC_Matrix.o $(OBJDIR)/file_io.o $(OBJDIR)/kallisto_util.o $(OBJDIR)/util.o 

all: main

main: $(OBJDIR)/main.o $(OBJS) 
	$(CXX) -o $(BINDIR)/$(EXENAME) $^

debug: $(OBJDIR)/debug_util.o $(OBJDIR)/util.o
	$(CXX) -o $(BINDIR)/debug $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/*
