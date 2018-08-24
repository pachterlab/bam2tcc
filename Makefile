SEQAN_PATH = /opt/local/bin/seqan

SRCDIR = src
OBJDIR = obj
BINDIR = bin

CXX = g++
CXXFLAGS = -g -I$(SEQAN_PATH)/include -Wall -ansi --std=c++14 

EXENAME = thing

OBJS = $(OBJDIR)/TCC_Matrix.o $(OBJDIR)/Read.o $(OBJDIR)/Transcript.o \
$(OBJDIR)/common.o $(OBJDIR)/FileUtil.o $(OBJDIR)/Mapper.o

all: thing

thing: $(OBJDIR)/main.o $(OBJS) 
	$(CXX) -o $(BINDIR)/$(EXENAME) $^ -pthread

debug: $(OBJDIR)/debugUtil.o $(OBJS)
	$(CXX) -o $(BINDIR)/debug $^ -pthread

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/*
