#CC = $(CC_ENV)
CC = g++
#CC = /usr/bin/g++-4.4
CCFLAGS = -Wall -g
DELPHES = -I.

# Delphes directory
DELPHES_DIR = ../
NSUB_DIR = ../tmp/external/fastjet/contribs/Nsubjettiness/

# root
ROOTCFLAGS = $(shell root-config --cflags)
ROOTGLIBS = $(shell root-config --glibs)


# File names
SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
EXEC = $(SOURCES:.cc=.exe)

# all
all: $(EXEC)

# To compile
%.exe: %.o ExRootTreeReader.o # -L$(NSUB_DIR) %.o Nsubjettiness.o $< -o $@
	$(CC) $(CCFLAGS) $(ROOTGLIBS)  -L$(DELPHES_DIR) -lDelphes ExRootTreeReader.o $< -o $@ # -L$(NSUB_DIR) -lNsubjettiness $< -o $@

# To obtain object files
%.o: %.cc
	$(CC) -c $(CCFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS) $(DELPHES) $< -o $@

# needed everywhere: ExRootTreeReader
ExRootTreeReader.o: ExRootAnalysis/ExRootTreeReader.cc ExRootAnalysis/ExRootTreeReader.h
	$(CC) $(CCFLAGS) $(ROOTCFLAGS) $(ROOTGLIBS) $(DELPHES) -c ExRootAnalysis/ExRootTreeReader.cc

# needed everywhere: ExRootTreeReader
# Nsubjettiness.o: fastjet/contribs/Nsubjettiness/Nsubjettiness.cc fastjet/contribs/Nsubjettiness/Nsubjettiness.h
#	$(CC) $(CCFLAGS) $(ROOTCFLAGS) $(ROOTGLIBS) $(DELPHES) -c fastjet/contribs/Nsubjettiness/Nsubjettiness.cc

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

