#    C++ makefile
#    Purpose: src/cpp makefile
#
#    @author J.M.G. Salmerón
#    @version 1.5 17/11/2018
#        
#    Created on: 20/09/2018

CC := g++
CVERSION := -std=c++17
#CFLAGS := -O3
CFLAGS := -g -Wall

SRCDIR := src/cpp
TESTDIR := src/cpp/test
BUILDDIR := build
TARGET := run

MKLROOT := ~/intel/mkl

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

LIB_MKL := -L ${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt# Single Dynamic Library
LIB_PY := -lpython3.6m
LIB_GEN := -lpthread -lm -ldl
LIB := $(LIB_MKL) $(LIB_MATLAB) $(LIB_PY) $(LIB_GEN)

INC_MKL := -m64 -I ${MKLROOT}/include# Single Dynamic Library
INC_PY := -I /usr/include/python3.6m
INC := $(INC_MKL) $(INC_PY)

$(TARGET): $(OBJECTS)
	@echo " $(OBJECTS) "
	@echo " Linking..."
	@echo " $(CC) $(CVERSION) $(CFLAGS) $^ -o $(TARGET) $(LIB)"; $(CC) $(CVERSION) $(CFLAGS) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CVERSION) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CVERSION) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

total:
	make clean
	make
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/intel/mkl/lib/intel64

	#./run -n 3 -D facet
	#./run -n 5 -D facet
	./run -n 14 -M dimacs/Brock14.txt -t 5 -D facet
	#./run -n 16 -M dimacs/1tc16_Clique.txt -t 8 -D facet
	#./run -n 28 -M dimacs/Johnson8-2-4.txt -t 4 -D facet
	#./run -n 32 -M dimacs/1tc32_Clique.txt -t 8 -D facet
	
	#./run -n 3 -D zbund
	#./run -n 5 -D zbund
	./run -n 14 -M dimacs/Brock14.txt -t 5 -D zbund
	#./run -n 16 -M dimacs/1tc16_Clique.txt -t 8 -D zbund
	#./run -n 28 -M dimacs/Johnson8-2-4.txt -t 4 -D zbund
	#./run -n 32 -M dimacs/1tc32_Clique.txt -t 8 -D zbund


.PHONY: clean
