CC = g++
OPT = -O3
FLAGS = -std=c++11 -ffast-math
LIBS = -fopenmp
BIN_PATH = bin/
SRC_PATH = src/

BINARY = reflow.out

FILES_N = main.cpp mesh.cpp variables.cpp solver.cpp thermodynamics.cpp initial_cond.cpp boundary_cond.cpp reflow.cpp particle.cpp lagrange_solver.cpp chem_solver.cpp geometry.cpp
OBJECTS_N = main.o mesh.o variables.o solver.o thermodynamics.o initial_cond.o boundary_cond.o reflow.o particle.o lagrange_solver.o chem_solver.o geometry.o

OBJECTS = $(foreach F,$(OBJECTS_N),$(BIN_PATH)$(F))
FILES = $(foreach F,$(FILES),$(SRC_PATH)$(F))

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CC) -o  $@ $^ $(LIBS)

bin/%.o: src/%.cpp
	$(CC) $(FLAGS) $(OPT) -c $^ -o $@ $(LIBS)

clean: 
	rm $(BINARY) $(OBJECTS)