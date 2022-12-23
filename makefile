CC = g++
OPT = -O3
FLAGS = -Wall
BIN_PATH = bin/
SRC_PATH = src/

BINARY = reflow.out

FILES_N = main.cpp mesh.cpp variables.cpp solver.cpp thermodynamics.cpp reflow.cpp
OBJECTS_N = main.o mesh.o variables.o solver.o thermodynamics.o reflow.o

OBJECTS = $(foreach F,$(OBJECTS_N),$(BIN_PATH)$(F))
FILES = $(foreach F,$(FILES),$(SRC_PATH)$(F))

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CC) -o $@ $^ 

bin/%.o: src/%.cpp
	$(CC) $(FLAGS) -c $^ -o $@

clean: 
	rm $(BINARY) $(OBJECTS)