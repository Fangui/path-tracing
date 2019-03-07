CXX ?= g++
VPATH=src/
CXXFLAGS += -Wall -Wextra -std=c++17 -pedantic -O3 -fopenmp

SRC = main.cc vector.cc kdtree.cc triangle.cc material.cc parse.cc light.cc \
	compute_light.cc
OBJS = ${SRC:.cc=.o}
BIN = main

all: $(BIN)

main: ${OBJS}
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(BIN)

check: CXXFLAGS = -g3 -O0 -fno-inline -fopenmp
check: $(BIN)

.PHONY: clean check
clean:
	${RM} ${OBJS}
	${RM} $(BIN)

test:
	./$(BIN) ./examples/scenes/iron.json && feh out.ppm
	./$(BIN) ./examples/scenes/car.json && feh out.ppm
	./$(BIN) ./examples/scenes/eye.json && feh out.ppm