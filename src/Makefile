CXX ?= g++
CXXFLAGS += -Wall -Wextra -std=c++17 -pedantic -O3 -fopenmp

SRC = main.cc vector.cc kdtree.cc triangle.cc parse.cc
OBJS = ${SRC:.cc=.o}
BIN = main

all: $(BIN)

main: ${OBJS}
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN)

check: CXXFLAGS = -g3 -fopenmp
check: $(BIN)

.PHONY: clean check
clean:
	${RM} ${OBJS}
	${RM} $(BIN)
