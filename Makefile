CXX ?= g++
VPATH=src/
CXXFLAGS += -Wall -Wextra -std=c++17 -pedantic -O3 -flto -fopenmp -march=native # for perf -g3 -fno-omit-frame-pointer
CXXLIBS += -lSDL2 -lSDL2_image

SRC = main.cc vector.cc kdtree.cc triangle.cc material.cc parse.cc light.cc \
	compute_light.cc matrix.cc texture.cc sphere_light.cc
OBJS = ${SRC:.cc=.o}
BIN = main

all: $(BIN)

main: ${OBJS}
	$(CXX) $(CXXFLAGS) $(OBJS) $(CXXLIBS) -o $(BIN)

check: CXXFLAGS += -g3 -O0 -fno-inline -fsanitize=address
check: $(BIN)

.PHONY: clean check
clean:
	${RM} ${OBJS}
	${RM} $(BIN)

test: $(BIN)
	./$(BIN) ./examples/scenes/test.json && feh out.ppm
	./$(BIN) ./examples/scenes/iron.json && feh out.ppm
	./$(BIN) ./examples/scenes/car.json && feh out.ppm
	./$(BIN) ./examples/scenes/eye.json && feh out.ppm

tar:
	tar -cvjf examples.tar.bz2 examples Textures
untar:
	tar -xvf examples.tar.bz2
