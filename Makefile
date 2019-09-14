CXX = clang++
CXX_FLAGS = -std=c++11 -fopenmp -O3
CXX_LINK_FLAGS = -lm -lfftw3

TARGET = bin/sample
MAIN = sample.cpp

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=obj/%.o)
INC = $(wildcard include/*.hpp)
INCLUDE_DIR = include

obj/%.o : src/%.cpp include/%.hpp
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE_DIR) -c -o $@ $<

$(TARGET) : $(MAIN) $(OBJ)
	$(CXX) $(CXX_FLAGS) $(CXX_LINK_FLAGS) -I$(INCLUDE_DIR) -o $@ $^

all : $(TARGET)

clean :
	rm -rf obj/*.o bin/*