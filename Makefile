all : integral

integral : Integral.cpp
	clang++ $^ -fopenmp -O3 -o bin/$@