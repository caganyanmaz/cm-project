
main: main.c matrix.c matrix.h
	gcc main.c matrix.c -lm -o main

matrix_tests: matrix_tests.c matrix.c matrix.h
	gcc matrix_tests.c matrix.c -lm -o matrix_tests && ./matrix_tests


clear:
	rm main matrix_tests

all: main



