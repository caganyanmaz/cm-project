
main: main.c matrix.c matrix.h
	gcc main.c matrix.c -o main

matrix_tests: matrix_tests.c matrix.c matrix.h
	gcc matrix_tests.c matrix.c -o matrix_tests && ./matrix_tests


clear:
	rm main matrix_tests

all: main



