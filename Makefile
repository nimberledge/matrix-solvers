main: main.o Matrix.o CSRMatrix.o Interface.o
	g++ -std=c++17 bin/main.o bin/Matrix.o bin/CSRMatrix.o bin/Interface.o -o bin/main

interactive: Matrix.o CSRMatrix.o Interface.o interactive.o
	g++ -std=c++17 bin/Interface.o bin/Matrix.o bin/CSRMatrix.o bin/Interactive.o -o bin/ladsSolve

interactive.o: src/ladsSolverInteractive.cpp
	g++ -std=c++17 -I src/ -c src/ladsSolverInteractive.cpp -o bin/Interactive.o

benchmark: Matrix.o CSRMatrix.o benchmark.o
	g++ -std=c++17 bin/Matrix.o bin/CSRMatrix.o bin/benchmark.o -o bin/benchmark

benchmark.o: test/benchmark.cpp
	g++ -std=c++17 -I src/ -c test/benchmark.cpp -o bin/benchmark.o

main.o: src/main.cpp
	g++ -std=c++17 -c src/main.cpp -o bin/main.o

test: test1 test2 test3
	echo "tests successfully built"

test1: test1.o Matrix.o
	g++ -std=c++17 bin/test1.o bin/Matrix.o -o bin/test1

test2: test2.o Matrix.o CSRMatrix.o
	g++ -std=c++17 bin/test2.o bin/Matrix.o bin/CSRMatrix.o -o bin/test2

test3: test3.o Matrix.o CSRMatrix.o Interface.o
	g++ -std=c++17 bin/test3.o bin/Matrix.o bin/CSRMatrix.o bin/Interface.o -o bin/test3

test2.o: src/Matrix.h src/CSRMatrix.h test/test2.cpp
	g++ -std=c++17 -I src/ -c test/test2.cpp -o bin/test2.o

test3.o : src/Matrix.h src/CSRMatrix.h src/Interface.h test/test3.cpp
	g++ -std=c++17 -I src/ -c test/test3.cpp -o bin/test3.o

CSRMatrix.o: Matrix.o src/CSRMatrix.h src/CSRMatrix.cpp
	g++ -std=c++17 -I src/ -c src/CSRMatrix.cpp -o bin/CSRMatrix.o

Matrix.o: src/Matrix.cpp src/Matrix.h
	g++ -std=c++17 -c src/Matrix.cpp -o bin/Matrix.o

Interface.o: src/Interface.cpp src/Interface.h
	g++ -std=c++17 -c src/Interface.cpp -o bin/Interface.o

test1.o: src/Matrix.h test/test1.cpp
	g++ -std=c++17 -I src/ -c test/test1.cpp -o bin/test1.o

clean:
	rm bin/*