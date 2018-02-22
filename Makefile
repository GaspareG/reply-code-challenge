all:
	g++ -Wall --pedantic -O2 -fopenmp -lpthread solution.cpp -o solution
	g++ -Wall --pedantic -O2 -fopenmp -lpthread solution_2d.cpp -o solution_2d

test: all
	time ./solution_2d < input_1.txt > output_1.txt
	time ./solution_2d < input_2.txt > output_2.txt
	time ./solution_2d < input_3.txt > output_3.txt
	time ./solution_2d < input_4.txt > output_4.txt
