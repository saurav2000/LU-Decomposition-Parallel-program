pthread:
	g++ LUD_pthread.cpp -pthread
omp:
	g++ LUD_omp.cpp -fopenmp
seq:
	g++ LUD.cpp