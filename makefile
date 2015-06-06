all:
	gcc  -L/usr/include -lgsl -lm omp_scaling.c gw.c -o omp_scaling -fopenmp  