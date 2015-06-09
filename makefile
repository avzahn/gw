all: scaling scaling2 gaussian

scaling:
	gcc -lgsl -lgslcblas -lm omp_scaling.c gw.c -o omp_scaling -fopenmp 
	
scaling2:
	gcc -lgsl -lgslcblas -lm omp_scaling2.c gw.c -o omp_scaling2 -fopenmp 
	
gaussian:
	gcc -lgsl -lgslcblas -lm gaussian.c gw.c -o gaussian -fopenmp -ggdb

clean:
	rm -rf omp_scaling omp_scaling2 gaussian gaussian.out
