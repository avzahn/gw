#include "gw.h"

double lnprob(double * walker, int ndim) {
	return 1.0;
}

int main(){

	int i;
	int j;
	int t;
	gw_ensemble * s;
	FILE * f = fopen("omp_scaling.out","w");
	double mu[] = {1,1,1,1};


	s = gw_alloc(4,1000,4);
	s->lnprob = lnprob;
	gw_gaussian(s,mu,1);

	for(j=1; j < 4; j++){

		s->nthreads = j;

		t = time(NULL);

		for(i=0;i<10000;i++){
			gw_metropolis(s);
		}

		printf("%f\n", time(NULL)-t);

	}



	fclose(f);
	gw_free(s);
}