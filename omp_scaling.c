#include "gw.h"

double lnprob(double * walker, int ndim) {
	int i;
	double ret = 3;
	for(i = 1;i <= 1; i++){
		ret *= i;
		ret /= ret + 1;
	}
	
	
	return -ret * walker[0]*walker[0] ;
}

int main(){

	int i;
	int j;
	int t;
	double t1;
	int nwalkers;
	int ndim;
	gw_ensemble * s;
	double accum;
	double mu[] = {1,1,1,1};
	
	ndim = 4;
	nwalkers = 1000;


	s = gw_alloc(4,nwalkers,ndim);
	s->lnprob = lnprob;
	gw_gaussian(s,mu,1);
	
	
	
	t = clock();

	for(i=0;i<=100;i++){
			accum += s->lnprob(s->w,ndim);
	}
	t1 = (double)(clock()-t)/((double)CLOCKS_PER_SEC);
	printf("(%f) 100x lnprob time: %f\n",accum, t1);
	fflush(stdout);
	
	

	for(j=1; j <= 4; j++){

		s->nthreads = j;

		t = clock();

		for(i=0;i<100;i++){
			gw_metropolis(s);
		}
		t1 = (clock()-t)/((double)CLOCKS_PER_SEC);
		printf("%f,", t1);
		fflush(stdout);

	}
	printf("\n");

	gw_free(s);
}
