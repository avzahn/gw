#include "gw.h"

/* gaussian centered at 0 */
double lnprob(double * walker, int ndim){
		return -walker[0] * walker[0];
}

int main(){
	
	int i;
	double mu[] = {1,1,1,1};

	gw_ensemble * s;
	s = gw_alloc(1,1000,4);
	s->lnprob = lnprob;
	
	/* initialize with some other guassian centered elsewhere */
	gw_gaussian(s,mu,1);
	
	
	for(i=0;i<1000;i++){
		gw_metropolis(s);
	}
	
	
	gw_write_text(s,"gaussian.out");
	
	gw_free(s);
	
	
}
