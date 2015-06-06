#include "pt.h"

/* Inverse temperature ladder. Must be even in length */
double ladder[] = {16,8,4,2};
int nnodes;
int nthreads;
int nwalkers;
int ndim;
int nsteps;
int here;
int there;
int size;

gw_ensemble * self;
gw_ensemble * partner;

double * lnprob(double * walker, int ndim){
	return -walker[0]*walker[0]
}


int main(int argc, char ** argv ) {

	int i;

	MPI_init(&argc,&argv);

	nnodes = 4;
	nthreads = 4;
	ndim = 1;
	nsteps = 100000;

	MPI_Comm_rank(MPI_COMM_WORLD, here);
	MPI_Comm_size(MPI_COMM_WORLD, size);

	if (here == 0){
		there = size-1;
	}
	else if (here == size -1){
		there = 0;
	}

	self = gw_alloc(nthreads,nwalkers,ndim);
	partner = gw_alloc(nthreads,nwalkers,ndim);

	self->B = ladder[here];
	partner->B = ladder[there];
	self->lnprob = lnprob;
	self->ndim = ndim;
	self->nwalkers = nwalkers;
	self->nnodes = nnodes;
	self->nthreads = nthreads;
	partner->ndim = ndim;
	partner->nwalkers = nwalkers;
	partner->nnodes = nnodes;
	partner->nthreads = nthreads;

	/* intialize the ensemble somewhere
	 * far away from our target distribution,
	 * which is centered at zero. */
	gw_gaussian(self,100,10);

	for(i=0; i<nsteps; i++){

		gw_metropolis(self);
		pt_mpi_temper(self,partner);

	}

	if (here == 0){
		gw_write_text(self,"mpi_test.out")
	}

	gw_free(self);
	gw_free(partner);

	MPI_Finalize();

	return 0;
}
