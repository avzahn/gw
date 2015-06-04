#include "gw.h"
#include <mpi.h>

void pt_recv_ensemble(gw_ensemble * s, int src) {
/* Fill a gw_ensemble with walkers and log probabilities
 * from an ensemble on another node. Filling in the
 * other fields however is the user's job.
 *
 * In particular, this won't work properly without
 * filling in the array size fields correctly first. */

	MPI_Recv(w,s->w_arr_size,MPI_BYTE,src,
		0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(lnp,s->lnp_arr_size,MPI_BYTE,src,
		1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void pt_send_ensemble(gw_ensemble * s, int dest) {

	MPI_Send(w,s->w_arr_size,MPI_BYTE,src,
		0,MPI_COMM_WORLD);
	MPI_Send(lnp,s->lnp_arr_size,MPI_BYTE,src,
		1,MPI_COMM_WORLD);

}

void pt_update(gw_ensemble * s0, gw_ensemble * s1) {
/* Exchange walker data between two gw_ensembles. As before,
 * don't do anything about any of the other fields. */

	double * tmp;

	tmp = s0->w;
	s0->w = s1->w;
	s1->w = tmp;

	tmp = s0->lnp;
	s0->lnp = s1->lnp;
	s1->lnp = tmp;
}

void pt_swap(gw_ensemble * s0, gw_ensemble * s1) {
/* Implement a single tempering step.

 * To perform a tempering step, we do one 
 * metropolis iteration for each walker in
 * s0, drawing a random position from s1
 * as the proposal position. If the proposal
 * is accepted, we swap walkers between
 * ensembles. 

 * Only updates s0's acceptance count */

 	int i;
 	int j;
 	int tid;
 	double seed;
 	double lnpdiff;
 	double * tmp;
	gsl_rng * r;

	accept = 0;
	seed = rand();

	omp_set_num_threads(s0->nthreads);

 	#pragma omp parallel \
 		share(s0,s1) \
 		private(r,tid,lnpdiff,i,j,tmp) \
 		firstprivate(seed,accept)
 	{
 		tid = omp_get_thread_num();
		r = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(r,seed + tid * time(NULL));

		/* scratch memory for swapping */
		tmp = (double*)malloc(s0->stride);

		#pragma omp for 
		for(i=0;i<s0->nwalkers;i++) {

			lnpdiff = *gw_get_lnp(s0,i) - *gw_get_lnp(s1,i);
			lnpdiff *= ((s0->B)-(s1->B));

			if(lnpdiff > 0 || log(gsl_rng_uniform(r)) < lnpdiff) {

				/* Accept the swap */
				memcpy(tmp,gw_get(s0,i),
					s0->ndim*sizeof(double));
				memcpy(gw_get(s0,i),gw_get(s1,i),
					s0->ndim*sizeof(double));
				memcpy(gw_get(s1,i),tmp,
					s0->ndim*sizeof(double));

				/* recycle lnpdiff as space to swap */
				lnpdiff = *gw_get_lnp(s0,i);
				*gw_get_lnp(s0,i) = *gw_get_lnp(s1,i);
				*gw_get_lnp(s1,i) = lnpdiff;

				accept += 1;
			}

			/* I think this is an ok way to do this */
			s0->accept += accept;
		}

		free(tmp);
		free(r);
	}
}

