#include "gw.h"

gw_ensemble * gw_alloc(int nthreads,int nwalkers,int ndim){

	gw_ensemble * s;
	int i;
	int alloc_fail = 0;

	size_t w_size;
	size_t w_arr_size;
	size_t lnp_arr_size;
	size_t pad;
	size_t stride;

	double * lnp;

	/* For Goodman-Weare:
	 *
	 * Every thread is going to need to simultaneously
	 * read from one group of walkers while writing 
	 * to the other. Since a double is only 8 bytes,
	 * we have a signficant risk of cache contention when
	 * ndim is near 8. The risk is much worsened for
	 * low values of nwalkers.
	 *
	 * For Metropolis-Hastings:
	 * 
	 * Every thread is going to need to potentially read
	 * and write every walker at every iteration. 
	 *
	 * We have two options to deal with this. The first
	 * is to give each thread a separate copy of each 
	 * walker, and then synchronize all the copies in
	 * serial. Synchronization would be extremely
	 * expensive, however, given that way may have 
	 * several gigabytes worth of walkers.

	 * The second option, which we implement here, is 
	 * to pad the walkers inside the walker arrays and
	 * enforce cache line alignment of each walker. The
	 * price of this is that we have have to do some
	 * extra arithmetic to access a walker in the array. 
	 *
	 * See gw_inc and gw_get */

	s = (gw_ensemble *) malloc(sizeof(gw_ensemble));

	/* size of a single walker in memory */
	w_size = ndim * sizeof(double);

	/* stride will be the number of bytes between the 
	 * start of two adjacent walkers in the array */
	pad = w_size % CACHE_LINE_SIZE;
	stride = w_size + pad;
	w_arr_size = stride * nwalkers;

	/* Do the same thing for the lnprob values */
	lnp_arr_size = CACHE_LINE_SIZE * nwalkers;

	alloc_fail += posix_memalign((void**)&(s->w),
		CACHE_LINE_SIZE, w_arr_size);

	alloc_fail += posix_memalign((void**)&(s->lnp),
		CACHE_LINE_SIZE, lnp_arr_size);

	assert(alloc_fail == 0);

	for(i=0;i<lnp_arr_size;i++){
		lnp = gw_get_lnp(s,i);
		*lnp = 0;

	}


	s->a = 2;
	s->B = 1;
	s->accept = 0;
	s->mh_sigma = 1;
	s->nwalkers = nwalkers;
	s->ndim = ndim;
	s->stride = stride;
	s->w_arr_size = w_arr_size;
	s->lnp_arr_size = lnp_arr_size;
	s->lnp_arr_stride = CACHE_LINE_SIZE;

	/* going to need rand() later */
	srand(time(NULL));

	return s;
	
}

double * gw_get(gw_ensemble * s, int i) {
	return (double*)(((char*)s->w) + (s->stride * i));
}

double * gw_inc(size_t stride, double * w) {
	return (double*)(((char*)w) + stride);
}

double * gw_get_lnp(gw_ensemble * s, int i) {
	return (double*)(((char*)s->lnp) + (s->lnp_arr_stride * i));
}

void gw_metropolis(gw_ensemble * s) {
/* Do one metrpolis update of all the walkers
 * in the ensemble in parallel. */

	int tid,i,j,seed,accept;
	double lnp,lnpp;
	double * w;
	double *prop;
	double *_lnp;
	gsl_rng * r;

	accept = 0;
	seed = rand();

	omp_set_num_threads(s->nthreads);

	#pragma omp parallel \
		private(tid,r,w,prop,i,j,lnpp,lnp,_lnp) \
		firstprivate(seed,accept) \
		shared(s)
	{
		tid = omp_get_thread_num();

		/* Setup an rng for every thread. I would be 
		 * interested to know if there's a better way
		 * to ensure each thread has a fast and
		 * independent rng, however.  */

		 r = gsl_rng_alloc(gsl_rng_mt19937);

		/* Add seed to tid * time(NULL) to
		 * avoid making thread 0 always have
		 * the same seed at every iteration */

		 gsl_rng_set(r,seed + tid * time(NULL));

		 /* Every thread is going to need space
		  * to store a proposal move */
		 prop = (double*)malloc(sizeof(double) * s->ndim);

		 #pragma omp for
		 for(i=0;i<s->nwalkers;i++) {

		 	/* pointer to the walker of interest */
		 	w = gw_get(s,i);

			/* pick a gaussian proposal */
			for(j=0;j<s->ndim;j++){
				prop[j] = gsl_ran_gaussian(r,s->mh_sigma) + w[j];
			}

			/* get a pointer to log probability value of
			 * this walker position */
			_lnp = gw_get_lnp(s,i);

			/* Avoid calculating s->lnprob unless
			 * this is the first iteration */
			if (*_lnp == 0) {
				lnp = s->lnprob(w,s->ndim);
				*_lnp = lnp;
			}
			else {
				lnp = *_lnp;
			}
			
			lnpp = s->lnprob(prop,s->ndim);

			if( s->B * (lnpp-lnp) > 0 ||
				 log(gsl_ran_uniform(r)) < s->B * (lnpp-lnp)) {

				/* Accept proposed move*/

				accept += 1;
				*_lnp = lnpp;/* Save the lnprob value for posterity*/

				/* Copy over the proposal */

				for(j=0;j<s->ndim;j++){
					w[j] = prop[j];
				}

			}
		 }

		 free(r);
		 free(prop);

		 /* false sharing s->accept here should be ok for now*/
		 s->accept += accept;

	}
}

void gw_free(gw_ensemble * s) {

	free(s->lnp);
	free(s->w);

	free(s);
}

void gw_gaussian(gw_ensemble * s, double * mu, double sigma){
/* initialize ensemble with a symmetric gaussian ball */

	gsl_rng * r;
	double * w;
	int nwalkers;
	int stride;
	int ndim;
	int i;
	int j;

	ndim = s->ndim;
	stride = s->stride;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,rand());

	w = gw_get(s,0);
	for(i=0;i<nwalkers;i++){

		/* Recall that the ith walker spans the first
		 * ndim * sizeof(double) bytes after offset
		 * i*stride bytes  */		

		for(j=0;j<ndim;j++){
			w[j] = gsl_ran_gaussian(r,sigma) + mu[j];
		}

		w = gw_inc(s->stride,w);
	}
	gsl_rng_free(r);
}

void gw_write_text(gw_ensemble * s, char * fname) {

	int i;
	int j;
	double * w;

	FILE * f = fopen(fname,"w");

	w = gw_get(s,0);

	for(i=0;i<s->nwalkers;i++){
		for(j=0;j<s->ndim;j++){
			fprintf(f,"%f,",s->w[j]);
			fprintf(f,"\n");
		}
		w = gw_inc(s->stride,w);
	}
	fclose(f);
}


