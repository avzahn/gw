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
	 * We have two options to deal with this. The first
	 * is to give each thread a separate copy of each 
	 * walker, and then synchronize all the copies in
	 * serial. Synchronization would be extremely
	 * expensive, however, given that we would have to do
	 * it every time we updated a walker.
	 *
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
	
	

	alloc_fail += posix_memalign((void**)&(s->w),
		CACHE_LINE_SIZE, w_arr_size);

	alloc_fail += posix_memalign((void**)&(s->lnp),
		CACHE_LINE_SIZE, lnp_arr_size);

	assert(alloc_fail == 0);

	for(i=0;i<nwalkers;i++){
		lnp = gw_get_lnp(s,i);
		*lnp = 0;

	}




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

void gw_copy(gw_ensemble * src, gw_ensemble * dst){

	memcpy((void*)dst->w,(void*)src->w,src->w_arr_size);
	memcpy((void*)dst->lnp,(void*)src->lnp,src->lnp_arr_size);
	dst->w_arr_size = src->w_arr_size;
	dst->lnp_arr_size = src->lnp_arr_size;
	dst->stride = src->stride;
	dst->lnp_arr_stride = src->lnp_arr_stride;
	
}

void gw_copy_partial(gw_ensemble * src, gw_ensemble * dst, int idx0, int idx1){

	double * start, end;
	int size,size_lnp;

	size = src->stride * (idx1-idx0);
	size_lnp = src->lnp_arr_stride * (idx1-idx0);

	memcpy((void*)gw_get(dst,idx0),(void*)gw_get(src,idx0),size);
	memcpy((void*)gw_get_lnp(dst,idx0),(void*)gw_get_lnp(src,idx0),size_lnp);
	dst->w_arr_size = src->w_arr_size;
	dst->lnp_arr_size = src->lnp_arr_size;
	dst->stride = src->stride;
	dst->lnp_arr_stride = src->lnp_arr_stride;
	
}

void gw_metropolis2(gw_ensemble * s, int nsteps) {
	
	/* Not a stand in for the stretch move. Divide up
	 * the ensemble over a team of threads and have each
	 * thread independently iterate its subensemble
	 * on a local copy that gets copied back at the end of the
	 * function. */

	int tid,i,j,k,seed,accept,nwalkers;
	double lnp,lnpp,B,mh_sigma;
	double * w;
	double *prop;
	double *_lnp;
	gsl_rng * r;
	gw_ensemble * local;
	int ndim;

	accept = 0;
	seed = rand();
	mh_sigma = s->mh_sigma;
	nwalkers = s->nwalkers;
	ndim = s->ndim;
	B = s->B;
	
	omp_set_num_threads(s->nthreads);	
	
	#pragma omp parallel \
	private(local,tid,r,w,prop,lnp,_lnp,i,j,k) \
	firstprivate(seed,ndim,B,nsteps,nwalkers,mh_sigma) \
	shared(s,accept)
	{
		tid = omp_get_thread_num();
		r = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(r,seed + tid * time(NULL));
		prop = (double*)malloc(sizeof(double) * s->ndim);
		local = gw_alloc(s->nthreads,nwalkers,ndim);
		local->lnprob = s->lnprob;
		gw_copy(s,local);
		
		#pragma omp for reduction(+:accept)
		for(j = 0; j < nwalkers; j++){
			for(i=0;i<nsteps;i++){
				
				fprintf(stderr,"%i: %i, %i\n",tid,i,j);
				
			
				/* pointer to the walker of interest */
			 	w = gw_get(local,i);
			 	

	
				/* pick a gaussian proposal */
				for(j=0;j<ndim;j++){
					prop[j] = gsl_ran_gaussian(r,mh_sigma) + w[j];
				}
	
				/* get a pointer to log probability value of
				 * this walker position */
				_lnp = gw_get_lnp(local,i);
				
																
	
				/* Avoid calculating lnprob unless
				 * this is the first iteration */
				if (*_lnp == 0) {
					
					lnp = local->lnprob(w,ndim);
					
					*_lnp = lnp;					
				}
				else {
					lnp = *_lnp;
				}
				
											
				lnpp = s->lnprob(prop,ndim);
		
				
				if( s->B * (lnpp-lnp) > 0 ||
					 log(gsl_ran_flat(r,0,1)) < B * (lnpp-lnp)) {
	
					/* Accept proposed move*/
	
					accept += 1;
					*_lnp = lnpp;/* Save the lnprob value for posterity*/
	
					/* Copy over the proposal */

	
					for(k=0;k<s->ndim;k++){
						w[k] = prop[k];
					}
				}

			}
			/* Bottlenck step--every thread copies its results back */
			gw_copy_partial(local,s,j,j+1);
		}

		
		gsl_rng_free(r);
		gw_free(local);
		free(prop);
		
		
	}
	
	
}


void gw_metropolis(gw_ensemble * s) {
/* Do one metrpolis update of all the walkers
 * in the ensemble in parallel. */
 
 
/* This functions is a stand-in for the affine invariant 
 * stretch move; it happens do to a metrpolis step, but it
 * reads and writes a shared ensemble at every step exacty
 * the way the stretch move would have to */

	int tid,i,j,seed,accept;
	double lnp,lnpp;
	double * w;
	double *prop;
	double *_lnp;
	gsl_rng * r;
	
	int ndim,nwalkers;
	double B, mh_sigma;
	
	B = s->B;
	mh_sigma = s->mh_sigma;
	ndim = s->ndim;
	nwalkers = s->nwalkers;

	accept = 0;
	seed = rand();

	omp_set_num_threads(s->nthreads);

	#pragma omp parallel \
		private(tid,r,w,prop,i,j,lnpp,lnp,_lnp) \
		firstprivate(seed,B,ndim,nwalkers,mh_sigma) \
		shared(s,accept)
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
		 prop = (double*)malloc(sizeof(double) * ndim);
		 


		 #pragma omp for reduction(+:accept)
		 for(i=0;i<nwalkers;i++) {

		 	/* pointer to the walker of interest */
		 	w = gw_get(s,i);

			/* pick a gaussian proposal */
			for(j=0;j<ndim;j++){
				prop[j] = gsl_ran_gaussian(r,mh_sigma) + w[j];
			}

			/* get a pointer to log probability value of
			 * this walker position */
			_lnp = gw_get_lnp(s,i);

			/* Avoid calculating s->lnprob unless
			 * this is the first iteration */
			if (*_lnp == 0) {
				lnp = s->lnprob(w,ndim);
				*_lnp = lnp;
			}
			else {
				lnp = *_lnp;
			}
			
			lnpp = s->lnprob(prop,ndim);

			if( s->B * (lnpp-lnp) > 0 ||
				 log(gsl_ran_flat(r,0,1)) < B * (lnpp-lnp)) {

				/* Accept proposed move*/

				accept += 1;
				*_lnp = lnpp;/* Save the lnprob value for posterity*/

				/* Copy over the proposal */

				for(j=0;j<ndim;j++){
					w[j] = prop[j];
				}

			}
		 }

		 gsl_rng_free(r);
		 free(prop);

	}
    s->accept = accept;
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

	

	for(i=0;i<s->nwalkers;i++){
		w = gw_get(s,i);
		for(j=0;j<s->ndim;j++){
			fprintf(f,"%f,",s->w[j]);
			fprintf(f,"\n");
		}
	}
	fclose(f);
}


