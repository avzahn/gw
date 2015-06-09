#ifndef __gw__
#define __gw__

#include <gsl/gsl_rng.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Cache line size in bytes. Make sure this
 * value is right before compiling... */
#define CACHE_LINE_SIZE 64

typedef struct gw_ensemble {

/* Hold all the data needed by a fairly arbitrary MCMC.
 * Must be allocated with gw_alloc and freed with
 * gw_free.
 *
 * Note that walker positions need to be accessed with
 * gw_inc and gw_get, if you want to avoid knowing how
 * the walker array is laid out in memory. */

	/* Note: nwalkers must be even */

	int ndim; /* length of a walker, or dimension of parameter space */
	int nwalkers; /* number of walkers managed per node */
	int accept;
	int nthreads; /* number of openmp threads*/
	int nnodes; /* number of mpi nodes */
	double a; /* Goodman-Weare stretch parameter, almost always 2*/
	double mh_sigma; /* metropolis proposal width */
	double B; /* inverse temperature */
	
	/* traverse these two with stride */
	size_t stride; /* bytes */
	double *  w; /* walkers in the ensemble */
	size_t w_arr_size; /* size in bytes of w */

	size_t lnp_arr_stride;
	size_t lnp_arr_size;
	double * lnp; /* save last values of lnprob */

	double (*lnprob)(double * walker, int ndim);

} gw_ensemble;

gw_ensemble * gw_alloc(int nthreads,int nwalkers,int ndim);

void gw_free(gw_ensemble * s);

double * gw_get(gw_ensemble * s, int i);

double * gw_inc(size_t stride, double * w);

double * gw_get_lnp(gw_ensemble * s, int i);

double * gw_inc_lnp(size_t * stride, double * lnp);

void gw_stretch(gw_ensemble * s);

void gw_metropolis(gw_ensemble * s);

void gw_gaussian(gw_ensemble * s, double * mu, double sigma);

void gw_write_text(gw_ensemble * s, char * fname);

void gw_write_hdf(gw_ensemble * s, char * fname);

#endif
