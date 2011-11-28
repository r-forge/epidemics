/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions handle simulation parameters.
*/




/* L: length of the genomes */
/* nsus: number of hosts (carrying capacity */
/* mu: mutation rate per site and generation */
/* muL = mu*L */
/* beta: transmission rate */
/* rng: random number generator */
/* t1: number of unit time before a pathogen become infectious */
/* t2: number of unit time before a pathogen stops being infectious (i.e. dies) */
/* nstart: number of infections with wild genotype initiating the epidemic */
/* t_sample: array of integers giving the times at which to sample each isolate (in time steps from most recent isolate) */
/* n_sample: sample size, in number of pathogens */
/* npop: number of populations in the metapopulation */
/* duration: maximum number of steps to run simulations for; implicitely the duration of the epidemic until most recent sample */
struct param{
	int L, t1, t2, nstart, *t_sample, n_sample, duration, npop, *popsizes;
	double mu, muL, beta, *pdisp;
	gsl_rng * rng;
};




/* Free param */
void free_param(struct param *in);

/* Check parameters */
void check_param(struct param *in);

/* Print parameters */
void print_param(struct param *in);
