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
/* R: R0, i.e. number of infections created by an 
   infectious individual in an entirely susceptible 
   population at each time step. */
/* rng: random number generator */
/* t1: number of unit time before a pathogen become infectious */
/* t2: number of unit time before a pathogen stops being infectious (i.e. dies) */
/* nstart: number of infections with wild genotype initiating the epidemic */
/* t_sample: time at which to sample the outbreak */
/* n_sample: sample size, in number of pathogens */
struct param{
	int L, nsus, t1, t2, nstart, t_sample, n_sample, npop;
	double mu, muL, R, *pdisp;
	gsl_rng * rng;
};




/* Free param */
void free_param(struct param *in);

/* Check parameters */
void check_param(struct param *in);

/* Print parameters */
void print_param(struct param *in);
