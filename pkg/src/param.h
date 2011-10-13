/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions handle simulation parameters.
*/




/* L: length of the genomes */
/* mu: mutation rate per site and generation */
/* muL = mu*L */
struct param{
	int L;
	double mu, muL;
	gsl_rng * rng;
};




/* Free param */
void free_param(struct param *in);
