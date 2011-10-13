/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions handle simulation parameters.
*/




/* Free param */
void free_param(struct param *in){
	gsl_rng_free(in->rng);
}
