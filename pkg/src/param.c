/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions handle simulation parameters.
*/

#include "common.h"
#include "param.h"


/* Free param */
void free_param(struct param *in){
	gsl_rng_free(in->rng);
	free(in);
}



void check_param(struct param *in){
	/* nstart & K */
	if(in->nstart > in->K){
		fprintf(stderr, "\nParameter error: initial number of infections greater than host population.\n");
		exit(1);
	}

	/* K */
	if(in->K < 1){
		fprintf(stderr, "\nParameter error: less than one host in population.\n");
		exit(1);
	}

	/* nstart */
	if(in->nstart < 1){
		fprintf(stderr, "\nParameter error: less than one initial infection.\n");
		exit(1);
	}

	/* L */
	if(in->L < 1){
		fprintf(stderr, "\nParameter error: genome length less than one.\n");
		exit(1);
	}

	if(in->L < 1){
		fprintf(stderr, "\nParameter error: genome length less than one.\n");
		exit(1);
	}

	/* t1 */
	if(in->t1 < 1){
		fprintf(stderr, "\nParameter error: time to infectiousness less than one.\n");
		exit(1);
	}

	/* t2 */
	if(in->t2 < 1){
		fprintf(stderr, "\nParameter error: infection duration less than one.\n");
		exit(1);
	}

	/* t1 & t2 */
	if(in->t1 > in->t1){
		fprintf(stderr, "\nParameter error: infection last less than time to infectiousness.\n");
		exit(1);
	}


	/* mu */
	if(in->mu < 0.0){
		fprintf(stderr, "\nParameter error: negative mutation rate.\n");
		exit(1);
	}

	/* R */
	if(in->R < 0.0){
		fprintf(stderr, "\nParameter error: negative reproductive number.\n");
		exit(1);
	}

}
