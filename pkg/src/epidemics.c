/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "param.h"
#include "seqEvol.h"
#include "populations.h"




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void make_new_infections(struct * pathogen pat, struct *population pop, struct param * par){
	int nbnewinf=0;
	/* determine the number of descendents */

	/* reallocate the pathogen vector */
	if(nbnewinf>0){
		pop->pathogens = realloc(pop->pathogens, sizeof(struct pathogen *))
	}
	for(i=0;i<nbnewinf,i++){/* for all new infection, add new pathogen */
	}
}

/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void run_epidemics(int seqLength, double mutRate, int nHost){
	/* Initialize random number generator */
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	/* simulation parameters */
	struct param * par;
	par = (struct param *) calloc(1, sizeof(struct param));
	par->L = seqLenth;
	par->mu = 0.mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->K = nHost;

	/* initiate population */
	struct population * pop;


}




int main(){

	


	struct population * pop;

	pop = create_population(1000,10,0);

	print_population(pop);

	/* free memory */
	free_population(pop);
	gsl_rng_free(rng);

	return 0;
}
