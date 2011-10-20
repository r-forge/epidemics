/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "param.h"
#include "pathogens.h"
#include "populations.h"




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void process_infection(struct pathogen * pat, struct population * pop, struct param * par){
	int i, nbnewinf=0, orinsus=get_nsus(pop), orininfcum=get_ninfcum(pop);

	if(pat != NULL){ /* if infection is not a gost */
		/* determine the number of descendents */
		if(get_age(pat) > par->t1){
			nbnewinf = gsl_ran_poisson(par->rng, par->R);
		}

		/* adjust number of new infections to number of susceptibles */
		if(nbnewinf > orinsus) nbnewinf = orinsus;

		/* reallocate the pathogen vector */
		if(nbnewinf>0){
			/*pop->pathogens = realloc(pop->pathogens, (orininfcum+nbnewinf) * sizeof(struct pathogen *));*/

			/* for each new infection, add new pathogen */
			for(i=orininfcum;i<(orininfcum+nbnewinf);i++){
				replicate(pat, get_pathogens(pop)[i], par);
			}
		}

		/* pathogen ages */
		pat->age = pat->age + 1;
		if(get_age(pat) >= par->t2){
			/* pathogen dies: pointer turned to NULL */
			pat = NULL;
			pop->nrec = pop->nrec + 1;
			pop->ninf = pop->ninf - 1;
		}

		/* update number of susceptibles and infected */
		pop->nsus = orinsus - nbnewinf;
		pop->ninfcum = orininfcum + nbnewinf;
		pop->ninf = pop->ninf + nbnewinf;
	}
}





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void run_epidemics(int seqLength, double mutRate, int nHost, double Rzero, int nStart, int t1, int t2){
	int i, nstep=0;

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
	par->L = seqLength;
	par->mu = mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->K = nHost;
	par->R = Rzero;
	par->nstart = nStart;
	par->t1 = t1;
	par->t2 = t2;

	/* initiate population */
	struct population * pop;
	pop = create_population(par->K, par->nstart, 0);

	/* make population evolve */
	/* while(get_nsus(pop)>0 && get_ninf(pop)>0){ */
	/* 	printf("\n-- population a step %d",++nstep); */
	/* 	print_population(pop); */
	/* 	for(i=0;i<get_orinsus(pop);i++){ */
	/* 		process_infection(get_pathogens(pop)[i], pop, par); */
	/* 	} */
	/* } */

	/* free memory */
	free_population(pop);
	free_param(par);
}




int main(){

	run_epidemics(100, 0.01, 50, 3.5, 10, 1,3);

	return 0;
}
