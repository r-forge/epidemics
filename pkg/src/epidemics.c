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

/* seed new infections from a single pathogen */
void process_infection(struct pathogen * pat, struct population * pop, struct param * par){
	int i, nbnewinf=0, orinsus=get_nsus(pop), orininfcum=get_ninfcum(pop);

	if(pat != NULL){ /* if infection is not a gost */
		/* determine the number of descendents */
		nbnewinf = 2;
		printf("\nnumber of new infections %d",nbnewinf);
		if(get_age(pat) >= par->t1){
			nbnewinf = gsl_ran_poisson(par->rng, par->R);
			printf("\nnumber of new infections %d",nbnewinf);
		}

		/* adjust number of new infections to number of susceptibles */
		if(nbnewinf > orinsus) nbnewinf = orinsus;

		printf("\nnumber of new infections %d",nbnewinf);

		if(nbnewinf>0){
			/* for each new infection, add new pathogen */
			for(i=orininfcum;i<(orininfcum+nbnewinf);i++){
				(get_pathogens(pop))[i] = create_pathogen();
				replicate(pat, (get_pathogens(pop))[i], par);
			}
		}


		/* update number of susceptibles and infected */
		pop->nsus = orinsus - nbnewinf;
		pop->ninfcum = orininfcum + nbnewinf;
		pop->ninf = pop->ninf + nbnewinf;
	}
}





void age_population(struct population * pop, struct param * par){
	struct pathogen *ppat;
	int i;

	/* pathogens ages */
	for(i=0;i<get_orinsus(pop);i++){
		ppat = (pop->pathogens)[i]; /* to make code more readable*/
		if(ppat !=NULL){ /* if pathogen exists */
			ppat->age = ppat->age+1; /* get older */
			if(get_age(ppat) > par->t2) { /* die if you must */
				free_pathogen(ppat);
				(pop->pathogens)[i] = NULL;
				pop->nrec = pop->nrec + 1;
				pop->ninf = pop->ninf - 1;
			}
		}
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
	while(get_nsus(pop)>0 && get_ninf(pop)>0){
		printf("\n-- population a step %d",++nstep);
		print_population(pop);
		for(i=0;i<get_orinsus(pop);i++){
			process_infection(get_pathogens(pop)[i], pop, par);
		}
	}

	/* free memory */
	free_population(pop);
	free_param(par);
}




int main(){

	run_epidemics(100, 0.05, 50, 2, 10, 0,2);

	return 0;
}
