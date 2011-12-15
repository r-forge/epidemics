/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "auxiliary.h"
#include "param.h"
#include "pathogens.h"
#include "populations.h"
#include "dispersal.h"
#include "infection.h"




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* seed new infection from a single pathogen */
/* void make_new_infection(struct pathogen * ances, struct population * pop, struct param * par){ */

/* 	/\* GENERATE ERROR IF PATHOGEN IS INACTIVATED *\/ */
/* 	if(!isNULL_pathogen(pat) && get_age(pat) >= par->t1){ /\* if ancestor is active *\/ */
/* 		/\* UPDATE NUMBER OF SUSCEPTIBLES AND INFECTED IN THE POPULATION *\/ */
/* 		pop->nexpcum = pop->nexpcum + 1; */
/* 		pop->ninf = pop->ninf + 1; */
/* 		pop->nsus = pop->nsus - 1; */

/* 		/\* HANDLE GENOME REPLICATION *\/ */
/* 		replicate(pat, get_pathogens(pop)[Nexpcum], par); */
/* 		} */
/* } /\* end make_new_infection *\/ */






/* PROCESS ALL INFECTIONS IN ONE GIVEN POP, FOR ONE GIVEN TIME STEP */
void process_infections(struct population * pop, struct metapopulation * metapop, struct network *cn, struct param * par){
	int i, k, count, popid=get_popid(pop), nbNb=cn->nbNb[popid], nbnewcases, *nbnewcasesvec;
	double *lambdavec, lambda=0, proba=0;
	struct pathogen * ances;
	struct population *curpop;

	/* COMPUTE \lambda_j = \beta w_{j->k} I_j/N_j for each neighbouring population j */
	/* \lambda = \sum_j \lambda_j */
	lambdavec = (double *) malloc(nbNb * sizeof(double));
	for(i=0;i<nbNb;i++){
		curpop = metapop->populations[cn->listNb[popid][i]];
		lambdavec[i] = par->beta * cn->weights[popid][i] * ((double) get_ninf(curpop))/get_popsize(curpop);
		lambda += lambdavec[i];
	}

	/* COMPUTE PROBABILITY OF INFECTION PER SUSCEPTIBLE */
	proba = 1 - exp(-lambda);

	/* FIND NB OF NEW INFECTIONS SEEDED IN POP BY EACH NEIGHBOURING POPULATION */
	nbnewcases = gsl_ran_binomial(par->rng, proba, get_nsus(pop));

	/* DRAW NB OF ANCESTORS IN EACH NEIGHBOURING POPULATION */
	nbnewcasesvec = malloc(nbNb * sizeof(int));
	gsl_ran_multinomial(par->rng, nbNb, nbnewcases, lambdavec, (unsigned int *) nbnewcasesvec);

	/* PRODUCE NEW PATHOGENS */
	count = 0;
	for(k=0;k<nbNb;k++){
		curpop = metapop->populations[cn->listNb[popid][k]];
		for(i=0;i<nbnewcasesvec[k];i++){
			/* determine ancestor */
			ances = select_random_infectious_pathogen(curpop, par);
			/* produce new pathogen */
			pop->pathogens[pop->nexpcum + count++] = replicate(ances, par);
		}
	}

	/* UPDATE GROUP SIZES */
	pop->nsus = pop->nsus - nbnewcases;
	pop->nexpcum = pop->nexpcum + nbnewcases;
	pop->nexp = pop->nexp + nbnewcases;


	/* FREE MEMORY AND RETURN */
	free(lambdavec);
	free(nbnewcasesvec);
} /* end  process_infections */






/* gcc line:

   gcc -o infection param.c auxiliary.c pathogens.c populations.c dispersal.c infection.c -Wall -O0 -lgsl -lgslcblas

   valgrind --leak-check=yes infection

*/


/* int main(){ */
/* 	/\* Initialize random number generator *\/ */
/* 	time_t t; */
/* 	t = time(NULL); // time in seconds, used to change the seed of the random generator */
/* 	gsl_rng * rng; */
/* 	const gsl_rng_type *typ; */
/* 	gsl_rng_env_setup(); */
/* 	typ=gsl_rng_default; */
/* 	rng=gsl_rng_alloc(typ); */
/* 	gsl_rng_set(rng,t); // changes the seed of the random generator */
/* 	int i, j; */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->rng = rng; */
/* 	par->npop = 2; */
/* 	int popsizes[2] = {1000,200}; */
/* 	par->popsizes = popsizes; */
/* 	par->nstart = 10; */
/* 	par->t1 = 1; */
/* 	par->t2 = 2; */
/* 	par->beta = 1.1; */
/* 	int nbnb[2] = {2,2}; */
/* 	par->cn_nb_nb = nbnb; */
/* 	int listnb[4] = {0,1,1,0}; */
/* 	par->cn_list_nb = listnb; */
/* 	double weights[4] = {0.9,0.1,0.99,0.11}; */
/* 	par->cn_weights = weights; */
/* 	struct network *cn = create_network(par); */
/* 	par->mu = 0.01; */
/* 	par->L = 100; */
/* 	par->muL = par->mu*par->L; */

/* 	/\* CREATE METAPOPULATION *\/ */
/* 	struct metapopulation * metapop = create_metapopulation(par); */
/* 	printf("\n## CREATED METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */


/* 	/\* SIMULATE OUTBREAK OVER A FEW TIMESTEPS *\/ */
/* 	for(i=0;i<100;i++){ */
/* 		age_metapopulation(metapop, par); */
/* 		for(j=0;j<get_npop(metapop);j++){ */
/* 			process_infections(get_populations(metapop)[j], metapop, cn, par); */
/* 		} */
/* 		printf("\n - METAPOPULATION @ step %d -", i); */
/* 		print_metapopulation(metapop, FALSE); */

/* 	} */

/* 	printf("\n## RESULTING METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */

/* 	/\* free memory *\/ */
/* 	free_metapopulation(metapop); */
/* 	free_network(cn); */
/* 	free(par); */
/* 	gsl_rng_free(rng); */

/* 	return 0; */
/* } */
