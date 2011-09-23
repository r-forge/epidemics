/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
    
#include "seqEvol.h"

const gsl_rng *gBaseRand;

/* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
  
srand(time(NULL));                    /* initialization for rand() */
randSeed = rand();                    /* returns a non-negative integer */
gsl_rng_set(gBaseRand, randSeed);    /* seed the PRNG */

randomize(); /* seeds the basic function random() from sdlib.h*/


/*
   =============================
   === STRUCTURES DEFINITION ===
   =============================
*/




/*
   =================================
   === LOCAL AUXILIARY FUNCTIONS ===
   =================================
*/

/* Function to generate one mutation */
int create_mutation(int *L){
	randomize(); /* seeds the basic function random() from sdlib.h*/
	int out;
	out = random(*L)+1;
	return out;
}


/* Function replicating a genome, with mutations and back-mutations */
void replication(struct pathogen *in, int nbmut, int nbbackmut, int *L){
	int i;
	randomize(); /* seeds the basic function random() from sdlib.h*/
		
	for(i=0;i<nbmut;i++){
		if(random(*L)<get_nb_snps(in)){ /* back mutation */
			
		} else {
		
		}
	}
}


/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Basic accessors for pathogen objects */
/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in){
	return pathogen->length;
}

/* Returns the ID of the host, i.e. in->host array */
long long int get_host(struct pathogen *in){
	return pathogen->host;
}



/* Create empty pathogen */
struct pathogen * create_initial_pathogen(){
	struct pathogen *out;
	out = (struct pathogen *) calloc(1, sizeof(struct pathogen));
	if(out == NULL){
		fprintf(stderr, "No memory left for creating initial pathogen. Exiting.\n");
		exit(1);
	}
	out->snps = NULL;
	out->length = 0;
	out->host = 1; /* new infection starts with host 1 */
	LAST_HOST = 1; /* update host pool */
	return out;
}



/* Create a pathogen from an existant isolate */
/* */
struct pathogen * create_new_pathogen(struct pathogen *in){
	int i, nbmut=0, nbbackmut=0, *newmut;
	struct pathogen *out;
	out = (struct pathogen *) calloc(1, sizeof(struct pathogen));
	if(out == NULL){
		fprintf(stderr, "No memory left for creating new pathogen. Exiting.\n");
		exit(1);
	}
	
	/* determine nb of mutations */
	/* has to take into account possible back-mutations */
	/* TODO */
	nbmut = gsl_ran_binomial(gBaseRand,  MU, L);
	newmut = calloc(nbmut, sizeof(int));


	/* allocate new snps vector */
	out->snps = calloc(get_nb_snps(in)+nbmut, sizeof(int));
	
	/* inherite SNPs of the ancestor, but skip the last nbbackmut */
	for(i=0; i<(get_nb_snps(in)-nbbackmut); i++){
		out->snps[i] = in->snps[i];
	}
	
	/* add new mutations*/
	for(i=0; i<nbmut;i++){
		out->snps[i+get_nb_snps(in)] = create_mutation();
	}
	out->length = in->length+nbmut;
	out->host = in->host;
}






/* free base rand */
gsl_rng_free(gBaseRand);



/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/



/* TESTING in R */

/*
## test raw conversion
.C("testRaw", raw(256), 256L, PACKAGE="adegenet")


*/

