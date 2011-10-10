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
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "seqEvol.h"

/*
   ========================
   === GLOBAL VARIABLES ===
   ========================
*/

/* Initialize random number generators */
/* const gsl_rng *gBaseRand; */

/* /\* specifying to use Mersenne twister MT-19937 as the uniform PRNG *\/ */
/* gBaseRand = gsl_rng_alloc(gsl_rng_mt19937); */

/* srand(time(NULL));                    /\* initialization for rand() *\/ */
/* randSeed = rand();                    /\* returns a non-negative integer *\/ */
/* gsl_rng_set(gBaseRand, randSeed);    /\* seed the PRNG *\/ */

/* randomize(); /\* seeds the basic function random() from sdlib.h*\/ */




/*
   =================================
   === LOCAL AUXILIARY FUNCTIONS ===
   =================================
*/

/* Function to generate one mutation */
/* int create_mutation(int L){ */
/* 	int out; */
/* 	out = random(L)+1; */
/* 	return out; */
/* } */


/* Function replicating a genome, with mutations and back-mutations */
/* void replication(struct pathogen *in, struct pathogen *out, int nbmut, int nbbackmut, int L){ */
/* 	int i; */

/* 	for(i=0;i<nbmut;i++){ */
/* 		if(random(L)+1 < get_nb_snps(in)){ /\* back mutation *\/ */

/* 		} else { */

/* 		} */
/* 	} */
/* } */


/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Basic accessors for pathogen objects */
/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in){
	return in->length;
}

/* Returns the ID of the host, i.e. in->host array */
long long unsigned int get_host(struct pathogen *in){
	return in->host;
}

/* Returns SNP vector */
unsigned int * get_snps(struct pathogen *in){
	return in->snps;
}


/* Create empty pathogen */
struct pathogen * create_pathogen(){
	struct pathogen *out;
	out = (struct pathogen *) calloc(1, sizeof(struct pathogen));
	if(out == NULL){
		fprintf(stderr, "No memory left for creating initial pathogen. Exiting.\n");
		exit(1);
	}
	out->snps = NULL;
	out->length = 0;
	/* out->host = 1; /\* new infection starts with host 1 *\/ */
	/* LAST_HOST = 1; /\* update host pool *\/ */
	return out;
}


/* Copy pathogen */
/*  - memory allocation made outside the function */
void copy_pathogen(struct pathogen *in, struct pathogen *out){
	int i, N;

	N=get_nb_snps(in);

	out->snps = (unsigned int *) calloc(N, sizeof(unsigned int)); /* allocate memory for snps vector*/
	if(out->snps == NULL){
		fprintf(stderr, "No memory left for copying pathogen genome. Exiting.\n");
		exit(1);
	}

	for(i=0;i<N;i++){ /* copy snps */
		(out->snps)[i] = (get_snps(in)[i]);
	}
	out->length = N;
	out->host = ++LAST_HOST; /* this will have to be replaced */
}


/* Free pathogen */
void free_pathogen(struct pathogen *in){
	free(in->snps);
	free(in);
}


/* Print pathogen content */
void print_pathogen(struct pathogen *in){
	int i, N=get_nb_snps(in);
	printf("\n%d snps: ", N);
	for(i=0;i<N;i++) printf("%d", get_snps(in)[i]);
	printf("\nhost: %llu \n", get_host(in));
}



/* Replicate an isolate - mutation possible */
/* */
/* struct pathogen * replicate_pathogen(struct pathogen *in){ */
/* 	int i, nbmut=0, nbbackmut=0, *newmut; */
/* 	struct pathogen *out; */
/* 	out = (struct pathogen *) calloc(1, sizeof(struct pathogen)); */
/* 	if(out == NULL){ */
/* 		fprintf(stderr, "No memory left for creating new pathogen. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	/\* determine nb of mutations *\/ */
/* 	/\* has to take into account possible back-mutations *\/ */
/* 	/\* TODO *\/ */
/* 	nbmut = gsl_ran_binomial(gBaseRand,  MU, L); */
/* 	newmut = calloc(nbmut, sizeof(int)); */


/* 	/\* allocate new snps vector *\/ */
/* 	out->snps = calloc(get_nb_snps(in)+nbmut, sizeof(int)); */

/* 	/\* inherite SNPs of the ancestor, but skip the last nbbackmut *\/ */
/* 	for(i=0; i<(get_nb_snps(in)-nbbackmut); i++){ */
/* 		out->snps[i] = in->snps[i]; */
/* 	} */

/* 	/\* add new mutations*\/ */
/* 	for(i=0; i<nbmut;i++){ */
/* 		out->snps[i+get_nb_snps(in)] = create_mutation(); */
/* 	} */
/* 	out->length = in->length+nbmut; */
/* 	out->host = in->host; */
/* } */






/* free base rand */
/* gsl_rng_free(gBaseRand); */



/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/

main(){
	struct pathogen * ppat1, * ppat2;
	ppat1 = create_pathogen();
	ppat2 = create_pathogen();
	copy_pathogen(ppat1,ppat2);
	copy_pathogen(ppat1,ppat2);

	printf("\nppat 1");
	print_pathogen(ppat1);

	printf("\nppat 2");
	print_pathogen(ppat2);

	free_pathogen(ppat1);
	free_pathogen(ppat2);
}


/* TESTING in R */

/*
## test raw conversion
.C("testRaw", raw(256), 256L, PACKAGE="adegenet")


*/

