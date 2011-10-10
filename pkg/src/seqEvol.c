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
/* #include "gsl/gsl_rng.h" */
/* #include "gsl/gsl_randist.h" */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "seqEvol.h"





/*
   =================================
   === LOCAL AUXILIARY FUNCTIONS ===
   =================================
*/




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
void copy_pathogen(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, N;

	N=get_nb_snps(in);

	out->snps = (unsigned int *) calloc(N, sizeof(unsigned int)); /* allocate memory for snps vector*/
	if(get_snps(out) == NULL){
		fprintf(stderr, "No memory left for copying pathogen genome. Exiting.\n");
		exit(1);
	}

	for(i=0;i<N;i++){ /* copy snps */
		(out->snps)[i] = get_snps(in)[i];
	}
	out->length = N;
	out->host = get_host(in);
}





/* Function replicating a genome, with mutations and back-mutations */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, nbmut=0, nbbackmut=0, newsize, N, checkback;
	
	nbmut = gsl_ran_poisson(par->rng, par->muL);
	printf("mutation rate %.2f\n", par->mu);
	printf("length %d\n", par->L);
	printf("number of mutations %d\n", nbmut);
	
	/* determine the number of reverse mutations */
	if(nbmut>0){
		for(i=0;i<nbmut;i++){
			checkback = gsl_rng_uniform_int(par->rng, par->L)+1;
			if(checkback <= get_nb_snps(in)) nbbackmut++;
		}
	}
	
	nbmut -= nbbackmut; /* remove back mutations from new mutations */
	newsize = get_nb_snps(in) + nbmut - nbbackmut;

	/* reallocate memory for new vector */
	out->snps = (unsigned int *) calloc(newsize, sizeof(unsigned int));
	if(get_snps(out) == NULL){
		fprintf(stderr, "No memory left for replicating pathogen genome. Exiting.\n");
		exit(1);
	}

	/* inherit parental snps vector (except back mutations) */
	N = get_nb_snps(in)-nbbackmut; /* indices<N, ancestral snps; indices>=N, new snps */
	for(i=0;i<N;i++){ /* copy snps */
		out->snps[i] = get_snps(in)[i];
	}

	/* add new mutations */
	for(i=0;i<nbmut;i++){
/*careful: in theory, here it is possible to get twice the same mutation (albeit highly unprobably)*/
		(out->snps)[N+i] = gsl_rng_uniform_int(par->rng,par->L)+1;
	}

	/* finish to update new patogen data */
	out->length = newsize;
	par->lasthost += 1;
	out->host = par->lasthost;

} /*end replicate*/





/* Free pathogen */
void free_pathogen(struct pathogen *in){
	free(in->snps);
	free(in);
}


/* Print pathogen content */
void print_pathogen(struct pathogen *in){
	int i, N=get_nb_snps(in);
	printf("\n%d snps: ", N);
	for(i=0;i<N;i++) printf("%d ", get_snps(in)[i]);
	printf("\nhost: %llu \n", get_host(in));
}






/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/

void main(){
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
	par->L = 50;
	par->mu = 0.1;
	par->muL = par->mu * par->L;
	par->rng = rng;

	/* Test pathogen creation, copy */
	struct pathogen * ppat1, * ppat2;
	ppat1 = create_pathogen();
	ppat2 = create_pathogen();
	copy_pathogen(ppat1,ppat2,par);

	printf("\nppat 1");
	print_pathogen(ppat1);

	printf("\nppat 2");
	print_pathogen(ppat2);

	replicate(ppat1,ppat2,par);

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

