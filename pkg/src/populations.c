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
   =================
   === ACCESSORS ===
   =================
*/

unsigned int get_nsus(struct population *in){
	return in->nsus;
}


unsigned int get_ninf(struct population *in){
	return in->ninf;
}


unsigned int get_nrec(struct population *in){
	return in->nrec;
}






/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create new population */
struct population * create_population(unsigned int ns, unsigned int ni, unsigned int nr){
	
	/* create pointer to population */
	struct population *out;
	out = (struct population *) calloc(1, sizeof(struct population));
	if(out == NULL){
		fprintf(stderr, "\nNo memory left for creating new population. Exiting.\n");
		exit(1);
	}

	/* create the content */
	out->nsus = ns;
	out->ninf = ni;
	out->nrec = nr;

	/* susceptibles */
	/*if(ns==0){
		out->sus = NULL;
	} else {
		out->sus = (struct host **) calloc(nsus, sizeof(struct host *));
		if(out->sus == NULL){
			fprintf(stderr, "\nNo memory left for creating susceptibles in new population. Exiting.\n");
			exit(1);
		}
		for(i=0;i<nsus;i++){
			(out->sus)[i] = (struct host*) calloc(1, sizeof(struct host));
			if((out->sus)[i]==NULL){
				fprintf(stderr, "\nNo memory left for creating susceptibles in new population. Exiting.\n");
				exit(1);
			}
		}

	}*/
	
	/* infected */
	/*if(ninf==0){
		out->inf = NULL;
	} else {
		out->inf = (struct host *) calloc(ninf, sizeof(struct host));
		if(out->inf == NULL){
			fprintf(stderr, "\nNo memory left for creating infected in new population. Exiting.\n");
			exit(1);
		}
		for(i=0;i<ninf;i++){
			(out->inf)[i] = (struct host*) calloc(1, sizeof(struct host));
			if((out->inf)[i]==NULL){
				fprintf(stderr, "\nNo memory left for creating infected in new population. Exiting.\n");
				exit(1);
			}
		}
	}*/

	/* recovered */
	/*if(nrec==0){
		out->rec = NULL;
	} else {
		out->rec = (struct host *) calloc(nrec, sizeof(struct host));
		if(out->rec == NULL){
			fprintf(stderr, "\nNo memory left for creating recovered in new population. Exiting.\n");
			exit(1);
		}
		for(i=0;i<nrec;i++){
			(out->rec)[i] = (struct host*) calloc(1, sizeof(struct host));
			if((out->rec)[i]==NULL){
				fprintf(stderr, "\nNo memory left for creating recovered in new population. Exiting.\n");
				exit(1);
			}
		}
	}*/

	return out;
}







/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_population(struct population *in){
	/*int i;
	for(i=0;i<get_nsus(in);i++){
		free_host(in->sus);
	}
	for(i=0;i<get_ninf(in);i++){
		free_host(in->inf);
	}
	for(i=0;i<get_nrec(in);i++){
		free_host(in->rec);
	}
	free(in->sus);
	free(in->inf);
	free(in->rec);*/
	free(in);
}









/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Print population content */
void print_population(struct population *in){
	printf("\nnb susceptible: %d", get_nsus(in));
	printf("\nnb infected: %d", get_ninf(in));
	printf("\nnb recovered: %d\n", get_nrec(in));
}





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
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

	
/* 	/\* simulation parameters *\/ */
/* 	/\* struct param * par; *\/ */
/* 	/\* par = (struct param *) calloc(1, sizeof(struct param)); *\/ */
/* 	/\* par->L = 100; *\/ */
/* 	/\* par->mu = 0.01; *\/ */
/* 	/\* par->muL = par->mu * par->L; *\/ */
/* 	/\* par->rng = rng; *\/ */


/* 	struct population * pop; */

/* 	pop = create_population(1000,10,0); */

/* 	print_population(pop); */

/* 	/\* free memory *\/ */
/* 	free_population(pop); */
/* 	gsl_rng_free(rng); */

/* 	return 0; */
/* } */
