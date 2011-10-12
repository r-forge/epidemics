/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "seqEvol.h"



/*
   =================
   === ACCESSORS ===
   =================
*/

unsigned int get_host_id(struct host *in){
	return in->id;
}




unsigned int get_host_pop(struct host *in){
	return in->pop;
}




struct host * get_s(struct pop *in){
	return in->s;
}




struct host * get_i(struct pop *in){
	return in->i;
}




struct host * get_r(struct pop *in){
	return in->r;
}




unsigned int * get_ns(struct pop *in){
	return in->ns;
}




unsigned int * get_ni(struct pop *in){
	return in->ni;
}




unsigned int * get_nr(struct pop *in){
	return in->nr;
}








/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty host */
struct host * create_host(){
	struct host *out;
	out = (struct host *) calloc(1, sizeof(struct host));
	if(out == NULL){
		fprintf(stderr, "\nNo memory left for creating initial host. Exiting.\n");
		exit(1);
	}
	out->id = 1;
	return out;
}




/* Create new population */
struct pop * create_pop(unsigned int ns, unsigned int ni, unsigned int nr){
	/* create pointer to pop */
	struct pop *out;
	out = (struct pop *) calloc(1, sizeof(struct pop));
	if(out == NULL){
		fprintf(stderr, "\nNo memory left for creating new population. Exiting.\n");
		exit(1);
	}

	/* create the content */
	/* susceptibles */
	if(ns==0){
		out->s = NULL;
	} else {
		out->s = (struct host *) calloc(ns, sizeof(struct host));
		if(out->s == NULL){
			fprintf(stderr, "\nNo memory left for creating susceptibles in new population. Exiting.\n");
			exit(1);
		}

	}
	/* infected */
	if(ni==0){
		out->i = NULL;
	} else {
		out->i = (struct host *) calloc(ni, sizeof(struct host));
		if(out->i == NULL){
			fprintf(stderr, "\nNo memory left for creating infected in new population. Exiting.\n");
			exit(1);
		}

	}
	/* recovered */
	if(nr==0){
		out->r = NULL;
	} else {
		out->r = (struct host *) calloc(nr, sizeof(struct host));
		if(out->r == NULL){
			fprintf(stderr, "\nNo memory left for creating recovered in new population. Exiting.\n");
			exit(1);
		}

	}

	return out;
}








/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free host */
void free_host(struct host *in){
	free(in);
}




/* Free population */
void free_pop(struct pop *in){
	free(in->s);
	free(in->i);
	free(in->r);
	free(in);
}









/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Print host content */
void print_host(struct host *in){
	printf("\nhost %d", get_host_id(in));
}




/* Print population content */
void print_pop(struct pop *in){
	printf("\nnb susceptible: %d", get_ns(in));
	printf("\nnb infected: %d", get_ni(in));
	printf("\nnb recovered: %d", get_nr(in));
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
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

	int i;

	/* simulation parameters */
	struct param * par;
	par = (struct param *) calloc(1, sizeof(struct param));
	par->L = 100;
	par->mu = 0.01;
	par->muL = par->mu * par->L;
	par->rng = rng;


	int NREPLI = 1e3;

	struct pathogen ** ppat;

	/* allocate memory */
	ppat = (struct pathogen **) calloc(NREPLI, sizeof(struct pathogen *));
	if(ppat==NULL){
			fprintf(stderr, "\nNo memory left for creating new array of pathogens. Exiting.\n");
			exit(1);
	}
	for(i=1;i<NREPLI;i++){
		ppat[i] = (struct pathogen *) calloc(1, sizeof(struct pathogen));
		if(ppat[i]==NULL){
			fprintf(stderr, "\nNo memory left for expanding the array of pathogens. Exiting.\n");
			exit(1);
		}
	}

	/* initiate array of pathogens */
	ppat[0] = create_pathogen();

	/* replications */
	for(i=0;i<(NREPLI-1);i++){
		replicate(ppat[i],ppat[i+1],par);
	}

	for(i=0;i<NREPLI;i++){
		printf("\npathogen %d",i);
		print_pathogen(ppat[i]);
	}

	/* free memory */
	for(i=0;i<NREPLI;i++) free_pathogen(ppat[i]);
	free(ppat);
	free_param(par);
	free(par);
}
