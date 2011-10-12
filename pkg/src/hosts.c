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
/* Calls to GNU Scientific Library */
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */

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

