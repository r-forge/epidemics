/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/


#include "common.h"
#include "param.h"
#include "auxiliary.h"
#include "pathogens.h"





/*
   =================
   === ACCESSORS ===
   =================
*/

/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in){
	return in->snps->n;
}



/* Returns SNP vector */
struct vec_int * get_snps_vec(struct pathogen *in){
	return in->snps;
}


/* Returns SNP integer pointer */
int * get_snps(struct pathogen *in){
	return in->snps->values;
}


/* Returns the age of the pathogen - 0 when created */
int get_age(struct pathogen *in){
	return in->age;
}



/* Returns the population of the pathogen (-1 for inactive pathogen) */
int get_popid(struct pathogen *in){
	return in->popid;
}



/* Returns the population of the pathogen (-1 for inactive pathogen) */
struct pathogen * get_ances(struct pathogen *in){
	return in->ances;
}






/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty pathogen */
struct pathogen * create_pathogen(){
	struct pathogen *out;
	out = (struct pathogen *) calloc(1, sizeof(struct pathogen));
	if(out == NULL){
		fprintf(stderr, "\n[in: pathogen.c->create_pathogen]\nNo memory left for creating initial pathogen. Exiting.\n");
		exit(1);
	}
	out->snps = create_vec_int(0);
	out->age = 0;
	out->popid = 0;
	out->ances = NULL;
	return out;
}








/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free pathogen */
void free_pathogen(struct pathogen *in){
	if(in != NULL){
		free_vec_int(in->snps);
		free(in);
	}
}








/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Copy pathogen */
/*  (memory allocation for in/out made outside the function) */
void copy_pathogen(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, N;

	N=get_nb_snps(in);

	/* out->snps = (int *) calloc(N, sizeof(int)); /\* allocate memory for snps vector*\/ */
	free_vec_int(out->snps);
	out->snps = create_vec_int(N);
	if(get_snps(out) == NULL){
		fprintf(stderr, "\n[in: pathogen.c->copy_pathogen]\nNo memory left for copying pathogen genome. Exiting.\n");
		exit(1);
	}

	for(i=0;i<N;i++){ /* copy snps */
		out->snps->values[i] = get_snps(in)[i];
	}

	out->age = get_age(in);
	out->popid = get_popid(in);
	out->ances = get_ances(in);
}






/* generate a new, unique mutation (i.e., not a reverse mutation) */
int make_unique_mutation(struct pathogen *in, struct param *par){
	int x, i, N=get_nb_snps(in);

	if(N==0) return gsl_rng_uniform_int(par->rng,par->L)+1;
	do{
		i = 0;
		x = gsl_rng_uniform_int(par->rng,par->L)+1; /*generate mutation*/
		/* printf("\nmutation: %d",x); */
		while(i<N && x!=get_snps(in)[i]){ /* check if it exists */
			/* printf("\ncheck against: %d",get_snps(in)[i]); */
			i++;
		}
		/* printf("final value of i: %d N:%d",i,N); */
	} while (i != N);

	return x;
}




/* generate a mutation (possibly an existing one) */
int make_mutation(struct param *par){
	return gsl_rng_uniform_int(par->rng,par->L)+1;
}




/* Print pathogen content */
void print_pathogen(struct pathogen *in){
	int i, N=get_nb_snps(in);
	printf("\nin pop %d, age %d, %d snps:\n", get_popid(in), get_age(in), N);
	if(N>0) {
		for(i=0;i<N;i++) printf("%d ", get_snps(in)[i]);
	}
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome */
/* Create a new pathogen */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, nbmut=gsl_ran_poisson(par->rng, par->muL);

	/* check that output is OK */
	if(out == NULL){
		fprintf(stderr, "\n[in: pathogen.c->replicate]\nTrying to create a new pathogen but pointer is NULL. Exiting.\n");
		exit(1);
	}

	/* allocate memory for new vector */
	free_vec_int(out->snps);
	out->snps = create_vec_int(nbmut);

	/* add new mutations */
	for(i=0;i<nbmut;i++){
		out->snps->values[i] = make_mutation(par);
	}

	out->age = 0;
	out->popid = get_popid(in);
	out->ances = in;
} /*end replicate*/






int isNULL_pathogen(struct pathogen *in){
	if(in==NULL) {
		fprintf(stderr, "\n[in: pathogen.c->isNULL_pathogen]\nPointer to a pathogen is NULL. Exiting.\n");
		exit(1);
	}
	if(get_age(in)<0) return(1);
	return 0;
}








/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
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

/* 	int i; */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->L = 100; */
/* 	par->mu = 0.01; */
/* 	par->muL = par->mu * par->L; */
/* 	par->rng = rng; */


/* 	int NREPLI = 1e3; */

/* 	struct pathogen ** ppat; */

/* 	/\* allocate memory *\/ */
/* 	ppat = (struct pathogen **) calloc(NREPLI, sizeof(struct pathogen *)); */
/* 	if(ppat==NULL){ */
/* 			fprintf(stderr, "\nNo memory left for creating new array of pathogens. Exiting.\n"); */
/* 			exit(1); */
/* 	} */
/* 	for(i=1;i<NREPLI;i++){ */
/* 		ppat[i] = (struct pathogen *) calloc(1, sizeof(struct pathogen)); */
/* 		if(ppat[i]==NULL){ */
/* 			fprintf(stderr, "\nNo memory left for expanding the array of pathogens. Exiting.\n"); */
/* 			exit(1); */
/* 		} */
/* 	} */

/* 	/\* initiate array of pathogens *\/ */
/* 	ppat[0] = create_pathogen(); */

/* 	/\* replications *\/ */
/* 	for(i=0;i<(NREPLI-1);i++){ */
/* 		replicate(ppat[i],ppat[i+1],par); */
/* 	} */

/* 	for(i=0;i<NREPLI;i++){ */
/* 		printf("\npathogen %d",i); */
/* 		print_pathogen(ppat[i]); */
/* 	} */

/* 	/\* free memory *\/ */
/* 	for(i=0;i<NREPLI;i++) free_pathogen(ppat[i]); */
/* 	free(ppat); */
/* 	free_param(par); */

/* 	return 0; */
/* } */



