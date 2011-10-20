/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/


#include "common.h"
#include "param.h"
#include "pathogens.h"





/*
   =================
   === ACCESSORS ===
   =================
*/

/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in){
	return in->length;
}




/* Returns SNP vector */
unsigned int * get_snps(struct pathogen *in){
	return in->snps;
}



/* Returns the age of the pathogen - 0 when created */
int get_age(struct pathogen *in){
	return in->age;
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
		fprintf(stderr, "\nNo memory left for creating initial pathogen. Exiting.\n");
		exit(1);
	}
	out->snps = NULL;
	out->length = 0;
	out->age = 0;
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
		free(in->snps);
		/*free(in->host);*/
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

	out->snps = (unsigned int *) calloc(N, sizeof(unsigned int)); /* allocate memory for snps vector*/
	if(get_snps(out) == NULL){
		fprintf(stderr, "\nNo memory left for copying pathogen genome. Exiting.\n");
		exit(1);
	}

	for(i=0;i<N;i++){ /* copy snps */
		(out->snps)[i] = get_snps(in)[i];
	}
	out->length = N;
	out->age = get_age(in);
	/*out->host = get_host(in);*/
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




/* Print pathogen content */
void print_pathogen(struct pathogen *in){
	int i, N=get_nb_snps(in);
	printf("\nage: %d \n%d snps: ", get_age(in), N);
	if(N>0) {
		for(i=0;i<N;i++) printf("%d ", get_snps(in)[i]);
	}
	/* printf("\nhost: %llu \n", get_host(in)); */
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome, with mutations and back-mutations */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, nbmut=0, nbbackmut=0, newsize, N;
	double p;

	nbmut = gsl_ran_poisson(par->rng, par->muL);

	/* determine the number of reverse mutations */
	if(nbmut>0 && get_nb_snps(in)>0){
		p = ((double) get_nb_snps(in)) / ((double) par->L);
		nbbackmut =  gsl_ran_binomial(par->rng,p,nbmut);
		if(nbbackmut > get_nb_snps(in)) nbbackmut = get_nb_snps(in); /* can revert more than to wild genotype */
	}

	nbmut -= nbbackmut; /* remove back mutations from new mutations */
	newsize = get_nb_snps(in) + nbmut - nbbackmut;

	/* check that new size is not negative */
	if(newsize < 1){ /* go back to wild type */
		out->snps = NULL;
		out->length = 0;
	} else { /* new genotype */
		/* reallocate memory for new vector */
		out->snps = (unsigned int *) calloc(newsize, sizeof(unsigned int));
		if(get_snps(out) == NULL){
			fprintf(stderr, "\nNo memory left for replicating pathogen genome. Exiting.\n");
			exit(1);
		}

		/* inherit parental snps vector (except back mutations) */
		N = get_nb_snps(in)-nbbackmut; /* indices<N, ancestral snps; indices>=N, new snps */
		for(i=0;i<N;i++){ /* copy snps */
			out->snps[i] = get_snps(in)[i];
		}

		out->length=N;

		/* add new mutations */
		for(i=0;i<nbmut;i++){
			(out->snps)[N+i] = make_unique_mutation(out, par);
			out->length += 1;
		}
	} /* the genotype has been handled at this point */

	out->age = 0;

} /*end replicate*/








/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


/* void main(){ */
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
/* 	free(par); */
/* } */








/* TESTING in R */

/*
## test raw conversion
.C("myCfunction", arg1, arg2, ..., PACKAGE="epidemics")


*/

