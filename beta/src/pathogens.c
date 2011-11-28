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
	out->ances = NULL;
	return out;
}





/* Create empty lineage */
struct lineage * create_lineage(int n){
	struct lineage *out;
	out = (struct lineage *) calloc(1, sizeof(struct lineage));
	if(out == NULL){
		fprintf(stderr, "\n[in: lineage.c->create_lineage]\nNo memory left for creating initial lineage. Exiting.\n");
		exit(1);
	}
	out->pathogens = (struct pathogen **) calloc(n, sizeof(struct pathogen *));
	out->n = n;
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



/* Free lineage */
/* Note: free only pointers to pathogens, does not free pathogens themselves. */
void free_lineage(struct lineage *in){
	if(in != NULL){
		free(in->pathogens);
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
	if(out->snps != NULL) free_vec_int(out->snps); /* erase previous content */
	out->snps = create_vec_int(N);
	if(N>0 && get_snps(out) == NULL){
		fprintf(stderr, "\n[in: pathogen.c->copy_pathogen]\nNo memory left for copying pathogen genome. Exiting.\n");
		exit(1);
	}

	for(i=0;i<N;i++){ /* copy snps */
		out->snps->values[i] = get_snps(in)[i];
	}

	out->age = get_age(in);
	out->ances = get_ances(in);
}





/* generate a mutation (possibly an existing one) */
int make_mutation(struct param *par){
	return gsl_rng_uniform_int(par->rng,par->L)+1;
}




/* Print pathogen content */
void print_pathogen(struct pathogen *in){
	int i, N=get_nb_snps(in);
	printf("\n age %d, %d snps:\n", get_age(in),N);
	if(N>0) {
		for(i=0;i<N;i++) printf("%d ", get_snps(in)[i]);
	}
}






/* retrieve the lineage of a pathogen */
struct lineage * get_lineage(struct pathogen *in){
	int i, lineagesize=0;
	struct pathogen *curIsolate;
	struct lineage *out;

	/* find out lineage size */
	curIsolate = in;
	while(curIsolate != NULL){
		lineagesize++;
		curIsolate = get_ances(curIsolate); /* browse ancestry backward */
	}

	/* allocate and fill in output */
	out = create_lineage(lineagesize);
	curIsolate = in;
	for(i=0;i<lineagesize;i++){
		out->pathogens[i] = curIsolate;
		curIsolate = get_ances(curIsolate);
	}

	return out;
}






/* RECONSTRUCT GENOME OF AN ISOLATE */
struct pathogen * reconstruct_genome(struct pathogen *in){
	int i;
	struct lineage *line = get_lineage(in);
	struct pathogen *out;
	struct vec_int ** listSnpVec, *temp, *genome;


	/* get all snps in the lineage */
	listSnpVec = (struct vec_int **) calloc(line->n, sizeof(struct vec_int *));

	for(i=0;i<line->n;i++){
		listSnpVec[i] = get_snps_vec(line->pathogens[i]);
	}

	/* merge snps */
	temp = merge_vec_int(listSnpVec, line->n);

	/* remove reverse mutations */
	genome = keep_odd_int(temp);

	/* create output and fill it in */
	out = create_pathogen();
	free_vec_int(out->snps);
	out->snps = genome;
	out->age = in->age;
	out->ances = in->ances;

	/* free temporary allocation & return */
	free_lineage(line);
	free(listSnpVec);
	free_vec_int(temp);
	return out;
}





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome */
/* Assignment is done outside the function */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par){
	int i, nbmut=gsl_ran_poisson(par->rng, par->muL);


	/* check that output is OK */
	if(out == NULL){
		fprintf(stderr, "\n[in: pathogen.c->replicate]\nTrying to create a new pathogen but pointer is NULL. Exiting.\n");
		exit(1);
	}

	/* allocate memory for new vector */
	if(out->snps != NULL) free_vec_int(out->snps);
	out->snps = create_vec_int(nbmut);

	/* add new mutations */
	for(i=0;i<nbmut;i++){
		out->snps->values[i] = make_mutation(par);
	}

	out->age = 0;
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


/* gcc line:

   gcc -o pathogens param.c auxiliary.c pathogens.c -Wall -O0 -lgsl -lgslcblas
   valgrind --leak-check=yes pathogens
*/

/* int main(){ */
/* 	/\* Initialize random number generator *\/ */
/* 	time_t t; */
/* 	t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
/* 	gsl_rng * rng; */
/* 	const gsl_rng_type *typ; */
/* 	gsl_rng_env_setup(); */
/* 	typ=gsl_rng_default; */
/* 	rng=gsl_rng_alloc(typ); */
/* 	gsl_rng_set(rng,t); /\* changes the seed of the random generator *\/ */

/* 	int i; */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->L = 10; */
/* 	par->mu = 0.1; */
/* 	par->muL = par->mu * par->L; */
/* 	par->rng = rng; */


/* 	int NREPLI = 100; */

/* 	struct pathogen ** ppat; */

/* 	/\* allocate memory *\/ */
/* 	ppat = (struct pathogen **) calloc(NREPLI, sizeof(struct pathogen *)); */
/* 	if(ppat==NULL){ */
/* 			fprintf(stderr, "\nNo memory left for creating new array of pathogens. Exiting.\n"); */
/* 			exit(1); */
/* 	} */

/* 	/\* INITIATE ARRAY OF PATHOGENS *\/ */
/* 	for(i=0;i<NREPLI;i++) ppat[i] = create_pathogen(); */

/* 	/\* REPLICATIONS *\/ */
/* 	for(i=0;i<(NREPLI-1);i++){ */
/* 		printf("\n replication %d", i); */
/* 		replicate(ppat[i],ppat[i+1],par); */
/* 	} */

/* 	for(i=0;i<NREPLI;i++){ */
/* 		printf("\npathogen %d",i); */
/* 		print_pathogen(ppat[i]); */
/* 	} */


/* 	/\* TEST LINEAGES *\/ */
/* 	struct lineage * myline = get_lineage(ppat[10]); */
/* 	printf("\n\nChosen pathogen has a lineage of length %d", myline->n); */
/* 	for(i=0;i<myline->n;i++) print_pathogen(myline->pathogens[i]); */

/* 	/\* TEST RECONSTRUCTION OF THE GENOME *\/ */
/* 	struct pathogen * mypat; */
/* 	for(i=0;i<NREPLI;i++){ */
/* 		mypat = reconstruct_genome(ppat[i]); */
/* 		printf("\n\n RECONSTRUCTED PATHOGEN %d: ", i); */
/* 		print_pathogen(mypat); */
/* 		free_pathogen(mypat); */
/* 	} */

/* 	/\* free memory *\/ */
/* 	for(i=0;i<NREPLI;i++) free_pathogen(ppat[i]); */
/* 	free(ppat); */
/* 	free_param(par); */
/* 	free_lineage(myline); */

/* 	return 0; */
/* } */



