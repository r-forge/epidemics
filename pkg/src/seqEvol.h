#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h> /* GSL random number generator & distributions */
#include <gsl/gsl_randist.h> /* GSL random number generator & distributions */

#include "hosts.h"

#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

typedef short bool;



/*
   =============================
   === STRUCTURES DEFINITION ===
   =============================
*/

/* The structure 'pathogen' stores a vector of mutated alleles.
   Each integer indicates the position of a mutated allele.
   The wild genotye is an empty vector.
   - 'snps' is an array of unsigned integers
   - 'length' is the length if this array
   - 'host' is an integer identifying the host
 */
struct pathogen{
	unsigned int *snps;
	int length;
	struct host *host;
};




/* L: length of the genomes */
/* mu: mutation rate per site and generation */
struct param{
	int L;
	double mu, muL;
	gsl_rng * rng;
};








/*
   =================
   === ACCESSORS ===
   =================
*/

/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in);


/* Returns the ID of the host, i.e. in->host array */
long long unsigned int get_host(struct pathogen *in);


/* Returns SNP vector */
unsigned int * get_snps(struct pathogen *in);








/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty pathogen */
struct pathogen * create_pathogen();








/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free pathogen */
void free_pathogen(struct pathogen *in);


/* Free param */
void free_param(struct param *in);








/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Copy pathogen */
/*  (memory allocation for in/out made outside the function) */
void copy_pathogen(struct pathogen *in, struct pathogen *out, struct param *par);


/* generate a new, unique mutation (i.e., not a reverse mutation) */
int make_unique_mutation(struct pathogen *in, struct param *par);


/* Print pathogen content */
void print_pathogen(struct pathogen *in);








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome, with mutations and back-mutations */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par);
