/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/



/*
   =============================
   === STRUCTURES DEFINITION ===
   =============================
*/

/* The structure 'pathogen' stores a vector of mutated alleles.
   Each integer indicates the position of a mutated allele.
   The wild genotye is an empty vector.
   - 'snps' is an array of integers
   - 'length' is the length if this array
   - 'age' gives the age of the pathogen
   - 'popid' gives the index of the population in which the pathogen is
   - 'ances' is a pointer to the ancestor
 */
struct pathogen{
	int *snps;
	int length, age, popid;
	struct pathogen *ances;
};








/*
   =================
   === ACCESSORS ===
   =================
*/

/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in);


/* Returns SNP vector */
int * get_snps(struct pathogen *in);


/* Returns the age of the pathogen - 0 when created */
int get_age(struct pathogen *in);


/* Returns the population of the pathogen */
/* (-1 for inactive pathogen) */
int get_popid(struct pathogen *in);


/* Returns the ancestor of the pathogen */
struct pathogen * get_ances(struct pathogen *in);




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


/* generate a mutation (possibly an existing one) */
int make_mutation(struct pathogen *in, struct param *par);


/* Print pathogen content */
void print_pathogen(struct pathogen *in);








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome, with mutations and back-mutations */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par);


/* return 1 if pathogen is neutralized, i.e. aged -1; issues error if pointer is NULL */
int isNULL_pathogen(struct pathogen *in);
