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
	struct vec_int *snps;
	struct pathogen *ances;
	int age;
};



struct lineage{
	struct pathogen ** pathogens;
	int n;
};





/*
   =================
   === ACCESSORS ===
   =================
*/

/* Returns the number of mutated SNPs, i.e. length of in->snps array */
int get_nb_snps(struct pathogen *in);

/* Returns SNP vector */
struct vec_int * get_snps_vec(struct pathogen *in);

/* Returns SNP integer pointer */
int * get_snps(struct pathogen *in);


/* Returns the age of the pathogen - 0 when created */
int get_age(struct pathogen *in);


/* Returns the ancestor of the pathogen */
struct pathogen * get_ances(struct pathogen *in);




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty pathogen */
struct pathogen * create_pathogen();

/* Create empty lineage */
struct lineage * create_lineage(int n);





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free pathogen */
void free_pathogen(struct pathogen *in);


/* Free lineage */
/* Note: free only pointers to pathogens, does not free pathogens themselves. */
void free_lineage(struct lineage *in);





/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Copy pathogen */
/*  (memory allocation for in/out made outside the function) */
void copy_pathogen(struct pathogen *in, struct pathogen *out, struct param *par);



/* generate a mutation (possibly an existing one) */
int make_mutation(struct param *par);


/* Print pathogen content */
void print_pathogen(struct pathogen *in);

/* Get the lineage of a pathogen */
struct lineage * get_lineage(struct pathogen *in);

/* Reconstruct genome of an isolate */
struct pathogen * reconstruct_genome(struct pathogen *in);







/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Function replicating a genome, with mutations and back-mutations */
void replicate(struct pathogen *in, struct pathogen *out, struct param *par);


/* TEST IF PATHOGEN IS NULL OR INACTIVATED */
bool isNULL_pathogen(struct pathogen *in);


/* TEST IF PATHOGEN IS INFECTIOUS */
bool is_infectious(struct pathogen *in, struct param *par);
