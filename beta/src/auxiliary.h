/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/


/*
   ==================
   === STRUCTURES ===
   ==================
*/

struct table_int{
	int *items, *times, n;
};



struct vec_int{
	int *values, n;
};


struct distmat_int{
	int *x, n, length;
};



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/



struct distmat_int * create_distmat_int(int n);

struct vec_int * create_vec_int(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/


void free_distmat_int(struct distmat_int *in);

void free_table_int(struct table_int *in);

void free_vec_int(struct vec_int *in);




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
int int_in_vec(int x, int *vec, int vecSize);

int max_int(int *vec, int length);

int min_int(int *vec, int length);


/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/
struct table_int * get_table_int(int *vec, int length);

struct vec_int * sample_int_unif(int N, int I, gsl_rng * rng);

struct vec_int * sample_int_multinom(int N, int I, double * proba, gsl_rng * rng);

void print_table_int(struct table_int *in);

void print_vec_int(struct vec_int *in);

void print_distmat_int(struct distmat_int *in);

