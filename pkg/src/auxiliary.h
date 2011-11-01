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



struct distmat_int{
	int *x, n, length;
};



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/



struct distmat_int * create_distmat_int(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/


void free_distmat_int(struct distmat_int *in);




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
int int_in_vec(int x, int *vec, int vecSize);


/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_distmat_int(struct distmat_int *in);

