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

struct snplist{
	int *snps, length;
};





/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct snplist * create_snplist(int n);





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_snplist(struct snplist *in);






/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
int int_in_vec(int x, int *vec, int vecSize);


struct snplist * list_snps(struct sample *in, struct param *par);

/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_snplist(struct snplist *in);
