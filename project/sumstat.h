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
	int *snps, length, N;
};


struct allfreq{
	double *freq;
	int length;
};





/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct snplist * create_snplist(int n);

struct allfreq * create_allfreq(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_snplist(struct snplist *in);

void free_allfreq(struct allfreq *in);





/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

int dist_a_b(int *a, int *b, int na, int nb);

struct snplist * list_snps(struct sample *in, struct param *par);


/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_snplist(struct snplist *in);

void print_allfreq(struct allfreq *in);

struct allfreq * get_frequencies(struct sample *in, struct param *par);

double hs(struct sample *in, struct param *par);

double hs_full_genome(struct sample *in, struct param *par);

int nb_snps(struct sample *in, struct param *par);

double mean_nb_snps(struct sample *in);

double var_nb_snps(struct sample *in);

struct distmat_int * pairwise_dist(struct sample *in, struct param *par);

double mean_pairwise_dist(struct sample *in, struct param *par);

double var_pairwise_dist(struct sample *in, struct param *par);

double fst(struct sample *in, struct param *par);
