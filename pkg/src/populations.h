/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/


/*
# status legend # 
- s: susceptible
- i: infected
- r: removed/recovered
*/


struct metapopulation{
	struct pathogen ** pathogens;
	struct population ** populations;
	int maxnpat, npop, *popid;
};




struct population{
	/* nb of items */
	int nsus, ninf, nrec, ninfcum, orinsus;
};




struct sample{
	struct pathogen ** pathogens;
	/* nb of item; ninf is the length of **pathogens */
	int n;
};



/*
   =================
   === ACCESSORS ===
   =================
*/

struct pathogen ** get_pathogens(struct population *in);


struct population * get_populations(struct metapopulation *in);


int get_maxnpat(struct metapopulation *in);


int * get_popid(struct metapopulation *in);


int get_npop(struct metapopulation *in);


int get_nsus(struct population *in);


int get_ninf(struct population *in);


int get_nrec(struct population *in);


int get_ninfcum(struct population *in);


int get_orinsus(struct population *in);


int get_n(struct sample *in);


int get_total_nrec(struct metapopulation *in);


int get_total_ninfcum(struct metapopulation *in);


int get_total_orinsus(struct metapopulation *in);




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/
struct population * create_metapopulation(int maxnpat, int nini, int npop, int nsus);

struct population * create_population(int ns, int ni, int nr);


/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_metapopulation(struct population *in);

void free_population(struct population *in);

void free_sample(struct sample *in);


/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void print_metapopulation(struct metapopulation *in, bool showGen);

void print_population(struct population *in);

void print_sample(struct sample *in, bool showGen);


/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

/* draw sample from a population */
struct sample * draw_sample(struct population *in, struct param *par);
