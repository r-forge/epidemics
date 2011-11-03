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
	int maxnpat, npop;
};




struct population{
	/* nb of items */
	int nsus, ninf, nrec, ninfcum, orinsus;
};



struct ts_groupsizes{
	int *nsus, *ninf, *nrec, *ninfcum, length;
};





/*
   =================
   === ACCESSORS ===
   =================
*/

struct pathogen ** get_pathogens(struct metapopulation *in);


struct population ** get_populations(struct metapopulation *in);


int get_maxnpat(struct metapopulation *in);


int get_npop(struct metapopulation *in);


int get_nsus(struct population *in);


int get_ninf(struct population *in);


int get_nrec(struct population *in);


int get_ninfcum(struct population *in);


int get_orinsus(struct population *in);


int get_popsize(struct population *in);


int get_total_nsus(struct metapopulation *in);


int get_total_ninf(struct metapopulation *in);


int get_total_nrec(struct metapopulation *in);


int get_total_ninfcum(struct metapopulation *in);


int get_total_orinsus(struct metapopulation *in);




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/
struct metapopulation * create_metapopulation(struct param *par);

struct population * create_population(int ns, int ni, int nr);

struct ts_groupsizes * create_ts_groupsizes(int nsteps);


/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_metapopulation(struct metapopulation *in);

void free_population(struct population *in);

void free_ts_groupsizes(struct ts_groupsizes *in);


/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void print_metapopulation(struct metapopulation *in, bool showGen);

void print_population(struct population *in);



/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/


/* age metapopulation */
void age_metapopulation(struct metapopulation * metapop, struct param * par);

/* keep track of group sizes */
void fill_ts_groupsizes(struct ts_groupsizes *in, struct metapopulation *metapop, int step);
