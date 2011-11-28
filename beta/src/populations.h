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


struct population{
	int nsus, ninf, nrec, ninfcum, popsize;
	struct pathogen **pathogens;
};




struct metapopulation{
	struct population ** populations;
	int npop, *popsizes;
};




struct ts_groupsizes{
	int *nsus, *ninf, *nrec, *ninfcum, length;
};





/*
   =================
   === ACCESSORS ===
   =================
*/
/* FOR POPULATIONS */
int get_nsus(struct population *in);

int get_ninf(struct population *in);

int get_nrec(struct population *in);

int get_ninfcum(struct population *in);

int get_popsize(struct population *in);

struct pathogen ** get_pathogens(struct population *in);




/* FOR METAPOPULATIONS */
struct population ** get_populations(struct metapopulation *in);

int get_npop(struct metapopulation *in);

int * get_popsizes(struct metapopulation *in);

int get_total_nsus(struct metapopulation *in);

int get_total_ninf(struct metapopulation *in);

int get_total_nrec(struct metapopulation *in);

int get_total_ninfcum(struct metapopulation *in);

int get_total_popsize(struct metapopulation *in);




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/
struct population * create_population(int popsize, int nini);

struct metapopulation * create_metapopulation(struct param *par);

struct ts_groupsizes * create_ts_groupsizes(struct param *par);


/*

   ===================
   === DESTRUCTORS ===
   ===================
*/
void free_population(struct population *in);

void free_metapopulation(struct metapopulation *in);

void free_ts_groupsizes(struct ts_groupsizes *in);



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void print_population(struct population *in, bool showPat);

void print_metapopulation(struct metapopulation *in, bool showPat);




/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/
/* age population */
void age_population(struct population * in, struct param *par);

/* age metapopulation */
void age_metapopulation(struct metapopulation * metapop, struct param * par);

/* keep track of group sizes */
void fill_ts_groupsizes(struct ts_groupsizes *in, struct metapopulation *metapop, int step);
