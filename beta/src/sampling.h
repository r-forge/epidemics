/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/



struct sample{
	struct pathogen ** pathogens;
	int n, *popid;
};



/*
   =================
   === ACCESSORS ===
   =================
*/

int get_n(struct sample *in);

int get_npop_samp(struct sample *in);


/*
   ===================
   === CONTRUCTORS ===
   ===================
*/

struct sample * create_sample(int n);

struct sample ** create_list_sample(struct param * par);


/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_sample(struct sample *in);

void free_list_sample(struct sample **in, struct param *par);


/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

void print_sample(struct sample *in, bool showGen);



/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

/* draw sample from a metapopulation */
struct sample * draw_sample(struct metapopulation *in, int n, struct param *par);

/* draw sample from a population */
struct sample * draw_sample_onepop(struct population *in, int n, struct param *par);

/* merge several samples together */
struct sample *merge_samples(struct sample **in, int n, struct param *par);



/* translate sampling dates into simulation timestep */
void translate_dates(struct param *par);

/* slit data of a sample by population */
struct sample ** seppop(struct sample *in, struct param *par);
