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
	struct pathogen ** pathogens;
	/* nb of item; ninf is the length of **pathogens */
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


int get_nsus(struct population *in);


int get_ninf(struct population *in);


int get_nrec(struct population *in);


int get_ninfcum(struct population *in);


int get_orinsus(struct population *in);


int get_n(struct sample *in);


/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/
struct population * create_population(int ns, int ni, int nr);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_population(struct population *in);

void free_sample(struct sample *in);


/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void print_population(struct population *in, bool showGen);

void print_sample(struct sample *in, bool showGen);


/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

/* draw sample from a population */
struct sample * draw_sample(struct population *in, struct param *par);
