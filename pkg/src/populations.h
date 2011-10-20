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
	unsigned int nsus, ninf, nrec, ninfcum;
};




/*
   =================
   === ACCESSORS ===
   =================
*/

struct pathogen ** get_pathogens(struct population *in);


unsigned int get_nsus(struct population *in);


unsigned int get_ninf(struct population *in);


unsigned int get_nrec(struct population *in);


unsigned int get_ninfcum(struct population *in);




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/
struct population * create_population(unsigned int ns, unsigned int ni, unsigned int nr);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_population(struct population *in);




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/
void print_population(struct population *in);

