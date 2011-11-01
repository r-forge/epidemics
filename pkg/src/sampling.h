/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/



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

int get_n(struct sample *in);





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_sample(struct sample *in);




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

/* draw sample from a population */
struct sample * draw_sample(struct metapopulation *in, struct param *par);
