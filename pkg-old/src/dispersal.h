/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/




/*
   =======================
   === DATA STRUCTURES ===
   =======================
*/

/* mat is a list of arrays, where the i^th item gives dispersal probabilities to other populations */

struct dispmat{
	double ** mat;
	int n;
};




/*
   =================
   === ACCESSORS ===
   =================
*/
double ** get_mat(struct dispmat * in);






/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct dispmat * create_dispmat(struct param *par);






/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free dispmat */
void free_dispmat(struct dispmat *in);






/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

void print_dispmat(struct dispmat *in);




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* return the id of a new population */
int disperse(struct pathogen * pathogen, struct dispmat *disp, struct param *par);
