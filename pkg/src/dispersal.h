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

/* n: nb of vertices */
/* nbNb: nb of neighbours of each vertice */
/* listNb: list of neighbours for each vertice */
/* weights: weights ~ proba migration */
struct network{
	int n, *nbNb, **listNb;
	double ** weights;
};






/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct network * create_network(struct param *par);






/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free network */
void free_network(struct network *in);






/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

void print_network(struct network *in, bool detail);




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
