#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h> /* GSL random number generator & distributions */
#include <gsl/gsl_randist.h> /* GSL random number generator & distributions */


#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

typedef short bool;



/*
   =============================
   === STRUCTURES DEFINITION ===
   =============================
*/

/* The structure 'pathogen' stores a vector of mutated alleles.
   Each integer indicates the position of a mutated allele.
   The wild genotye is an empty vector.
   - 'snps' is an array of unsigned integers
   - 'length' is the length if this array
   - 'host' is an integer identifying the host
 */
struct pathogen{
	unsigned int *snps;
	int length;
	long long unsigned int host;
};


/* L: length of the genomes */
/* mu: mutation rate per site and generation */
struct param{
	int L;
	double mu, muL;
	gsl_rng * rng;
	long long unsigned int lasthost;
};


/*
   =================================
   === LOCAL AUXILIARY FUNCTIONS ===
   =================================
*/




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* Update once sorted out*/







/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/

