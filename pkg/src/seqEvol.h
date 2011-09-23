#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

typedef short bool;

long long int LAST_HOST=0;

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
	long long int host;
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
struct pathogen * create_new_pathogen(unsigned int *snps, int *length, int *host);
void create_initial_pathogen(struct pathogen *out, unsigned int *snps, int *length, int *host);
int get_nb_snps(struct pathogen *in);
long long int get_host(struct pathogen *in);







/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/

