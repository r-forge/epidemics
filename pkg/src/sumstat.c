/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "param.h"
#include "pathogens.h"
#include "populations.h"
#include "sumstat.h"



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct snplist * create_snplist(int n){
	struct snplist *out;
	out = calloc(1, sizeof(struct snplist));
	if(out == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_snplist]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->snps = calloc(n, sizeof(int));
	if(out->snps == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_snplist]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->length = n;

	return out;
}





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_snplist(struct snplist *in){
	if(in->snps != NULL) free(in->snps);
	if(in != NULL) free(in);
}





/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* check if an integer is in a vector of integers, and return the matching position */
int int_in_vec(int x, int *vec, int vecSize){
	int i=0;
	while(x!=vec[i] && i<vecSize) i++;
	if(i==vecSize || vecSize<1) return -1; /* -1 will mean: no match*/
	return i;
}



/* count and list number of snps in a sample */
struct snplist * list_snps(struct sample *in, struct param *par){
	int i=0, j=0, N=get_n(in), *pool, poolsize, curNbSnps;
	struct snplist *out;

	/* create pool of snps */
	pool = calloc(par->L, sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: sumstat.c->list_snps]\nNo memory left for creating pool of SNPs. Exiting.\n");
		exit(1);
	}

	/* list and count all SNPs */
	poolsize = 0;
	for(i=0;i<N;i++){
		curNbSnps = get_nb_snps(in->pathogens[i]);
		for(j=0;j<curNbSnps;j++){
			if(int_in_vec(get_snps(in->pathogens[i])[j], pool, poolsize) < 0){
				pool[poolsize++] = get_snps(in->pathogens[i])[j];
			}
		}
	}

	/* make output */
	out = create_snplist(poolsize);

	for(i=0;i<poolsize;i++) out->snps[i] = pool[i];

	/* free local pointers and return result */
	free(pool);
	return out;
}




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void print_snplist(struct snplist *in){
	int i;
	printf("\nlist of %d SNPs:\n", in->length);
	for(i=0;i<in->length;i++) printf("%d ",in->snps[i]);
	printf("\n");
}


/* double hs(struct sample *in, struct param *par){ */
/* 	int i; */
/* 	double out; */

/* 	return out; */
/* 	/\* *\/ */
/* } */




/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


/* int main(){ */
/* 	int  i, vec[5]={1,2,3,4,5}; */

/* 	for(i=0;i<10;i++) printf("\ni=%d, result:%d", i, int_in_vec(i,vec,5)); */
	
/* 	return 0; */
/* } */
