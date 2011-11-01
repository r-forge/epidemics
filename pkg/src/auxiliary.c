/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "auxiliary.h"




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* create distmat_int between n objects */
struct distmat_int * create_distmat_int(int n){
	struct distmat_int *out;
	int length=n*(n-1)/2;

	out = calloc(1, sizeof(struct distmat_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	out->x = calloc(length, sizeof(int));
	if(out->x == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	out->length = length;
	out->n = n;

	return out;
}




/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_distmat_int(struct distmat_int *in){
	if(in->x != NULL) free(in->x);
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
	while(i<vecSize && x!=vec[i]) i++; /* note: condition needs to be in this order */
	if(i==vecSize || vecSize<1) return -1; /* -1 will mean: no match*/
	return i;
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void print_distmat_int(struct distmat_int *in){
	int i,j, counter=0, N=in->n;

	printf("\npairwise distances between %d individuals:\n", in->n);
	for(i=0;i<N;i++) printf("\t'%d'",i);
	for(i=0;i<N-1;i++){
		printf("\n'%d'\t", i);
		for(j=0;j<i+1;j++){
			printf("- \t");
		}
		for(j=i+1;j<N;j++){
			printf("%d \t", in->x[counter++]);
		}
	}

	printf("\n");
}









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
