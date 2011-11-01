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

/* create empty distmat_int between n objects */
struct distmat_int * create_distmat_int(int n){
	struct distmat_int *out;
	int length=n*(n-1)/2;

	out = calloc(1, sizeof(struct distmat_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	out->x = calloc(length, sizeof(int));
	if(out->x == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
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



void free_table_int(struct table_int *in){
	if(in->items != NULL) free(in->items);
	if(in->times != NULL) free(in->times);
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

/* find max value in a vector of integers */
int max_int(int *vec, int length){
	int i, out=vec[0];
	for(i=0;i<length;i++) if(out<vec[i]) out=vec[i];
	return out;
}


/* find min value in a vector of integers */
int min_int(int *vec, int length){
	int i, out=vec[0];
	for(i=0;i<length;i++) if(out>vec[i]) out=vec[i];
	return out;
}





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

struct table_int * get_table_int(int *vec, int length){
	int i, j, nbitems=0, *pool, poolsize;
	struct table_int *out = calloc(1, sizeof(struct table_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_table_int]\nNo memory left for creating table of integers. Exiting.\n");
		exit(1);
	}

	/* enumerate nb of unique items */
	/* create pool of unique items */
	pool = calloc(length, sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->get_table_int]\nNo memory left for listing unique integers. Exiting.\n");
		exit(1);
	}

	/* list and count all SNPs */
	poolsize = 0;
	for(i=0;i<length;i++){
		if(int_in_vec(vec[i], pool, poolsize) < 0){
			pool[poolsize++] = vec[i];
		}
	}

	/* copy list of items to output */
	out->items = calloc(poolsize, sizeof(int));
	if(out->items == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_table_int]\nNo memory left for creating table of integers. Exiting.\n");
		exit(1);
	}

	for(i=0;i<poolsize;i++) out->items[i] = pool[i];
	out->n = poolsize;

	/* count number of occurences of each item */
	out->times = calloc(poolsize, sizeof(int));
	if(out->times == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_table_int]\nNo memory left for creating table of integers. Exiting.\n");
		exit(1);
	}

	for(i=0;i<poolsize;i++){
		for(j=0;j<length;j++){
			if(vec[j]==out->items[i]) out->times[i] = out->times[i]+1;
		}
	}

	/* free local pointers and return result */
	free(pool);
	return out;
}




void print_table_int(struct table_int *in){
	int i;
	printf("\nItems: ");
	for(i=0;i<in->n;i++) printf("%d\t", in->items[i]);
	printf("\nTimes: ");
	for(i=0;i<in->n;i++) printf("%d\t", in->times[i]);
	printf("\n");
}



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
/* 	int  i, vec[10]={1,2,1,4,3,2,2,2,1,5}, n=10; */
/* 	struct table_int *out; */

/* 	printf("\ninput: "); */
/* 	for(i=0;i<n;i++) printf("%d\t", vec[i]); */

/* 	out = get_table_int(vec,n); */

/* 	printf("\noutput"); */
/* 	print_table_int(out); */

/* 	free_table_int(out); */
/* 	return 0; */
/* } */
