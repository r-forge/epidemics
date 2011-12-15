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

	out = (struct distmat_int *) malloc(sizeof(struct distmat_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	out->x = malloc(length * sizeof(int));
	if(out->x == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_distmat]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	out->length = length;
	out->n = n;

	return out;
}




/* create a vector of integers of size n */
struct vec_int * create_vec_int(int n){
	struct vec_int *out = (struct vec_int *) malloc(sizeof(struct vec_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
		exit(1);
	}

	/* NOTE out->values is not allocated when n=0 */
	if(n>0){
		out->values = (int *) malloc(n * sizeof(int));
		if(out->values == NULL){
			fprintf(stderr, "\n[in: auxiliary.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
			exit(1);
		}
	}

	out->n = n;

	return(out);
}





/* create a vector of integers of size n initialized to zero */
struct vec_int * create_vec_int_zero(int n){
	struct vec_int *out = (struct vec_int *) malloc(sizeof(struct vec_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
		exit(1);
	}

	out->values = (int *) calloc(n, sizeof(int));
	if(out->values == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
		exit(1);
	}

	out->n = n;

	return(out);
}





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_distmat_int(struct distmat_int *in){
	free(in->x);
	free(in);
}


void free_table_int(struct table_int *in){
	free(in->items);
	free(in->times);
	free(in);
}

void free_vec_int(struct vec_int *in){
	if(in->n > 0) free(in->values);
	free(in);
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
/* compute the number of occurence of items in a vect of integers */
struct table_int * get_table_int(int *vec, int length){
	int i, j, *pool, poolsize;
	struct table_int *out = (struct table_int *) malloc(sizeof(struct table_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->get_table_int]\nNo memory left for creating table of integers. Exiting.\n");
		exit(1);
	}

	/* enumerate nb of unique items */
	/* create pool of unique items */
	pool = (int *) malloc(length * sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->get_table_int]\nNo memory left for listing unique integers. Exiting.\n");
		exit(1);
	}

	/* list and count all items */
	poolsize = 0;
	for(i=0;i<length;i++){
		if(int_in_vec(vec[i], pool, poolsize) < 0){
			pool[poolsize++] = vec[i];
		}
	}

	/* copy list of items to output */
	out->items = (int *) malloc(poolsize * sizeof(int));
	if(out->items == NULL){
		fprintf(stderr, "\n[in: auxiliary.c->create_table_int]\nNo memory left for creating table of integers. Exiting.\n");
		exit(1);
	}

	for(i=0;i<poolsize;i++) out->items[i] = pool[i];
	out->n = poolsize;

	/* count number of occurences of each item */
	out->times = (int *) calloc(poolsize, sizeof(int)); /* important to use calloc here */
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





/* sample N times from I items - with replacement, uniform proba */
/* returns a vector of size I with number of times each item was sampled */
/* the sum of all values being N.*/
struct vec_int * sample_int_unif(int N, int I, gsl_rng * rng){
	int i, temp;
	struct vec_int * out = create_vec_int(I);

	/* draw values */
	for(i=0;i<N;i++){
		temp=gsl_rng_uniform_int(rng, I);
		out->values[temp] = out->values[temp] + 1;
	}

	/* free local pointers and return result */
	return out;
}




/* sample N times from I items - with replacement, specified proba */
/* returns a vector of size I with number of times each item was sampled */
/* the sum of all values being N.*/
struct vec_int * sample_int_multinom(int N, int I, double * proba, gsl_rng * rng){
	struct vec_int * out = create_vec_int(I);

	gsl_ran_multinomial (rng, I, N, proba, (unsigned int *) out->values);

	/* free local pointers and return result */
	return out;
}




/* merge K vectors together */
struct vec_int * merge_vec_int(struct vec_int ** in, int nbvec){
	int i, j, newsize=0, count=0;
	/* find size of new vector */
	for(i=0;i<nbvec;i++){
		newsize += in[i]->n;
	}

	/* allocate output and fill it in */
	struct vec_int * out = create_vec_int(newsize);
	for(i=0;i<nbvec;i++){
		for(j=0;j<in[i]->n;j++){
			out->values[count++] = in[i]->values[j];
		}
	}

	return out;
}



/* keep only integers which occur an odd number of times */
struct vec_int * keep_odd_int(struct vec_int *in){
	int i, nbOdd=0, count=0;
	struct vec_int * out;

	/* find number of elements to retain */
	struct table_int * tab = get_table_int(in->values, in->n);
	for(i=0;i<tab->n;i++){
		if(tab->times[i] % 2 > 0) nbOdd++;
	}

	/* fill in the result */
	out = create_vec_int(nbOdd);
	for(i=0;i<tab->n;i++){
		if(tab->times[i] % 2 > 0) out->values[count++] = tab->items[i];
	}

	/* free local alloc and return */
	free_table_int(tab);
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




void print_vec_int(struct vec_int *in){
	int i;
	printf("\nVector of %d values: ", in->n);
	for(i=0;i<in->n;i++) printf("%d ", in->values[i]);
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
/* 	/\* Initialize random number generator *\/ */
/* 	time_t t; */
/* 	t = time(NULL); // time in seconds, used to change the seed of the random generator */
/* 	gsl_rng * rng; */
/* 	const gsl_rng_type *typ; */
/* 	gsl_rng_env_setup(); */
/* 	typ=gsl_rng_default; */
/* 	rng=gsl_rng_alloc(typ); */
/* 	gsl_rng_set(rng,t); // changes the seed of the random generator */


/* 	int  i, vec[10]={1,2,1,4,3,2,2,2,1,5}, n=10; */
/* 	struct table_int *out; */
/* 	struct vec_int *out2, *out3; */

/* 	printf("\ninput: "); */
/* 	for(i=0;i<n;i++) printf("%d\t", vec[i]); */

/* 	out = get_table_int(vec,n); */

/* 	printf("\noutput"); */
/* 	print_table_int(out); */

/* 	printf("\ndrawing 1000 times amongst 4 items, uniform proba\n"); */
/* 	out2 = sample_int_unif(1000, 4, rng); */
/* 	print_vec_int(out2); */

/* 	double proba[4] = {1.0, 2.0, 3.0, 4.0}; */
/* 	printf("\ndrawing 1000 times amongst 4 items, weights 1,2,3,4\n"); */
/* 	out3 = sample_int_multinom(1000, 4 , proba, rng); */
/* 	print_vec_int(out2); */


/* 	/\* /\\* test binomial vs poisson *\\/ *\/ */
/* 	/\* /\\* convergence ok *\\/ *\/ */
/* 	/\* time_t t1,t2; *\/ */
/* 	/\* n=1e7; *\/ */
/* 	/\* double p=1e-8, lambda=p*n; *\/ */
/* 	/\* time(&t1); *\/ */
/* 	/\* for(i=0;i<1e8;i++) gsl_ran_binomial (rng, p, n); *\/ */
/* 	/\* time(&t2); *\/ */

/* 	/\* printf("\n tirage binomial: %d seconds\n ", (int) t2-t1); *\/ */

/* 	/\* time(&t1); *\/ */
/* 	/\* for(i=0;i<1e8;i++) gsl_ran_poisson (rng, lambda); *\/ */
/* 	/\* time(&t2); *\/ */

/* 	/\* printf("\n tirage poisson: %d seconds\n ", (int) t2-t1); *\/ */


/* 	/\* /\\* test binomial vs poisson *\\/ *\/ */
/* 	/\* /\\* convergence not ok *\\/ *\/ */
/* 	/\* n=1e6; *\/ */
/* 	/\* p=1e-4; *\/ */
/* 	/\* lambda=p*n; *\/ */
/* 	/\* time(&t1); *\/ */
/* 	/\* for(i=0;i<1e7;i++) gsl_ran_binomial (rng, p, n); *\/ */
/* 	/\* time(&t2); *\/ */

/* 	/\* printf("\n tirage binomial: %d seconds\n ", (int) t2-t1); *\/ */

/* 	/\* time(&t1); *\/ */
/* 	/\* for(i=0;i<1e7;i++) gsl_ran_poisson (rng, lambda); *\/ */
/* 	/\* time(&t2); *\/ */

/* 	/\* printf("\n tirage poisson: %d seconds\n ", (int) t2-t1); *\/ */


/* 	/\* test vector merging *\/ */
/* 	printf("\nmerging vectors:"); */
/* 	print_vec_int(out2); */
/* 	print_vec_int(out2); */
/* 	print_vec_int(out3); */

/* 	struct vec_int ** in = calloc(3, sizeof(struct vec_int *)); */
/* 	in[0] = out2; */
/* 	in[1] = out2; */
/* 	in[2] = out3; */

/* 	struct vec_int *out4; */
/* 	out4 = merge_vec_int(in, 3); */
/* 	print_vec_int(out4); */


/* 	/\* test retain only items appearing an odd number of times *\/ */
/* 	printf("\nkeeping elements appearing odd nb of times:"); */
/* 	struct vec_int *out5; */
/* 	print_vec_int(out4); */
/* 	out5 = keep_odd_int(out4); */
/* 	print_vec_int(out5); */


/* 	/\* free & return *\/ */
/* 	free_table_int(out); */
/* 	free_vec_int(out2); */
/* 	free_vec_int(out3); */
/* 	free_vec_int(out4); */
/* 	free_vec_int(out5); */

/* 	gsl_rng_free(rng); */
/* 	free(in); */
/* 	return 0; */
/* } */
