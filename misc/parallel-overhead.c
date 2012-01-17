#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */


void waste_time(int x){
	int i=0, a=0;
	for(i=0;i<x;i++){
		//rand();
		a++;
		a--;
		a++;
	}
}

void waste_time_gsl(int x, gsl_rng *rng){
	int i=0, a=0;
	for(i=0;i<x;i++){
		a = gsl_rng_uniform_int(rng,100);

	}
}



/*

gcc -o para3 parallel-overhead.c  -lgsl -lgslcblas -fopenmp -O3


*/

int main(){
	const int I=1000, J=100000, WASTE_ARG=10;
	int i,j,count=0;
	time_t time1, time2;

	/* RANDOM NUMBER GENERATION USING GSL */
	/* SETUP */
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


/* 	/\* SERIAL CODE *\/ */
	printf("\nProcessing array, serial code");
	time(&time1);
	for(i=0;i<I;i++){
		for(j=0;j<J;j++){
			waste_time_gsl(WASTE_ARG,rng);
		}
	}
	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));


	/* PARALLEL CODE */
	printf("\nProcessing array, parallel code");
	time(&time1);
//#pragma omp parallel for private(j)
	for(i=0;i<I;i++){
#pragma omp parallel for
		for(j=0;j<J;j++){
			waste_time_gsl(WASTE_ARG,rng);
		}
	}

	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));





/* 	/\* SANITY CHECK *\/ */
/* 	count = I*J*WASTE_ARG; */

/* 	/\* SERIAL CODE *\/ */
/* 	printf("\n\n == Sanity check =="); */

/* 	printf("\nProcessing array, serial code"); */
/* 	time(&time1); */
/* 	for(i=0;i<count;i++){ */
/* 		//waste_time(1); */
/* 		waste_time_gsl(1, rng); */
/* 	} */
/* 	time(&time2); */

/* 	printf("\nthe operation took %d seconds\n", (int) (time2-time1)); */


/* 	/\* PARALLEL CODE *\/ */
/* 	printf("\nProcessing array, parallel code"); */
/* 	time(&time1); */
/* #pragma omp parallel for */
/* 	for(i=0;i<count;i++){ */
/* 		//waste_time(1); */
/* 		waste_time_gsl(1, rng); */
/* 	} */

/* 	time(&time2); */

/* 	printf("\nthe operation took %d seconds\n", (int) (time2-time1)); */

	return 0;

}
