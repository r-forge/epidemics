#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */

#define FAC 2



void do_smthg(gsl_rng * rng, int ntimes){
	int i;
	for(i=0;i<ntimes*FAC;i++){
		gsl_rng_uniform_int(rng,100);
	}
}



/*

  gcc -o para-pb-age parallel-AgeLike.c  -lgsl -lgslcblas -fopenmp -O3 && para-pb-age

*/
int main(){
	const int nstep = 20, nreps=200;
	int i, k, N;
	time_t time1, time2;

	/* GSL RNG INITIALIZATION */
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	/* =============================== */
	/* TEST USING FIXED-SIZE LOOPS */
	/* =============================== */
	N = 100000;
	printf("\n == Test using fixed-size loops ==\n ");

	/* SERIAL CODE */
	printf("\nSerial code");
	time(&time1);
	for(i=0;i<nreps;i++){
		for(k=0;k<nstep;k++){
			do_smthg(rng,N);
		}
	}
	time(&time2);
	printf("\nthe operation took %d seconds\n", (int) (time2-time1));

	/* PARALLEL CODE */
	printf("\nParallel code");
	time(&time1);
	for(i=0;i<nreps;i++){
#pragma omp parallel for schedule(static,1)
		for(k=0;k<nstep;k++){
			do_smthg(rng,N);
		}
	}
	time(&time2);
	printf("\nthe operation took %d seconds\n", (int) (time2-time1));





	/* =============================== */
	/* TEST USING VARIABLE-SIZE LOOPS */
	/* =============================== */
	printf("\n == Test using variable-size loops ==\n ");

	/* SERIAL CODE */
	printf("\nSerial code");
	time(&time1);
	for(i=0;i<nreps;i++){
		for(k=0;k<nstep;k++){
			N = (int) (pow(2.0,k));
			do_smthg(rng,N);
		}
	}
	time(&time2);
	printf("\nthe operation took %d seconds\n", (int) (time2-time1));

	/* PARALLEL CODE */
	printf("\nParallel code");
	time(&time1);
	for(i=0;i<nreps;i++){
#pragma omp parallel for private(N)
		for(k=0;k<nstep;k++){
			N = (int) (pow(2.0,k));
			do_smthg(rng,N);
		}
	}
	time(&time2);
	printf("\nthe operation took %d seconds\n", (int) (time2-time1));

	return 0;
}
