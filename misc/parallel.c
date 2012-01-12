#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */

int main(int argc, char *argv[]){
/* #pragma omp parallel */
/* 	{ */
/* 		printf("Hello, world.\n"); */
/* 	} */

/* 	const int N = 100000; */
/* 	int i, a[N], th_id, nthreads; */

/* #pragma omp parallel private(i) */
/* 	{ */
/* 		printf("\nthread number: %d \n", omp_get_thread_num()); */
/* 	} */

/* 	omp_set_num_threads(5); */

/* #pragma omp parallel for */
/* 	for (i = 0; i < N; i++){ */
/* 		a[i] = 2 * i; */
/* 	} */

/* #pragma omp parallel private(th_id) */
/* 	{ */
/* 		th_id = omp_get_thread_num(); */
/* 		printf("Hello World from thread %d\n", th_id); */
/* #pragma omp barrier */
/* 		if ( th_id == 0 ) { */
/* 			nthreads = omp_get_num_threads(); */
/* 			printf("There are %d threads\n",nthreads); */
/* 		} */
/* 	} */
	
/* 	printf("\ntest conditional compiling (yes):\n"); */
/* #define USE_OMP 1 */
/* #if USE_OMP */
/* #pragma omp parallel */
/* #endif */
/* 	{ */
/* 		printf("Hello from thread %d\n", omp_get_thread_num()); */
/* 	} */
	
/* 	printf("\ntest conditional compiling (no):\n"); */
/* #ifdef USE_OMP */
/* #undef USE_OMP */
/* #endif */
/* #define USE_OMP 0 */
/* #if USE_OMP */
/* #pragma omp parallel */
/* #endif */
/* 	{ */
/* 		printf("Hello from thread %d\n", omp_get_thread_num()); */
/* 	} */



/* 	printf("\ntest for loop:\n"); */

/* /\* #undef USE_OMP *\/ */
/* /\* #define USE_OMP 1 *\/ */
/* /\* #if USE_OMP *\/ */
/* /\* #pragma omp for *\/ */
/* /\* #endif *\/ */
/* #pragma omp parallel for */
/* 	for(i=0;i<10;i++) */
/* 		printf("Hello from thread %d\n", omp_get_thread_num()); */





/* ================================================= */
/* ================================================= */
/* ================================================= */

/* /\* TRY PARALLEL WITH COMMON OUTPUT ARRAY *\/ */
/* 	#define VECSIZE 100 */
/* 	int i, j, k, temp, counter, *vec1, *vec2, I=2, K=2; */
/* 	vec1 = (int*) calloc(VECSIZE, sizeof(int)); */
/* 	vec2 = (int*) calloc(VECSIZE, sizeof(int)); */

/* 	time_t time1,time2; */

/* 	printf("\n== Serial for loop ==\n"); */
/* 	time(&time1); */

/* 	counter = 0; */
/* 	for(k=0;k<K;k++){ */
/* 		for(i=0;i<I;i++){ */
/* 			vec1[counter] = counter; */
/* 			counter++; */
/* 			j=0; */
/* 			while(j<100e6) j++; */
/* 		} */
/* 	} */

/* 	time(&time2); */
/* 	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1)); */


/* 	printf("\n== Parallel for loop ==\n"); */
/* 	time(&time1); */

/* 	counter = 0; */
/* //#pragma omp parallel for private(i,temp,j) */
/* 	for(k=0;k<K;k++){ */
/* #pragma omp parallel for private(temp,j) */
/* 		for(i=0;i<I;i++){ */
/* 			temp = (k*I) + i; */
/* 			vec2[temp] = temp; */
/* 			j=0; */
/* 			while(j<100e6) j++; */
/* 		} */
/* 	} */

/* 	time(&time2); */
/* 	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1)); */


/* 	printf("\nArray (serial):\n"); */
/* 	for(i=0;i<(I*K);i++){ */
/* 		printf("%d ", vec1[i]); */
/* 	} */
/* 	printf("\n"); */

/* 	printf("\nArray (parallel):\n"); */
/* 	for(i=0;i<(I*K);i++){ */
/* 		printf("%d ", vec2[i]); */
/* 	} */
/* 	printf("\n"); */

/* 	for(i=1;i<0;i++) printf("!%d ",i); */

/* 	free(vec1); */
/* 	free(vec2); */


/* /\* TRY PARALLEL WITH COUNTER AND REDUCTION *\/ */
/* 	int count=0, count2=0; */

/* 	printf("\nserial loop\n"); */
/* 	time(&time1); */

/* 	for(i=0;i<50;i++){ */
/* 		count2 = 0; */
/* 		while(count2<50e6) count2++; */
/* 		{ */
/* 			count++; */
/* 		} */
/* 		/\* printf("\nvalue of count: %d", count); *\/ */
/* 	} */

/* 	/\* printf("\nvalue of count after loop: %d\n", count); *\/ */

/* 	time(&time2); */
/* 	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1)); */



/* 	printf("\nparallel loop with counter in critical\n"); */
/* 	time(&time1); */

/* 	count=0; */
/* #pragma omp parallel for shared(count) private(count2) */
/* 	for(i=0;i<50;i++){ */
/* 		count2 = 0; */
/* 		while(count2<50e6) count2++; */
/* 		#pragma omp critical */
/* 		{ */
/* 			count++; */
/* 		/\* printf("\nvalue of count: %d", count); *\/ */
/* 		} */
/* 	} */

/* 	/\* printf("\nvalue of count after loop: %d\n", count); *\/ */

/* 	time(&time2); */
/* 	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1)); */





/* ================================================= */
/* ================================================= */
/* ================================================= */

	/* RANDOM NUMBER GENERATION USING GSL */
	/* SETUP */
	int i;
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	printf("\nRandom number generation, serial loop:\n");
	for(i=0;i<30;i++){
		printf("%d\n",gsl_rng_uniform_int(rng,100));
	}


	printf("\nRandom number generation, parallel loop:\n");
#pragma omp parallel for schedule(static,5) ordered
	for(i=0;i<30;i++){
		printf("%d\n",gsl_rng_uniform_int(rng,100));
	}
	

	return 0;
}
