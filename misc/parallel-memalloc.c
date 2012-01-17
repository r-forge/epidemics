#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define B_SIZE 1000



struct foo{
	int b[B_SIZE];
};






struct foo * create_foo(){
	struct foo *out = (struct foo *) calloc(1,sizeof(struct foo));
	return(out);
}




void free_foo(struct foo *in){
	free(in);
}




/*

gcc -o para4 parallel-memalloc.c  -lgsl -lgslcblas -fopenmp -O3


*/
int main(){
	const int N=1000000;
	int i, k;
	time_t time1, time2;
	struct foo **x = (struct foo **) malloc(N*sizeof(struct foo *));

	/* SERIAL CODE */
	printf("\nAllocating memory, serial code");
	time(&time1);
	for(i=0;i<N;i++){
		x[i] = create_foo();
	}
	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));


	printf("\nFreeing memory, serial code");
	time(&time1);
	for(i=0;i<N;i++){
		free_foo(x[i]);
	}
	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));



	/* PARALLEL CODE */
	printf("\nAllocating memory, parallel code");
	time(&time1);
#pragma omp parallel for schedule(static,100000)
	for(i=0;i<N;i++){
		x[i] = create_foo();
	}
	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));


	printf("\nFreeing memory, parallel code");
	time(&time1);
#pragma omp parallel for
	for(i=0;i<N;i++){
		free_foo(x[i]);
	}
	time(&time2);

	printf("\nthe operation took %d seconds\n", (int) (time2-time1));

	return 0;
}
