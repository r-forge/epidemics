#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define B_SIZE 10000000



struct foo{
	int a, *b, count;
};






struct foo ** create_foo_ar(int N){
	int i;
	struct foo **out = (struct foo **) calloc(N,sizeof(struct foo*));
	for(i=0;i<N;i++){
		out[i] = (struct foo *) calloc(1,sizeof(struct foo));
		out[i]->a = 0;
		out[i]->b = (int *) calloc(B_SIZE,sizeof(int));
		out[i]->count = 0;
	}
	return(out);
}




void free_foo_ar(struct foo **in, int N){
	int i;
	for(i=0;i<N;i++){
		free(in[i]->b);
		free(in[i]);
	}

	free(in);
}




/* void f1(){ */
/* 	int i,j; */
/* 	for(i=0;i<5000;i++){ */
/* 		j++; */
/* 	} */
/* } */




/* void ini_foo(struct foo *in, int val){ */
/* 	int i; */
/* 	for(i=0;i<B_SIZE;i++){ */
/* 		in->b[i] = val + i; */
/* 	} */
/* 	in->a = in->a + val; */
/* } */



void incr_b(struct foo *in){
	int i;
	for(i=0;i<B_SIZE;i++){
		in->b[i] = in->b[i] + 1;
		if(in->b[i] > i) in->b[i] = 0;
	}
	in->count = in->count + 1;
}

void incr_b_until(struct foo *in, int until){
	int i, x = (until > B_SIZE) ? B_SIZE : until;
	for(i=0;i<until;i++){
		in->b[i] = in->b[i] + 1;
		if(in->b[i] > i) in->b[i] = 0;
	}
	in->count = in->count + 1;
}


void incr_b_paral(struct foo *in){
	int i;
#pragma omp parallel for
	for(i=0;i<B_SIZE;i++){
		in->b[i] = in->b[i] + 1;
		if(in->b[i] > i) in->b[i] = 0;
	}
	in->count = in->count + 1;
}



/*

gcc -o para2 parallel-bigarray.c  -lgsl -lgslcblas -fopenmp -O3


*/
int main(){
	const int N=25, nstep = 100;
	int i, k;
	time_t time1, time2;
	struct foo **x = create_foo_ar(N);

	/* SERIAL CODE */
	printf("\nProcessing array, serial code");
	time(&time1);
	for(i=0;i<nstep;i++){
		for(k=0;k<N;k++){
			//incr_b(x[k]);
			//incr_b_until(x[k], i * B_SIZE/nstep);
			incr_b_paral(x[k]);
		}
	}
	time(&time2);

	printf("\nthe operation took %d-%d seconds\n", (int) (time2-time1)-1, (int) (time2-time1)+1);


	/* PARALLEL CODE */
	printf("\nProcessing array, parallel code");
	time(&time1);
/* #pragma omp parallel for private(k) */
	for(i=0;i<nstep;i++){
//#pragma omp parallel for schedule(static,1)
		for(k=0;k<N;k++){
			//incr_b(x[k]);
			//incr_b_until(x[k], i * B_SIZE/nstep);
			incr_b_paral(x[k]);
		}
	}
	time(&time2);

	printf("\nthe operation took %d-%d seconds\n", (int) (time2-time1)-1, (int) (time2-time1)+1);


	/* FREE STUFF AND RETURN */
	free_foo_ar(x,N);
	return 0;
}
