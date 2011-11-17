/* timer.c */
 
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
 
main()
{  
	const int DIM=1000;
	int i,j, k;
	time_t t1,t2;
	int A[DIM][DIM], B[DIM][DIM];

/* first dimension first*/
	time(&t1);
	for(k=0;k<200;k++){
		for(i=0;i<DIM;i++){
			for(j=0;j<DIM;j++){
				A[i][j] = i+j;
				B[i][j] = 0;
			}
		}

		for(i=0;i<DIM;i++){
			for(j=0;j<DIM;j++){
				B[i][j] = A[i][j]*2*3*4*5*6*7*8*9/2/3/4/5/6/7/8/9;
			}
		}
	}
	time(&t2);

	printf("\nMethod 1: %d ", (int) t2-t1);


/* last dimension first*/
	time(&t1);
	
	for(k=0;k<200;k++){
		for(j=0;j<DIM;j++){
			for(i=0;i<DIM;i++){
				A[i][j] = i+j;
				B[i][j] = 0;
			}
		}

		for(j=0;j<DIM;j++){
			for(i=0;i<DIM;i++){
				B[i][j] = A[i][j]*2*3*4*5*6*7*8*9/2/3/4/5/6/7/8/9;
			}
		}
	}
	time(&t2);

	printf("\nMethod 2: %d\n ", (int) t2-t1);


       

}
