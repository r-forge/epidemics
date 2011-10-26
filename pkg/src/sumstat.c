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




/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* check if an integer is in a vector of integers */
bool int_in_vec(int x, int *vec, int vecSize){
	int i=0;
	while(x!=vec[i] && i<vecSize) i++;
	if(i==vecSize || vecSize<1) return FALSE;
	return TRUE;
}






/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

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
int main(){
	int  i, vec[5]={1,2,3,4,5};

	for(i=0;i<10;i++) printf("\ni=%d, result:%d", i, int_in_vec(i,vec,5));
	
	return 0;
}
