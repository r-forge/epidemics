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
#include "dispersal.h"




/*
   =================
   === ACCESSORS ===
   =================
*/

double ** get_mat(struct dispmat * in){
	return(in->mat);
}



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty dispmat */
/* 'in' is a vector of k*k doubles, k being the number of populations */
struct dispmat * create_dispmat(struct param *par){
	int i, j, k=par->npop, counter;
	double temprowsum;
	struct dispmat *out;

	/* allocate memory */
	out = (struct dispmat *) calloc(1, sizeof(struct dispmat));
	if(out == NULL){
		fprintf(stderr, "\n[in: dispersal.c->create_dispmat]\nNo memory left for creating matrix of dispersal. Exiting.\n");
		exit(1);
	}

	out->mat = (double **) calloc(k, sizeof(double *));
	if(out->mat == NULL){
		fprintf(stderr, "\n[in: dispersal.c->create_dispmat]\nNo memory left for creating matrix of dispersal. Exiting.\n");
		exit(1);
	}

	for(i=0;i<k;i++){
		out->mat[i] = (double *) calloc(k, sizeof(double));
		if(out->mat[i] == NULL){
			fprintf(stderr, "\n[in: dispersal.c->create_dispmat]\nNo memory left for creating matrix of dispersal. Exiting.\n");
			exit(1);
		}
	}

	/* fill in the matrix */
	if(par->pdisp == NULL) { /* no dispersal */
		for(i=0;i<k;i++){
			for(j=0;j<k;j++){
				if(i==j) {
					out->mat[i][j] = 1.0;
				} else {
					out->mat[i][j] = 0.0;
				}
			}
		}
	} else { /* there is dispersal */
		counter = 0;
		for(i=0;i<k;i++){
			for(j=0;j<k;j++){
				out->mat[i][j] = par->pdisp[counter++];
			}
		}
	}


	/* row-standardize the matrix */
	for(i=0;i<k;i++){
		temprowsum = 0.0;
		for(j=0;j<k;j++){
			temprowsum += out->mat[i][j];
		}
		for(j=0;j<k;j++){
			out->mat[i][j] = out->mat[i][j] / temprowsum;
		}
	}


	out->n = k;
	return out;
}






/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free dispmat */
void free_dispmat(struct dispmat *in){
	int i, k=in->n;
	if(in != NULL){
		for(i=0;i<k;i++){
			free(in->mat[i]);
		}
	}
	free(in->mat);
	free(in);
}






/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

void print_dispmat(struct dispmat *in){
	int i, j, k=in->n;
	for(i=0;i<k;i++){
		printf("\npopulation %d\n",i);
		for(j=0;j<k;j++){
			printf("%.2f \t",in->mat[i][j]);
		}
	}
}







/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* return the id of a new population */
int disperse(struct pathogen * pathogen, struct dispmat *disp, struct param *par){
	int i=1, k=disp->n, popid = get_popid(pathogen), cumprob=get_mat(disp)[popid][0];
	double x=gsl_rng_uniform(par->rng); /* nb between 0 and 1 */

	while(x > cumprob && i<k){
		cumprob += get_mat(disp)[popid][i++]; /* to check */
	}
	printf("new pop: %d",i);
	return i;
}
