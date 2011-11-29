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
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct network * create_network(struct param *par){
	int i, j, counter;
	struct network *out;

	/* ALLOCATE MEMORY FOR OUTPUT */
	out = (struct network *) calloc(1, sizeof(struct network));
	if(out == NULL){
		fprintf(stderr, "\n[in: dispersal.c->create_network]\nNo memory left for creating connection network. Exiting.\n");
		exit(1);
	}

	/* FILL CONTENT */
	/* allocate memory */
	out->n = par->npop;
	out->nbNb = (int *) calloc(out->n, sizeof(int));
	out->listNb = (int **) calloc(out->n, sizeof(int *));
	out->weights = (double **) calloc(out->n, sizeof(double *));

	if(out->nbNb == NULL || out->listNb == NULL || out->weights == NULL){
		fprintf(stderr, "\n[in: dispersal.c->create_network]\nNo memory left for creating connection network. Exiting.\n");
		exit(1);
	}

	/* nb of neighbours*/
	for(i=0;i<par->npop;i++){
		out->nbNb[i] = par->cn_nb_nb[i];
	}

	/* list of neighbours and weights */
	counter = 0;
	for(i=0;i<par->npop;i++){
		out->listNb[i] = (int *) calloc(out->nbNb[i], sizeof(int));
		out->weights[i] = (double *) calloc(out->nbNb[i], sizeof(double));
		if(out->listNb[i] == NULL || out->weights[i] == NULL){
			fprintf(stderr, "\n[in: dispersal.c->create_network]\nNo memory left for creating connection network. Exiting.\n");
			exit(1);
		}

		for(j=0;j<out->nbNb[i];j++){
			out->listNb[i][j] = par->cn_list_nb[counter];
			out->weights[i][j] = par->cn_weights[counter++];
		}
	}

	return out;
}






/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free network */
void free_network(struct network *in){
	int i;
	if(in != NULL){
		for(i=0;i<in->n;i++){
			free(in->listNb[i]);
			free(in->weights[i]);
		}
		free(in->nbNb);
		free(in->listNb);
		free(in->weights);
		free(in);
	}
}







/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

void print_network(struct network *in, bool detail){
	int i, j, temp=0;
	printf("\nconnection network with %d vertices", in->n);
	printf("\ndegrees:");
	for(i=0;i<in->n;i++) {
		printf("%d ", in->nbNb[i]);
		temp += in->nbNb[i];
	}
	printf("\ntotal nb of connections: %d\n", temp);
	if(detail){
	printf("\nconnections (weights): %d\n", temp);
		for(i=0;i<in->n;i++) {
			printf("item %d: ", i);
			for(j=0;j<in->nbNb[i];j++){
				printf("%d(%.1f) ", in->listNb[i][j], in->weights[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/




/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


/* gcc line:

   gcc -o dispersal param.c auxiliary.c dispersal.c -Wall -O0 -lgsl -lgslcblas
  
   valgrind --leak-check=yes dispersal

*/

/* int main(){ */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->npop=3; */

/* 	int nbNb[3] = {1,2,1}; */
/* 	int listNb[4] = {1,0,2,1}; */
/* 	double weights[4] = {1,0.2,0.5,0.1}; */
/* 	par->cn_nb_nb= nbNb; */
/* 	par->cn_list_nb = listNb; */
/* 	par->cn_weights = weights; */

/* 	struct network * cn = create_network(par); */

/* 	print_network(cn, TRUE); */

/* 	/\* free memory *\/ */
/* 	free(par); */
/* 	free_network(cn); */

/* 	return 0; */
/* } */



