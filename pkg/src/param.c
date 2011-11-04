/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions handle simulation parameters.
*/

#include "common.h"
#include "auxiliary.h"
#include "param.h"

/* Free param */
void free_param(struct param *in){
	gsl_rng_free(in->rng);
	free(in);
}



void check_param(struct param *in){
	int i, checkOK=0;

	/* nstart & K */
	if(in->nstart > in->nsus[0]){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: initial number of infections greater than host population.\n");
		exit(1);
	}

	/* nsus */
	for(i=0;i<in->npop;i++){
		if(in->nsus[i] < 1){
			fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: less than one host in population %d.\n",i);
			exit(1);
		}
	}

	/* nstart */
	if(in->nstart < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: less than one initial infection.\n");
		exit(1);
	}

	/* L */
	if(in->L < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: genome length less than one.\n");
		exit(1);
	}

	if(in->L < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: genome length less than one.\n");
		exit(1);
	}

	/* t1 */
	if(in->t1 < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: time to infectiousness less than one.\n");
		exit(1);
	}

	/* t2 */
	if(in->t2 < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: infection duration less than one.\n");
		exit(1);
	}

	/* t1 & t2 */
	if(in->t1 > in->t2){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: infections last less than time to infectiousness.\n");
		exit(1);
	}


	/* mu */
	if(in->mu < 0.0){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: negative mutation rate.\n");
		exit(1);
	}

	/* beta */
	if(in->beta < 0.0){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: negative reproductive number.\n");
		exit(1);
	}

	/* n_sample */
	if(in->n_sample < 0){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: negative sample size.\n");
		exit(1);
	}

	/* t_sample */
	if(min_int(in->t_sample, in->n_sample) < 0){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: negative sampling time detected.\n");
		exit(1);
	}

	if(max_int(in->t_sample, in->n_sample) > in->duration){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: sampling time span (%d) longer than epidemic duration (%d).\n", max_int(in->t_sample, in->n_sample), in->duration);
		exit(1);
	}

	/* npop */
	if(in->npop < 1){
		fprintf(stderr, "\n[in: param.c->check_param]\nParameter error: number of populations < 1.\n");
		exit(1);
	}

}




/* print parameters */
void print_param(struct param *in){
	int i, totnsus=0;
	printf("\n-- simulation parameters --");

	/* genetic parameters */
	printf("\nmutation rate: %f   genome length: %d", in->mu, in->L);

	/* epidemiological parameters*/
	printf("\nnb of populations: %d", in->npop);
	printf("\nduration of the epidemic: %d   ", in->duration);
	printf("\nnb susceptible per populations:");
	for(i=0;i<in->npop;i++) printf("%d\t", in->nsus[i]);
	printf("\ntransmission rate (beta): %.2f", in->beta);
	for(i=0;i<in->npop;i++) totnsus += in->nsus[i];
	printf("\ntotal nb of susceptible: %d", totnsus);
	printf("\nnb initial infections: %d", in->nstart);
	printf("\nstart infectious period: %d   ", in->t1);
	printf("\nend infectious period: %d   ", in->t2);
	/* printf("\nsampling time: %d   ", in->t_sample); */
	printf("\nsample size: %d   ", in->n_sample);
	printf("\n");
}
