/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "auxiliary.h"
#include "param.h"
#include "pathogens.h"
#include "populations.h"
#include "sampling.h"
#include "dispersal.h"
#include "infection.h"
#include "sumstat.h"



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* create snplist object */
struct snplist * create_snplist(int n){
	struct snplist *out = (struct snplist *) malloc(sizeof(struct snplist));
	if(out == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_snplist]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->snps = calloc(n, sizeof(int));
	if(out->snps == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_snplist]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->length = n;

	return out;
}


/* create allfreq object */
struct allfreq * create_allfreq(int n){
	struct allfreq *out = (struct allfreq *) malloc(sizeof(struct allfreq));
	if(out == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_allfreq]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->freq = (double *) calloc(n, sizeof(double));
	if(out->freq == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_allfreq]\nNo memory left for listing SNPs. Exiting.\n");
		exit(1);
	}

	out->length = n;

	return out;
}




/* Create ts_sumstat */
struct ts_sumstat * create_ts_sumstat(struct param *par){
	int nsteps = par->duration;
	struct ts_sumstat * out = (struct ts_sumstat *) malloc(sizeof(struct ts_sumstat));
	if(out == NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_ts_sumstat]\nNo memory left for storing summary statistics. Exiting.\n");
		exit(1);
	}

	/* use calloc here - values default to 0*/
	out->steps = (int *) calloc(nsteps, sizeof(int));
	out->nbSnps = (int *) calloc(nsteps, sizeof(int));
	out->Hs = (double *) calloc(nsteps, sizeof(double));
	out->meanNbSnps = (double *) calloc(nsteps, sizeof(double));
	out->varNbSnps = (double *) calloc(nsteps, sizeof(double));
	out->meanPairwiseDist = (double *) calloc(nsteps, sizeof(double));
	out->varPairwiseDist = (double *) calloc(nsteps, sizeof(double));
	out->meanPairwiseDistStd = (double *) calloc(nsteps, sizeof(double));
	out->varPairwiseDistStd = (double *) calloc(nsteps, sizeof(double));
	out->Fst = (double *) calloc(nsteps, sizeof(double));

	if(out->nbSnps==NULL || out->Hs==NULL || out->meanNbSnps==NULL || out->varNbSnps==NULL || out->meanPairwiseDist==NULL || out->varPairwiseDist==NULL || out->Fst==NULL){
		fprintf(stderr, "\n[in: sumstat.c->create_ts_sumstat]\nNo memory left for storing summary statistics. Exiting.\n");
		exit(1);
	}
	
	out->maxlength=par->duration;
	out->length=0;
	return out;
}








/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_snplist(struct snplist *in){
	if(in->snps != NULL) free(in->snps);
	if(in != NULL) free(in);
}

void free_allfreq(struct allfreq *in){
	if(in->freq != NULL) free(in->freq);
	if(in != NULL) free(in);
}

void free_ts_sumstat(struct ts_sumstat *in){
	if(in!=NULL){
		free(in->steps);
		free(in->nbSnps);
		free(in->Hs);
		free(in->meanNbSnps);
		free(in->varNbSnps);
		free(in->meanPairwiseDist);
		free(in->varPairwiseDist);
		free(in->meanPairwiseDistStd);
		free(in->varPairwiseDistStd);
		free(in->Fst);
	}
	free(in);
}



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* compute the number of integers differing between two sets */
int dist_a_b(int *a, int *b, int na, int nb){
	int i, out=0;
	for(i=0;i<na;i++){
		if(int_in_vec(a[i], b, nb)<0) out++;
	}
	for(i=0;i<nb;i++){
		if(int_in_vec(b[i], a, na)<0) out++;
	}
	return out;
}




/* count and list number of snps in a sample */
struct snplist * list_snps(struct sample *in, struct param *par){
	int i=0, j=0, N=get_n(in), *pool, poolsize, curNbSnps;
	struct snplist *out;

	/* create pool of snps */
	pool = malloc(par->L * sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: sumstat.c->list_snps]\nNo memory left for creating pool of SNPs. Exiting.\n");
		exit(1);
	}

	/* list and count all SNPs */
	poolsize = 0;
	for(i=0;i<N;i++){
		curNbSnps = get_nb_snps(in->pathogens[i]);
		for(j=0;j<curNbSnps;j++){
			if(int_in_vec(get_snps(in->pathogens[i])[j], pool, poolsize) < 0){
				pool[poolsize++] = get_snps(in->pathogens[i])[j];
			}
		}
	}

	/* make output */
	out = create_snplist(poolsize);

	for(i=0;i<poolsize;i++) out->snps[i] = pool[i];

	/* free local pointers and return result */
	free(pool);
	return out;
}







/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void print_snplist(struct snplist *in){
	int i;
	printf("\nlist of %d SNPs:\n", in->length);
	for(i=0;i<in->length;i++) printf("%d ",in->snps[i]);
	printf("\n");
}



void print_allfreq(struct allfreq *in){
	int i;
	printf("\nallele frequencies for %d SNPs:\n", in->length);
	for(i=0;i<in->length;i++) printf("%.2f ",in->freq[i]);
	printf("\n");
}







struct allfreq * get_frequencies(struct sample *in, struct param *par){
	int i, j, N=get_n(in);
	struct snplist *alleles;
	struct allfreq *out;

	/* list and count alleles */
	alleles = list_snps(in, par);

	/* allocate output */
	out = create_allfreq(alleles->length);
	out->length = alleles->length;

	/* compute frequencies */
	for(i=0;i<N;i++){
		for(j=0;j<alleles->length;j++){
			if(int_in_vec(alleles->snps[j], get_snps(in->pathogens[i]), get_nb_snps(in->pathogens[i])) > -1) 
				out->freq[j] = out->freq[j] + 1.0;
		}
	}

	for(j=0;j<alleles->length;j++) out->freq[j] = out->freq[j]/((double) N);

	/* free memory and return results */
	free_snplist(alleles);
	return out;
}





double hs(struct sample *in, struct param *par){
	int i;
	double out;
	struct allfreq *freq;

	/* get allele frequencies */
	freq = get_frequencies(in, par);

	/* compute Hs */
	out = 0.0;
	for(i=0;i<freq->length;i++){
		out = out + freq->freq[i] * freq->freq[i];
	}
	out = out / freq->length;
	out = 1.0 - out;

	/* free local pointers and return */
	free_allfreq(freq);
	return out;
}





double hs_full_genome(struct sample *in, struct param *par){
	int i;
	double out;
	struct allfreq *freq;

	/* get allele frequencies */
	freq = get_frequencies(in, par);

	/* compute Hs */
	out = 0.0;
	for(i=0;i<freq->length;i++){
		out = out + freq->freq[i] * freq->freq[i];
	}
	out = out + (par->L - freq->length); /* fixed loci */
	out = out / par->L;
	out = 1.0 - out;

	/* free local pointers and return */
	free_allfreq(freq);
	return out;
}




int nb_snps(struct sample *in, struct param *par){
	struct snplist *alleles;
	int out;
	alleles = list_snps(in, par);
	out = alleles->length;
	free_snplist(alleles);
	return out;
}




double mean_nb_snps(struct sample *in){
	int i, N=get_n(in);
	double out=0;
	for(i=0;i<N;i++){
		out += get_nb_snps(in->pathogens[i]);
	}
	out = (double) out / (double) N;
	return out;
}




double var_nb_snps(struct sample *in){
	int i, N=get_n(in);
	double out=0, mu=mean_nb_snps(in);
	for(i=0;i<N;i++){
		out += pow(get_nb_snps(in->pathogens[i]) - mu,2); /* \sum x_i - \bar{x}*/
	}
	out = (double) out / (double) (N-1.0);
	return out;
}




struct distmat_int * pairwise_dist(struct sample *in, struct param *par){
	int i, j, N=get_n(in), counter=0, length=N*(N-1)/2;
	struct distmat_int * out;
	out = create_distmat_int(length);

	/* computations */
	for(i=0;i<N-1;i++){
		for(j=i+1;j<N;j++){
			out->x[counter++] = dist_a_b(get_snps(in->pathogens[i]),get_snps(in->pathogens[j]),get_nb_snps(in->pathogens[i]),get_nb_snps(in->pathogens[j]));
		}
	}

	out->n=N;
	out->length=length;

	return out;
}





double mean_pairwise_dist(struct sample *in, struct param *par){
	struct distmat_int * mat = pairwise_dist(in, par);
	int i, n=mat->length;
	double out=0.0;

	/* computations */
	for(i=0;i<n;i++) out += mat->x[i];
	out = out / (double) n;

	/* free memory and return */
	free_distmat_int(mat);
	return out;
}




double mean_pairwise_dist_std(struct sample *in, struct param *par){
	double out = mean_pairwise_dist(in, par);
	int n = nb_snps(in, par);
	out= out/(double) n;
	return out;
}





double var_pairwise_dist(struct sample *in, struct param *par){
	struct distmat_int * mat = pairwise_dist(in, par);
	int i, n=mat->length;
	double mu=0.0, out=0.0;

	/* computations */
	for(i=0;i<n;i++) mu += mat->x[i];
	mu = mu/n;
	for(i=0;i<n;i++) out += pow((double) mat->x[i] - mu, 2);
	out = out / (double) (n-1);

	/* free memory and return */
	free_distmat_int(mat);
	return out;
}





double var_pairwise_dist_std(struct sample *in, struct param *par){
	double out = var_pairwise_dist(in, par);
	int n = nb_snps(in, par);
	out=out/(double) (n*n);
	return out;
}





double fst(struct sample *in, struct param *par){
	/* Fst = 1 - Hsbar/Ht with
	Hsbar: expected H averaged over groups
	Ht: expected H over all data */

	int i, npop=get_npop_samp(in), sumweights=0;
	double Ht, Hsbar=0, out;
	struct sample ** listsamp;

	/* global exp heteroz */
	Ht = hs(in, par);

	/* get Hs per population*/
	listsamp = seppop(in, par);
	for(i=0;i<npop;i++){
		Hsbar += hs(listsamp[i], par)*get_n(listsamp[i]);
		sumweights += get_n(listsamp[i]);
	}

	Hsbar = Hsbar/(double) sumweights;

	/* Fst */
	out = 1.0 - (Hsbar/Ht);

	/* free local pointers and return */
	for(i=0;i<npop;i++) free_sample(listsamp[i]);
	free(listsamp);
	return out;
}






void fill_ts_sumstat(struct ts_sumstat *in, struct sample *samp, int step, struct param *par){
	int idx = in->length;

	if(idx > in->maxlength){
		fprintf(stderr, "\n[in: sumstat.c->fill_ts_sumstat]\n. ts_sumstat object is not long enough to store output of step %d. Exiting.\n", step);
		exit(1);
	}

	in->steps[idx] = step;
	in->nbSnps[idx] = nb_snps(samp, par);
	in->Hs[idx] = hs(samp, par);
	in->meanNbSnps[idx] = mean_nb_snps(samp);
	in->varNbSnps[idx] = var_nb_snps(samp);
	in->meanPairwiseDist[idx] = mean_pairwise_dist(samp, par);
	in->varPairwiseDist[idx] = var_pairwise_dist(samp, par);
	in->meanPairwiseDistStd[idx] = mean_pairwise_dist_std(samp, par);
	in->varPairwiseDistStd[idx] = var_pairwise_dist_std(samp, par);
	in->Fst[idx] = fst(samp, par);
	in->length = in->length + 1;
}







/*
  =========================
   === TESTING FUNCTIONS ===
   =========================
*/



/* gcc line:

   gcc -o sumstat param.c auxiliary.c pathogens.c populations.c dispersal.c infection.c sampling.c sumstat.c -Wall -O0 -lgsl -lgslcblas

   valgrind --leak-check=yes sumstat
*/




/* int main(){ */
/* 	/\* Initialize random number generator *\/ */
/* 	time_t t; */
/* 	t = time(NULL); // time in seconds, used to change the seed of the random generator */
/* 	gsl_rng * rng; */
/* 	const gsl_rng_type *typ; */
/* 	gsl_rng_env_setup(); */
/* 	typ=gsl_rng_default; */
/* 	rng=gsl_rng_alloc(typ); */
/* 	gsl_rng_set(rng,t); // changes the seed of the random generator */
/* 	int i, j; */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->rng = rng; */
/* 	par->npop = 2; */
/* 	int popsizes[2] = {1000,200}; */
/* 	par->popsizes = popsizes; */
/* 	par->nstart = 10; */
/* 	par->t1 = 1; */
/* 	par->t2 = 2; */
/* 	par->beta = 1.1; */
/* 	int nbnb[2] = {2,2}; */
/* 	par->cn_nb_nb = nbnb; */
/* 	int listnb[4] = {0,1,1,0}; */
/* 	par->cn_list_nb = listnb; */
/* 	double weights[4] = {0.9,0.1,0.99,0.11}; */
/* 	par->cn_weights = weights; */
/* 	struct network *cn = create_network(par); */
/* 	par->mu = 0.01; */
/* 	par->L = 100; */
/* 	par->muL = par->mu*par->L; */

/* 	/\* CREATE METAPOPULATION *\/ */
/* 	struct metapopulation * metapop = create_metapopulation(par); */
/* 	printf("\n## CREATED METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */


/* 	/\* SIMULATE OUTBREAK OVER A FEW TIMESTEPS *\/ */
/* 	for(i=0;i<3;i++){ */
/* 		age_metapopulation(metapop, par); */
/* 		for(j=0;j<get_npop(metapop);j++){ */
/* 			process_infections(get_populations(metapop)[j], metapop, cn, par); */
/* 		} */
/* 		printf("\n - METAPOPULATION @ step %d -", i); */
/* 		print_metapopulation(metapop, FALSE); */

/* 	} */

/* 	printf("\n## RESULTING METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */


/* 	printf("\n## RESULTING SAMPLE ##"); */
/* 	struct sample *samp; */
/* 	samp = draw_sample(metapop,20,par); */
/* 	print_sample(samp, TRUE); */

/* 	/\* TEST SUMMARY STATISTICS *\/ */
/* 	/\* test allele listing *\/ */
/* 	struct snplist *snpbilan; */
/* 	snpbilan = list_snps(samp, par); */
/* 	print_snplist(snpbilan); */

/* 	/\* test allele frequencies *\/ */
/* 	struct allfreq *freq; */
/* 	freq = get_frequencies(samp, par); */
/* 	print_allfreq(freq); */

/* 	/\* test Hs*\/ */
/* 	double Hs = hs(samp,par); */
/* 	printf("\nHs = %0.3f\n", Hs); */

/* 	/\* test Hs full genome *\/ */
/* 	Hs = hs_full_genome(samp,par); */
/* 	printf("\nHs (full genome) = %0.5f\n", Hs); */

/* 	/\* test nb of snps *\/ */
/* 	int nball = nb_snps(samp,par); */
/* 	printf("\nnumber of SNPs = %d\n", nball); */

/* 	/\* test mean nb of snps *\/ */
/* 	double temp = mean_nb_snps(samp); */
/* 	printf("\nmean number of SNPs = %.2f\n", temp); */

/* 	/\* test var nb of snps *\/ */
/* 	temp = var_nb_snps(samp); */
/* 	printf("\nvariance of number of alleles = %.2f\n", temp); */

/* 	/\* test pairwise distances *\/ */
/* 	struct distmat_int *mat = pairwise_dist(samp, par); */
/* 	print_distmat_int(mat); */

/* 	/\* test mean pairwise distances *\/ */
/* 	temp = mean_pairwise_dist(samp,par); */
/* 	printf("\nmean pairwise distance: %.2f", temp); */

/* 	/\* test variance of pairwise distances *\/ */
/* 	temp = var_pairwise_dist(samp,par); */
/* 	printf("\nvar pairwise distance: %.2f", temp); */

/* 	/\* test Fst *\/ */
/* 	temp = fst(samp,par); */
/* 	printf("\nfst: %.2f", temp); */

/* 	printf("\n\n"); */


/* 	/\* free memory *\/ */
/* 	free_metapopulation(metapop); */
/* 	free_sample(samp); */
/* 	free_network(cn); */
/* 	free(par); */
/* 	gsl_rng_free(rng); */
/* 	free_snplist(snpbilan); */
/* 	free_allfreq(freq); */
/* 	free_distmat_int(mat); */

/* 	return 0; */
/* } */
