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
#include "sumstat.h"
#include "dispersal.h"
#include "infection.h"
#include "inout.h"









/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Function to be called from R */
void R_epidemics(int *seqLength, double *mutRate, int *npop, int *nHostPerPop, double *beta, int *nStart, int *t1, int *t2, int *Nsample, int *Tsample, int *duration, int *nbnb, int *listnb, double *pdisp){
	int i, nstep, counter_sample = 0, tabidx;

	/* Initialize random number generator */
	int j;
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	/* transfer simulation parameters */
	struct param * par;
	par = (struct param *) malloc(sizeof(struct param));
	par->L = *seqLength;
	par->mu = *mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = *npop;
	par->popsizes = nHostPerPop;
	par->beta = *beta;
	par->nstart = *nStart;
	par->t1 = *t1;
	par->t2 = *t2;
	par->t_sample = Tsample;
	par->n_sample = *Nsample;
	par->duration = *duration;
	par->cn_nb_nb = nbnb;
	par->cn_list_nb = listnb;
	par->cn_weights = pdisp;
	par->samp_min_size = -1; /* no minimum sample size */

	/* check/print parameters */
	check_param(par);
	printf("\n ! Beta version used !\n");
	print_param(par);

	/* dispersal matrix */
	struct network *cn = create_network(par);
	/* print_network(cn, TRUE); */

	/* group sizes */
	struct ts_groupsizes ** grpsizes = create_list_ts_groupsizes(par);

	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);

	/* get sampling schemes (timestep+effectives) */
	translate_dates(par);
	struct table_int *tabdates = get_table_int(par->t_sample, par->n_sample);
	printf("\n\nsampling at timesteps:");
	print_table_int(tabdates);

	/* create sample */
	struct sample ** samplist = (struct sample **) malloc(tabdates->n * sizeof(struct sample *));
	struct sample *samp;


	/* MAKE METAPOPULATION EVOLVE */
	nstep = 0;
	while(get_total_nsus(metapop)>0 && (get_total_ninf(metapop)+get_total_nexp(metapop))>0 && nstep<par->duration){
		nstep++;

		/* age metapopulation */
		age_metapopulation(metapop, par);

		/* process infections */
		for(j=0;j<get_npop(metapop);j++){
			process_infections(get_populations(metapop)[j], metapop, cn, par);
		}

		/* draw samples */
		if((tabidx = int_in_vec(nstep, tabdates->items, tabdates->n)) > -1){ /* TRUE if step must be sampled */
			samplist[counter_sample++] = draw_sample(metapop, tabdates->times[tabidx], par);
		}

		fill_list_ts_groupsizes(grpsizes, metapop, nstep);

	}

	/* we stopped after 'nstep' steps */
	if(nstep < par->duration){
		printf("\nEpidemics ended at time %d, before last sampling time (%d).\n", nstep, par->duration);
	} else {

		/* printf("\n\n-- FINAL METAPOPULATION --"); */
		/* print_metapopulation(metapop, FALSE); */

		/* merge samples */
		samp = merge_samples(samplist, tabdates->n, par);

		/* write sample to file */
		printf("\n\nWriting sample to file 'out-sample.txt'\n");
		write_sample(samp);

		/* free memory */
		free_sample(samp);

	}

	/* write group sizes to file */
	printf("\n\nPrinting group sizes to file 'out-popsize.txt'\n");
	write_list_ts_groupsizes(grpsizes, par);


	/* free memory */
	free_metapopulation(metapop);
	free_param(par);
	for(i=0;i<counter_sample;i++) free_sample(samplist[i]);
	free(samplist);
	free_table_int(tabdates);
	free_network(cn);
	free_list_ts_groupsizes(grpsizes, par);
}







/* Function to be called from R */
void R_monitor_epidemics(int *seqLength, double *mutRate, int *npop, int *nHostPerPop, double *beta, int *nStart, int *t1, int *t2, int *Nsample, int *Tsample, int *duration, int *nbnb, int *listnb, double *pdisp, int *minSize){
		int nstep;

	/* Initialize random number generator */
	int j;
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	/* transfer simulation parameters */
	struct param * par;
	par = (struct param *) malloc(sizeof(struct param));
	par->L = *seqLength;
	par->mu = *mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = *npop;
	par->popsizes = nHostPerPop;
	par->beta = *beta;
	par->nstart = *nStart;
	par->t1 = *t1;
	par->t2 = *t2;
	par->t_sample = Tsample;
	par->n_sample = *Nsample;
	par->duration = *duration;
	par->cn_nb_nb = nbnb;
	par->cn_list_nb = listnb;
	par->cn_weights = pdisp;
	par->samp_min_size = *minSize;

	/* check/print parameters */
	check_param(par);
	printf("\n ! Beta version used !\n");
	print_param(par);

	/* dispersal matrix */
	struct network *cn = create_network(par);
	/* print_network(cn, TRUE); */

	/* group sizes */
	struct ts_groupsizes ** grpsizes = create_list_ts_groupsizes(par);
	struct ts_sumstat ** sumstats = create_list_ts_sumstat(par);

	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);

	/* create sample */
	struct sample **samples = create_list_sample(par);


	/* MAKE METAPOPULATION EVOLVE */
	nstep = 0;
	while(get_total_nsus(metapop)>0 && (get_total_ninf(metapop)+get_total_nexp(metapop))>0 && nstep<par->duration){
		nstep++;

		/* age metapopulation */
		age_metapopulation(metapop, par);

		/* process infections */
		for(j=0;j<get_npop(metapop);j++){
			process_infections(get_populations(metapop)[j], metapop, cn, par);
		}


		/* draw samples */
		samples[0] = draw_sample(metapop, par->n_sample, par); /* sample from the metapop */
		for(j=0;j<get_npop(metapop);j++){
			samples[j+1] = draw_sample_onepop(get_populations(metapop)[j],par->n_sample, par);
		}

		/* compute statistics */
		fill_list_ts_sumstat(sumstats, samples, nstep, par);

		/* get group sizes */
		fill_list_ts_groupsizes(grpsizes, metapop, nstep);
	}

	/* write group sizes to file */
	printf("\n\nWriting results to file...");
	write_list_ts_groupsizes(grpsizes, par);
	write_list_ts_sumstat(sumstats,par);
	printf("done.\n\n");

	/* free memory */
	free_metapopulation(metapop);
	free_network(cn);
	free_list_sample(samples, par);
	free_list_ts_groupsizes(grpsizes, par);
	free_list_ts_sumstat(sumstats,par);
	free_param(par);
} /* end R_monitor_epidemics */



















/* all-in-one function testing epidemics growth, summary statistics, etc. */
void test_epidemics(int seqLength, double mutRate, int npop, int *nHostPerPop, double beta, int nStart, int t1, int t2, int Nsample, int *Tsample, int duration, int *nbnb, int *listnb, double *pdisp){
	int i, j, nstep=0, tabidx, counter_sample = 0;

	/* Initialize random number generator */
	time_t t;
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	/* transfer simulation parameters */
	struct param * par;
	par = (struct param *) malloc(sizeof(struct param));
	par->L = seqLength;
	par->mu = mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = npop;
	par->npop = npop;
	par->popsizes = nHostPerPop;
	par->beta = beta;
	par->nstart = nStart;
	par->t1 = t1;
	par->t2 = t2;
	par->t_sample = Tsample;
	par->n_sample = Nsample;
	par->duration = duration;
	par->cn_nb_nb = nbnb;
	par->cn_list_nb = listnb;
	par->cn_weights = pdisp;

	/* check/print parameters */
	check_param(par);
	/* print_param(par); */

	/* dispersal matrix */
	struct network *cn = create_network(par);

	/* group sizes */
	struct ts_groupsizes ** grpsizes = create_list_ts_groupsizes(par);


	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);

	/* /\* check memory representation *\/ */
	/* printf("\n"); */
	/* printf("\nAddress of metapop: %p (dec: %u)", metapop, (unsigned int) metapop); */
	/* printf("\nAddresses of populations:"); */
	/* for(i=0;i<par->npop;i++){ */
	/* 	printf("\nAddress of population %d: %p (dec: %u)", i, metapop->populations[i], (unsigned int) metapop->populations[i]); */
	/* } */
	/* printf("\n"); */


	/* get sampling schemes (timestep+effectives) */
	translate_dates(par);
	struct table_int *tabdates = get_table_int(par->t_sample, par->n_sample);
	/* printf("\n\nsampling at timesteps:"); */
	/* print_table_int(tabdates); */

	/* create sample */
	struct sample ** samplist = (struct sample **) malloc(tabdates->n * sizeof(struct sample *));
	struct sample *samp;


/* 	/\* SANITY CHECK: TEST PARALLELISM *\/ */
/* THIS ONE IS OK 25-26s VS 16-20s */
/* #pragma omp parallel for schedule(static,1) */
/* 	for(i=0;i<10;i++){ */
/* 		for(j=0;j<100000000;j++){ */
/* 			gsl_rng_uniform_int(par->rng,par->L); */
/* 		} */
/* 	} */




	/* MAKE METAPOPULATION EVOLVE */
	nstep = 0;
	while(get_total_nsus(metapop)>0 && (get_total_ninf(metapop)+get_total_nexp(metapop))>0 && nstep<par->duration){
		nstep++;


/* #pragma omp parallel for schedule(static,1) */
/* 	for(j=0;j<get_npop(metapop);j++){ */
/* 		//testpop1(get_populations(metapop)[j], 5e7, rng); /\* WORKS! *\/ */
/* 		//testpop2(get_populations(metapop)[j], 5e7, rng, par); /\* WORKS! *\/ */
/* 		//testpop3(get_populations(metapop)[j], 5e7, rng, metapop, par); /\* WORKS! *\/ */
/* 		//testpop4(get_populations(metapop)[j], 5e5, rng, metapop, par); /\* WORKS! *\/ */
/* 		testpop5(get_populations(metapop)[j], 5e5, rng, metapop, par); /\* DOES NOT WORK! *\/ */
/* 	} */

		/* age metapopulation */
		age_metapopulation(metapop, par);
		//age_metapopulation2(metapop, t1,t2);

		/* process infections */
/* #if USE_OMP */
/* #pragma omp parallel for */
/* #endif */
		for(j=0;j<get_npop(metapop);j++){
			process_infections(get_populations(metapop)[j], metapop, cn, par);
		}


		/* draw samples */
		if((tabidx = int_in_vec(nstep, tabdates->items, tabdates->n)) > -1){ /* TRUE if step must be sampled */
			samplist[counter_sample++] = draw_sample(metapop, tabdates->times[tabidx], par);
		}

		fill_list_ts_groupsizes(grpsizes, metapop, nstep);

	}


	/* /\* we stopped after 'nstep' steps *\/ */
	/* if(nstep < par->duration){ */
	/* 	printf("\nEpidemics ended at time %d, before last sampling time (%d).\n", nstep, par->duration); */
	/* } else { */

	/* 	printf("\n\n-- FINAL METAPOPULATION --"); */
	/* 	print_metapopulation(metapop, TRUE); */

	/* 	/\* test samples *\/ */
	/* 	for(i=0;i<tabdates->n;i++) { */
	/* 		printf("\nsample %d\n", i); */
	/* 		print_sample(samplist[i], TRUE); */
	/* 	} */
	/* 	samp = merge_samples(samplist, tabdates->n, par) ; */
	/* 	print_sample(samp, TRUE); */

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
	/* 	free_sample(samp); */
	/* 	free_snplist(snpbilan); */
	/* 	free_allfreq(freq); */
	/* 	free_distmat_int(mat); */

	/* } */

	/* /\* write group sizes to file *\/ */
	/* printf("\n\nPrinting group sizes to file 'out-popsize.txt'"); */
	/* write_ts_groupsizes(grpsizes); */

	/* free memory */
	free_metapopulation(metapop);
	free_param(par);
	for(i=0;i<counter_sample;i++) free_sample(samplist[i]);
	free(samplist);
	free_table_int(tabdates);
	free_network(cn);
	free_list_ts_groupsizes(grpsizes, par);
}




int main(){
/* args: (int seqLength, double mutRate, int npop, int nHostPerPop, double beta, int nStart, int t1, int t2,int Tsample, int Nsample) */
	/* ONE SINGLE PATCH */
	/* double mu=1e-6, beta=2, pdisp[1]={1.0}; */
	/* time_t time1,time2; */
	/* int genoL=1e5, duration=100, npop=1, nstart=10, t1=1, t2=2, nsamp=10; */
	/* int tsamp[10] = {1,1,1,1,1,1,1,0,0,0}, popsize[1] = {50e6}; */
	/* int nbnb[1] = {1}, listnb[1] = {0}; */


	/* 4x4 LATTICE */
	double mu=1e-6, beta=2;
	time_t time1,time2;
	int genoL=1e5, duration=100, npop=25, nstart=10, t1=1, t2=2, nsamp=100;
	int tsamp[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int popsize[25] = {2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000,2000000};
	int nbnb[25] = {3,4,4,4,3,4,5,5,5,4,4,5,5,5,4,4,5,5,5,4,3,4,4,4,3};
	int listnb[105] = {0,1,5,1,0,2,6,2,1,3,7,3,2,4,8,4,3,9,5,0,6,10,6,1,5,7,11,7,2,6,8,12,8,3,7,9,13,9,4,8,14,10,5,11,15,11,6,10,12,16,12,7,11,13,17,13,8,12,14,18,14,9,13,19,15,10,16,20,16,11,15,17,21,17,12,16,18,22,18,13,17,19,23,19,14,18,24,20,15,21,21,16,20,22,22,17,21,23,23,18,22,24,24,19,23};
	double pdisp[105] = {9.999e-1,5e-5,5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,5e-5,5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,2.5e-5,2.5e-5,2.5e-5,2.5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,5e-5,5e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,3.333333e-5,3.333333e-5,3.333333e-5,9.999e-1,5e-5,5e-5};

	time(&time1);
	test_epidemics(genoL, mu, npop, popsize, beta, nstart, t1, t2, nsamp, tsamp, duration, nbnb, listnb, pdisp);
	time(&time2);
	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1));
	return 0;
}




/* gcc line:
## OPTIMIZED COMPILE - CHECK TIME ##

   gcc -o epidemics param.c auxiliary.c pathogens.c populations.c dispersal.c infection.c sampling.c sumstat.c inout.c epidemics.c -Wall -O3 -lgsl -lgslcblas -fopenmp

   ./epidemics

 
## FOR MEMORY LEAKS ## 
 
   gcc -o epidemics param.c auxiliary.c pathogens.c populations.c dispersal.c infection.c sampling.c sumstat.c inout.c epidemics.c -Wall -O0 -lgsl -lgslcblas 

 
   valgrind --leak-check=yes epidemics 
 
 
## FOR PROFILING ## 
   gcc -o epidemics param.c auxiliary.c pathogens.c populations.c dispersal.c infection.c sampling.c sumstat.c inout.c epidemics.c -Wall -O3 -pg -lgsl -lgslcblas 
 
   ./epidemics 
 
   gprof epidemics > epidemics-prof.txt 
 
   emacs epidemics-prof.txt & 
 
*/ 
 
