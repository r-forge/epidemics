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
#include "inout.h"


/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* seed new infections from a single pathogen */
void process_infection(struct pathogen * pat, struct metapopulation * metapop, struct dispmat *D, struct param * par){
	struct population * pop;
	int i, Nsus, Npop, Ninfcum, newpopid, nbnewinf;

	/* GENERATE ERROR IF PATHOGEN IS INACTIVATED */
	if(!isNULL_pathogen(pat) && get_age(pat) >= par->t1){ /* if infection is not a gost and age is OK */

		/* GET POPULATION INFORMATION */
		/* determine the pathogen's original population , Nsus, Ninfcum*/
		newpopid = disperse(pat, D, par);
		/* printf("\n new pop id %d", newpopid); */
		pop = get_populations(metapop)[newpopid]; /* with dispersal */
		/* pop = get_populations(metapop)[get_popid(pat)]; */ /* no dispersal*/
		Nsus=get_nsus(pop);
		Npop=get_popsize(pop);
		Ninfcum=get_total_ninfcum(metapop);

		/* DETERMINE NUMBER OF NEW INFECTIONS */
		nbnewinf = gsl_ran_binomial(par->rng, par->beta/((double) Npop), Nsus); /* ~ B(beta/N_t, S_t) */

		/* HANDLE GENOME REPLICATION */
		for(i=0;i<nbnewinf;i++){
			/* printf("\n## trying to write on pathogen %d", Ninfcum+i); */
			replicate(pat, (get_pathogens(metapop))[Ninfcum+i], par);
			(metapop->pathogens[Ninfcum+i])->popid = newpopid;
		}

		/* UPDATE NUMBER OF SUSCEPTIBLES AND INFECTED IN THE POPULATION */
		pop->nsus = pop->nsus - nbnewinf;
		pop->ninfcum = pop->ninfcum + nbnewinf;
		pop->ninf = pop->ninf + nbnewinf;
	}
} /* end process_infection */








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Function to be called from R */
void R_epidemics(int *seqLength, double *mutRate, int *npop, int *nHostPerPop, double *beta, int *nStart, int *t1, int *t2, int *Nsample, int *Tsample, int *duration, double *pdisp){
	int i, nstep=0, maxnpat;

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
	par = (struct param *) calloc(1, sizeof(struct param));
	par->L = *seqLength;
	par->mu = *mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = *npop;
	par->nsus = nHostPerPop;
	par->beta = *beta;
	par->nstart = *nStart;
	par->t1 = *t1;
	par->t2 = *t2;
	par->t_sample = Tsample;
	par->n_sample = *Nsample;
	par->duration = *duration;
	par->pdisp = pdisp;

	/* check/print parameters */
	check_param(par);
	/* printf("\n ! Beta version used !\n"); */
	print_param(par);

	/* dispersal matrix */
	struct dispmat *D;
	D = create_dispmat(par);
	/* printf("\ndispersal matrix:"); */
	/* print_dispmat(D); */

	/* group sizes */
	struct ts_groupsizes * grpsizes = create_ts_groupsizes(par);

	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);
	maxnpat = get_maxnpat(metapop);

	/* get sampling schemes (timestep+effectives) */
	translate_dates(par);
	struct table_int *tabdates = get_table_int(par->t_sample, par->n_sample);
	printf("\n\nsampling at timesteps:");
	print_table_int(tabdates);

	/* create sample */
	struct sample ** samplist = (struct sample **) calloc(tabdates->n, sizeof(struct sample *));
	struct sample *samp;
	int counter_sample = 0, tabidx;
	int firstActiveIdx=0; /* idx of first active pathogens */
	bool activeMark; /* TRUE if active pathogens have been reached */

	/* make metapopulation evolve */
	while(get_total_nsus(metapop)>0 && get_total_ninf(metapop)>0 && nstep<par->duration){
		nstep++;

		/* age metapopulation */
		age_metapopulation(metapop, par);

		/* handle replication for each infection */
		activeMark = FALSE;
		for(i=firstActiveIdx;i<get_total_ninfcum(metapop);i++){
			if(!activeMark && !isNULL_pathogen(get_pathogens(metapop)[i])) {
				firstActiveIdx = i;
				activeMark=TRUE;
			}
			process_infection(get_pathogens(metapop)[i], metapop, D, par);
		}

		/* for(i=0;i<maxnpat;i++){ */
		/* 	process_infection(get_pathogens(metapop)[i], metapop, D, par); */
		/* } */

		/* draw samples */
		if((tabidx = int_in_vec(nstep, tabdates->items, tabdates->n)) > -1){
			samplist[counter_sample++] = draw_sample(metapop, tabdates->times[tabidx], par);
		}

		fill_ts_groupsizes(grpsizes, metapop, nstep);

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
		printf("\n\nPrinting sample to file 'out-sample.txt'\n");
		write_sample(samp);

		/* free memory */
		free_sample(samp);

	}

	/* write group sizes to file */
	printf("\n\nPrinting group sizes to file 'out-popsize.txt'\n");
	write_ts_groupsizes(grpsizes);


	/* free memory */
	free_metapopulation(metapop);
	free_param(par);
	for(i=0;i<counter_sample;i++) free_sample(samplist[i]);
	free(samplist);
	free_table_int(tabdates);
	free_dispmat(D);
	free_ts_groupsizes(grpsizes);
}







/* Function to be called from R */
void R_monitor_epidemics(int *seqLength, double *mutRate, int *npop, int *nHostPerPop, double *beta, int *nStart, int *t1, int *t2, int *Nsample, int *duration, double *pdisp, int *minSize){
	int i, nstep=0, maxnpat;

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
	par = (struct param *) calloc(1, sizeof(struct param));
	par->L = *seqLength;
	par->mu = *mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = *npop;
	par->nsus = nHostPerPop;
	par->beta = *beta;
	par->nstart = *nStart;
	par->t1 = *t1;
	par->t2 = *t2;
	par->t_sample = NULL;
	par->n_sample = *Nsample;
	par->duration = *duration;
	par->pdisp = pdisp;

	/* check/print parameters */
	check_param(par);
	printf("\n ! Beta version used !\n");
	print_param(par);

	/* dispersal matrix */
	struct dispmat *D;
	D = create_dispmat(par);
	printf("\ndispersal matrix:");
	print_dispmat(D);

	/* group sizes */
	struct ts_groupsizes * grpsizes = create_ts_groupsizes(par);
	struct ts_sumstat * sumstats = create_ts_sumstat(par);

	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);
	maxnpat = get_maxnpat(metapop);

	/* memory allocations for sample and results */
	struct sample *samp;
	int firstActiveIdx=0; /* idx of first active pathogens */
	bool activeMark=FALSE; /* TRUE if active pathogens have been reached */

	/* make metapopulation evolve */
	while(get_total_nsus(metapop)>0 && get_total_ninf(metapop)>0 && nstep<par->duration){
		nstep++;

		/* age metapopulation */
		age_metapopulation(metapop, par);

		/* handle replication for each infection */
		for(i=firstActiveIdx;i<get_total_ninfcum(metapop);i++){
			if(!activeMark && !isNULL_pathogen(get_pathogens(metapop)[i])) {
				firstActiveIdx = i;
				activeMark=TRUE;
			}
			process_infection(get_pathogens(metapop)[i], metapop, D, par);
		}

		/* for(i=0;i<maxnpat;i++){ */
		/* 	process_infection(get_pathogens(metapop)[i], metapop, D, par); */
		/* } */

		/* draw sample */
		samp = draw_sample(metapop, par->n_sample, par);

		/* compute statistics */
		if(get_total_ninf(metapop)> *minSize) fill_ts_sumstat(sumstats, samp, nstep, par);

		/* get group sizes */
		fill_ts_groupsizes(grpsizes, metapop, nstep);
	}

	/* write group sizes to file */
	printf("\n\nPrinting results to file...");
	write_ts_groupsizes(grpsizes);
	write_ts_sumstat(sumstats);
	printf("done.\n\n");

	/* free memory */
	free_sample(samp);
	free_metapopulation(metapop);
	free_param(par);
	free_dispmat(D);
	free_ts_groupsizes(grpsizes);
	free_ts_sumstat(sumstats);
}








/* all-in-one function testing epidemics growth, summary statistics, etc. */
void test_epidemics(int seqLength, double mutRate, int npop, int *nHostPerPop, double beta, int nStart, int t1, int t2, int Nsample, int *Tsample, int duration, double *pdisp){
	int i, nstep=0, maxnpat;

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
	par = (struct param *) calloc(1, sizeof(struct param));
	par->L = seqLength;
	par->mu = mutRate;
	par->muL = par->mu * par->L;
	par->rng = rng;
	par->npop = npop;
	par->nsus = nHostPerPop;
	par->beta = beta;
	par->nstart = nStart;
	par->t1 = t1;
	par->t2 = t2;
	par->t_sample = Tsample;
	par->n_sample = Nsample;
	par->duration = duration;
	par->pdisp = pdisp;

	/* check/print parameters */
	check_param(par);
	print_param(par);

	/* dispersal matrix */
	struct dispmat *D;
	D = create_dispmat(par);
	printf("\ndispersal matrix:");
	print_dispmat(D);

	/* group sizes */
	struct ts_groupsizes * grpsizes = create_ts_groupsizes(par);


	/* initiate population */
	struct metapopulation * metapop;
	metapop = create_metapopulation(par);
	maxnpat = get_maxnpat(metapop);


	/* get sampling schemes (timestep+effectives) */
	translate_dates(par);
	struct table_int *tabdates = get_table_int(par->t_sample, par->n_sample);
	printf("\n\nsampling at timesteps:");
	print_table_int(tabdates);

	/* create sample */
	struct sample ** samplist = (struct sample **) calloc(tabdates->n, sizeof(struct sample *));
	struct sample *samp;
	int counter_sample = 0, tabidx, firstActiveIdx=0;
	bool activeMark=FALSE;

	/* make metapopulation evolve */
	while(get_total_nsus(metapop)>0 && get_total_ninf(metapop)>0 && nstep<par->duration){
		nstep++;

		/* printf("\nmetapop check, time step %d", nstep); */
		/* for(i=0;i<maxnpat;i++) if(metapop->pathogens[i]==NULL) printf("\npathogen %d is NULL", i); */
		/* printf("...ok"); */

		/* age metapopulation */
		age_metapopulation(metapop, par);

		/* handle replication for each infection */
		for(i=firstActiveIdx;i<get_total_ninfcum(metapop);i++){
			if(!activeMark && !isNULL_pathogen(get_pathogens(metapop)[i])) {
				firstActiveIdx = i;
				activeMark=TRUE;
			}
			process_infection(get_pathogens(metapop)[i], metapop, D, par);
		}

		/* draw samples */
		if((tabidx = int_in_vec(nstep, tabdates->items, tabdates->n)) > -1){
			samplist[counter_sample++] = draw_sample(metapop, tabdates->times[tabidx], par);
		}

		fill_ts_groupsizes(grpsizes, metapop, nstep);

	}

	/* we stopped after 'nstep' steps */
	if(nstep < par->duration){
		printf("\nEpidemics ended at time %d, before last sampling time (%d).\n", nstep, par->duration);
	} else {

		printf("\n\n-- FINAL METAPOPULATION --");
		print_metapopulation(metapop, TRUE);

		/* test samples */
		for(i=0;i<tabdates->n;i++) {
			printf("\nsample %d\n", i);
			print_sample(samplist[i], TRUE);
		}
		samp = merge_samples(samplist, tabdates->n, par) ;
		print_sample(samp, TRUE);

		/* test allele listing */
		struct snplist *snpbilan;
		snpbilan = list_snps(samp, par);
		print_snplist(snpbilan);

		/* test allele frequencies */
		struct allfreq *freq;
		freq = get_frequencies(samp, par);
		print_allfreq(freq);

		/* test Hs*/
		double Hs = hs(samp,par);
		printf("\nHs = %0.3f\n", Hs);

		/* test Hs full genome */
		Hs = hs_full_genome(samp,par);
		printf("\nHs (full genome) = %0.5f\n", Hs);

		/* test nb of snps */
		int nball = nb_snps(samp,par);
		printf("\nnumber of SNPs = %d\n", nball);

		/* test mean nb of snps */
		double temp = mean_nb_snps(samp);
		printf("\nmean number of SNPs = %.2f\n", temp);

		/* test var nb of snps */
		temp = var_nb_snps(samp);
		printf("\nvariance of number of alleles = %.2f\n", temp);

		/* test pairwise distances */
		struct distmat_int *mat = pairwise_dist(samp, par);
		print_distmat_int(mat);

		/* test mean pairwise distances */
		temp = mean_pairwise_dist(samp,par);
		printf("\nmean pairwise distance: %.2f", temp);

		/* test variance of pairwise distances */
		temp = var_pairwise_dist(samp,par);
		printf("\nvar pairwise distance: %.2f", temp);

		/* test Fst */
		temp = fst(samp,par);
		printf("\nfst: %.2f", temp);


		printf("\n\n");

		/* free memory */
		free_sample(samp);
		free_snplist(snpbilan);
		free_allfreq(freq);
		free_distmat_int(mat);

	}

	/* write group sizes to file */
	printf("\n\nPrinting group sizes to file 'out-popsize.txt'");
	write_ts_groupsizes(grpsizes);

	/* free memory */
	free_metapopulation(metapop);
	free_param(par);
	for(i=0;i<counter_sample;i++) free_sample(samplist[i]);
	free(samplist);
	free_table_int(tabdates);
	free_dispmat(D);
	free_ts_groupsizes(grpsizes);
}




int main(){
/* args: (int seqLength, double mutRate, int npop, int nHostPerPop, double beta, int nStart, int t1, int t2,int Tsample, int Nsample) */
	double mu=1e-4, beta=1.1, pdisp[1]={1.0}; //pdisp[9] = {0.5,0.25,0.25,0.0,0.5,0.5,0.0,0.0,1.0};
	time_t time1,time2;
	int genoL=1e4, duration=50, npop=1, nstart=10, t1=1, t2=2, nsamp=10;
	int tsamp[10] = {1,1,1,1,1,1,1,0,0,0}, popsize[1]={1e6}; //popsize[3]={10,1,1};

	time(&time1);
	test_epidemics(genoL, mu, npop, popsize, beta, nstart, t1, t2, nsamp, tsamp, duration, pdisp);
	time(&time2);
	printf("\ntime ellapsed: %d seconds \n", (int) (time2-time1));
	return 0;
}
