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
#include "sumstat.h"



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* seed new infections from a single pathogen */
void process_infection(struct pathogen * pat, struct metapopulation * metapop, struct population *pop, int popid, struct param * par){
	int i, nbnewinf=0,  Nsus=get_nsus(pop), Ninfcum=get_ninfcum(pop);

	if(!isNULL_pathogen(pat)){ /* if infection is not a gost */
		/* determine the number of descendents */
		if(get_age(pat) >= par->t1){
			nbnewinf = gsl_ran_poisson(par->rng, par->R);

			/* adjust number of new infections to number of susceptibles */
			if(nbnewinf > Nsus) nbnewinf =  Nsus;
		}

		if(nbnewinf>0){
			/* for each new infection, add new pathogen */
			for(i=Ninfcum;i<(Ninfcum+nbnewinf);i++){
 				/* printf("\n## trying to write on pathogen %d", i); */
				replicate(pat, (get_pathogens(metapop))[i], par);
				metapop->popid[i] = popid;
			}

			/* update number of susceptibles and infected */
			pop->nsus = pop->nsus - nbnewinf;
			pop->ninfcum = pop->ninfcum + nbnewinf;
			pop->ninf = pop->ninf + nbnewinf;
		}
	}
} /* end process_infection */








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void run_epidemics(int seqLength, double mutRate, int npop, int nHostPerPop, double incid, int nStart, int t1, int t2,int Tsample, int Nsample){
	int i, nstep=0, curPopId;

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
	par->R = incid;
	par->nstart = nStart;
	par->t1 = t1;
	par->t2 = t2;
	par->t_sample = Tsample;
	par->n_sample = Nsample;


	/* check/print parameters */
	check_param(par);
	print_param(par);

	/* initiate population */
	struct metapopulation * metapop;
	struct sample * samp;

	metapop = create_metapopulation(par);
	print_metapopulation(metapop, FALSE);

	/* make metapopulation evolve */
	printf("\ntotal nsus: %d", get_total_nsus(metapop));
	printf("\ntotal ninf: %d", get_total_ninf(metapop));
	/* while(get_total_nsus(metapop)>0 && get_total_ninf(metapop)>0 && nstep<par->t_sample){ */
	/* 	nstep++; */

	/* 	/\* handle replication for each infection *\/ */
	/* 	for(i=0;i<get_maxnpat(metapop);i++){ */
	/* 		curPopId = get_popid(metapop)[i]; */
	/* 		process_infection(get_pathogens(metapop)[i], metapop, get_populations(metapop)[curPopId], curPopId, par); */
	/* 	} */

	/* 	/\* age metapopulation *\/ */
	/* 	age_metapopulation(metapop, par); */
	/* } */

	/* /\* we stopped after 'nstep' steps *\/ */
	/* if(nstep < par->t_sample){ */
	/* 	printf("\nEpidemics ended at time %d, before sampling time (%d).\n", nstep, par->t_sample); */
	/* } */


	printf("\n\n-- FINAL POPULATION --");
	print_metapopulation(metapop, FALSE);

	/* test sampling */
	samp = draw_sample(metapop, par);
	printf("\n\n-- SAMPLE --");
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
	printf("\nnumber of alleles = %d\n", nball);

	/* test mean nb of snps */
	double temp = mean_nb_snps(samp);
	printf("\nmean number of alleles = %.2f\n", temp);

	/* test var nb of snps */
	temp = var_nb_snps(samp);
	printf("\nvariance of number of alleles = %.2f\n", temp);

	printf("\n");

	/* free memory */
	free_metapopulation(metapop);
	free_param(par);
	free_sample(samp);
	free_snplist(snpbilan);
	free_allfreq(freq);
}




int main(){
/* args: (int seqLength, double mutRate, int npop, int nHostPerPop, double incid, int nStart, int t1, int t2,int Tsample, int Nsample) */
	run_epidemics(1e4, 1e-4, 3, 1000, 1.05, 10, 1,1, 5, 10);

	return 0;
}
