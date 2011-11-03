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
   =================
   === ACCESSORS ===
   =================
*/

struct pathogen ** get_pathogens(struct metapopulation *in){
	return in->pathogens;
}


struct population ** get_populations(struct metapopulation *in){
	return in->populations;
}


int get_maxnpat(struct metapopulation *in){
	return in->maxnpat;
}


int get_npop(struct metapopulation *in){
	return in->npop;
}


int get_nsus(struct population *in){
	return in->nsus;
}


int get_ninf(struct population *in){
	return in->ninf;
}


int get_nrec(struct population *in){
	return in->nrec;
}


int get_ninfcum(struct population *in){
	return in->ninfcum;
}


int get_orinsus(struct population *in){
	return in->orinsus;
}


int get_popsize(struct population *in){
	return in->nsus + in->ninf + in->nrec;
}


int get_total_nsus(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_nsus(get_populations(in)[i]);
	}
	return out;
}


int get_total_ninf(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_ninf(get_populations(in)[i]);
	}
	return out;
}


int get_total_nrec(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_nrec(get_populations(in)[i]);
	}
	return out;
}


int get_total_ninfcum(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_ninfcum(get_populations(in)[i]);
	}
	return out;
}


int get_total_orinsus(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_orinsus(get_populations(in)[i]);
	}
	return out;
}




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create new metapopulation */
struct metapopulation * create_metapopulation(struct param *par){
	int i, npop, nsus, maxnpat, nini ;
	npop = par->npop;
	nsus = par->nsus;
	maxnpat = npop*nsus;
	nini = par->nstart;

	/* create pointer to metapopulation */
	struct metapopulation *out;
	out = (struct metapopulation *) calloc(1, sizeof(struct metapopulation));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating new metapopulation. Exiting.\n");
		exit(1);
	}

	/* set maxnpat and npop */
	out->maxnpat = maxnpat;
	out->npop = npop;

	/* allocate pathogen array */
	out->pathogens = (struct pathogen **) calloc(maxnpat, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating pathogen array in the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate population array */
	out->populations = (struct population **) calloc(npop, sizeof(struct population *));
	if(out->populations == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating populations array in the metapopulation. Exiting.\n");
		exit(1);
	}


	/* fill in the pathogens and popid arrays */
	for(i=0;i<maxnpat;i++){
		(out->pathogens)[i] = create_pathogen();
		if(i<nini){ /* there are nini intial pathogens in the metapopulation */
			(out->pathogens[i])->age = 1; /* 'active' pathogen */
			(out->pathogens[i])->popid = 0;

		} else {
			(out->pathogens[i])->age = -1; /* 'neutralised' pathogen */
			(out->pathogens[i])->popid = -1;
		}
	}

	/* fill in the populations arrays */
	out->populations[0] = create_population(nsus, nini, 0);
	for(i=1;i<npop;i++) {
		out->populations[i] = create_population(nsus, 0, 0);
	}

	return out;
}






/* Create new population */
struct population * create_population(int ns, int ni, int nr){
	struct population *out;
	out = (struct population *) calloc(1, sizeof(struct population));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_population]\nNo memory left for creating new population. Exiting.\n");
		exit(1);
	}

	/* fill the content */
	out->orinsus = ns;
	out->nsus = ns-ni; /* remove susc. because of initial infections */
	out->ninf = ni;
	out->nrec = nr;
	out->ninfcum = ni;

	return out;
}




/* Create ts_groupsizes */
struct ts_groupsizes * create_ts_groupsizes(int nsteps){
	struct ts_groupsizes * out = (struct ts_groupsizes *) calloc(1, sizeof(struct ts_groupsizes));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_ts_groupsizes]\nNo memory left for storing group sizes. Exiting.\n");
		exit(1);
	}

	out->nsus = (int *) calloc(nsteps, sizeof(int));
	out->ninf = (int *) calloc(nsteps, sizeof(int));
	out->nrec = (int *) calloc(nsteps, sizeof(int));
	out->ninfcum = (int *) calloc(nsteps, sizeof(int));
	out->length = nsteps;

	if(out->nsus==NULL || out->ninf==NULL || out->nrec==NULL || out->ninfcum==NULL){
		fprintf(stderr, "\n[in: population.c->create_ts_groupsizes]\nNo memory left for storing group sizes. Exiting.\n");
		exit(1);
	}

	return out;
}




/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free metapopulation */
void free_metapopulation(struct metapopulation *in){
	int i, npat=get_maxnpat(in), npop=get_npop(in);
	for(i=0;i<npat;i++){
		if(in->pathogens[i] != NULL) free_pathogen((in->pathogens)[i]);
	}

	for(i=0;i<npop;i++){
		if(in->populations[i] != NULL) free_population(in->populations[i]);
	}

	free(in->pathogens);
	free(in->populations);
	free(in);
}




/* Free metapopulation */
void free_population(struct population *in){
	free(in);
}



/* Free ts_groupsizes */
void free_ts_groupsizes(struct ts_groupsizes *in){
	if(in!=NULL){
		free(in->nsus);
		free(in->ninf);
		free(in->nrec);
		free(in->ninfcum);
		free(in);
	}
}





/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Print metapopulation content */
void print_metapopulation(struct metapopulation *in, bool showGen){
	int i, k, K=get_npop(in);
	struct population *curPop;
	struct pathogen * curpat;

	/* display general info */
	printf("\nnb of populations: %d", K);
	printf("\ntotal nb susceptible: %d", get_total_nsus(in));
	printf("\ntotal nb infected: %d", get_total_ninf(in));
	printf("\ntotal nb recovered: %d", get_total_nrec(in));
	printf("\ntotal population size: %d\n", get_total_nsus(in)+get_total_ninf(in)+get_total_nrec(in));

	/* display populations */
	for(k=0;k<K;k++){
		curPop = get_populations(in)[k];
		printf("\npopulation %d", k);
		print_population(curPop);
		if(showGen){
			for(i=0;i<get_maxnpat(in);i++){
				curpat = get_pathogens(in)[i];
				if(!isNULL_pathogen(curpat) && get_popid(curpat)==k) print_pathogen(curpat);
			}
			printf("\n");
		}
	}
}



/* Print metapopulation content */
void print_population(struct population *in){
	printf("\nnb susceptible: %d", get_nsus(in));
	printf("\nnb infected: %d", get_ninf(in));
	printf("\nnb recovered: %d\n", get_nrec(in));
}






/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void age_metapopulation(struct metapopulation * metapop, struct param * par){
	struct pathogen *ppat;
	int i, maxnpat = get_maxnpat(metapop);

	/* pathogens ages */
	for(i=0;i<maxnpat;i++){
		ppat = (metapop->pathogens)[i]; /* to make code more readable*/
		if(!isNULL_pathogen(ppat)){ /* if pathogen exists */
			ppat->age = ppat->age+1; /* get older */
			if(get_age(ppat) > par->t2) { /* die if you must */
				ppat->age = -1; /* inactivate pathogen */

				/* update nrec and ninf in corresponding population */
				(metapop->populations[ppat->popid])->nrec = (metapop->populations[ppat->popid])->nrec + 1;
				(metapop->populations[ppat->popid])->ninf = (metapop->populations[ppat->popid])->ninf - 1;
				ppat->popid = -1; /* inactivate pathogen */
			}
		}
	}
} /* end age_metapopulation */




/* keep track of group sizes */
void fill_ts_groupsizes(struct ts_groupsizes *in, struct metapopulation *metapop, int step){
	if(step>in->length){
		fprintf(stderr, "\n[in: population.c->fill_ts_groupsizes]\n. ts_groupsizes object is not long enough to store output of step %d. Exiting.\n", step);
		exit(1);
	}

	in->nsus[step-1] = get_total_nsus(metapop);
	in->ninf[step-1] = get_total_ninf(metapop);
	in->nrec[step-1] = get_total_nrec(metapop);
	in->ninfcum[step-1] = get_total_ninfcum(metapop);
}








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


/* 	/\* simulation parameters *\/ */
/* 	/\* struct param * par; *\/ */
/* 	/\* par = (struct param *) calloc(1, sizeof(struct param)); *\/ */
/* 	/\* par->L = 100; *\/ */
/* 	/\* par->mu = 0.01; *\/ */
/* 	/\* par->muL = par->mu * par->L; *\/ */
/* 	/\* par->rng = rng; *\/ */


/* 	struct population * pop; */

/* 	pop = create_population(1000,10,0); */

/* 	print_population(pop); */

/* 	/\* free memory *\/ */
/* 	free_population(pop); */
/* 	gsl_rng_free(rng); */

/* 	return 0; */
/* } */
