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


struct population * get_populations(struct metapopulation *in){
	return in->populations;
}


int get_maxnpat(struct metapopulation *in){
	return in->maxnpat;
}


int get_npop(struct metapopulation *in){
	return in->npop;
}


int * get_popid(struct metapopulation *in){
	return in->popid;
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


int get_n(struct sample *in){
	return in->n;
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
struct metapopulation * create_metapopulation(int maxnpat, int nini, int npop, int nsus){
	int i;
	/* create pointer to metapopulation */
	struct metapopulation *out;
	out = (struct metapopulation *) calloc(1, sizeof(struct metapopulation));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating new metapopulation. Exiting.\n");
		exit(1);
	}


	/* allocate pathogen array */
	out->pathogens = (struct pathogen **) calloc(maxnpat, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating pathogen array in the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate population array */
	out->populations = (struct population **) calloc(1, sizeof(struct population *));
	if(out->populations == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating populations array in the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate population identifier array */
	out->popid = (int *) calloc(maxnpat, sizeof(int));
	if(out->popid == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating population identifier array in the metapopulation. Exiting.\n");
		exit(1);
	}


	/* fill in the pathogens and popid arrays */
	for(i=0;i<maxnpat;i++){
		(out->pathogens)[i] = create_pathogen();
		if(i<nini){ /* there are nini intial pathogens in the metapopulation */
			(out->pathogens)[i]->age = 1; /* 'active' pathogen */
			out->popid[i] = 0;

		} else {
			(out->pathogens)[i]->age = -1; /* 'neutralised' pathogen */
			out->popid[i] = -1;

		}
	}

	/* fill in the populations arrays */
	out->population[0] = create_population(nsus, nini, 0);
	for(i=1;i<npop;i++) {
		out->population[i] = create_population(nsus, 0, 0);
	}

	return out;
}






/* Create new population */
struct population create_population(int ns, int ni, int nr){
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
	free(in->popid);
	free(in);
}




/* Free metapopulation */
void free_population(struct population *in){
	free(in);
}




/* Free sample */
void free_sample(struct sample *in){
	if(in->pathogens != NULL) free(in->pathogens);
	free(in);
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

	/* general display */
	printf("\nnb of populations: %d", K);
	printf("\ntotal nb infected: %d", get_total_ninf(in));
	printf("\ntotal nb recovered: %d\n", get_total_nrec(in));

	/* display populations */
	for(k=0;k<K;k++){
		curPop = get_populations(in)[k];
		print_population(curPop);
		if(showGen){
			for(i=0;i<get_maxnpat(in);i++){
				if(!isNULL_pathogen(get_pathogens(curPop)[i]) && popid[i]==k) print_pathogen(get_pathogens(in)[i]);
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



/* Print sample content */
void print_sample(struct sample *in, bool showGen){
	int i;
	printf("\n%d pathogens", in->n);
	if(showGen){
		for(i=0;i<in->n;i++){
			//if(!isNULL_pathogen((in->pathogens)[i])) print_pathogen((in->pathogens)[i]);
			print_pathogen((in->pathogens)[i]);
		}
		printf("\n");
	}
}





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Get sample of isolates */
struct sample * draw_sample(struct metapopulation *in, struct param *par){
	int i, j, popSize=get_orinsus(in), n=par->n_sample, id, nIsolates=0;
	int *availIsolates; 

	/* create pointer to pathogens */
	struct sample *out;

	out = (struct sample *) calloc(1, sizeof(struct sample));

	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left to sample the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate memory for pathogens */
	out->pathogens = (struct pathogen **) calloc(n, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left sample pathogens from the metapopulation. Exiting.\n");
		exit(1);
	}

	/* get the number of isolates that can be sampled */
	for(i=0;i<popSize;i++){
		if(!isNULL_pathogen((in->pathogens)[i])) nIsolates++;
	}

	if(nIsolates != get_ninf(in)){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNumber of available isolates (%d) does not match number of infected (%d). Exiting.\n", nIsolates, get_ninf(in));
		exit(1);
	}

	/* escape if no isolate available */
	if(nIsolates < 1){
		printf("\nMetapopulation without infections - sample will be empty.\n");
		out->n = 0;
		out->pathogens = NULL;
		return out;
	}

	/* make vector of indices of available isolates */
	availIsolates = (int *) calloc(nIsolates, sizeof(int));
	if(availIsolates == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left to isolate available pathogens. Exiting.\n");
		exit(1);
	}

	j=0;
	for(i=0;i<nIsolates;i++){
		while(isNULL_pathogen((in->pathogens)[j])) j++;
		availIsolates[i] = j++;
	}

	/* choose from available pathogens */
	for(i=0;i<n;i++){
		id=gsl_rng_uniform_int(par->rng,nIsolates);
		(out->pathogens)[i] = (in->pathogens)[availIsolates[id]];
	}

	out->n = n;

	/* free local pointers */
	free(availIsolates);

	return out;
} /* end draw_sample */









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
