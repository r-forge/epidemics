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

struct pathogen ** get_pathogens(struct population *in){
	return in->pathogens;
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




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create new population */
struct population * create_population(int ns, int ni, int nr){
	int i;
	/* create pointer to population */
	struct population *out;
	out = (struct population *) calloc(1, sizeof(struct population));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_population]\nNo memory left for creating new population. Exiting.\n");
		exit(1);
	}

	/* create the content */
	out->orinsus = ns;
	out->nsus = ns-ni; /* remove susc. because of initial infections */
	out->ninf = ni;
	out->nrec = nr;
	out->ninfcum = ni;

	/* infected */

	out->pathogens = (struct pathogen **) calloc(ns, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_population]\nNo memory left for creating pathogen array in the population. Exiting.\n");
		exit(1);
	}

	for(i=0;i<ns;i++){
		(out->pathogens)[i] = create_pathogen();
		if(i<ni){
			(out->pathogens)[i]->age = 1; /* 'active' pathogen */
		} else {
			(out->pathogens)[i]->age = -1; /* 'neutralised' pathogen */
		}
	}

	return out;
}








/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free population */
void free_population(struct population *in){
	int i;
	for(i=0;i<get_orinsus(in);i++){
		if((in->pathogens)[i] != NULL) free_pathogen((in->pathogens)[i]);
	}

	free(in->pathogens);
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

/* Print population content */
void print_population(struct population *in, bool showGen){
	int i;
	printf("\nnb susceptible: %d", get_nsus(in));
	printf("\nnb infected: %d", get_ninf(in));
	printf("\nnb recovered: %d\n", get_nrec(in));
	if(showGen){
		for(i=0;i<get_orinsus(in);i++){
			if(!isNULL_pathogen(get_pathogens(in)[i])) print_pathogen(get_pathogens(in)[i]);
		}
		printf("\n");
	}
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
struct sample * draw_sample(struct population *in, struct param *par){
	int i, j, popSize=get_orinsus(in), n=par->n_sample, id, nIsolates=0;
	int *availIsolates; 

	/* create pointer to pathogens */
	struct sample *out;

	out = (struct sample *) calloc(1, sizeof(struct sample));

	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left to sample the population. Exiting.\n");
		exit(1);
	}

	/* allocate memory for pathogens */
	out->pathogens = (struct pathogen **) calloc(n, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left sample pathogens from the population. Exiting.\n");
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
		printf("\nPopulation without infections - sample will be empty.\n");
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
