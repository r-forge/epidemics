
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

/* FOR POPULATIONS */
int get_nsus(struct population *in){
	return in->nsus;
}


int get_nexp(struct population *in){
	return in->nexp;
}


int get_ninf(struct population *in){
	return in->ninf;
}


int get_nrec(struct population *in){
	return in->nrec;
}


int get_nexpcum(struct population *in){
	return in->nexpcum;
}


int get_popsize(struct population *in){
	return in->popsize;
}


int get_popid(struct population *in){
	return in->popid;
}


int get_idfirstinfectious(struct population *in){
	return in->idfirstinfectious;
}


int get_idlastinfectious(struct population *in){
	return in->idlastinfectious;
}


struct pathogen ** get_pathogens(struct population *in){
	return in->pathogens;
}




/* FOR METAPOPULATIONS */
struct population ** get_populations(struct metapopulation *in){
	return in->populations;
}


int get_npop(struct metapopulation *in){
	return in->npop;
}


int * get_popsizes(struct metapopulation *in){
	return in->popsizes;
}


int get_total_nsus(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_nsus(get_populations(in)[i]);
	}
	return out;
}


int get_total_nexp(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_nexp(get_populations(in)[i]);
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


int get_total_nexpcum(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_nexpcum(get_populations(in)[i]);
	}
	return out;
}


int get_total_popsize(struct metapopulation *in){
	int i, k=get_npop(in), out=0;
	for(i=0;i<k;i++) {
		out += get_popsize(get_populations(in)[i]);
	}
	return out;
}




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create new population */
struct population * create_population(int popsize, int nini, int popid){
	int i;

	/* allocate output */
	struct population *out;
	out = (struct population *) malloc(sizeof(struct population));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_population]\nNo memory left for creating new population. Exiting.\n");
		exit(1);
	}

	/* fill the content */
	out->popsize = popsize;
	out->nexp = nini;
	out->ninf = 0;
	out->nsus = popsize-nini; /* remove susc. because of initial infections */
	out->nrec = 0;
	out->nexpcum = nini;
	out->popid = popid;

	/* allocate pathogen array */
	out->pathogens = (struct pathogen **) malloc(popsize * sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_population]\nNo memory left for creating pathogen array in the population. Exiting.\n");
		exit(1);
	}

	/* fill in the pathogens array */
	for(i=0;i<popsize;i++){
		(out->pathogens)[i] = create_pathogen();
		if(i<nini){ /* there are nini intial pathogens in the metapopulation */
			(out->pathogens[i])->age = 0; /* 'active' pathogen */
		} else {
			(out->pathogens[i])->age = -1; /* 'neutralised' pathogen */
		}
	}

	return out;
}




/* Create new metapopulation */
struct metapopulation * create_metapopulation(struct param *par){
	int i, nini = par->nstart;

	/* allocate output */
	struct metapopulation *out;
	out = (struct metapopulation *) malloc(sizeof(struct metapopulation));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating new metapopulation. Exiting.\n");
		exit(1);
	}

	/* set content */
	out->npop = par->npop;
	out->popsizes = par->popsizes;

	/* allocate population array */
	out->populations = (struct population **) malloc(out->npop * sizeof(struct population *));
	if(out->populations == NULL){
		fprintf(stderr, "\n[in: population.c->create_metapopulation]\nNo memory left for creating populations array in the metapopulation. Exiting.\n");
		exit(1);
	}

	out->populations[0] = create_population(out->popsizes[0], nini, 0); /* pop 0 has some active pathogens */
	for(i=1;i<out->npop;i++) {
		out->populations[i] = create_population(out->popsizes[i], 0, i);
	}

	return out;
}










/* Create ts_groupsizes */
struct ts_groupsizes * create_ts_groupsizes(struct param * par){
	int nsteps = par->duration;
	struct ts_groupsizes * out = (struct ts_groupsizes *) malloc(sizeof(struct ts_groupsizes));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_ts_groupsizes]\nNo memory left for storing group sizes. Exiting.\n");
		exit(1);
	}

	/* calloc is needed there */
	out->nsus = (int *) calloc(nsteps, sizeof(int));
	out->nexp = (int *) calloc(nsteps, sizeof(int));
	out->ninf = (int *) calloc(nsteps, sizeof(int));
	out->nrec = (int *) calloc(nsteps, sizeof(int));
	out->nexpcum = (int *) calloc(nsteps, sizeof(int));
	out->length = nsteps;

	if(out->nsus==NULL || out->nexp==NULL || out->ninf==NULL || out->nrec==NULL || out->nexpcum==NULL){
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
void free_population(struct population *in){
	int i;

	for(i=0;i<in->popsize;i++){
		if(in->pathogens[i] != NULL) free_pathogen((in->pathogens)[i]);
	}

	free(in->pathogens);
	free(in);
}



/* Free metapopulation */
void free_metapopulation(struct metapopulation *in){
	int i, npop=get_npop(in);

	for(i=0;i<npop;i++){
		if(in->populations[i] != NULL) free_population(in->populations[i]);
	}

	free(in->populations);
	free(in);
}






/* Free ts_groupsizes */
void free_ts_groupsizes(struct ts_groupsizes *in){
	if(in!=NULL){
		free(in->nsus);
		free(in->nexp);
		free(in->ninf);
		free(in->nrec);
		free(in->nexpcum);
		free(in);
	}
}





/*
  ===========================
  === AUXILIARY FUNCTIONS ===
  ===========================
*/
/* PRINT POPULATION CONTENT */
void print_population(struct population *in, bool showPat){
	int i, nrec=get_nrec(in), nexpcum=get_nexpcum(in);

	printf("\npopulation %d", get_popid(in));
	printf("\nnb susceptible: %d", get_nsus(in));
	printf("\nnb exposed: %d", get_nexp(in));
	printf("\nnb infectious: %d", get_ninf(in));
	printf("\nnb recovered: %d", get_nrec(in));
	printf("\npathogens:");
	if(showPat){
		for(i=nrec;i<nexpcum;i++){
			print_pathogen(get_pathogens(in)[i]);
		}
	}
	printf("\n");
}



/* PRINT METAPOPULATION CONTENT */
void print_metapopulation(struct metapopulation *in, bool showPat){
	int i, npop=get_npop(in);
	struct population *curPop;

	/* display general info */
	printf("\nnb of populations: %d", npop);
	printf("\npopulation sizes: ");
	for(i=0;i<get_npop(in);i++) printf("%d ", get_popsizes(in)[i]);
	printf("\ntotal nb susceptible: %d", get_total_nsus(in));
	printf("\ntotal nb exposed: %d", get_total_nexp(in));
	printf("\ntotal nb infectious: %d", get_total_ninf(in));
	printf("\ntotal nb recovered: %d", get_total_nrec(in));
	printf("\ntotal population size: %d\n", get_total_popsize(in));

	/* display populations */
	for(i=0;i<npop;i++){
		curPop = get_populations(in)[i];
		print_population(curPop,showPat);
	}
	printf("\n");
}









/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
void age_population(struct population * in, struct param *par){
	int i, nrec=get_nrec(in), nexpcum=get_nexpcum(in), nbnewinf=0, nbnewrec=0;
	struct pathogen *ppat;

	for(i=nrec;i<nexpcum;i++){
		ppat = get_pathogens(in)[i];
		if(!isNULL_pathogen(get_pathogens(in)[i])){ /* if pathogen is active */
			ppat->age = ppat->age + 1; /* get older */
			if(get_age(ppat) == par->t1){ /* becomes infectious this time step */
				nbnewinf++;
			}
			if(get_age(ppat) >= par->t2) { /* die if you must! */
				ppat->age = -1; /* inactivate pathogen */
				nbnewrec++;
			}
		}
	}

	/* update nexp, ninf, nrec in corresponding population */
	in->nexp = in->nexp - nbnewinf;
	in->ninf = in->ninf + nbnewinf - nbnewrec;
	in->nrec = in->nrec + nbnewrec;

} /* end age_population */





/* AGE METATPOPULATION */
void age_metapopulation(struct metapopulation * in, struct param * par){
	int i, npop=get_npop(in);

	/* age each population */
	for(i=0;i<npop;i++){
		age_population(get_populations(in)[i], par);
	}
} /* end age_metapopulation */





/* keep track of group sizes */
void fill_ts_groupsizes(struct ts_groupsizes *in, struct metapopulation *metapop, int step){
	if(step>in->length){
		fprintf(stderr, "\n[in: population.c->fill_ts_groupsizes]\n. ts_groupsizes object is not long enough to store output of step %d. Exiting.\n", step);
		exit(1);
	}

	in->nsus[step-1] = get_total_nsus(metapop);
	in->nexp[step-1] = get_total_nexp(metapop);
	in->ninf[step-1] = get_total_ninf(metapop);
	in->nrec[step-1] = get_total_nrec(metapop);
	in->nexpcum[step-1] = get_total_nexpcum(metapop);
}




/* /\* FIND INDEX OF THE FIRST ACTIVE PATHOGEN IN THE PATHOGEN ARRAY *\/ */
/* void update_first_active_pathogen(struct population *in, struct param *par){ */
/* 	int out = get_nrec(in), max=get_popsize(in); */
/* 	while(out < max && !is_infectious(get_pathogens(in)[out], par)) out++; */
/* 	if(out == max) return -1; */
/* 	in->idfirstinfectious = out; */
/* } */




/* /\* FIND INDEX OF THE LAST ACTIVE PATHOGEN IN THE PATHOGEN ARRAY *\/ */
/* void update_last_active_pathogen(struct population *in, struct param *par){ */
/* 	int max=get_popsize(in), out = find_id_first_active_pathogen(in, par); */
/* 	if(out < 0) return -1; */
/* 	while(out < max && is_infectious(get_pathogens(in)[out], par)) out++; */
/* 	in->idlastinfectious = out -1; */
/* } */




/* SELECT A RANDOM INFECTIOUS PATHOGEN FROM THE POPULATION */
struct pathogen * select_random_infectious_pathogen(struct population *in, struct param *par){
	/* int id; */
	/* if(in->ninf < 1) return NULL; */
	/* if(in->ninf==1) return get_pathogens(in)[in->nrec]; /\* gsl_rng_unif does not like a range of 0 *\/ */
	/* id = in->nrec + gsl_rng_uniform_int(par->rng, in->ninf); */
	/* return get_pathogens(in)[id]; */
	
	/* if(in->ninf < 1) return NULL; /\* should never happen *\/ */
	if(in->ninf==1) return get_pathogens(in)[in->nrec]; /* gsl_rng_unif does not like a range of 0 */
	return  get_pathogens(in)[in->nrec + gsl_rng_uniform_int(par->rng, in->ninf)];
}





/* SELECT A RANDOM PATHOGEN (host exp or inf) FROM THE POPULATION */
struct pathogen * select_random_pathogen(struct population *in, struct param *par){
	int id, nbavail=in->ninf+in->nexp;
	/* printf("\nfirst pathogen: %d     last pathogen: %d ", in->nrec, in->nrec + nbavail - 1); */
	if(nbavail < 1) return NULL;
	if(nbavail == 1) return get_pathogens(in)[in->nrec]; /* gsl_rng_unif does not like a range of 0 */
	id = in->nrec + gsl_rng_uniform_int(par->rng, nbavail);
	return get_pathogens(in)[id];
}



/* gcc line:

   gcc -o populations param.c auxiliary.c pathogens.c populations.c -Wall -O0 -lgsl -lgslcblas

   valgrind --leak-check=yes populations

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
/* 	int i; */

/* 	/\* simulation parameters *\/ */
/* 	struct param * par; */
/* 	par = (struct param *) calloc(1, sizeof(struct param)); */
/* 	par->rng = rng; */
/* 	par->npop = 3; */
/* 	int popsizes[3] = {1000,200,300}; */
/* 	par->popsizes = popsizes; */
/* 	par->nstart = 10; */
/* 	par->t1 = 1; */
/* 	par->t2 = 3; */

/* 	/\* TRY POPULATION *\/ */
/* 	struct population * pop = create_population(1000,10,69); */
/* 	printf("\nPOPULATION"); */
/* 	print_population(pop, TRUE); */

/* 	/\* TRY METAPOPULATION *\/ */
/* 	struct metapopulation * metapop = create_metapopulation(par); */
/* 	printf("\n## METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */

/* 	/\* TRY SELECT RANDOM INFECTIOUS PATHOGEN *\/ */
/* 	struct pathogen *ppat; */
/* 	ppat = select_random_infectious_pathogen(get_populations(metapop)[0], par); */
/* 	printf("\nselected pathogen:"); */
/* 	if(ppat != NULL) print_pathogen(ppat); else printf(" NULL\n"); */

/* 	/\* TRY SELECT RANDOM PATHOGEN *\/ */
/* 	ppat = select_random_pathogen(get_populations(metapop)[0], par); */
/* 	printf("\nselected pathogen:"); */
/* 	if(ppat != NULL) print_pathogen(ppat); else printf(" NULL\n"); */

/* 	/\* TRY AGEING *\/ */
/* 	age_metapopulation(metapop, par); */
/* 	printf("\n## AGED METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */

/* 	printf("\n10 random infectious pathogens"); */
/* 	for(i=0;i<10;i++){ */
/* 		ppat = select_random_infectious_pathogen(get_populations(metapop)[0], par); */
/* 		printf("\nselected pathogen:"); */
/* 		if(ppat != NULL) print_pathogen(ppat); else printf(" NULL\n"); */
/* 		printf("\npathogen address: %d", ppat); */
/*       	} */

/* 	printf("\n10 random pathogens"); */
/* 	for(i=0;i<10;i++){ */
/* 		ppat = select_random_pathogen(get_populations(metapop)[0], par); */
/* 		printf("\nselected pathogen:"); */
/* 		if(ppat != NULL) print_pathogen(ppat); else printf(" NULL\n"); */
/* 		printf("\npathogen address: %d", ppat); */
/* 	} */


/* 	age_metapopulation(metapop, par); */
/* 	printf("\n## AGEDx2 METAPOPULATION ##"); */
/* 	print_metapopulation(metapop, TRUE); */

/* 	/\* free memory *\/ */
/* 	free_population(pop); */
/* 	free_metapopulation(metapop); */
/* 	free(par); */
/* 	gsl_rng_free(rng); */

/* 	return 0; */
/* } */
