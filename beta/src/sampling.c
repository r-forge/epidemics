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
#include "dispersal.h"
#include "infection.h"
#include "sampling.h"


/*
   =================
   === ACCESSORS ===
   =================
*/


int get_n(struct sample *in){
	return in->n;
}



/* get nb of populations in a sample */
int get_npop_samp(struct sample *in){
	int i, n=get_n(in), *pool, npop;

	/* ENUMERATE NB OF UNIQUE ITEMS */
	/* create pool of unique items */
	pool = (int *) malloc(n * sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: sampling.c->get_npop]\nNo memory left to sample per population. Exiting.\n");
		exit(1);
	}

	/* list and count pop occurences */
	npop = 0;
	for(i=0;i<n;i++){
		if(int_in_vec(in->popid[i], pool, npop) < 0){
			pool[npop++] = in->popid[i];
		}
	}

	/* free memory and return */
	free(pool);
	return(npop);
}






/*
   ===================
   === CONTRUCTORS ===
   ===================
*/

struct sample * create_sample(int n){
	struct sample *out = calloc(1, sizeof(struct sample));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_sample]\nNo memory left to sample the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate memory for pathogens */
	out->pathogens = (struct pathogen **) calloc(n, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_sample]\nNo memory left to sample pathogens from the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate memory for popid */
	out->popid = (int *) calloc(n, sizeof(int));
	if(out->popid == NULL){
		fprintf(stderr, "\n[in: population.c->create_sample]\nNo memory left to sample pathogens from the metapopulation. Exiting.\n");
		exit(1);
	}

	out->n = n;
	return out;

}



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free sample */
void free_sample(struct sample *in){
	int i, n=get_n(in);
	if(in->pathogens != NULL) {
		for(i=0;i<n;i++) {
			if(in->pathogens[i] != NULL) free_pathogen(in->pathogens[i]);
		}
		free(in->pathogens);
	}
	if(in->popid != NULL) free(in->popid);
	free(in);
}







/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Print sample content */
void print_sample(struct sample *in, bool showGen){
	int i;
	printf("\n - sample of pathogens -");
	printf("\n populations:");
	for(i=0;i<in->n;i++) printf("%d ", in->popid[i]);
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

/* GET SAMPLE OF ISOLATES */
/* Isolates are COPIED, so that any modification of the sample does not alter */
/* the metapopulation. */
struct sample * draw_sample(struct metapopulation *in, int n, struct param *par){
	int i, j, *nIsolatesPerPop, count;
	double *nAvailPerPop;
	struct pathogen *ppat;

	/* create pointer to pathogens */
	struct sample *out=create_sample(n);

	/* escape if no isolate available */
	if(get_total_ninf(in) + get_total_nexp(in) < 1){
		printf("\nMetapopulation without infections - sample will be empty.\n");
		out->n = 0;
		out->pathogens = NULL;
		out->popid = NULL;
		return out;
	}

	/* get nb of isolates available in each population */
	nAvailPerPop = (double *) calloc(get_npop(in), sizeof(double));
	if(nAvailPerPop == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left to isolate available pathogens. Exiting.\n");
		exit(1);
	}

	for(j=0;j<get_npop(in);j++){
		nAvailPerPop[j] = (double) get_nexp(get_populations(in)[j]) + get_ninf(get_populations(in)[j]);
	}

	/* get nb of isolates sampled in each population*/
	nIsolatesPerPop = (int *) calloc(get_npop(in), sizeof(int));
	if(nIsolatesPerPop == NULL){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNo memory left to isolate available pathogens. Exiting.\n");
		exit(1);
	}

	gsl_ran_multinomial(par->rng, get_npop(in), n, nAvailPerPop, (unsigned int *) nIsolatesPerPop);

	/* fill in the sample pathogens */
	count = 0;
	for(j=0;j<get_npop(in);j++){ /* for each population */
		for(i=0;i<nIsolatesPerPop[j];i++){
			ppat = select_random_pathogen(get_populations(in)[j], par);
			/* free_pathogen(out->pathogens[count]); */
			out->pathogens[count] = reconstruct_genome(ppat);
			out->popid[count++] = j;
		}
	}

	/* free local pointers */
	free(nAvailPerPop);
	free(nIsolatesPerPop);

	return out;
} /* end draw_sample */






/* merge several samples together */
struct sample * merge_samples(struct sample **in, int nsamp, struct param *par){
	int i, j, newsize=0, counter=0;


	/* create output */
	for(i=0;i<nsamp;i++) newsize += get_n(in[i]);
	struct sample * out = create_sample(newsize);

	/* fill in output */
	for(i=0;i<nsamp;i++){
		for(j=0;j<get_n(in[i]);j++){
			out->pathogens[counter] = copy_pathogen(in[i]->pathogens[j]);
			out->popid[counter++] = in[i]->popid[j];
		}
	}

	out->n = newsize;
	return out;
}





/* translate sampling dates into simulation timestep */
void translate_dates(struct param *par){
	int i;
	for(i=0;i<par->n_sample;i++) par->t_sample[i] = par->duration - par->t_sample[i]; /* minimum date must be 1 */
}




/* SPLIT DATA OF A SAMPLE BY POPULATION */
struct sample ** seppop(struct sample *in, struct param *par){
	int i, j, counter, n=get_n(in), npop;
	struct table_int * tabpop;
	struct sample ** out;

	/* get table of population sizes */
	tabpop = get_table_int(in->popid, n);
	npop = tabpop->n;

	/* allocate memory */
	out = (struct sample **) calloc(npop, sizeof(struct sample *));
	if(out==NULL){
		fprintf(stderr, "\n[in: sampling.c->seppop]\nNo memory left to separate isolates per population. Exiting.\n");
		exit(1);
	}

	for(i=0;i<npop;i++){
		out[i] = create_sample(tabpop->times[i]);
		if(out[i]==NULL){
			fprintf(stderr, "\n[in: sampling.c->seppop]\nNo memory left to separate isolates per population. Exiting.\n");
			exit(1);
		}
	}

	/* copy pathogens */
	for(i=0;i<npop;i++){
		counter=0;
		for(j=0;j<n;j++){
			if(in->popid[j]==tabpop->items[i]) {
				out[i]->pathogens[counter] = copy_pathogen(in->pathogens[j]);
				out[i]->popid[counter++] = in->popid[j];
			}
		}
	}

	/* free memory and return */
	free_table_int(tabpop);
	return out;
}











/* gcc line:

   gcc -o sampling param.c auxiliary.c pathogens.c populations.c dispersal.c sampling.c -Wall -O0 -lgsl -lgslcblas

   valgrind --leak-check=yes sampling

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

/* 	/\* free memory *\/ */
/* 	free_metapopulation(metapop); */
/* 	free_sample(samp); */
/* 	free_network(cn); */
/* 	free(par); */
/* 	gsl_rng_free(rng); */

/* 	return 0; */
/* } */
