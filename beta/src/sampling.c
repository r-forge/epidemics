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
	int i, n=get_n(in), *popid, *pool, npop;

	/* get vector of pop id */
	popid = (int *) calloc(n, sizeof(int));
	if(popid == NULL){
		fprintf(stderr, "\n[in: sampling.c->get_npop]\nNo memory left to sample per population. Exiting.\n");
		exit(1);
	}
	for(i=0;i<n;i++) popid[i] = get_popid(in->pathogens[i]);


	/* enumerate nb of unique items */
	/* create pool of unique items */
	pool = (int *) calloc(n, sizeof(int));
	if(pool == NULL){
		fprintf(stderr, "\n[in: sampling.c->get_npop]\nNo memory left to sample per population. Exiting.\n");
		exit(1);
	}

	/* list and count all SNPs */
	npop = 0;
	for(i=0;i<n;i++){
		if(int_in_vec(popid[i], pool, npop) < 0){
			pool[npop++] = popid[i];
		}
	}

	/* free memory and return */
	free(popid);
	free(pool);
	return(npop);
}






/*
   ===================
   === CONTRUCTORS ===
   ===================
*/

struct sample * create_sample(int n){
	int i;
	struct sample *out = calloc(1, sizeof(struct sample));
	if(out == NULL){
		fprintf(stderr, "\n[in: population.c->create_sample]\nNo memory left to sample the metapopulation. Exiting.\n");
		exit(1);
	}

	/* allocate memory for pathogens */
	out->pathogens = (struct pathogen **) calloc(n, sizeof(struct pathogen *));
	if(out->pathogens == NULL){
		fprintf(stderr, "\n[in: population.c->create_sample]\nNo memory left sample pathogens from the metapopulation. Exiting.\n");
		exit(1);
	}

	for(i=0;i<n;i++) out->pathogens[i] = create_pathogen();
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
	printf("\n%d pathogens", in->n);
	if(showGen){
		for(i=0;i<in->n;i++){
			//if(!isNULL_pathogen((in->pathogens)[i])) print_pathogen((in->pathogens)[i]);
			print_pathogen((in->pathogens)[i]);
		}
		printf("\n");
	}
}




/* reconstruct genome of an isolate */
struct vec_int * reconstruct_genome(struct pathogen *in, struct metapopulation * metapop){
	int i, lineagesize=1;
	struct pathogen *curAnces = get_ances(in);
	struct vec_int ** lineage, *temp, *genome;

	/* identify lineage */
	while(curAnces != NULL){
		lineagesize++;
		curAnces = get_ances(curAnces);
	}

	/* get all snps in the lineage */
	lineage = (struct vec_int **) calloc(lineagesize, sizeof(struct vec_int *));

	lineage[0] = get_snps_vec(in);
	for(i=1;i<lineagesize;i++){
		lineage[i] = get_snps_vec(curAnces);
		curAnces = get_ances(curAnces);
	}

	/* merge snps */
	temp = merge_vec_int(lineage, lineagesize);

	/* remove reverse mutations */
	genome = keep_odd_int(temp);

	/* free temporary allocation & return */
	free(lineage);
	free_vec_int(temp);
	return genome;
}




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* Get sample of isolates */
/* Isolates are COPIED, so that any modification of the sample does not alter */
/* the metapopulation. */
struct sample * draw_sample(struct metapopulation *in, int n, struct param *par){
	int i, j, id, nIsolates=0, maxnpat=get_maxnpat(in);
	int *availIsolates;
	struct vec_int *temp;

	/* create pointer to pathogens */
	struct sample *out=create_sample(n);

	/* get the number of isolates that can be sampled */
	for(i=0;i<maxnpat;i++){
		if(!isNULL_pathogen((get_pathogens(in)[i]))) nIsolates++;
	}

	if(nIsolates != get_total_ninf(in)){
		fprintf(stderr, "\n[in: population.c->draw_sample]\nNumber of available isolates (%d) does not match total number of infected (%d). Exiting.\n", nIsolates, get_total_ninf(in));
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
		while(isNULL_pathogen(get_pathogens(in)[j])) j++;
		availIsolates[i] = j++;
	}

	/* choose from available pathogens */
	for(i=0;i<n;i++){
		id=gsl_rng_uniform_int(par->rng,nIsolates);
		/* (out->pathogens)[i] = (in->pathogens)[availIsolates[id]]; */ /* this just copies addresses; need to copy content! */
		copy_pathogen(in->pathogens[availIsolates[id]], out->pathogens[i], par);

		/* reconstruct genomes */
		printf("\npathogen %d before", i);
		print_pathogen(out->pathogens[i]);
		temp = reconstruct_genome(out->pathogens[i], in);
		free_vec_int(out->pathogens[i]->snps); /* free old snps */
		out->pathogens[i]->snps = temp; /* replace with reconstructed genome */
		printf("\npathogen %d after", i);
		print_pathogen(out->pathogens[i]);
	}

	/* free local pointers */
	free(availIsolates);

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
			copy_pathogen(in[i]->pathogens[j], out->pathogens[counter++], par);
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




/* slit data of a sample by population */
struct sample ** seppop(struct sample *in, struct param *par){
	int i, j, counter, *popid, n=get_n(in), npop;
	struct table_int * tabpop;
	struct sample ** out;

	/* get table of population sizes */
	popid = (int *) calloc(n, sizeof(int));
	if(popid==NULL){
		fprintf(stderr, "\n[in: sampling.c->seppop]\nNo memory left to separate isolates per population. Exiting.\n");
		exit(1);
	}

	for(i=0;i<n;i++) popid[i] = get_popid(in->pathogens[i]);
	tabpop = get_table_int(popid, n);
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
			if(popid[j]==tabpop->items[i]) {
				copy_pathogen(in->pathogens[j], out[i]->pathogens[counter++], par);
			}
		}
	}

	/* free memory and return */
	free(popid);
	free_table_int(tabpop);
	return out;
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
