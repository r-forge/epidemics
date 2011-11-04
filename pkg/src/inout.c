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
#include "sampling.h"
#include "sumstat.h"
#include "inout.h"



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/







/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

/* write group compositions */
void write_ts_groupsizes(struct ts_groupsizes *in){
	int i;
	FILE *outfile = fopen( "out-popsize.txt", "w");
	if(outfile==NULL){
		fprintf(stderr, "\n[in: inout.c->write_popsize]\nUnable to open file 'out-popsize.txt'.\n");
		exit(1);
	}

	fprintf(outfile, "step\tnsus\tninf\tnrec\tninfcum\n");
	for(i=0;i<in->length;i++)
		fprintf(outfile, "%d\t%d\t%d\t%d\t%d\n", i+1, in->nsus[i], in->ninf[i], in->nrec[i], in->ninfcum[i]);
	fclose(outfile);
}




/* write summary statistics */
void write_ts_sumstat(struct ts_sumstat *in){
	int i;
	FILE *outfile = fopen( "out-sumstat.txt", "w");
	if(outfile==NULL){
		fprintf(stderr, "\n[in: inout.c->write_ts_sumstat]\nUnable to open file 'out-sumstat.txt'.\n");
		exit(1);
	}

	fprintf(outfile, "step\tnbSnps\tHs\tmeanNbSnps\tvarNbSnps\tmeanPairwiseDist\tvarPairwiseDist\tFst\n");
	for(i=0;i<in->length;i++){
		fprintf(outfile, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", in->steps[i], in->nbSnps[i], in->Hs[i], in->meanNbSnps[i], in->varNbSnps[i], in->meanPairwiseDist[i], in->varPairwiseDist[i], in->Fst[i]);
	}
	fclose(outfile);
}





/* write a sample */
void write_sample(struct sample *in){
	int i, j, nbSnps;
	FILE *outfile = fopen( "out-sample.txt", "w");
	if(outfile==NULL){
		fprintf(stderr, "\n[in: inout.c->write_popsize]\nUnable to open file 'out-popsize.txt'.\n");
		exit(1);
	}

	for(i=0;i<in->n;i++){
		nbSnps = get_nb_snps(in->pathogens[i]);
		fprintf(outfile, "pop %d\n", get_popid(in->pathogens[i]));
		for(j=0;j<nbSnps;j++){
			fprintf(outfile, "%d ", get_snps(in->pathogens[i])[j]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}

