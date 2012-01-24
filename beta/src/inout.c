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
void write_list_ts_groupsizes(struct ts_groupsizes **in, struct param *par){
	int i, j, listsize = par->npop+1;
	FILE *outfile = fopen( "out-popsize.txt", "w");
	if(outfile==NULL){
		fprintf(stderr, "\n[in: inout.c->write_popsize]\nUnable to open file 'out-popsize.txt'.\n");
		exit(1);
	}

	fprintf(outfile, "patch\tstep\tnsus\tnexp\tninf\tnrec\tnexpcum\n");
	for(j=0;j<listsize;j++){
		for(i=0;i<in[j]->length;i++)
			fprintf(outfile, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", j, i+1, in[j]->nsus[i], in[j]->nexp[i], in[j]->ninf[i], in[j]->nrec[i], in[j]->nexpcum[i]);
	}
	fclose(outfile);
}




/* write summary statistics */
void write_list_ts_sumstat(struct ts_sumstat **in, struct param *par){
	int i, j, listsize = par->npop+1;
	FILE *outfile = fopen( "out-sumstat.txt", "w");
	if(outfile==NULL){
		fprintf(stderr, "\n[in: inout.c->write_ts_sumstat]\nUnable to open file 'out-sumstat.txt'.\n");
		exit(1);
	}

	fprintf(outfile, "patch\tstep\tnbSnps\tHs\tmeanNbSnps\tvarNbSnps\tmeanPairwiseDist\tvarPairwiseDist\tmeanPairwiseDistStd\tvarPairwiseDistStd\tFst\n");
	for(j=0;j<listsize;j++){
		for(i=0;i<in[j]->length;i++){
			fprintf(outfile, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", j, in[j]->steps[i], in[j]->nbSnps[i], in[j]->Hs[i], in[j]->meanNbSnps[i], in[j]->varNbSnps[i], in[j]->meanPairwiseDist[i], in[j]->varPairwiseDist[i], in[j]->meanPairwiseDistStd[i], in[j]->varPairwiseDistStd[i], in[j]->Fst[i]);
		}
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
		fprintf(outfile, "pop %d\n", in->popid[i]);
		for(j=0;j<nbSnps;j++){
			fprintf(outfile, "%d ", get_snps(in->pathogens[i])[j]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}

