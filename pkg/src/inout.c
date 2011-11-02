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
