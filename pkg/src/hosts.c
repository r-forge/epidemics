/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/

#include "common.h"
#include "hosts.h"
#include "seqEvol.h"



/*
   =================
   === ACCESSORS ===
   =================
*/

unsigned int get_host_id(struct host *in){
	return in->id;
}




short int get_host_ninf(struct host *in){
	return in->ninf;
}




struct pathogen ** get_host_inf(struct host *in){
	return in->infections;
}





/*
struct host * get_sus(struct population *in){
	return in->sus;
}

struct host * get_inf(struct population *in){
	return in->inf;
}

struct host * get_rec(struct population *in){
	return in->rec;
}
*/








/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

/* Create empty host */
struct host * create_host(){
	struct host *out;
	out = (struct host *) calloc(1, sizeof(struct host));
	if(out == NULL){
		fprintf(stderr, "\nNo memory left for creating new host. Exiting.\n");
		exit(1);
	}
	out->id = 1;
	out->infections = NULL; /* new host created without infections */
	return out;
}













/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

/* Free host */
void free_host(struct host *in){
	int i;
	for(i=0;i<get_host_ninf(in);i++){
		free_pathogen((in->infections)[i]);
	}
	free(in->infections);
	free(in);
}






/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/

/* Print host content */
void print_host(struct host *in){
	printf("\nhost %d", get_host_id(in));
	printf("\n%d infections",get_host_ninf(in));
}





bool is_infected(struct host * in){
	return get_host_ninf(in)==0;
}








/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/
/* void infect_new_host(struct host *host1, struct host *host2){ */
	/* create new infection vector in host */
	/* make pathogen replication */
/*} */



