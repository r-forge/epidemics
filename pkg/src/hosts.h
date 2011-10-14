/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/



/*
# status legend # 
- s: susceptible
- i: infected
- r: removed
*/
struct host{
	unsigned int id;
	unsigned short int ninf;
	struct pathogen **infections;
};


struct population{
	struct host **sus, **inf, **rec; /* arrays */
	unsigned int nsus, ninf, nrec; /* size of each array */
}

