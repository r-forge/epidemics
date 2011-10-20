/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/


/*
# status legend # 
- s: susceptible
- i: infected
- r: removed/recovered
*/


struct population{
	struct pathogen ** pathogens;
	/* nb of item; ninf is the length of **pathogens */
	unsigned int nsus, ninf, nrec;
};

