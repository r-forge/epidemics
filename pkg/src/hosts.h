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
};


struct pop{
	struct host *s, *i, *r; /* arrays */
	unsigned int ns, ni, nr; /* size of each array */
}

