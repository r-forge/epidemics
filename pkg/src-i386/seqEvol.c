/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating sequence evolution.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "seqEvol.h"





/*
   =============================
   === STRUCTURES DEFINITION ===
   =============================
*/




/*
   =================================
   === LOCAL AUXILIARY FUNCTIONS ===
   =================================
*/




/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/








/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/



/* TESTING in R */

/*
## test raw conversion
.C("testRaw", raw(256), 256L, PACKAGE="adegenet")
.C("testSizePointer", integer(1), integer(1), integer(1), PACKAGE="adegenet")

## test raw->int conversion
x <- sample(0:1,800,replace=TRUE)
toto <- .bin2raw(x)$snp
all(.C("bytesToBinInt", toto, length(toto), integer(length(toto)*8))[[3]]==x)

## test raw vec -> binary integers
.C("bytesToBinInt",as.raw(c(12,11)), 2L, integer(16), PACKAGE="adegenet")

## test several raw vec -> int (allele counts, any ploidy)
.C("bytesToInt",as.raw(c(12,11)), 1L, 2L, integer(8), integer(16), PACKAGE="adegenet")


*/

