/// \file 
/// \brief simple example of using a function from the gsl library

#include <cstdio>
#include <gsl_rng.h>
#include <gsl_randist.h>

int main (int argc, char *argv[])
{
  /* set up GSL RNG */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  /* end of GSL setup */

  int i,n=4;   // generate 4 random pairs by default 
  double gauss,gamma;  

  if (argc>1) n=atoi(argv[1]);  // generate # of pairs requested on cmd line
  for (i=0;i<n;i++)
    {
      gauss=gsl_ran_gaussian(r,2.0);
      gamma=gsl_ran_gamma(r,2.0,3.0);
      printf("%2.4f %2.4f\n", gauss,gamma);
    }
  return(0);
}

