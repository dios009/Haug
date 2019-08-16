//Differential cross section of electron-positron Bremsstrahlung - Eberhard Huag

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "huag.h"
#include "cuba.h"

#define k 10.

double E2, p2;


static int Integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata) {


  double x2= (1.-x2min(E2, k, p2))*xx[0]+x2min(E2, k, p2);

 
    ff[0] = p2/(k*ro(E2,k,p2,x2))*(M0(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +MX1(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +MX2(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +MX3(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +MX4(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +MX5(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +ML1(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +ML2(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +ML3(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +ML4(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
    +M0(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +MX1(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +MX2(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +MX3(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +MX4(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +MX5(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +ML1(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +ML2(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +ML3(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2))
    +ML4(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2)))*(1.-x2min(E2,k,p2))*(-0.57946727);
    
    

  return 0;
}

/*********************************************************************/

#define NDIM 1
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 5000 /*era 0 */
#define MAXEVAL 5000000 /*era 50000*/

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main() {


  int nregions, neval, fail;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
	double suma =0;

/*Do loop */



	E2=emin(k);

	FILE *fptr;
	fptr = fopen("datos-Vegas-Huag2.dat","w");

	do
	{
	  p2 = sqrt(E2*E2-1.);

	printf("\n E2=%f \n",E2);
	
			
	  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	    EPSREL, EPSABS, VERBOSE, SEED,
	    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
	    GRIDNO, STATEFILE, SPIN,
	    &neval, &fail, integral, error, prob);

	  printf("VEGAS RESULT:\tneval %d\tfail %d\n",neval, fail);

    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral[0], error[0], prob[0]);
	
	fprintf(fptr,"%f \t %f \n",E2,integral[0]);
	
		E2 += (emax(k)-emin(k))/100;

	suma+= integral[0]*(emax(k)-emin(k))/100;
	}
	while(E2 < emax(k));
	
	printf("\n Suma =%e \n", suma);

	fclose(fptr);

	return 0;
}
