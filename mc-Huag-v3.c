// MonteCarlo
//Differential cross section of electron-positron Bremsstrahlung - Eberhard Huag
// version con E2 y x2 unitarios para calcular la seccion total!!!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "huag.h"

#define K 100.
#define EMIN emin(K)
#define EMAX emax(K)

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float rand2(long *idum)
/*Long period (> 2  1018) random number generator of L'Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.*/
{
	int j;
	long kk; //Cambie k por kk , para evitar problema con el k de la energia definido arriba
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0) 
	{ //Initialize.
		if (-(*idum) < 1) 
			*idum=1; //Be sure to prevent idum = 0.
		else *idum = -(*idum);
			idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) { //Load the shuffle table (after 8 warm-ups).
			kk=(*idum)/IQ1;
			*idum=IA1*(*idum-kk*IQ1)-kk*IR1;
			if (*idum < 0) 
				*idum += IM1;
			if (j < NTAB) 
				iv[j] = *idum;
		}
		iy=iv[0];
	}
	kk=(*idum)/IQ1; //Start here when not initializing.
	*idum=IA1*(*idum-kk*IQ1)-kk*IR1; //Compute idum=(IA1*idum) % IM1 without
	if (*idum < 0) 
		*idum += IM1; //overflows by Schrage's method.
	kk=idum2/IQ2;
	idum2=IA2*(idum2-kk*IQ2)-kk*IR2; //Compute idum2=(IA2*idum) % IM2 likewise.
	if (idum2 < 0) 
		idum2 += IM2;
	j=iy/NDIV; //Will be in the range 0..NTAB-1.
	iy=iv[j]-idum2; //Here idum is shuffled, idum and idum2 are
	iv[j] = *idum; //combined to generate output.
	if (iy < 1) 
		iy += IMM1;
	if ((temp=AM*iy) > RNMX) 
		return RNMX; //Because users don't expect endpoint values.
	else return temp;
}

// Defino la matriz densidad con x2 y E2 que varíen entre 0 y 1 

double MM (double E2, double k, double p2, double x2){
return 
p2/(k*ro(E2,k,p2,x2))*(M0(k1(k),k2(E2,k,p2,x2),ro(E2,k,p2,x2),tau(E2))
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
+ML4(k2(E2,k,p2,x2),k1(k),ro(E2,k,p2,x2),tau(E2)))*(-0.57946727);
}
/*********************************************************************/


int main() {


int i;
long init = -1;

double NEVAL = 1000000000.;

double count = 0.;
for (i = 1; i < NEVAL+1 ; ++i)
  {
    double E2 = rand2(&init);
	double x2 = rand2(&init);
    double p2 = sqrt((EMIN+(EMAX-EMIN)*E2)*(EMIN+(EMAX-EMIN)*E2)-1.); //definido para E2 entre 0 y 1
	count += MM(EMIN+(EMAX-EMIN)*E2,K,p2,(1.-x2min(EMIN+(EMAX-EMIN)*E2,K,p2))*x2+x2min(EMIN+(EMAX-EMIN)*E2,K,p2))*(1.-x2min(EMIN+(EMAX-EMIN)*E2,K,p2))*(EMAX-EMIN);//definida en el hypercubo unitario
  }
  printf("TOTAL = %f \n", count/NEVAL);

	return 0;
}