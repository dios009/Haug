#include <stdio.h>
#include <math.h>
#include "huag.h"

double k = 100.;
double E2 = 1.01;
double x2 = 1.;

double p2x, k1x, k2x, rox, taux ;


int main()
{

p2x = sqrt(E2*E2-1.);
k1x = -k;
k2x = k*(-p2x*x2+E2);
rox = sqrt(2.*(k+1-E2-k*(-p2x*x2+E2)));
taux = -E2;

    printf("\n  p2= %lf , \t k2= %lf , \t k1= %lf , \t ro=%lf \t tau=%f \n", p2x, k2x, k1x, rox, taux);
/*

*/
    printf("\n  M0= %lf \n", M0(k1x,k2x,rox,taux));
    printf("\n  MX1= %lf \n", MX1(k1x,k2x,rox,taux));
    printf("\n  MX2= %lf \n", MX2(k1x,k2x,rox,taux));
    printf("\n  MX3= %lf \n", MX3(k1x,k2x,rox,taux));
    printf("\n  MX4= %lf \n", MX4(k1x,k2x,rox,taux));
    printf("\n  MX5= %lf \n", MX5(k1x,k2x,rox,taux));
    printf("\n  ML1= %lf \n", ML1(k1x,k2x,rox,taux));
    printf("\n  ML2= %lf \n", ML2(k1x,k2x,rox,taux));
    printf("\n  ML3= %lf \n", ML3(k1x,k2x,rox,taux));
    printf("\n  ML4= %lf \n \n", ML4(k1x,k2x,rox,taux));

    printf("\n  M0ex= %lf \n", M0(k2x,k1x,rox,taux));
    printf("\n  MX1ex= %lf \n", MX1(k2x,k1x,rox,taux));
    printf("\n  MX2ex= %lf \n", MX2(k2x,k1x,rox,taux));
    printf("\n  MX3ex= %lf \n", MX3(k2x,k1x,rox,taux));
    printf("\n  MX4ex= %lf \n", MX4(k2x,k1x,rox,taux));
    printf("\n  MX5ex= %lf \n", MX5(k2x,k1x,rox,taux));
    printf("\n  ML1ex= %lf \n", ML1(k2x,k1x,rox,taux));
    printf("\n  ML2ex= %lf \n", ML2(k2x,k1x,rox,taux));
    printf("\n  ML3ex= %lf \n", ML3(k2x,k1x,rox,taux));
    printf("\n  ML4ex= %lf \n \n", ML4(k2x,k1x,rox,taux));

printf("\n  Coef= %lf \n", (-0.57946727)*p2x/(k*rox));

   return 0;
}
