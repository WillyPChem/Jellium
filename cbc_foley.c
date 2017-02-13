#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>

int i,j;
double MAX_I;
double x1, x2;
double x, xpdx, dx, fx, fxpdx, fp, m, fxp;
double psi1, psi2, psi3, psi4, ee, rep;
double pi;
double sum, total;

int main()

{

	pi = 4.*atan(1.0);
	MAX_I = 10.;
	dx = 10./MAX_I;
	sum = 0.;
	i = 0;
	j = 0;

	for(i=0; i<=MAX_I; i++)
	{
		x1 = dx*i;

		for(j=0; j<=MAX_I; j++)
		{

		x2 = dx*j;

		psi1 = sqrt(2./10)*sin((pi*x1)/10.);
		psi2 = sqrt(2./10)*sin((pi*2.*x1)/10.);
		psi3 = sqrt(2./10)*sin((pi*x2)/10.);
		psi4 = sqrt(2./10)*sin((pi*2.*x2)/10.);

		ee = 1./fabs(x1-x2);
		sum += psi1*psi2*ee*psi3*psi4;
		printf("%i %i %f %f %f %f %f %f \n",i,j,psi1,psi2,psi3,psi4,ee,sum);
		
		}

	total +=sum;
	

	}

	printf(" %f\n",total);
	
}




