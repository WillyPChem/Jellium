#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>

// Coding Boot Camp

int main() 
{
	int i, MAX_I;
	double x, dx, phi, phi2;
	double pi = 4.*atan(1.0);

// 2B. Evaluate sqrt 2/L * sin (pi*x / L) from O to L when L = 10

	MAX_I = 6.;
	dx = 10./MAX_I;

	for(i=0; i<=MAX_I; i++) 
	{
		x = dx*i;
		phi = sqrt(2./10)*sin(pi*x/10.);
		double phi2 = pow(phi,2);
		printf(" %i %f %f \n",i,x,phi2);
	}

}




