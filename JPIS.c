//
//  JPIS.c
//  
//
//  Physical Chemistry II Project
//
//

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>

int i,j;
double MAX_I;
double x1, x2,dx;
double psi1, psi2, psi3, psi4, ee;
double pi;
double sum, total;

//Basic Parameters
int nmax, nelec, ntotal; // nmax = highest eigenfunction value for n, nelec = total # of electrons in system, ntotal = total # of orbitals.
int nocc, nuno; // nocc = number of occupied orbitals, which is total # of electrons / 2, with 2 electrons per orbital.
                // nuno is remaining unoccupied orbitals.
int ncis, nstates; // ncis is total # of single excited configurations, nstates is total number of ncis plus ground state configuration.


// Required energy functions...


int main()

{
    
    // DEFINE BASIC PARAMETERS (try for r only first...)
    pi = 4.*atan(1.0);
    
    nmax = 3; // highest eigenfunction value for n.
    nelec = 10; // varies for given system. 10 for water.
    
    // NEED EXPLANATION ON HOW THIS LIMITS TO S,P,D ORBITALS
    //
    // total number of orbitals... note we are limiting l<=2 in this case (s, p, and d orbs)
    ntotal=0;
    for (i=1; i<=nmax; i++) {
        
        if (i<=3) {
            ntotal+= i*i;
        }
        else {
            ntotal += 9;
        }
    }
    // -------------------------------------------------------
    
    
    // CONTINUE DEFINING BASIC PAREMETERS
    nocc = nelec / 2.;
    nuno = ntotal - nocc;
    
    printf(" total # of orbitals is %i\n",ntotal);
    printf(" total # of occupied orbitals is %i\n",nocc);
    printf(" virtual unocc orbitals is %i\n",nuno);

    ncis = nocc * nuno; // NUMBER OF SINGLE EXCITED CONFIGURATIONS
    nstates = ncis + 1; // NUMBER OF SINGLE EXCITED CONFIGURATIONS + GROUND STATE = TOTAL STATES
    
    printf(" # of single excited states is %i\n",ncis);
    printf(" # of total states is %i\n",nstates);
    
    
    MAX_I = 10.; // nmax => The length of the box.
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
            //printf("%i %i %f %f %f %f %f %f \n",i,j,psi1,psi2,psi3,psi4,ee,sum);
            
        }
        
        total +=sum;
        
        
    }
    
    //printf(" %f\n",total);
    
}

