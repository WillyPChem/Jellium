//
//  JPIS.c
//  
//
//  Physical Chemistry II Project
//
//

#include<stdio.h>
#include<iostream>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<malloc.h>
#include<complex.h>
#include<time.h>

int i,j;
double MAX_I;
double x1, x2,dx;
double psi1, psi2, psi3, psi4, ee;
double pi;
double sum, total;

double L;

double mass = 1.;
double hbar = 1.;

void Phi(double R, int nmin, int nmax, int lmin, int lmax, double *T);
void Spherical_Y(int l, int m, double theta, double phi, double *YR, double *YI);
double Legendre(int l, int m, double theta);
double Bessel(double R, double r, int n, int l);
void OrderPsis(int norbs, int *E, int **MO);
double prefac(int m, int l);
double factorial(int n);
double plgndr(int l, int m, double x);

//Basic Parameters
int nelec, ntotal; // nmax = highest eigenfunction value for n, nelec = total # of electrons in system, ntotal = total # of orbitals.
int nocc, nuno; // nocc = number of occupied orbitals, which is total # of electrons / 2, with 2 electrons per orbital.
                // nuno is remaining unoccupied orbitals.
int ncis, nstates; // ncis is total # of single excited configurations, nstates is total number of ncis plus ground state configuration.

// Relevant Hartree Fock Matrices
double enuc, *S, *AO, *V, *Hcore, *E;

int nmin, nmax, lmin, lmax;
int n,l,m;
int *nval, *lval, *mval;


// Required energy functions.

int main()

{

    
    // DEFINE BASIC PARAMETERS (try for r only first...)
    pi = 4.*atan(1.0);
    
    nmin = 0; // lowest eigenfunction value for n.
    nmax = 3; // highest eigenfunction value for n.
    
    lmin = 0; // lowest l
    lmin = 3; // highest l 

    nelec = 4; // varies for given system. 10 for water.
    

	nval = (int *)malloc(nmax*nmax*nmax*sizeof(int));
	lval = (int *)malloc(nmax*nmax*nmax*sizeof(int));
	mval = (int *)malloc(nmax*nmax*nmax*sizeof(int));

     AO = (double *)malloc(nmax*nmax*nmax*sizeof(double)); // Atomic Orbital Energy Matrix
   // E = (double *)malloc(nmax*nmax*sizeof(double));


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

    MAX_I = 100.; // nmax => The length of the box. Er, radius of particle. Err...
    L = 1;

    Phi(L, nmin, nmax, lmin, lmax, AO);
    
    // HARTREE FOCK RELEVANT MATRICIES
    S = (double *)malloc(nmax*nmax*sizeof(double)); // A-O Overlap Matrix
    Hcore = (double *)malloc(nmax*nmax*sizeof(double)); // Hamiltonian Core for HF

    // Calculate orbital energies for HF.
    
}


/* 
/ Function to Determine Energy Calculation
/ Take function and loop through n to keep all atomic orbital energy levels.
*/

void Phi(double R, int nmin, int nmax, int lmin, int lmax, double *T)
{
    // n,l,m
    int idx = 0;
    for ( n=nmin; n<=nmax; n++)
    {
        for ( l=lmin; l<=lmax; l++)
        {
            for ( m = -l; m<=l; m++)
            {
                nval[idx] = n;
                lval[idx] = l;
                mval[idx] = m;

                T[idx] = ((pow(hbar,2)*pow(pi,2))/(8*mass*pow(R,2)))*pow((2.*n+l+1),2); 
                printf(" for phi=%i, n=%i l=%i m=%i energy is %f\n",idx, nval[idx], lval[idx], mval[idx], T[idx]);
                idx++;
        
            }
        }
    }

    for(int jdx=0.; jdx<=nmax; jdx++)
    {

   // printf(" for phi=%i, n=%i l=%i m=%i energy is %f\n",jdx, nval[jdx], lval[jdx], mval[jdx], T[jdx]); 
    
    }
}


//  Evaluates Spherical Harmonics at values of theta and phi given
//  quantum numbers l and m... real part of the function goes to YR, imag goes to YI
void Spherical_Y(int l, int m, double theta, double phi, double *YR, double *YI) {
    
    int mp;
    double ctheta, pfac, P;
    double complex y;
    
    mp = m;
    // Legendre Polynomial function will only take positive values of m
    if (m<0) {
        
        mp = abs(m);
    }
    
    // Prefactor for Y
    pfac = prefac( mp, l );
    
    // cosine of theta
    ctheta = cos(theta);
    
    // Legendre Polynomial P_l^m (cos(theta))
    P = plgndr( l, mp, ctheta);
    
    // Spherical Harmonic = prefac*P_l^m(cos(theta))*exp(i*m*phi)
    y = pfac*P*cexp(I*m*phi);
    
    *YR = creal(y);
    *YI = cimag(y);
    
}

//  Uses asymptotic approximation to Spherical Bessel Function of
//  quantum numnber n and l (see The Kraus and Schatz, JCP 79, 6130 (1983); doi: 10.1063/1.445794)
//  Big R is radius of the particle, little r is current value of r variable
//  Function returns value of the function at r
double Bessel(double R, double r, int n, int l) {
    
    double cterm1, cterm2, cterm3, pre, val;
    
    val = 0;
    
    if (r>R) val=0;
    
    else {
        
        pre = 2./( sqrt(R)*r);
        cterm1 = (2*n+l+2)*r/R;
        cterm2 = l+1;
        cterm3 = (pi/2)*(cterm1-cterm2);
        
        //val = pre*cos(cterm3)/sqrt(2);
        val = pre*cos(cterm3);
    }
    
    return val;
}

//  This is from Numerical Recipes in C!
double plgndr(int l, int m, double x) {
    //Computes the associated Legendre polynomial P m
    //l (x). Here m and l are integers satisfying
    //0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.
    //void nrerror(char error_text[]);
    
    double fact,pll,pmm,pmmp1,somx2;
    
    int i,ll;
    
    if (m < 0 || m > l || fabs(x) > 1.0) {
        //if (m>l || fabs(x) > 1.0) {
        printf("Bad arguments in routine plgndr\n");
        exit(0);
    }
    
    pmm=1.0;   //Compute P^m_m .
    
    if (m > 0) {
        
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
        
        for (i=1;i<=m;i++) {
            
            pmm *= -fact*somx2;
            fact += 2.0;
            
        }
    }
    
    if (l == m)
        return pmm;
    
    else {    //Compute P^m_m+1
        
        pmmp1=x*(2*m+1)*pmm;
        
        if (l == (m+1))
            
            return pmmp1;
        
        else {   //Compute P^m_l, l>m+1
            
            for (ll=m+2;ll<=l;ll++) {
                
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm=pmmp1;
                pmmp1=pll;
                
            }
            
            return pll;
        }
    }
}
//  Computes factorials!
double factorial(int n) {
    int c;
    double result = 1;
    
    for (c = 1; c <= n; c++)
        result = result * c;
    
    return result;
}
//  Computes normalization constant for Spherical Harmonics for a given m and l
double prefac(int m, int l) {
    
    
    double p, num1, num2, denom1, denom2;
    
    num1 = 2*l+1;
    num2 = factorial( (l-m) );
    
    denom1 = 4*pi;
    denom2 = factorial( (l+m) );
    
    
    p = sqrt((num1/denom1)*(num2/denom2));
    
    return p;
    
}
//  Returns value of normalized Legendre polynomial at value of theta
//  given quantum number l and m
double Legendre(int l, int m, double theta) {
    
    
    int mp;
    double ctheta, pfac, P;
    double y;
    
    mp = m;
    // Legendre Polynomial function will only take positive values of m
    if (m<0) {
        
        mp = abs(m);
    }
    // Prefactor 
    pfac = prefac( mp, l );
    
    // cosine of theta
    ctheta = cos(theta);
    
    // Legendre Polynomial P_l^m (cos(theta))
    P = plgndr( l, mp, ctheta);
    
    // Spherical Harmonic = prefac*P_l^m(cos(theta))*exp(i*m*phi)
    y = pfac*P;
    
    return y;
}

