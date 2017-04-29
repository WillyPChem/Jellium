#include<stdio.h>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include"blas.h"
#include<malloc.h>
#include<complex.h>
#include<time.h>
#include<string.h>

int pi;
int i,j;
int dim;
int nmax;
double L, mass, hbar;

// Relevant HF functions
void BuildDensity(int dim, int occ, double *C, double *D);
int DIAG_N(int dim, int number, double *mat, double *en, double *wfn);
void Diagonalize(double*M,long int dim, double*eigval,double*eigvec);
void print_matrix( char* desc, int m, int n, double* a, int lna);
void LoopMM(int dim, double *a, char *transa, double *b, char *transb, double *c);
double E_Total(int dim, double *D, double *HCore, double *F, double Enuc);
void ReadEI(int dim, FILE *fp, double *EE);
int FourDIndx(int i, int j, int k, int l, int dim);
double DensityDiff(int dim, double *D, double *Dnew);
void UpdateF(int dim, double *D, double *Hcore, double *EI, double *Fnew);

// Custom cubic HF functions
void CubicPhi();
void AtomicOrbitalOverlap();
void CrawdadFormat();

//Basic Parameters
int nelec, ntotal; // nmax = highest eigenfunction value for n, nelec = total # of electrons in system, ntotal = total # of orbitals.
int nocc, nuno; // nocc = number of occupied orbitals, which is total # of electrons / 2, with 2 electrons per orbital. // nuno is remaining unoccupied orbitals.
int ncis, nstates; // ncis is total # of single excited configurations, nstates is total number of ncis plus ground state configuration.

// Relevant Hartree Fock Variables
double *S, *sqrtS;
double *T;
double *Hcore;
double *E, *E1, *E2;
double *A;
double *lambdasqrt, *temp, *Fock;
double Enuc;

// Phi Variables
int *NPOrbE, *NPOrb_x, *NPOrb_y, *NPOrb_z;

int n;

int main()
{


	// Definition of pi
    pi = 4.*atan(1.0);
    
    // Atomic Units
    // ------------
    L = 1;
    mass = 1;
    hbar = 1;

    // Highest eigenfunction value for N.
    nmax = 2;

    // Define dimensions (nmax*nmax*nmax)
    //  dim = nmax*nmax*nmax;
    dim = 7;

    // Number of electrons in system.
    nelec = 2; 
   
    // total number of orbitals... note we are limiting l<=2 in this case (s, p, and d orbs)

    ntotal=0;
  /*  for (i=1; i<=nmax; i++) {
        
        if (i<=3) {
            ntotal+= i*i;
        }
        else {
            ntotal += 9;
        }
    } */
    // ------------------------------------------------------- 
    
    // Continue definition of basic parameters
    nocc = nelec / 2.;
    nuno = ntotal - nocc;
    
    printf(" total # of orbitals is %i\n",ntotal);
    printf(" total # of occupied orbitals is %i\n",nocc);
    printf(" virtual unocc orbitals is %i\n",nuno);

    ncis = nocc * nuno; // NUMBER OF SINGLE EXCITED CONFIGURATIONS
    nstates = ncis + 1; // NUMBER OF SINGLE EXCITED CONFIGURATIONS + GROUND STATE = TOTAL STATES
    
    printf(" # of single excited states is %i\n",ncis);
    printf(" # of total states is %i\n",nstates);

    // Hartree Fock Code
	FILE *nucfp, *nAttract, *kinEnergy, *eeRep;




}




