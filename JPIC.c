//
//  JPIS.c
//  
//
//  Physical Chemistry II Project
//
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

// Hartree Fock H20 Files



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
void KineticEnergyIntegrals();
void CrawdadFormat();

// Peter Gill Two Electron Repulsion Integral 
double ERI(int dim, double *xa, double *w, double *a, double *b, double *c, double *d);
// Revised Version to fit HF code.
void TwoERICalc();
double pq_int(double px, double py, double pz, double qx, double qy, double qz);
double g_pq(double p, double q, double r);

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

// Two Electron Repulsion Integral Variables

double *adim, *bdim, *cdim, *ddim;
double *ERIa, *ERIb, *ERIc, *ERId, *teri;
double *xb, *w;
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
    nelec = 5; 
   
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

    // HARTREE FOCK CODE
    // ------------------------------------

    FILE *enucfp, *overlap, *nucatt, *ekin, *EEfp;
    double val, enuc, *Sc, *Vc, *Tc, *Hcorec, *lambda, *lambdasquareroot, *Ls, *Fockc, *squarerootS, *temporary;
    double *eps, *Cp, *C, *D, *Dn, sum, Eelec, *Fnew, *EE;
    int ij,kl;
    double *Svals, *Svecs, *SqrtSvals, *SqrtS;
    double ESCF, ESCF_i, deltaE, deltaDD, tolE, tolDD;
    int iter, itermax;

    itermax = 100;

    // Initialize HF relevant matricies
    //
    S = (double *)malloc(dim*dim*sizeof(double)); // A-O Overlap Matrix
    Hcore = (double *)malloc(dim*dim*sizeof(double)); // Hamiltonian Core for HF

    NPOrb_x = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // X corresponding to phi
    NPOrb_y = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // Y corresponding to phi
    NPOrb_z = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // Z corresponding to phi
    NPOrbE = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // energy corresponding

    E1 = (double *)malloc(1+nmax*nmax*nmax*sizeof(double));
    E2 = (double *)malloc(1+nmax*nmax*nmax*sizeof(double));
    A = (double *)malloc(nmax*nmax*nmax*sizeof(double)); // Atomic Orbital Integrals
    T = (double *)malloc(nmax*nmax*nmax*sizeof(double)); // Kinetic Energy Integrals

    lambdasqrt = (double *)malloc(dim*dim*sizeof(double));
    sqrtS = (double *)malloc(dim*dim*sizeof(double)); // Sqrt S matrix
    temp  = (double *)malloc(dim*dim*sizeof(double)); // Temp matrix
    Fock  = (double *)malloc(dim*dim*sizeof(double)); // Fock matrix

    // Initialize 2eri arrays
    adim = (double *)malloc(3*sizeof(double));
    bdim = (double *)malloc(3*sizeof(double));
    cdim = (double *)malloc(3*sizeof(double));
    ddim = (double *)malloc(3*sizeof(double));

    // Storage of 2eri values after ERI function.
    ERIa = (double *)malloc(dim*dim*sizeof(double));
    ERIb = (double *)malloc(dim*dim*sizeof(double));
    ERIc = (double *)malloc(dim*dim*sizeof(double));
    ERId = (double *)malloc(dim*dim*sizeof(double));
    teri = (double *)malloc(dim*dim*sizeof(double));

    n = 320;
    xb = (double *)malloc(n*sizeof(double)); 
    w = (double *)malloc(n*sizeof(double));


    // HF H2O INFO
    Sc = (double *)malloc(dim*dim*sizeof(double));
    Tc = (double *)malloc(dim*dim*sizeof(double));
    Vc = (double *)malloc(dim*dim*sizeof(double));
    Hcorec = (double *)malloc(dim*dim*sizeof(double));
    lambda = (double *)malloc(dim*sizeof(double));
    lambdasquareroot = (double *)malloc(dim*dim*sizeof(double));
    Ls = (double *)malloc(dim*dim*sizeof(double));
    Fockc = (double *)malloc(dim*dim*sizeof(double));
    Fnew = (double *)malloc(dim*dim*sizeof(double));
    squarerootS = (double *)malloc(dim*dim*sizeof(double));
    temporary = (double *)malloc(dim*dim*sizeof(double));
    EE = (double *)malloc(dim*dim*dim*dim*sizeof(double));
    
    eps = (double *)malloc(dim*sizeof(double));
    Cp = (double *)malloc(dim*dim*sizeof(double));
    C = (double *)malloc(dim*dim*sizeof(double));

    D = (double *)malloc(dim*dim*sizeof(double));
    Dn = (double *)malloc(dim*dim*sizeof(double));

    Svals = (double *)malloc(dim*sizeof(double));
    SqrtSvals = (double *)malloc(dim*dim*sizeof(double));
    Svecs = (double *)malloc(dim*dim*sizeof(double));
    SqrtS = (double *)malloc(dim*dim*sizeof(double));	

    //---------------------------------------------------------------------------
    // Step #1: Nuclear Repulsion Energy
    //---------------------------------------------------------------------------

    // Considered a constant. Can skip this for now?
    // HF H20, enuc.dat
    // Can put reading into a function.

    //---------------------------------------------------------------------------
    // Step #2: AO-Overlap, KE Integrals, Building of Hcore.
    //--------------------------------------------------------------------------- 

    // CUSTOM:____________________________________________________________________________
    //------------------------------------------------------------------------------------

    // Calculate AO energies
    // ---------------------
     CubicPhi();
     AtomicOrbitalOverlap();

    // Calculate KE energy integrals (T matrix)
    // ----------------------------------------
     KineticEnergyIntegrals();
     CrawdadFormat();

    // Print Hamiltonian Core
    // ----------------------
    print_matrix(" Hcore ", dim, dim, Hcore, dim); 

    // Define S matrix
    // ---------------
    // print_matrix(" S ", dim, dim, S, dim);

    
    // WATER: __________________________________________________________________________
    //----------------------------------------------------------------------------------    

    enucfp = fopen("./enuc.dat", "r");
    overlap = fopen("./s.dat", "r");
    nucatt = fopen("./v.dat", "r");
    ekin = fopen("./t.dat", "r");
    EEfp = fopen("./eri.dat", "r");
    fscanf(enucfp,"%lf",&Enuc);


    for(i=0; i<dim; i++) {
        for(j=0; j<=i; j++) {

        // Overlap (S)

        fscanf(overlap,"%i",&ij);
        fscanf(overlap,"%i",&kl);
        fscanf(overlap,"%lf",&val);
        Sc[i*dim+j] = val;
        Sc[j*dim+i] = val;

        // Nuclear Attraction (V)

        fscanf(nucatt,"%i",&ij);
        fscanf(nucatt,"%i",&kl);
        fscanf(nucatt,"%lf",&val);
        Vc[i*dim+j] = val;
        Vc[j*dim+i] = val;

        // Kinetic Energy (T)   

        fscanf(ekin, "%lf", &ij);
        fscanf(ekin, "%lf", &kl);
        fscanf(ekin, "%lf", &val);
        Tc[i*dim+j] = val;
        Tc[j*dim+i] = val;

        Hcorec[i*dim+j] = Tc[i*dim+j] + Vc[i*dim+j];
        Hcorec[j*dim+i] = Tc[j*dim+i] + Vc[j*dim+i];

    }
}

    //---------------------------------------------------------------------------
    // Step #3: Two-Electron Repulsion Integrals
    //---------------------------------------------------------------------------

    // WATER:____________________________________________________________________________
    //-----------------------------------------------------------------------------------

    // Read 2-electron integrals
    ReadEI(dim, EEfp, EE);

    // CUSTOM:____________________________________________________________________________
    //------------------------------------------------------------------------------------

    // Peter Gill Gaussian Legendre formatted for our code. Needs to be tested..
    TwoERICalc();

    //---------------------------------------------------------------------------
    // Step #4: Build the Orthogonalization Matrix
    //---------------------------------------------------------------------------

    // WATER:____________________________________________________________________________
    //-----------------------------------------------------------------------------------

    // Diagonalize the overlap matrix.
     	DIAG_N(dim, dim, S, Svals, Svecs);
        
	for (i=0; i<dim; i++) 
	{
    	SqrtSvals[i*dim + i] = pow(Svals[i],-1./2);
	}

    for (i=0; i<dim*dim; i++)
    {
        Dn[i] = 0.;
    }

    // Build the symmetric orthoigonalization matrix: Form S^{-1/2} = L_S s^{-1/2} L_S^t
	 LoopMM(dim, SqrtSvals, "n", Svecs, "t", temp);
	 LoopMM(dim, Svecs, "n", temp, "n", SqrtS);
	// print_matrix(" S^-1/2 ", dim, dim, SqrtS, dim);

    // CUSTOM:____________________________________________________________________________
    //------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------
    // Step #5: Build the initial guess density matrix
    //---------------------------------------------------------------------------


    // WATER:____________________________________________________________________________
    //-----------------------------------------------------------------------------------
      
    // Form an initial (guess) Fock matrix in orthonormal basis using core Hamiltonian as a guess: Fock matrix F = S^{-1/2}^t H_core S^{1/2}
	LoopMM(dim, Hcore, "n", SqrtS, "n", temp);
	LoopMM(dim, SqrtS, "t", temp, "n", Fockc);
	//print_matrix("  Fock", dim, dim, Fockc, dim);

	// Get Guess MO matrix from diagnoalizing Fock matrix:
	DIAG_N(dim, dim, Fock, eps, Cp);
//	print_matrix(" Initial Coeff ", dim, dim, Cp, dim);
	
    // Transform the eigenvectors into the original (non-orthogonal) AO basis.
    LoopMM(dim, SqrtS, "n", Cp, "n", C);

	// Build initial density matrix using the occupied MO's: D = sum (m to occ) C * C.
	BuildDensity(dim, nelec, C, D);
//	print_matrix("  Coefficients", dim, dim, C, dim);
	print_matrix("  Density Matrix", dim, dim, D, dim);

    // CUSTOM:____________________________________________________________________________
    //------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #6: Compute the Initial SCF Energy
    //---------------------------------------------------------------------------
	ESCF_i = E_Total(dim, D, Hcore, Fock, Enuc);
	printf("  Initial E_SCF is %12.10f\n",ESCF);

	int die = 1;
  	iter=0;
	printf("  ITERATION 0:  RHF ENERGY IS %18.14f\n",ESCF_i);

    //---------------------------------------------------------------------------
    // Step #7: Compute the New Fock Matrix
    //---------------------------------------------------------------------------    

    do {

    // Update Fock matrix
    UpdateF(dim, D, Hcorec, EE, Fnew);

    // Form Fock matrix F = S^{-1/2}^t H_core S^{1/2}
    LoopMM(dim, Fnew, "n", SqrtS, "n", temp);
    LoopMM(dim, SqrtS, "t", temp, "n", Fockc);
 //   print_matrix(" Fnew ", dim, dim, Fnew, dim);
  
    // Diagonalize new Fock matrix
    DIAG_N(dim, dim, Fnew, eps, Cp);

    // Get new MO coefficients
    LoopMM(dim, SqrtS, "n", Cp, "n", C);
//    print_matrix("  Coefficients", dim, dim, C, dim);
//   print_matrix(" Initial Coeff ", dim, dim, Cp, dim);


//    print_matrix(" Dnew ", dim, dim, Dn, dim);

    // Build new Density matrix
    BuildDensity(dim, 5, C, Dn);
//    print_matrix("  New Density Matrix ", dim, dim, Dn, dim);

    // Compute new Energy
    ESCF = E_Total(dim, Dn, Hcorec, Fnew, Enuc);

    // Get RMS_D for density matrix, copy new density matrix to D array
    deltaDD = DensityDiff(dim, D, Dn);

    // get change in energy
    deltaE = ESCF - ESCF_i;
    // call current energy ESCF_i for next iteration
    ESCF_i = ESCF;
    
    if (fabs(deltaE)<tolE && deltaDD<tolDD) die=0;
    else if (iter>itermax) die=0;
    
    iter++;
    printf("  ITERATION %5i:  RHF ENERGY IS %18.14f  DeltaE is %18.14f  DeltaD is %18.14f\n",iter, ESCF, deltaE, deltaDD);

  }while(die);

    return 0;


    /*

    // Diagonalize overlap matrix
    DIAG_N(dim, dim, S, Svals, Svecs);

    // Cannot diagonalize here, fortran code dependency???
    // Works in linux...

    for (i = 0; i<dim; i++)
    {
        lambdasqrt[i*dim+i] = pow(Svals[i],-0.5);
        printf(" %12.10f\n",lambdasqrt[i*dim+i]);
    }

    print_matrix(" lambdasqrt ", dim, dim, lambdasqrt, dim);

    LoopMM(dim, lambdasqrt, "n", Svecs, "t", temp);

    LoopMM(dim, Svecs, "n", temp, "n", sqrtS);

    print_matrix( "S^1/2", dim, dim, sqrtS, dim);

    // Continue from here... 

    */
    
    /*
  
    LoopMM(dim, Hcore, "n", sqrtS, "n", temp);

    LoopMM(dim, sqrtS, "n", temp, "n", Fock);

    print_matrix(" Fock ", dim, dim, Fock, dim);

    DIAG_N(dim, dim, Fock, Fvals, Fvecs);

  //LoopMM(dim, sqrtS, "n", )

    //---------------------------------------------------------------------------
    // Step #6: Compute the Initial SCF Energy
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #7: Compute the New Fock Matrix
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #8: Build the New Density Matrix
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #9: Compute the New SCF Energy
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #10: Test for Convergence
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Step #11: Building of Hamiltonian from HF MO's for PIS
    //---------------------------------------------------------------------------

*/

}




void CrawdadFormat()
{
    int i,j,idx;
    idx = 0;

do {

    for(i=0; i<=dim*idx; i++) 
    {
        for(j=0; j<=i; j++) 
        {

            if ( i == j)
            {
                Hcore[idx*dim+idx] = T[idx];
                S[idx*dim+idx] = A[idx];
		printf(" %f %f\n",Hcore[idx*dim+idx],T[idx]);
            }
        }        
    
    }

    idx++;

}   while (idx < dim);

}


// Kinetic energy operator following crawdad labeling.
void KineticEnergyIntegrals()
{
    int i, j, imax, z, m;
    double factor;
    factor = (hbar*hbar*pi*pi) / (2*L*L);

    imax = 0;
    z = 0;

     do 
    {
        for (i=imax; i<=z; i++)
        {
            for(j=0; j<=i; j++)
            {
                if ( i == j)
                {
                    
                     T[i] = factor * (pow(NPOrb_x[i],2) + pow(NPOrb_y[i],2) + pow(NPOrb_z[i],2));
                     //printf("%i %i %f\n",i,j,T[i]);
                //     printf("for nx=%i ny=%i nz=%i phi = %i %f\n",NPOrb_x[i],NPOrb_y[i],NPOrb_z[i],i,T[i]);
                }
                else if ( i != j) 
                {   
                    
                     T[i] = 0;
                //    printf("%i %i %f\n",i,j,T[i]);
                }
            }
        }

        imax++;
        z++;

    } while ( z <= nmax*nmax*nmax-1 );
}


// Complex conjugate of phi(x) * phi(y) integrated over -infty to infty = 1 when x == y and 0 otherwise.
// Anyway, this is confusing. Thought phi(x) = sqrt(2./L)sin(pi*n*x/L) ?? Well.. thats the energy eigenfunction. Right?
void AtomicOrbitalOverlap()
{
    int i,j,imax,z;

    imax = 0;
    z=0;

    do 
    {
        for (i=imax; i<=z; i++)
        {
            for(j=0; j<=i; j++)
            {
                if ( i == j)
                {
                    E1[i] = i;
                    E2[i] = j;
                     A[i] = 1;
                    //printf("%i %i %f\n",i,j,A[i]);
                }
                else if ( i != j) 
                
                {   
                    E1[i] = i;
                    E2[i] = j;
                     A[i] = 0;
                    //printf("%i %i %f\n",i,j,A[i]);
                }
            }
        }

        imax++;
        z++;

    } while ( z <= nmax*nmax*nmax-1 );
}


void CubicPhi() 
{
    int nx, ny, nz;
    int idx, l;

    // variables to use for ordering orbitals in increasing energy.
    int cond, Ecur, swap, c, d;

    idx = 0;

    for(nx=0; nx<nmax; nx++)
    {
        for(ny=0; ny<nmax; ny++)
        {
            for(nz=0; nz<nmax; nz++)
            {
                idx = nx*nmax*nmax + ny*nmax + nz;
                l = (nx+1)*(nx+1) + (ny+1)*(ny+1) + (nz+1)*(nz+1);
                NPOrbE[idx] = l;

            }
        }
    }

    for(c=0; c < (nmax*nmax*nmax-1); c++)
    {
        for(d=0; d < (nmax*nmax*nmax-c-1); d++)
        {
            if (NPOrbE[d] > NPOrbE[d+1]) 
            {
                swap = NPOrbE[d];
                NPOrbE[d] = NPOrbE[d+1];
                NPOrbE[d+1] = swap;
            }
        }
    }

    c=0;
    do 
    {
        Ecur = NPOrbE[c];
        nx = 0;
        do 
        {
            nx++;
            ny=0;
            do 
            {
                ny++;
                nz=0;
                do 
                {
                    nz++;
                    cond = Ecur-(nx*nx + ny*ny + nz*nz);

                    if (cond == 0)
                    {
                        NPOrb_x[c] = nx;
                        NPOrb_y[c] = ny;
                        NPOrb_z[c] = nz;

                        printf(" for phi=%i, x=%i y=%i z=%i energy is %d\n",c,NPOrb_x[c],NPOrb_y[c],NPOrb_z[c],NPOrbE[c]);

                        c++;
                    }
                } while (Ecur == NPOrbE[c] && nz<nmax);
            } while (Ecur == NPOrbE[c] && ny<nmax);
        } while (Ecur == NPOrbE[c] && nx<nmax);
    } while (c<nmax*nmax*nmax);
}

void print_matrix( char* desc, int m, int n, double* a, int lna ) {
        int i, j;
        printf("\n\n----------------------");
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %12.9f", a[i*lna+j] );
                printf( "\n" );
        }
        printf("----------------------\n\n");
}



int DIAG_N(int dim, int number, double *mat, double *en, double *wfn) {
  int i,j,ind, state_max, count;
  double *pacMat, *eigval, *eigvec;

  pacMat = (double *)malloc((dim*(dim+1)/2)*sizeof(double));
  eigval = (double *)malloc(dim*sizeof(double));
  eigvec = (double *)malloc(dim*dim*sizeof(double));

  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      if (i<=j) {
        ind =  j*(j+1)/2 + i; // Position(i,j);
        pacMat[ind] = mat[i*dim+j];
      }
    }

  }

 Diagonalize(pacMat,dim,eigval,eigvec);

  count=0;
  for (i=0; i<number; i++) {
    en[i] = eigval[i];
    if (en[i]<=0) count++;
    for (j=0; j<dim; j++) {
      wfn[j*dim+i] = eigvec[i*dim+j];
    }
  }

  return count;
  free(pacMat);
  free(eigval);
  free(eigvec);


}

 void Diagonalize(double*M,long int dim, double*eigval,double*eigvec){
  integer one,info,edim,fdim,*ifail,tdim,i,il,iu,m;
  doublereal zero,vl,vu;
  doublereal tol;
  char N, U;
  doublereal *work;
  integer*iwork;
  edim = 8*dim;
  fdim = 5*dim;
  tdim = 3*dim;
  N    = 'V'; // 'N' for eigenvalues only, 'V' for eigenvectors, too
  U    = 'U';
  one  = dim;   // if N='N', one=1; otherwise, one=dim;
  work  = (doublereal*)malloc(edim*sizeof(doublereal));
  DSPEV(N,U,dim,M,eigval,eigvec,one,work,info);

  //for (i=0; i<dim; i++) printf("  Eig %i  %12.10f\n",i+1,eigval[i]);
  free(work);
} 

void LoopMM(int dim, double *a, char *transa, double *b, char *transb, double *c) {
  int i, j, k; 
  double sum;

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
  
      sum = 0.;
      for (k=0; k<dim; k++) {

        if (strcmp(transa,"t")==0 && strcmp(transb,"t")==0) {
          sum += a[k*dim+i]*b[j*dim+k];  
        }

        else if (strcmp(transa,"n")==0 && strcmp(transb,"t")==0) {
          sum += a[i*dim+k]*b[j*dim+k];
        }
        else if (strcmp(transa,"t")==0 && strcmp(transb,"n")==0) {
          sum += a[k*dim+i]*b[k*dim+j];
        }
        else {
          sum += a[i*dim+k]*b[k*dim+j];
        }
        
      }
      c[i*dim+j] = sum;
    }
  }
} 

double DensityDiff(int dim, double *D, double *Dnew) {
  int m, n;
  double sum;

  sum = 0;

  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {

      sum += (Dnew[m*dim+n]-D[m*dim+n])*(Dnew[m*dim+n]-D[m*dim+n]);
      D[m*dim+n] = Dnew[m*dim+n];

    }
  }   

  return sqrt(sum);

}

void ReadEI(int dim, FILE *fp, double *EE) {
  int i, j, k, l, ij, kl, ijkl;
  double val;

  while(fscanf(fp, "%d %d %d %d %lf",&i,&j,&k,&l,&val) !=EOF) {
    i--;
    j--;
    k--;
    l--;
    // ijkl
    //ij = i*(i+1)/2 + j;
    //kl = k*(k+1)/2 + l;
    //ijkl = ij*(ij+1)/2 + kl;
    EE[FourDIndx(i,j,k,l,dim)] = val;
    EE[FourDIndx(j,i,k,l,dim)] = val;
    EE[FourDIndx(i,j,l,k,dim)] = val;
    EE[FourDIndx(j,i,l,k,dim)] = val;
    EE[FourDIndx(k,l,i,j,dim)] = val;
    EE[FourDIndx(l,k,i,j,dim)] = val;
    EE[FourDIndx(k,l,j,i,dim)] = val;
    EE[FourDIndx(l,k,j,i,dim)] = val;
  }

}

int FourDIndx(int i, int j, int k, int l, int dim) {

  return i*dim*dim*dim+j*dim*dim+k*dim+l;

}

void UpdateF(int dim, double *D, double *Hcore, double *EI, double *Fnew) {

  int m, n, l, s, mnls, mlns;
  double sum;

  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {

      sum = 0.;
      for (l=0; l<dim; l++) {
        for (s=0; s<dim; s++) {

          mnls = FourDIndx(m, n, l, s, dim);
          mlns = FourDIndx(m, l, n, s, dim); 
          
          sum += D[l*dim+s]*(2*EI[mnls]-EI[mlns]);

        }
      }
      Fnew[m*dim+n] = Hcore[m*dim+n] +  sum;
      //Fnew[m*dim+n] = sum;
    }
  }
}



// Need to loop through all phi's. 4 orbitals, so 4 for loops.
void TwoERICalc()
{
    int a,b,c,d;
    int index;
    

    index = 0;

    for(a=0; a<NPOrbE[dim]; a++)
    {
        for(b=0; b<NPOrbE[dim]; b++)
        {
            for(c=0; c<NPOrbE[dim]; c++)
            {
                for(d=0; d<NPOrbE[dim]; d++)
                {

                    // Assign phi nx, ny & nz values for ERI format.
                    adim[0] = NPOrb_x[a];
                    adim[1] = NPOrb_y[a];
                    adim[2] = NPOrb_z[a];

                    bdim[0] = NPOrb_x[b];
                    bdim[1] = NPOrb_y[b];
                    bdim[2] = NPOrb_z[b];

                    cdim[0] = NPOrb_x[c];
                    cdim[1] = NPOrb_y[c];
                    cdim[2] = NPOrb_z[c];

                    ddim[0] = NPOrb_x[d];
                    ddim[1] = NPOrb_y[d];
                    ddim[2] = NPOrb_z[d];

                  
                    // Store ERI for values in teri array & coords.
                    ERIa[index] = a;
                    ERIb[index] = b;
                    ERIc[index] = c;
                    ERId[index] = d;
                    teri[index] = ERI(n, *xb, *w, *adim, *bdim, *cdim, *ddim);
                  
                    index++;

                }
            }
        }
    }
}



//  Arguments:  dim = number of points for gauss-legendre grid
//              xa[]  = points on gauss-legendre grid
//              w[]   = weights from gauss-legendre grid
//              a[]   = array of nx, ny, and nz for orbital a
//              b[]   = array of nx, ny, and nz for orbital b
//              c[]   = array of nx, ny, and nz for orbital c
//              d[]   = array of nx, ny, and nz for orbital d
//  This function computes the ERI (a b | c d) where a, b, c, d are
//  all associated with three unique quantum numbers (nx, ny, nz)
//  According to Gill paper, each ERI can be written as a linear combination of (p|q) 
//  integrals where p is related to (a-b) or (a+b) and q is related to (c-d) or (c+d)
//  This function automatically enumerates all the appropriate (p|q), computes them, and
//  accumulates the total... Hopefully it works!
double ERI(int dim, double *xa, double *w, double *a, double *b, double *c, double *d) {

  int i, j, k, l, m, n;
  double *x1, *x2, *y1, *y2, *z1, *z2;
  double faci, facj, fack, facl, facm, facn, fac;;
  //char *cp, *cq, *cr, *cs;
  double eri_val;
  static const char *cx1[] = {"px x1", "qx x1"};
  static const char *cx2[] = {"rx x2", "sx x2"};
  static const char *cy1[] = {"py y1", "qy y1"};
  static const char *cy2[] = {"ry y2", "sy y2"};
  static const char *cz1[] = {"pz z1", "qz z1"};
  static const char *cz2[] = {"rz z2", "sz z2"};  


  x1 = (double *)malloc(3*sizeof(double));
  x2 = (double *)malloc(3*sizeof(double));
  y1 = (double *)malloc(3*sizeof(double));
  y2 = (double *)malloc(3*sizeof(double));
  z1 = (double *)malloc(3*sizeof(double));
  z2 = (double *)malloc(3*sizeof(double));

  //x1[0] = ax-bx, x1[1] = ax+bx
  x1[0] = a[0] - b[0];
  x1[1] = a[0] + b[0];
  y1[0] = a[1] - b[1];
  y1[1] = a[1] + b[1];
  z1[0] = a[2] - b[2];
  z1[1] = a[2] + b[2];

  //x1[0] = cx-dx, x1[1] = cx+dx
  x2[0] = c[0] - d[0];
  x2[1] = c[0] + d[0];
  y2[0] = c[1] - d[1];
  y2[1] = c[1] + d[1];
  z2[0] = c[2] - d[2];
  z2[1] = c[2] + d[2];

  double tempval = 0.;
  eri_val = 0.;
  // Generate all combinations of phi_a phi_b phi_c phi_d in expanded cosine form
  for (i=0; i<2; i++) {
    faci = pow(-1,i);
    for (j=0; j<2; j++) {
      facj = pow(-1,j);
      for (k=0; k<2; k++) {
        fack = pow(-1,k);
        for (l=0; l<2; l++) { 
          facl = pow(-1,l);
          for (m=0; m<2; m++) {
            facm = pow(-1,m);
            for (n=0; n<2; n++) {
              facn = pow(-1,n);
   
              fac = faci*facj*fack*facl*facm*facn;          
             
              // Uncomment to see the functions being integrated in each call to pq_int 
              //printf(" + %f Cos[%s] Cos[%s] Cos[%s] Cos[%s] Cos[%s] Cos[%s] \n",
              //fac,cx1[n],cx2[m],cy1[l],cy2[k],cz1[j],cz2[i]);
              // recall pq_int args are -> dim, *xa, *w, px, py, pz, qx, qy, qz
              // order of indices to get these values is a bit strange, see print statement
              // for example of ordering!
              tempval = pq_int( x1[n], y1[l], z1[j], x2[m], y2[k], z2[i]);
              printf("  (%f %f %f | %f %f %f) -> %17.14f\n",x1[n], y1[l], z1[j], x2[m], y2[k], z2[i],tempval);
              eri_val += fac*tempval;
              

            }
          } 
        }
      }
    }
  }

 
  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(z1);
  free(z2);

  return eri_val;

}     



void BuildDensity(int dim, int occ, double *C, double *D) {
  int i, j, m;
  double sum;

	sum = 0.;

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      sum = 0.;
      for (m=0; m<occ; m++) {
        sum += C[i*dim+m]*C[j*dim+m];
      }
      D[i*dim+j] = sum;
    }
  }
}

double E_Total(int dim, double *D, double *Hc, double *F, double Enuc) 
{

  int m, n;  
  double sum;

  sum = 0.;
  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {
     // sum += D[m*dim+n]*(Hc[m*dim+n] + F[m*dim+n]);
            sum += D[m*dim+n]*(Hc[m*dim+n] + F[m*dim+n]);
                }
                  }
                    return sum + Enuc;
}



double pq_int(double px, double py, double pz, double qx, double qy, double qz) {

  double dr = 0.005;
  int max = 200;
  double sum = 0.;
  double num, denom;
  double x, y, z;

  for (int i=0; i<max; i++) {
    x = dr*i + 0.005;
    for (int j=0; j<max; j++) {
      y = dr*j + 0.005;
      for (int k=0; k<max; k++) {
        z = dr*k + 0.005;
        num = g_pq(px, qx, x)*g_pq(py, qy, y)*g_pq(pz, qz, z);
        denom = sqrt(x*x+y*y+z*z);
        sum += (num/denom)*dr*dr*dr;
        printf("  sum %f  x %f  y %f  z %f\n",sum, x, y, z);
      }
    }
  }       

  return (8./pi)*sum;

}

double g_pq(double p, double q, double x) {

  int d;
  d = (int)(fabs(p-q));
  double g;
  g = 0.;
  if (p == q && p == 0) {

    g = 1 - x;

  }
  else if ( p == q && p > 0 ) {

    g = (1 - x)*cos(p*pi*x)/2. - sin(p*pi*x)/(2*p*pi);

  }
  else if ( (d % 2)==0) {

    g = (q*sin(q*pi*x) - p*sin(p*pi*x))/((p*p-q*q)*pi);
  } 
  else g = 0.;

  return g;
}


