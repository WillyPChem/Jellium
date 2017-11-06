#include<stdio.h>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>
#include<time.h>
#include<string.h>

void CubicPhi();


// Peter Gill Two Electron Repulsion Integral 
void TwoERICalc(int number, double *xa, double *wa);
double ERI(int dim, double *xa, double *w, double *a, double *b, double *c, double *d);
double g_pq(double p, double q, double r);
double pq_int(int dim, double *x, double *w, double px, double py, double pz, double qx, double qy, double qz);
double E0_Int(int dim, double *xa, double *w);
double Vab_Int(int dim, double *xa, double *w, double *a, double *b);
// Gauss-Legendre quadrature functions
void legendre_compute_glr ( int n, double x[], double w[] );
void legendre_compute_glr0 ( int n, double *p, double *pp );
void legendre_compute_glr1 ( int n, double *roots, double *ders );
void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
void legendre_handle ( int n, double a, double b );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
double rk2_leg ( double t, double tn, double x, int n );
void timestamp ( void );
double ts_mult ( double *u, double h, int n );
double wtime ( );

double *adim, *bdim, *cdim, *ddim;
double *ERIa, *ERIb, *ERIc, *ERId, *teri;
int nmax, dim;


// Phi Variables

int pi;
int *NPOrbE, *NPOrb_x, *NPOrb_y, *NPOrb_z;
int n;

int main() 
{
    dim = 7;

    NPOrb_x = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // X corresponding to phi
    NPOrb_y = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // Y corresponding to phi
    NPOrb_z = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // Z corresponding to phi
    NPOrbE = (int *)malloc(nmax*nmax*nmax*sizeof(int)); // energy corresponding

    // Initialize 2eri arrays
    adim = (double *)malloc(3*sizeof(double));
    bdim = (double *)malloc(3*sizeof(double));
    cdim = (double *)malloc(3*sizeof(double));
    ddim = (double *)malloc(3*sizeof(double));

    // Storage of 2eri values after ERI function.
    ERIa = (double *)malloc(dim*dim*dim*dim*sizeof(double));
    ERIb = (double *)malloc(dim*dim*dim*dim*sizeof(double));
    ERIc = (double *)malloc(dim*dim*dim*dim*sizeof(double));
    ERId = (double *)malloc(dim*dim*dim*dim*sizeof(double));
    teri = (double *)malloc(dim*dim*dim*dim*sizeof(double));

    n = 50;


     double *x, *w;

    x = (double *)malloc(n*sizeof(double));
    w = (double *)malloc(n*sizeof(double));
   

    legendre_compute_glr(n, x, w);
    rescale (0, 1, n, x, w); 

     // Definition of pi
    pi = 4.*atan(1.0);

    nmax = 4;

    CubicPhi();
    TwoERICalc(n, x, w);


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

void TwoERICalc(int number, double *xa, double *wa)
{
    int a,b,c,d;
    int index;
    
    a = 0;
    b = 1;

    index = 0;


    for(a=0; a<nmax; a++)
    {
        for(b=0; b<nmax; b++) 
        {

           //

        }
    }

    for(a=0; a<=nmax; a++)
    {
        for(b=0; b<=nmax; b++)
        {
            for(c=0; c<=nmax; c++)
            {
                for(d=0; d<=nmax; d++)
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

                /*

                    adim[0] = 1;
                    adim[1] = 1;
                    adim[2] = 1;

                    bdim[0] = 1;
                    bdim[1] = 1;
                    bdim[2] = 1;

                    cdim[0] = 1;
                    cdim[1] = 1;
                    cdim[2] = 1;

                    ddim[0] = 1;
                    ddim[1] = 1;
                    ddim[2] = 1;

                 */

                    // printf(" %i - %i - %i - %i : %f \n",a,b,c,d,teri[index]);

                  
                    // Store ERI for values in teri array & coords.
                    bool dup = 0;


                    if(a != b && c != d && a == c && b == d)
                    {
                        for(int i=0; i<index; i++)
                        {
                            // Determines if 1,3,1,3 == 3,1,3,1 b/c same integral. 

                            if (a != ERIb[i] && b != ERIa[i] && c != ERId[i] && d != ERIc[i])
                            {

                               dup = 1;
                            }

                        else 
                        {
                              dup = 0;
                              i = index+1;
                        }

                        }

                    if (dup = 1)
                    {

                        ERIa[index] = a;
                        ERIb[index] = b;
                        ERIc[index] = c;
                        ERId[index] = d;
                        teri[index] = ERI(number, xa, wa, adim, bdim, cdim, ddim);
                        printf(" %i - %i - %i - %i : %f \n",a,b,c,d,teri[index]);
                        dup = 0;
                        index++;
                    }


                    } 
                     else if ( a == b && a == c && a == d) 
                        {

                        ERIa[index] = a;
                        ERIb[index] = b;
                        ERIc[index] = c;
                        ERId[index] = d;
                        teri[index] = ERI(number, xa, wa, adim, bdim, cdim, ddim);
                        printf(" %i - %i - %i - %i : %f \n",a,b,c,d,teri[index]);
                        index++;
                               
                        } else 
                        {
                        
                        ERIa[index] = a;
                        ERIb[index] = b;
                        ERIc[index] = c;
                        ERId[index] = d;
                        teri[index] = ERI(number, xa, wa, adim, bdim, cdim, ddim);
                        printf(" %i - %i - %i - %i : %f \n",a,b,c,d,teri[index]);
                        index++;
                        
                        }

                }
            }
        }
    }
}



/******************************************************************************/

void legendre_compute_glr ( int n, double x[], double w[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   19 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, int N, the order.
 *
 *                                                         Output, double X[N], the abscissas.
 *
 *                                                             Output, double W[N], the weights.
 *                                                             */
{
  int i;
  double p;
  double pp;
  double w_sum;
/*
 *   Get the value and derivative of the N-th Legendre polynomial at 0.
 *   */
  legendre_compute_glr0 ( n, &p, &pp );
/*
 *   Either zero is a root, or we have to call a function to find the first root.
 *   */  
  if ( n % 2 == 1 )
  {
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }
  else
  {
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
/*
 *   Get the complete set of roots and derivatives.
 *   */
  legendre_compute_glr1 ( n, x, w );
/*
 *   Compute the weights.
 *   */
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr0 ( int n, double *p, double *pp )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   19 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, int N, the order of the Legendre polynomial.
 *
 *                                                         Output, double *P, *PP, the value of the N-th Legendre polynomial
 *                                                             and its derivative at 0.
 *                                                             */
{
  double dk;
  int k;
  double pm1;
  double pm2;
  double ppm1;
  double ppm2;

  pm2 = 0.0;
  pm1 = 1.0;
  ppm2 = 0.0;
  ppm1 = 0.0;

  for ( k = 0; k < n; k++ )
  {
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr1 ( int n, double *x, double *ders )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
 *
 *         Discussion:
 *
 *             This routine requires that a starting estimate be provided for one
 *                 root and its derivative.  This information will be stored in entry
 *                     (N+1)/2 if N is odd, or N/2 if N is even, of ROOTS and DERS.
 *
 *                       Licensing:
 *
 *                           This code is distributed under the GNU LGPL license. 
 *
 *                             Modified:
 *
 *                                 19 October 2009
 *
 *                                   Author:
 *
 *                                       Original C version by Nick Hale.
 *                                           This C version by John Burkardt.
 *
 *                                             Reference:
 *
 *                                                 Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                                     A fast algorithm for the calculation of the roots of special functions, 
 *                                                         SIAM Journal on Scientific Computing,
 *                                                             Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                               Parameters:
 *
 *                                                                   Input, int N, the order of the Legendre polynomial.
 *
 *                                                                       Input/output, double X[N].  On input, a starting value
 *                                                                           has been set in one entry.  On output, the roots of the Legendre 
 *                                                                               polynomial.
 *
 *                                                                                   Input/output, double DERS[N].  On input, a starting value
 *                                                                                       has been set in one entry.  On output, the derivatives of the Legendre 
 *                                                                                           polynomial at the zeros.
 *
 *                                                                                             Local Parameters:
 *
 *                                                                                                 Local, int M, the number of terms in the Taylor expansion.
 *                                                                                                 */
{
  double dk;
  double dn;
  double h;
  int j;
  int k;
  int l;
  int m = 30;
  int n2;
  const double pi = 3.141592653589793;
  int s;
  double *u;
  double *up;
  double xp;

  if ( n % 2 == 1 )
  {
    n2 = ( n - 1 ) / 2;
    s = 1;
  }
  else
  {
    n2 = n / 2;
    s = 0;
  }

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;

  for ( j = n2; j < n - 1; j++ )
  {
    xp = x[j];

    h = rk2_leg ( pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = ders[j];

    up[0] = 0.0;
    up[1] = u[2];

    for ( k = 0; k <= m - 2; k++ )
    {
      dk = ( double ) k;

      u[k+3] = 
      ( 
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 )
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }

    for ( l = 0; l < 5; l++ )
    { 
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    ders[j+1] = ts_mult ( up, h, m-1 );
  }

  free ( u );
  free ( up );

  for ( k = 0; k < n2 + s; k++ )
  {
    x[k] = - x[n-k-1];
    ders[k] = ders[n-k-1];
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr2 ( double pn0, int n, double *x1,  double *d1 )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR2 finds the first real root.
 *
 *         Discussion:
 *
 *             This routine is only called if N is even.
 *
 *                 Thanks to Morten Welinder, for pointing out a typographical error
 *                     in indexing, 17 May 2013.
 *
 *                       Licensing:
 *
 *                           This code is distributed under the GNU LGPL license. 
 *
 *                             Modified:
 *
 *                                 17 May 2013
 *
 *                                   Author:
 *
 *                                       Original C version by Nick Hale.
 *                                           This C version by John Burkardt.
 *
 *                                             Reference:
 *
 *                                                 Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                                     A fast algorithm for the calculation of the roots of special functions, 
 *                                                         SIAM Journal on Scientific Computing,
 *                                                             Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                               Parameters:
 *
 *                                                                   Input, double PN0, the value of the N-th Legendre polynomial at 0.
 *
 *                                                                       Input, int N, the order of the Legendre polynomial.
 *
 *                                                                           Output, double *X1, the first real root.
 *
 *                                                                               Output, double *D1, the derivative at X1.
 *
 *                                                                                 Local Parameters:
 *
 *                                                                                     Local, int M, the number of terms in the Taylor expansion.
 *                                                                                     */
{
  double dk;
  double dn;
  int k;
  int l;
  int m = 30;
  const double pi = 3.141592653589793;
  double t;
  double *u;
  double *up;

  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;
/*
 *   U[0] and UP[0] are never used.
 *     U[M+1] is set, but not used, and UP[M] is set and not used.
 *       What gives?
 *       */
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;
 
  for ( k = 0; k <= m - 2; k = k + 2 )
  {
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / ( dk + 1.0 ) / ( dk + 2.0 );
 
    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }
  
  for ( l = 0; l < 5; l++ )
  {
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  free ( u );
  free ( up) ;

  return;
}
/******************************************************************************/

void legendre_handle ( int n, double a, double b )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Input, int N, the order of the rule.
 *
 *                                   Input, double A, B, the left and right endpoints.
 *                                   */ 
{
  int i;
  char output_r[255];
  char output_w[255];
  char output_x[255];
  double *r;
  double t;
  double *w;
  double *x;

  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );

  r[0] = a;
  r[1] = b;
/*
 *   Compute the rule.
 *   */
  t = wtime ( );
  legendre_compute_glr ( n, x, w );
  t = wtime ( ) - t;

  printf ( "\n" );
  printf ( "  Elapsed time during computation was %g seconds.\n", t );
/*
 *   Rescale the rule to [A,B].
 *   */
  rescale ( a, b, n, x, w );
/*
 *   Write the rule to 3 files.
 *   */
  sprintf ( output_w, "leg_o%d_w.txt", n );
  sprintf ( output_x, "leg_o%d_x.txt", n );
  sprintf ( output_r, "leg_o%d_r.txt", n );

  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", output_w );
  printf ( "  Abscissa file will be \"%s\".\n", output_x );
  printf ( "  Region file will be   \"%s\".\n", output_r );
            
  r8mat_write ( output_w, 1, n, w );
  r8mat_write ( output_x, 1, n, x );
  r8mat_write ( output_r, 1, 2, r );

  free ( r );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       R8MAT_WRITE writes an R8MAT file.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   01 June 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Input, char *OUTPUT_FILENAME, the output filename.
 *
 *                                   Input, int M, the spatial dimension.
 *
 *                                       Input, int N, the number of points.
 *
 *                                           Input, double TABLE[M*N], the table data.
 *                                           */
{
  int i;
  int j;
  FILE *output;
/*
 *   Open the file.
 *   */
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
 *   Write the data.
 *   */
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
 *   Close the file.
 *   */
  fclose ( output );

  return;
}
/******************************************************************************/

void rescale ( double a, double b, int n, double x[], double w[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         Original MATLAB version by Nick Hale.
 *                             C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, double A, B, the endpoints of the new interval.
 *
 *                                                         Input, int N, the order.
 *
 *                                                             Input/output, double X[N], on input, the abscissas for [-1,+1].
 *                                                                 On output, the abscissas for [A,B].
 *
 *                                                                     Input/output, double W[N], on input, the weights for [-1,+1].
 *                                                                         On output, the weights for [A,B].
 *                                                                         */
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
/******************************************************************************/

double rk2_leg ( double t1, double t2, double x, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       RK2_LEG advances the value of X(T) using a Runge-Kutta method.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Parameters:
 *
 *                                   Input, double T1, T2, the range of the integration interval.
 *
 *                                       Input, double X, the value of X at T1.
 *
 *                                           Input, int N, the number of steps to take.
 *
 *                                               Output, double RK2_LEG, the value of X at T2.
 *                                               */
{
  double f;
  double h;
  int j;
  double k1;
  double k2;
  int m = 10;
  double snn1;
  double t;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );

  t = t1;

  for ( j = 0; j < m; j++ )
  {
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;

    t = t + h;

    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );   
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       TIMESTAMP prints the current YMDHMS date as a time stamp.
 *
 *         Example:
 *
 *             31 May 2001 09:45:54 AM
 *
 *               Licensing:
 *
 *                   This code is distributed under the GNU LGPL license. 
 *
 *                     Modified:
 *
 *                         24 September 2003
 *
 *                           Author:
 *
 *                               John Burkardt
 *
 *                                 Parameters:
 *
 *                                     None
 *                                     */
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double ts_mult ( double *u, double h, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       TS_MULT evaluates a polynomial.
 *
 *         Discussion:
 *
 *             TS_MULT = U[1] + U[2] * H + ... + U[N] * H^(N-1).
 *
 *               Licensing:
 *
 *                   This code is distributed under the GNU LGPL license. 
 *
 *                     Modified:
 *
 *                         17 May 2013
 *
 *                           Author:
 *
 *                               Original C version by Nick Hale.
 *                                   This C version by John Burkardt.
 *
 *                                     Parameters:
 *
 *                                         Input, double U[N+1], the polynomial coefficients.
 *                                             U[0] is ignored.
 *
 *                                                 Input, double H, the polynomial argument.
 *
 *                                                     Input, int N, the number of terms to compute.
 *
 *                                                         Output, double TS_MULT, the value of the polynomial.
 *                                                         */
{
  double hk;
  int k;
  double ts;
  
  ts = 0.0;
  hk = 1.0;
  for ( k = 1; k<= n; k++ )
  {
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}
/******************************************************************************/

double wtime ( void )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       WTIME estimates the elapsed wall clock time.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   21 October 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Output, double WTIME, the current elapsed wall clock time.
 *                               */
{
  double now;

  now = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC; 

  return now;
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


// From Eq. 3.6 in the Peter Gill paper, 
// the Vab integrals are -1/pi^3 \int (phi_a phi_b )/(|r1 - r2|) dr1 dr2
// This essentially leads to (p|q) integrals with 
// px = a_x - b_x
// py = a_y - b_y
// pz = a_z - b_z
// qx = a_x + b_x
// qy = a_y + b_y
// qz = a_z + b_z
// Note the expansion of the trigonetric identity:
// Cos[px x1] Cos[py y1] Cos[pz z1] - Cos[qx x1] Cos[py y1] Cos[pz z1] - 
// Cos[px x1] Cos[qy y1] Cos[pz z1] + Cos[qx x1] Cos[qy y1] Cos[pz z1] -
// Cos[px x1] Cos[py y1] Cos[qz z1] + 
// Cos[qx x1] Cos[py y1] Cos[qz z1] + Cos[px x1] Cos[qy y1] Cos[qz z1] -
// Cos[qx x1] Cos[qy y1] Cos[qz z1]
// In order to be consistent with the defintiion of the (p|q) integrals, 
// the term Cos[px x1] Cos[py y1] Cos[pz z1] -> Cos[px x1] Cos[py y1] Cos[pz z1] Cos[0 x2] Cos[0 y2] Cos[0 z2]
// In terms of how the pq_int function is called for the above integral, it should be
// pq_int(dim, xa, w, px, py, pz, 0, 0, 0)

double Vab_Int(int dim, double *xa, double *w, double *a, double *b){
  double px, py, pz, qx, qy, qz;
  double Vab;
  px = a[0] - b[0];
  py = a[1] - b[1];
  pz = a[2] - b[2];
  qx = a[0] + b[0];
  qy = a[1] + b[1];
  qz = a[2] + b[2];

  Vab = 0.;
  // Cos[px x1] Cos[py y1] Cos[pz z1]
  Vab += pq_int(dim, xa, w, px, py, pz, 0, 0, 0);
  // - Cos[qx x1] Cos[py y1] Cos[pz z1]
  Vab -= pq_int(dim, xa, w, 0,  py, pz, qx,0, 0);
  // - Cos[px x1] Cos[qy y1] Cos[pz z1]
  Vab -= pq_int(dim, xa, w, px, 0, pz, 0, qy, 0);
  // + Cos[qx x1] Cos[qy y1] Cos[pz z1]
  Vab += pq_int(dim, xa, w, 0, 0, pz, qx, qy, 0);   
  // -Cos[px x1] Cos[py y1] Cos[qz z1]  
  Vab -= pq_int(dim, xa, w, px, py, 0, 0, 0, qz);
  // +Cos[qx x1] Cos[py y1] Cos[qz z1] 
  Vab += pq_int(dim, xa, w, 0, py, 0, qx, 0, qz);
  // + Cos[px x1] Cos[qy y1] Cos[qz z1]
  Vab += pq_int(dim, xa, w, px, 0, 0, 0, qy, qz);
  // -Cos[qx x1] Cos[qy y1] Cos[qz z1]
  Vab -= pq_int(dim, xa, w, 0, 0, 0, qx, qy, qz);

  return -Vab/(pi*pi*pi);

}

// the E integral is 1/pi^6 \int 1/(|r1 - r2|) dr1 dr2
// which is equivalent to 
// 1/pi^6 \int cos(0 x1) cos(0 y1) cos(0 z1) cos(0 x2) cos(0 y2) cos(0 z2)/|r1-r2| dr1 dr2
//
double E0_Int(int dim, double *xa, double *w) {
  double Eint;

  Eint = pq_int(dim, xa, w, 0, 0, 0, 0, 0, 0);
  return Eint/(pow(pi,6));

}

// TWO ELECTRON REPULSION CALCULATIONS

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
             // printf(" + %f Cos[%s] Cos[%s] Cos[%s] Cos[%s] Cos[%s] Cos[%s] \n",
            //  fac,cx1[n],cx2[m],cy1[l],cy2[k],cz1[j],cz2[i]);
              // recall pq_int args are -> dim, *xa, *w, px, py, pz, qx, qy, qz
              // order of indices to get these values is a bit strange, see print statement
              // for example of ordering!
              tempval = pq_int(dim, xa, w, x1[n], y1[l], z1[j], x2[m], y2[k], z2[i]);
            //  printf("  (%f %f %f | %f %f %f) -> %17.14f\n",x1[n], y1[l], z1[j], x2[m], y2[k], z2[i],tempval);
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

// This function implements Eq. 4.7 and 4.8 in Peter Gills paper on 2-electrons in a cube
// Gauss-Legendre quadrature is used for the 3d integral on the range 0->1 for x, y, and z
// int dim is the number of points on this grid, double *xa is a vector containing the actual points on this grid, and
// double *w is a vector containing the weights associated with this grid (analogous to differential length elements
// in rectangle rule integration).
// double px, py, pz, qx, qy, qz has the same interpretation as it does in the Gill paper.
double pq_int(int dim, double *xa, double *w, double px, double py, double pz, double qx, double qy, double qz) {

  double sum = 0.;
  double num, denom;
  double x, y, z, dx, dy, dz;
  double gx, gy, gz;
  for (int i=0; i<dim; i++) {
    x = xa[i];
    dx = w[i];
    gx = g_pq(px, qx, x);
    for (int j=0; j<dim; j++) {
      y = xa[j];
      dy = w[j];
      gy = g_pq(py, qy, y);
      for (int k=0; k<dim; k++) {
        z = xa[k];
        dz = w[k];
        gz = g_pq(pz, qz, z);
        num = gx*gy*gz;
        denom = sqrt(x*x+y*y+z*z);
        sum += (num/denom)*dx*dy*dz;
        //printf("  sum %f  x %f  y %f  z %f\n",sum, x, y, z);
      }
    }
  }

  return (8./pi)*sum;

}  
