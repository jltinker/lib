/* VOID FFTW3D(int n, fftw_real *scr, int iopt, char *wisdomfile) 

   --- void fftw3d(int n, fftw_real *scr, int iopt, char *wisdomfile) 
   --- Fourier transform the 3-d array scr using FFTW

      * n = number of cells per cube side
      * scr = array to be Fourier transformed (in src3ft format)
      * iopt = 1 --> forward transform (real to complex)
             = -1 --> invsrse transform (complex to real)
      * wisdomfile = name wisdom file

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

/* #include  "rfftw.h" */
#include  "fftw3.h"
#define indf(a,b,c) (a)*(n+2)*n+(b)*(n+2)+(c)
#define ind(a,b,c) (a)*n*n+(b)*n+(c)

void fftw3d(int n, float *scr, int iopt, char *wisdomfile)
{
  static fftw_plan p,pinv; /* keep plans from one call to next */
  static int icall=1,n_plan=0,n_planinv=0;
  FILE *fp1;
  float vol;
  int i,j,k,indx;
  fftw_complex *A;
  fftw_real *B;
  vol = (float)n*n*n;

  if(icall==1)
    { 
      icall++;

      /* Read Wisdom from wisdomfile, if given */

      if(strcmp(wisdomfile,"")!=0)
	{/* read wisdom */
	  fp1=fopen(wisdomfile,"r");
	  if(FFTW_FAILURE==fftw_import_wisdom_from_file(fp1))
	    fprintf(stderr,"Error reading wisdom\n");
	  fclose(fp1);
	}
      
    }      

  if(iopt==1) /* Forward FFT */
    {
      /* create plan if it is not present*/
      
      if(p==NULL || n != n_plan) /*check if n changed? if yes make new plan*/
	{
	  fprintf(stderr,"fftw3d> Creating plan for fwd. FFT n=%d\n",n);
	  p = fftw_plan_dft_r2c_3d(n,n,n,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
	  n_plan=n; /* keeps track of value of n for current plan*/
	}
      B=(fftw_real*) scr;
      rfftwnd_one_real_to_complex(p,B,NULL);
      scr=(float*) B;
    }

  if(iopt==-1) /* Inverse FFT */
    {
      /* create plan if it is not present*/

      if(pinv==NULL || n != n_planinv)
	{
	  fprintf(stderr,"fftw3d> Creating plan for inverse FFT n=%d\n",n);
	  pinv = fftw_plan_dft_c2r_3d(n,n,n,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
	  n_planinv=n;
	}
      A=(fftw_complex*) scr;
      rfftwnd_one_complex_to_real(pinv,A,NULL);
      scr=(float*) A;

      for(i=0;i<n;i++)  /* this normalizes the inverse transform */
	  for(j=0;j<n;j++)
	    for(k=0;k<n;k++)
	      {
	      indx=indf(i,j,k);
	      *(scr+indx) = *(scr+indx)/vol;
	      }
    }

  return;
}




