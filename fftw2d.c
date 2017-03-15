/* VOID FFTW2D(int n, fftw_real *scr, int iopt, char *wisdomfile) 

   --- void fftw2d(int n, fftw_real *scr, int iopt, char *wisdomfile) 
   --- Fourier transform the 2-d array scr using FFTW

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

#include  "rfftw.h"
#include  "fftw.h"
#define ind2d(a,b) (a)*n+(b)
#define ind2df(a,b) (a)*(n+2)+(b)

void fftw2d(int n, float *scr, int iopt, char *wisdomfile)
{
  static rfftwnd_plan p,pinv; /* keep plans from one call to next */
  static int icall=1,n_plan=0,n_planinv=0;
  FILE *fp1;
  float vol;
  int i,j,k,indx;
  fftw_complex *A;
  fftw_real *B;
  vol = (float)n*n;

  if(icall==1)
    { 
      icall++;

      /* Read Wisdom from wisdomfile, if given */

      fprintf(stderr,"fftw2d> wisdomfile = %s\n",wisdomfile);

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
	  fprintf(stderr,"fftw2d> Creating plan for fwd. FFT n=%d\n",n);
	  p = rfftw2d_create_plan(n,n,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
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
	  fprintf(stderr,"fftw2d> Creating plan for inverse FFT n=%d\n",n);
	  pinv = rfftw2d_create_plan(n,n,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
	  n_planinv=n;
	}
      A=(fftw_complex*) scr;
      rfftwnd_one_complex_to_real(pinv,A,NULL);
      scr=(float*) A;

      for(i=0;i<n;i++)  /* this normalizes the inverse transform */
	  for(j=0;j<n;j++)
	      {
	      indx=ind2df(i,j);
	      *(scr+indx) = *(scr+indx)/vol;
	      }
    }

  return;
}




