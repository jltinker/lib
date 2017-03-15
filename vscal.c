/* PROGRAM VSCAL

--- vscal(npr,x,a,incx)
--- scale the elements of a vector by a scalar a
   * npr = nuber of cells per cube side
   * x = array x
   * a = scaling constant.
   * incx = increment of x()
   Note : for a 3-d array, npr=n*n*n

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

void vscal(npr,x,a,incx)
int npr,incx;
float *x,a;

{
   int i;

   for(i=0;i<npr;i+=incx)
     (*(x+i))*=a;
   return;
}
