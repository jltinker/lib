/* cwrite_() -- write an unformatted fortran record to stdout
   cwrite_(ptr,obytes)
      * ptr = memory location for start of write (e.g. name of fortran
	      variable or array)
      * obytes = (pointer to) number of bytes to be written, on return
		 contains number of bytes successfully written
*/
#include <stdio.h>

cwrite_(ptr,obytes) 
char *ptr ;
int *obytes ;

{
    int nbytes, nitem1 ;
    int errno ;
    unsigned nitems,
	     size = 1 ;

    nitems = *obytes ;
    nbytes = size*nitems ;
    if ( fwrite(&nbytes,sizeof(int),1,stdout) != 1 )  {
	fprintf(stderr,"write error, is the file open ? \n") ;
	}
    nitem1 = fwrite(ptr,size,nitems,stdout) ;
    if ( nitem1 != nitems ) {
	fprintf(stderr,"write error, %d items requested, %d items written. \n",
		nitems,nitem1) ;
    }
    if ( fwrite(&nbytes,sizeof(int),1,stdout) != 1 )  {
	fprintf(stderr,"write error on second byte label \n") ;
    }

   *obytes = nitems ;
}
