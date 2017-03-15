#include <stdio.h>
#include <stdlib.h>

/* This takes a file and reads the number of lines in it,
 * rewinds the file and returns the lines.
 */

int filesize(FILE *fp)
{
  int i=-1;
  char a[1000];

  while(!feof(fp))
    {
      i++;
      fgets(a,1000,fp);
    }
  rewind(fp);
  return(i);
}

