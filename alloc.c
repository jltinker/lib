/*
 * Arrayalloc.c -- routines to provide dynamic memory allocation for 2 and 3
 * dimensional arrays.  The idea is to implement vectored arrays in such a
 * way that they are compatible with normal multi-dimension arrays as far
 * as syntax goes.
 *
 * $Log:	arrayalloc.c,v $
 * Revision 1.4  85/08/12  12:36:27  roy
 * *** empty log message ***
 * 
 * Revision 1.4  85/08/12  12:36:27  roy
 * Random commenting and de-lintifying.
 * 
 * Revision 1.3  85/08/12  12:12:51  roy
 * Added random comments in preperation for distribution.
 * 
 * Revision 1.2  85/08/06  18:30:27  roy
 * Added debugging statments in arrayalloc().
 * 
 * Revision 1.1  85/08/05  18:45:20  roy
 * Initial revision
 * 
 */

# include <stdio.h>

static char *rcsid = "$Header: arrayalloc.c,v 1.4 85/08/12 12:36:27 roy Rel $";

/*
 * Arrayalloc () -- allocate an imax by jmax vectored array of "size" byte
 * elements.  If memory can't be allocated, either for the main array or for
 * the row address vector, we return NULL.  See accompanying documentation
 * for more details.
 */
char **arrayalloc (imax, jmax, size)
unsigned imax, jmax, size;
{
	char *calloc(), *malloc ();
	register char **vector, *array;
	register int k, stride;

	/*
	 * Get memory for main array.
	 */
	if ((array = calloc (imax * jmax , size)) == NULL)
		return (NULL);

# ifdef DEBUG
	printf ("array = %x\n", array);
# endif

	/*
	 * Get memory for intermediate row address vector.
	 */
	if ((vector = (char **) malloc (imax * sizeof (char *))) == NULL) {
		free(array) ;
		return (NULL);
	}

	/*
	 * Initialize the address vector so each element points to the
	 * first element in the corresponding row in the main array.
	 */
	stride = jmax * size;
	for (k = 0; k < imax; k++)
	{
		vector [k] = &array [k*stride];
# ifdef DEBUG
		printf ("vector [%d] = %x\n", k, vector[k]);
# endif
	}

	return (vector);
}

/*
 * Arrayfree () -- free the memory acquired from arrayalloc ().  No checks
 * are made to make sure things are as they should be, so it is the user's
 * responsibility to make sure that you don't arrayfree() anything that you
 * didn't arrayalloc() in the first place.  Eventually, checks will be added
 * to make sure the user hasn't screwed things up.  We have to first free the
 * real array memory, and then free the intermediate vector.  This sounds
 * more complicated than it really is.
 */
void arrayfree (v)
char **v;
{
	free (v[0]);
	free ((char *) v);
}
/*
 * A3alloc () -- allocate an imax by jmax by kmax vectored array of "size" byte
 * elements.  If memory can't be allocated, either for the main array or for
 * the row address vector, we return NULL.  See accompanying documentation
 * for more details.
 */
char ***a3alloc (imax, jmax, kmax, size)
unsigned imax, jmax, kmax, size;
{
	char *malloc(), *calloc ();
	register char ***column, **vector, *array;
	register int k, j, stride;

	/*
	 * Get memory for main array.
	 */
	if ((array = calloc (imax * jmax * kmax * size,  sizeof (char)))
	    							== NULL)
		return (NULL);

# ifdef DEBUG
	printf ("array = %x\n", array);
# endif

	/*
	 * Get memory for intermediate row address vector.
	 */
	if ((vector = (char **) malloc (imax*jmax*sizeof(char *))) == NULL) {
		free(array) ;
		return (NULL);
	}

	/*
	 * Get memory for intermediate column address vector.
	 */
	if ((column = (char ***) malloc (imax * sizeof (char *))) == NULL) {
		free(array) ;
		free (vector);
		return (NULL);
	}

	/*
	 * Initialize the address vector so each element points to the
	 * first element in the corresponding row in the main array.
	 */
	stride = kmax * size;
	for (k = 0; k < imax; k++)
	{
	    column [k] = &vector [k*jmax];
	    for (j = 0;  j < jmax;  j++) {
	    	column[k][j] = &array[(k*jmax + j) * stride];
# ifdef DEBUG
		printf ("vector [%d] = %x\n", k, vector[k]);
# endif
	    }
	}

	return (column);
}

/*
 * A3free () -- free the memory acquired from a3alloc ().  No checks
 * are made to make sure things are as they should be, so it is the user's
 * responsibility to make sure that you don't a3free() anything that you
 * didn't a3alloc() in the first place.  Eventually, checks will be added
 * to make sure the user hasn't screwed things up.  We have to first free the
 * real array memory, and then free the intermediate vector.  This sounds
 * more complicated than it really is.
 */
void a3free (v)
char ***v;
{
        free (v[0][0]);
	free ((char *) v[0]);
	free ((char *) v);
}
