c--- sscal -- scale a vector by a scalar constant
c--- sscal(n,a,x,incx)
c---   * n = number of elements
c---   * a = scalar constant
c---   * x() = vector
c---   * incx = increment of x()

      subroutine sscal(n,a,x,incx)
      real x(*)

      do 10 i=0,n-1
	j=1+i*incx
	x(j) = a*x(j)
10    continue
      return
      end
