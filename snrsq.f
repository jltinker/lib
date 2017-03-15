c--- snrsq -- sum the squares of the elements of a vector
c--- snrsq(n,x,incx)
c---   * n = number of elements
c---   * x() = vector
c---   * incx = increment between vector elements

      real function snrsq(n,x,incx)
      real x(*)
      snrsq = 0.0
      do 10 i=0,n-1
	snrsq = snrsq + x(1+i*incx)**2
10    continue
      return
      end

