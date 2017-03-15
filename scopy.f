c--- scopy -- copy vector x into vector y
c--- scopy(n,x,incx,y,incy)
c---   * n = number of elements to copy
c---   * x = x vector
c---   * incx = increment of x vector
c---   * y = y vector
c---   * incy = increment of y vector

      subroutine scopy(n,x,incx,y,incy)
      real x(*),y(*)

      do 10 i=0,n-1
	y(1+i*incy)=x(1+i*incx)
10    continue
      return
      end
