c--- src3ft -- real-to-complex 3-d FFT; new version, allows factors of 3, 5
c--- src3ft(x,l1,l2,l3,ldx,mdx,iopt,ier)
c---   * x = data array:  real *4 x(ldx,mdx,l3)
c---   * l1 = number of rows of data, power of 2
c---   * l2 = number of columns of data, power of 2
c---   * l3 = number of planes of data, power of 2
c---   * ldx = leading dimension of x, ldx >= l1+2
c---   * mdx = middle dimension of x, mdx >= l2
c---   * iopt = option flag:
c---      >= 0   compute forward transform (real to half complex)
c---      < 0    compute inverse transform (half complex to real)
c---   * ier = error flag, returns -1 if incorrectly called

      subroutine src3ft(x,l1,l2,l3,ldx,mdx,iopt,ier)
      parameter (nmax=400)
      real x(*),speq(2*nmax*nmax)

c--- check requirements on dimensions
      ier=0
      if ((ldx.lt.l1+2).or.(mdx.lt.l2).or.(l2.gt.nmax).or.(l3.gt.nmax))
     >then
	ier=-1
	write(0,*) 'FFT: error in array dimensions'
	stop
	return
      end if

c--- realfft() assumes that data completely fill the arrays and it takes 
c--- separately a full 3-d array and a 2*2-d array for the Nyquist frequency
c--- components.  Some reordering before and after the call to realfft()
c--- is therefore needed, so that the i/o matches src3ft.  Unfortunately, there
c--- seems to be no practical way to carry out this reordering without the
c--- supplementary matrix speq(), even though there is enough actual space in
c--- the input array if it is correctly dimensioned.

c--- forward transform, starting with a real matrix
      if (iopt.ge.0) then
	isign=1
        do 20 i3=0,l3-1
	do 20 i2=0,l2-1
	do 20 i1=0,l1-1
 	  x(1+i1+i2*l1+i3*l1*l2) = x(1+i1+i2*ldx+i3*ldx*mdx)
20      continue
	call realfft(x,l1,l2,l3,speq,isign)
	do 30 i3=l3-1,0,-1
	do 30 i2=l2-1,0,-1
	do 30 i1=l1-1,0,-1
	  x(1+i1+i2*ldx+i3*ldx*mdx) = x(1+i1+i2*l1+i3*l1*l2)
30      continue
	do 40 i3=0,l3-1
	do 40 i2=0,l2-1
	  x(1+l1+i2*ldx+i3*ldx*mdx) = speq(1+i2*2+i3*2*l2)
	  x(2+l1+i2*ldx+i3*ldx*mdx) = speq(2+i2*2+i3*2*l2)
40      continue

c--- inverse transform, starting with a half-complex matrix
      else
	isign=-1
	do 50 i3=0,l3-1
	do 50 i2=0,l2-1
	  speq(1+i2*2+i3*2*l2) = x(1+l1+i2*ldx+i3*ldx*mdx)
	  speq(2+i2*2+i3*2*l2) = x(2+l1+i2*ldx+i3*ldx*mdx)
50      continue
	do 60 i3=0,l3-1
	do 60 i2=0,l2-1
	do 60 i1=0,l1-1
	  x(1+i1+i2*l1+i3*l1*l2) = x(1+i1+i2*ldx+i3*ldx*mdx)
60      continue
	call realfft(x,l1,l2,l3,speq,isign)
	xfac = 2./float(l1*l2*l3)
	do 70 i3=l3-1,0,-1
	do 70 i2=l2-1,0,-1
	do 70 i1=l1-1,0,-1
	  x(1+i1+i2*ldx+i3*ldx*mdx) = x(1+i1+i2*l1+i3*l1*l2)
70      continue

      end if

      return
      end


      subroutine realfft(data,n1,n2,n3,speq,isign)
c---  FFT code for a 3D real data den(n1,n2,n3) & work array speq(2,n2,n3)
c---  The dimensions n1, n2, n3 should be integer powers of 2, 3 & 5.

c---  [ Forward FFT ] isign = 1
c---  Input  : DATA contains real values at grid points (n1,n2,n3)
c---           SPEQ does not affects the result
c---  Output : DATA contains Fourier components at grid points
c---                (n1/2,n2,n3); i.e. DATA(1,1,1) and DATA(2,1,1)
c---                are the real and complex components of the 1st
c---                (kx=ky=kz=0) Fourier component and so on.
c---                This means we use only the half plane kx >= 0
c---                for FFTs of a real array to save time and memory.
c---           SPEQ contains the Nyquist frequency values of kx 
c---                frequency component; i.e. SPEQ(1,1,1) and SPEQ(2,1,1)
c---                are the real and complex components of the component
c---                with kx=n1/2, ky=kz=0 and so on.

c---  [ Inverse FFT ] isign = -1
c---  Input  : DATA contains Fourier components at grid points
c---                (n1/2,n2,n3); i.e. same as the output of forward FFT
c---           SPEQ contains the Nyquist frequency values of kx 
c---                frequency component; i.e. same as the output of the
c---                forward FFT
c---  Output : DATA contains real values at grid points (n1,n2,n3)
c---           SPEQ can be ignored
c---  Ref:Press,W.H. and Teukolsky, S.A. 1989, Computers in Physics.

c---  Modification (Jan. 30, 1991) : To endow more flexability and
c---  efficiency, cfft99 routine is used instead of the inferior FOURN.
c---  Cfft99 is originally written by Clive Temperton at ecmwf in 1978.
c---  This routine is modified so that it doesn't use twice more
c---  memory than the array size, which is crucial to n-body work.

      double precision wr,wi,wpr,wtemp,theta
      real*4 data(n1,n2,n3),speq(2,n2,n3)
      ip1= 1
      im1= 1
      c1 = 0.5
      c2 = -0.5*isign
      theta = 6.28318530717959d0/dble(isign*n1)
      wpr = -2.d0*dsin(0.5d0*theta)**2
      wpi = dsin(theta)
c---  Forward Transform
      if(isign.eq.1)then
      call cxfft3(data,n1/2,n2,n3,isign,speq)
      do 12 i3=1,n3
	 do 11 i2=1,n2
	 speq(1,i2,i3) = data(1,i2,i3)
	 speq(2,i2,i3) = data(2,i2,i3)
   11    continue
   12 continue
      ip1 = -1
      else
      im1 = -1
c---  Do-loop 10 can be omitted if you do this in your main program
c---  input is divided by (n1/2)*n2*n3
      factor = 2./float(n1*n2*n3)
      do 10 k=1,n3
      do 10 j=1,n2
      speq(1,j,k) = speq(1,j,k)*factor
      speq(2,j,k) = speq(2,j,k)*factor
      do 10 i=1,n1
   10 data(i,j,k) = data(i,j,k)*factor
      end if

      do 16 i3=1,n3
      j3 = 1
      if(i3.ne.1)j3=n3-i3+2
      wr=1.d0
      wi=0.d0
	 do 13 i2=1,n2
         j2 = 1
	 if(i2.ne.1)j2 = n2-i2+2
	   d1ps1 = data(1,i2,i3)+speq(1,j2,j3)
	   d1ms1 = data(1,i2,i3)-speq(1,j2,j3)
	   d2ps2 =(data(2,i2,i3)+speq(2,j2,j3))*im1
	   d2ms2 =(data(2,i2,i3)-speq(2,j2,j3))*im1
	   data(1,i2,i3) =  c1*d1ps1-c2*d2ps2
	   data(2,i2,i3) = (c1*d2ms2+c2*d1ms1)*ip1
	   speq(1,j2,j3) =  c1*d1ps1+c2*d2ps2
	   speq(2,j2,j3) =(-c1*d2ms2+c2*d1ms1)*ip1
   13    continue
      wtemp = wr
      wr = wr*wpr-wi*wpi+wr
      wi = wi*wpr+wtemp*wpi+wi

      do 15 i1=2,n1/4+1
      j1 = n1/2-i1+2
      i1t = 2*i1
      j1t = 2*j1
	   d1ps1 = data(i1t-1,1,i3)+data(j1t-1,1,j3)
	   d1ms1 = data(i1t-1,1,i3)-data(j1t-1,1,j3)
	   d2ps2 =(data(i1t,1,i3)+data(j1t,1,j3))*im1
	   d2ms2 =(data(i1t,1,i3)-data(j1t,1,j3))*im1
	   data(i1t-1,1,i3) =  c1*d1ps1-wr*c2*d2ps2-wi*c2*d1ms1
	   data(i1t,1,i3)   = (c1*d2ms2-wi*c2*d2ps2+wr*c2*d1ms1)*ip1
	   data(j1t-1,1,j3) =  c1*d1ps1+wr*c2*d2ps2+wi*c2*d1ms1
	   data(j1t,1,j3)   =(-c1*d2ms2-wi*c2*d2ps2+wr*c2*d1ms1)*ip1
	 do 14 i2=2,n2
	 j2 = n2-i2+2
	   d1ps1 = data(i1t-1,i2,i3)+data(j1t-1,j2,j3)
	   d1ms1 = data(i1t-1,i2,i3)-data(j1t-1,j2,j3)
	   d2ps2 =(data(i1t,i2,i3)+data(j1t,j2,j3))*im1
	   d2ms2 =(data(i1t,i2,i3)-data(j1t,j2,j3))*im1
	   data(i1t-1,i2,i3) =  c1*d1ps1-wr*c2*d2ps2-wi*c2*d1ms1
	   data(i1t,i2,i3)   = (c1*d2ms2-wi*c2*d2ps2+wr*c2*d1ms1)*ip1
	   data(j1t-1,j2,j3) =  c1*d1ps1+wr*c2*d2ps2+wi*c2*d1ms1
	   data(j1t,j2,j3)   =(-c1*d2ms2-wi*c2*d2ps2+wr*c2*d1ms1)*ip1
   14    continue
      wtemp = wr
      wr = wr*wpr-wi*wpi+wr
      wi = wi*wpr+wtemp*wpi+wi
   15 continue
   16 continue

c---  Inverse Transform
   20 continue
      if(isign.eq.-1)then
      call cxfft3(data,n1/2,n2,n3,isign,speq)
      end if
      return
      end


      subroutine cxfft3(cx,mx,my,mz,isign,cw)
c---  Modified version of cxfft3 : the size of work array is 2*my*mz
      parameter (mxtrig=2400)
      integer ifax(13),ifay(13),ifaz(13)
      complex cx(mx,my,mz), cw(1,my,mz)
      character cmsgnm*70
      real trig(mxtrig)
      external cfft99, cftfax
      save icall
      data icall/1/
      if (icall .eq. 1) then
          ntrig = 2* (mx+my+mz)
          if (ntrig .gt. mxtrig) then
              write (0,fmt='(2i5,a)') ntrig, mxtrig,
     $          ' 2*nx+2*ny+2*nz .gt. mxtrig. increase mxtrig '
              stop
          end if
          ltrigx = 1
          ltrigy = 1 + 2*mx
          ltrigz = 1 + 2*mx + 2*my
          call cftfax(mx,ifax,trig(ltrigx))
          call cftfax(my,ifay,trig(ltrigy))
          call cftfax(mz,ifaz,trig(ltrigz))
          icall = 2
      end if

c --- 1st ,2nd & 3rd dimension by planes
      do 10 k=1,mz
   10 call cfft99(cx(1,1,k),cw,trig(ltrigx),ifax,1,mx,mx,my,isign)
      do 20 k=1,mz
   20 call cfft99(cx(1,1,k),cw,trig(ltrigy),ifay,mx,1,my,mx,isign)
      do 30 j=1,my
   30 call cfft99(cx(1,j,1),cw,trig(ltrigz),ifaz,mx*my,1,mz,mx,isign)
      return
      end

      subroutine cfft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC
C              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.
C              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE
C              PACKAGE RETURNS A SET OF COMPLEX FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS
C              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.
C
C              THE PACKAGE CFFT99 CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE CFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO CFFT99
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINE CFFT99
C                THE ACTUAL TRANSFORM ROUTINE ROUTINE, CABABLE OF
C                PERFORMING BOTH THE TRANSFORM AND ITS INVERSE.
C                HOWEVER, AS THE TRANSFORMS ARE NOT NORMALIZED,
C                THE APPLICATION OF A TRANSFORM FOLLOWED BY ITS
C                INVERSE WILL YIELD THE ORIGINAL VALUES MULTIPLIED
C                BY N.
C
C
C ACCESS       *FORTRAN,P=XLIB,SN=CFFT99
C
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF
C              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER
C              COEFFICIENT VECTORS OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(2*N)
C                   COMPLEX A(...), WORK(...)
C
C                   CALL CFTFAX (N, IFAX, TRIGS)
C                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND
C              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR
C              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED
C              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.
C              SEE THE INDIVIDUAL WRITE-UPS FOR CFTFAX AND
C              CFFT99 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.  IT WAS
C              FURTHER MODIFIED FOR THE FULLY COMPLEX CASE BY DAVE
C              FULKER IN NOVEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,
C              PROVIDED THAT N IS NOT CHANGED.
C
C ARGUMENT     IFAX(13),TRIGS(2*N)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF
C               THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.
C
C              TRIGS
C               A REAL ARRAY OF DIMENSION 2*N
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS
C               GREATER THAN 5, IFAX(1) IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE CFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS (UNNORMALIZED) COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS.  GIVEN A SET OF COMPLEX GRIDPOINT
C              VECTORS, THE PACKAGE RETURNS A SET OF
C              COMPLEX FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE A
C              NUMBER HAVING NO PRIME FACTORS OTHER THAN
C              2, 3, AND 5.  THIS ROUTINE IS
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     COMPLEX A(N*INC+(LOT-1)*JUMP), WORK(N*LOT)
C DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13)
C
C ARGUMENTS
C
C ON INPUT     A
C               A COMPLEX ARRAY OF LENGTH N*INC+(LOT-1)*JUMP CONTAINING
C               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY
C               OVERWRITTEN BY THE RESULTS.
C
C               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF
C               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS
C               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT
C               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER"S
C               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS
C               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE
C               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A
C               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF
C               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE
C               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND
C               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND
C               IMAGINARY PORTIONS EXPLICITLY.
C
C              WORK
C               A COMPLEX WORK ARRAY OF LENGTH N*LOT OR A REAL ARRAY
C               OF LENGTH 2*N*LOT.  SEE N.B. ABOVE.
C
C              TRIGS
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C
C               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED
C               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED
C               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.
C
C              INC
C               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENT
C               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR
C               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS
C               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              LOT
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C
C ON OUTPUT    A
C               IF ISIGN = -1, AND LOT GRIDPOINT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N COEFFICIENT VALUES:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 C(K) = SUM(J=0,...,N-1)( G(J)*EXP(-2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               IF ISIGN = +1, AND LOT COEFFICIENT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N GRIDPOINT VALUES:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 G(J) = SUM(K=0,...,N-1)( G(K)*EXP(+2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               A CALL WITH ISIGN=-1 FOLLOWED BY A CALL WITH ISIGN=+1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED
C               BY THE FACTOR N.
C
C
C EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN
C               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL
C               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,
C               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE
C               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE
C               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE
C               ORTHOGONAL DIRECTION.
C
C               COMPLEX A(66,9), WORK(64,9)
C               REAL TRIGS1(128), TRIGS2(18)
C               INTEGER IFAX1(13), IFAX2(13)
C
C               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:
C
C               CALL CFTFAX(64, IFAX1, TRIGS1)
C               CALL CFTFAX( 9, IFAX2, TRIGS2)
C
C               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE
C               STORED IN MEMORY AS FOLLOWS (USING U AND V TO
C               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND
C               ASSUMING CONVENTIONAL FORTRAN STORAGE):
C
C   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS,
C
C   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS,
C
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C
C   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS.
C
C               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE
C               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 9, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST
C               DIMENSION OF A).
C
C               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)
C
C               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT
C               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFIN
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 64, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 1
C
C               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)
C
C               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONA
C               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT
C               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN
C               WITH ISIGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF
C               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).
C
C
C-----------------------------------------------------------------------

      dimension a(*),work(*),trigs(*),ifax(*)

      nn = n+n
      ink=inc+inc
      jum = jump+jump
      nfax=ifax(1)
      jnk = 2
      jst = 2
      if (isign.ge.0) go to 30
      jnk = -2
      jst = nn-2
      if (mod(nfax,2).eq.1) goto 40
      ibase = 1
      ilast = (n-1)*ink
      nh = n/2
      do 20 l=1,lot
      i1 = ibase+ink
      i2 = ibase+ilast

      do 10 m=1,nh
      hreal = a(i1)
      himag = a(i1+1)
      a(i1) = a(i2)
      a(i1+1) = a(i2+1)
      a(i2) = hreal
      a(i2+1) = himag
      i1 = i1+ink
      i2 = i2-ink
   10 continue
      ibase = ibase+jum
   20 continue
      goto 100
   30 continue
      if (mod(nfax,2).eq.0) goto 100
   40 continue
      ibase=1
      jbase=1
      do 60 l=1,lot
      work(jbase) = a(ibase)
      work(jbase+1) = a(ibase+1)
      i=ibase+ink
      j=jbase+jst

      do 50 m=2,n
      work(j) = a(i)
      work(j+1) = a(i+1)
      i=i+ink
      j=j+jnk
   50 continue
      ibase=ibase+jum
      jbase=jbase+nn
   60 continue
  100 continue

      igo = 110
      if (mod(nfax,2).eq.1) igo = 120
      la=1
      do 140 k=1,nfax
      if (igo.eq.120) go to 120
  110 continue
      call vpassm(a(1),a(2),work(1),work(2),trigs,
     *   ink,2,jum,nn,lot,n,ifax(k+1),la)
      igo=120
      go to 130
  120 continue
      call vpassm(work(1),work(2),a(1),a(2),trigs,
     *    2,ink,nn,jum,lot,n,ifax(k+1),la)
      igo=110
  130 continue
      la=la*ifax(k+1)
  140 continue
      return
      end


      subroutine cftfax(n,ifax,trigs)
      dimension ifax(13),trigs(*)

      call fact(n,ifax)
      k = ifax(1)
      if (k .lt. 1 .or. ifax(k+1) .gt. 5) ifax(1) = -99
      if (ifax(1) .le. 0)  then
	  write(0,*) 'FFT: error in array dimensions'
	  stop
      end if
cc      if (ifax(1) .le. 0 )  call errmsg('fatal','cfft99f',
cc     $ '  fftfax - invalid n')
      call cftrig (n, trigs)
      return
      end
      subroutine fact(n,ifax)
      dimension ifax(13)
      if (n.gt.1) go to 10
      ifax(1) = 0
      if (n.lt.1) ifax(1) = -99
      return
   10 nn=n
      k=1
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
c     now find remaining factors
   50 l=5
      max = sqrt(float(nn))
      inc=2

   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 if (l.gt.max) go to 75
      l=l+inc
      inc=6-inc
      go to 60
   75 k = k+1
      ifax(k) = nn
   80 ifax(1)=k-1
      return
      end


      subroutine cftrig(n,trigs)
      dimension trigs(*)
      pi=2.0*asin(1.0)
      del=(pi+pi)/float(n)
      l=n+n
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      return
      end


      subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
      dimension a(n),b(n),c(n),d(n),trigs(n)
      data sin36/0.587785252292473/,cos36/0.809016994374947/,
     *     sin72/0.951056516295154/,cos72/0.309016994374947/,
     *     sin60/0.866025403784437/

      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo

c     coding for factor 2
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase

      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase

      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return

c     coding for factor 3
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase

      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase

      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return

c     coding for factor 4
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase

      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase

      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return

c     coding for factor 5
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase

      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase

      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
      end

