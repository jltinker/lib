      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE)
      PARAMETER (EPS=1.E-5,JMAX=16,JMAXP=JMAX+1,KM=4,K=KM+1)
      DIMENSION S(JMAXP),H(JMAXP)
      EXTERNAL FUNC,CHOOSE
      H(1)=1.
      DO 11 J=1,JMAX
        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
c           print*,J,H(J-KM),S(J-KM),K,SS,DSS,'//'
          CALL POLINT(H(J-KM),S(J-KM),K,0.0,SS,DSS)
c          print*,SS,DSS,'//'
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      PAUSE 'Too many steps.'
      END

      SUBROUTINE MIDPNT(FUNC,A,B,S,N)
      EXTERNAL FUNC
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END

      SUBROUTINE MIDINF(FUNK,AA,BB,S,N)
      EXTERNAL FUNK
      FUNC(X)=FUNK(1./X)/X**2
      B=1./AA
      A=1./BB
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

