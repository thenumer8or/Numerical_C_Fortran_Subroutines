      FUNCTION ROFUNC(B)
      PARAMETER (NMAX=1000)
      COMMON /ARRAYS/ NDATA,X(NMAX),Y(NMAX),ARR(NMAX),AA,ABDEV
      N1=NDATA+1
      NML=N1/2
      NMH=N1-NML
      DO 11 J=1,NDATA
        ARR(J)=Y(J)-B*X(J)
11    CONTINUE
      CALL SORT(NDATA,ARR)
      AA=0.5*(ARR(NML)+ARR(NMH))
      SUM=0.
      ABDEV=0.
      DO 12 J=1,NDATA
        D=Y(J)-(B*X(J)+AA)
        ABDEV=ABDEV+ABS(D)
        SUM=SUM+X(J)*SIGN(1.0,D)
12    CONTINUE
      ROFUNC=SUM
      RETURN
      END
