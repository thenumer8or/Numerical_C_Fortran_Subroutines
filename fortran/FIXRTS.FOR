      SUBROUTINE FIXRTS(D,NPOLES)
      PARAMETER (NPMAX=100)
      DIMENSION D(NPOLES)
      LOGICAL POLISH
      COMPLEX A(NPMAX),ROOTS(NPMAX)
      A(NPOLES+1)=CMPLX(1.,0.)
      DO 11 J=NPOLES,1,-1
        A(J)=CMPLX(-D(NPOLES+1-J),0.)
11    CONTINUE
      POLISH=.TRUE.
      CALL ZROOTS(A,NPOLES,ROOTS,POLISH)
      DO 12 J=1,NPOLES
        IF(CABS(ROOTS(J)).GT.1.)THEN
          ROOTS(J)=1./CONJG(ROOTS(J))
        ENDIF
12    CONTINUE
      A(1)=-ROOTS(1)
      A(2)=CMPLX(1.,0.)
      DO 14 J=2,NPOLES
        A(J+1)=CMPLX(1.,0.)
        DO 13 I=J,2,-1
          A(I)=A(I-1)-ROOTS(J)*A(I)
13      CONTINUE
        A(1)=-ROOTS(J)*A(1)
14    CONTINUE
      DO 15 J=1,NPOLES
        D(NPOLES+1-J)=-REAL(A(J))
15    CONTINUE
      RETURN
      END
