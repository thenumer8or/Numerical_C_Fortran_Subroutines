C***********************************************************************
C                                                                      *
C                      BROYDEN'S METHOD                                *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION OF THE NONLINEAR SYSTEM F(X)=0
C     GIVEN AN INITIAL APPROXIMATION X.
C
C     INPUT:   NUMBER N OF EQUATIONS AND UNKNOWNS; INITIAL
C              APPROXIMATION X=(X(1),...,X(N)); TOLERANCE TOL;
C              MAXIMUM NUMBER OF ITERATIONS N.
C
C     OUTPUT:  APPROXIMATE SOLUTION X=(X(1),...,X(N)) OR A MESSAGE
C              THAT THE NUMBER OF ITERATIONS WAS EXCEEDED.
C
      DIMENSION A(3,3),C(3,3),X(3),Y(3),S(3),Z(3),V(3)
C     USE NN INSTEAD OF N FOR THE MAXIMUM NUMBER OF ITERATIONS
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Broyden Method for Nonlinear Systems.'
      WRITE(6,*) 'Has the function F been defined and have the'
      WRITE(6,*) 'partial derivatives been defined as follows:'
      WRITE(6,*) '1.  function F(N,I,X) where I is the number of the'
      WRITE(6,*) '    component function, X is the vector at which'
      WRITE(6,*) '    to evaluate and N is the dimension.'
      WRITE(6,*) '2.  function FP(N,I,J,X) where I is the number of'
      WRITE(6,*) '    the component function and J is the number'
      WRITE(6,*) '    of the variable with respect to which'
      WRITE(6,*) '    partial differentiation is performed,'
      WRITE(6,*) '    X is the vector at which to evaluate and N'
      WRITE(6,*) '    is the dimension.'
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AAA
      IF(( AAA .EQ. 'Y' ) .OR. ( AAA .EQ. 'y' )) THEN
         OK = .FALSE.
101      IF (OK) GOTO 11
            WRITE(6,*) 'Input the number n of equations.'
            WRITE(6,*) ' '
            READ(5,*) N
            IF (N.GE.2) THEN
               OK = .TRUE.
            ELSE
               WRITE(6,*) 'N must be greater than 1.'
            ENDIF
         GOTO 101
11       OK = .FALSE.
14       IF (OK) GOTO 15
            WRITE(6,*) 'Input tolerance '
            WRITE(6,*) ' '
            READ(5,*) TOL
            IF (TOL.LE.0.0) THEN
               WRITE(6,*) 'Tolerance must be positive '
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 14
15       OK=.FALSE.
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input the maximum number of iterations.'
            WRITE(6,*) ' '
            READ(5,*) NN
            IF ( NN .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       DO 16 I=1,N
            WRITE(6,*) 'Input initial approximation X(',I,')'
            WRITE(6,*) ' '
            READ(5,*) X(I)
16       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'F(N,I,X) and FP(N,I,J,X) can be created '
         OK = .FALSE.
      ENDIF
      IF(.NOT.OK) GOTO 400
      WRITE(6,*) 'Select output destinations: '
      WRITE(6,*) '1. Screen '
      WRITE(6,*) '2. Text file '
      WRITE(6,*) 'Enter 1 or 2 '
      WRITE(6,*) ' '
      READ(5,*) FLAG
      IF ( FLAG .EQ. 2 ) THEN
         WRITE(6,*) 'Input the file name in the form - '
         WRITE(6,*) 'drive:name.ext'
         WRITE(6,*) 'with the name contained within quotes'
         WRITE(6,*) 'as example:   ''A:OUTPUT.DTA'' '
         WRITE(6,*) ' '
         READ(5,*) NAME1
         OUP = 3
         OPEN(UNIT=OUP,FILE=NAME1,STATUS='NEW')
      ELSE
         OUP = 6
      ENDIF
      WRITE(OUP,*) 'BROYDEN METHOD FOR NONLINEAR SYSTEMS'
      WRITE(OUP,*) 'Iteration, approximation x(1),...,x(n) and error'
      PI = 4*ATAN(1.)
C     SN WILL REPRESENT THE L2 NORM OF THE ERROR
      SN=0
      K=0
      WRITE(OUP,3) K,(X(I),I=1,N),SN
C     STEP 1
C     A WILL HOLD THE JACOBIAN FOR THE INITIAL APPROXIMATION
C     THE SUBPROGRAM FP COMPUTES THE ENTRIES OF THE JACOBIAN
      DO 10 I=1,N
           DO 10 J=1,N
10         A(I,J)=FP(N,I,J,X)
C     COMPUTE V = F(X(0))
C     THE SUBPROGRAM F(N,I,X) COMPUTES THE ITH COMPONENT OF F AT X
      DO 20 I=1,N
20         V(I)=F(N,I,X)
C     STEP 2
      CALL INVER(N,A,C)
C     INVER FINDS THE INVERSE OF THE N BY N MATRIX A AND RETURNS IT IN A
C     C IS A DUMMY PARAMETER USED TO RESERVE STORAGE IN THE SUBROUTINE
C     STEP 3
      K=1
C     NOTE: S=S(1)
      CALL MULT(N,A,V,S,SN)
C     MULT COMPUTES THE PRODUCT S=-A*V AND THE L2-NORM SN OF S
      DO 30 I=1,N
30         X(I)=X(I)+S(I)
C     NOTE: X=X(1)
      WRITE(OUP,3) K,(X(I),I=1,N),SN
C     STEP 4
200   IF ( K.GT.NN ) GOTO 300
C          STEP 5
C          THE VECTOR W IS NOT USED SINCE THE FUNCTION F IS EVALUATED
C          COMPONENT BY COMPONENT
           DO 40 I=1,N
                VV=F(N,I,X)
                Y(I)=VV-V(I)
40         V(I)=VV
C          NOTE: V=F(X(K)) AND Y=Y(K)
C          STEP 6
           CALL MULT(N,A,Y,Z,ZN)
C          NOTE: Z=-A(K-1)**-1 * Y(K)
C          STEP 7
           P=0
C          P WILL BE S(K)**T * A(K-1)**-1 * Y(K)
           DO 60 I=1,N
                P=P-S(I)*Z(I)
C               STEP 8
                DO 60 J=1,N
60                   C(I,J)=(S(I)+Z(I))*S(J)
           DO 100 I=1,N
100             C(I,I)=C(I,I)+P
C          C=S(K)**T * A(K-1)**-1 * Y(K)*I + (S(K)+A(K-1)**-1 * Y(K))*
C                S(K)**T )
C          STEPS 9 AND 10
           DO 70 I=1,N
                DO 80 J=1,N
                     ACC=0
                     DO 90 L=1,N
C                         ACC IS THE (I,J) ENTRY OF THE INVERSE OF A
90                        ACC=ACC+C(I,L)*A(L,J)/P
C                    Z WILL HOLD THE J-TH COLUMN OF THE INVERSE OF A
80                   Z(J)=ACC
                DO 70 J=1,N
70                   A(I,J)=Z(J)
           CALL MULT(N,A,V,S,SN)
C          NOTE: A=A(K)**-1 AND S=-A(K)**-1 * F(X(K))
C          STEP 11
           DO 110 I=1,N
110             X(I)=X(I)+S(I)
C          NOTE: X = X(K+1)
           WRITE(OUP,3) K+1,(X(I),I=1,N),SN
C          STEP 12
           IF (SN.LE.TOL) THEN
C               PROCEDURE COMPLETED SUCCESSFULLY
                WRITE(OUP,6)
                GOTO 400
           END IF
C          STEP 13
           K=K+1
      GOTO 200
C     STEP 14
300   WRITE(OUP,4)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,I2,3(1X,E15.8,1X),E15.8)
4     FORMAT(1X,'MAXIMUM NUMBER OF ITERATIONS EXCEEDED')
6     FORMAT(1X,'PROCEDURE COMPLETED SUCCESSFULLY')
      END
C
      FUNCTION FP(N,I,J,X)
      DIMENSION X(N)
      IF (I.EQ.1) THEN
                IF(J.EQ.1) FP=3
                IF(J.EQ.2) FP=X(3)*SIN(X(2)*X(3))
                IF(J.EQ.3) FP=X(2)*SIN(X(2)*X(3))
      ELSE
         IF (I.EQ.2) THEN
                IF(J.EQ.1) FP=2*X(1)
                IF(J.EQ.2) FP=-162*(X(2)+.1)
                IF(J.EQ.3) FP= COS(X(3))
         ELSE
                IF(J.EQ.1) FP=-X(2)*EXP(-X(1)*X(2))
                IF(J.EQ.2) FP= -X(1)*EXP(-X(1)*X(2))
                IF(J.EQ.3) FP=20
         ENDIF
      ENDIF
      RETURN
      END
C
      FUNCTION F(N,I,X)
      DIMENSION X(N)
      IF(I.EQ.1) F=3*X(1)-COS(X(2)*X(3))-.5
      IF(I.EQ.2) F=X(1)**2-81*(X(2)+.1)**2+SIN(X(3))+1.06
      IF(I.EQ.3) F=EXP(-X(1)*X(2))+20*X(3)+(10*3.141593-3)/3
      RETURN
      END
C
      SUBROUTINE MULT(N,A,Y,Z,ZN)
      DIMENSION A(N,N),Y(N),Z(N)
      ZN=0
      DO 10 I=1,N
           Z(I)=0
           DO 15 J=1,N
15         Z(I)=Z(I)-A(I,J)*Y(J)
10    ZN=ZN+Z(I)*Z(I)
      ZN=SQRT(ZN)
      RETURN
      END
C
      SUBROUTINE INVER(N,A,B)
      DIMENSION A(N,N),B(N,N)
      DO 10 I=1,N
           DO 10 J=1,N
           B(I,J)=0
10    IF(J.EQ.I) B(I,J)=1
      DO 20 I=1,N
           I1=I+1
           I2=I
           IF(I.NE.N) THEN
                DO 30 J=I1,N
30              IF(ABS(A(J,I)).GT.ABS(A(I2,I))) I2=J
                IF(I2.NE.I) THEN
                        DO 40 J=1,N
                             C=A(I,J)
                             A(I,J)=A(I2,J)
                             A(I2,J)=C
                             C=B(I,J)
                             B(I,J)=B(I2,J)
                             B(I2,J)=C
40                      CONTINUE
                END IF
           END IF
           DO 50 J=1,N
                IF (J.NE.I) THEN
                     C=A(J,I)/A(I,I)
                     DO 60 K=1,N
                          A(J,K)=A(J,K)-C*A(I,K)
                          B(J,K)=B(J,K)-C*B(I,K)
60                   CONTINUE
                END IF
50         CONTINUE
20    CONTINUE
      DO 70 I=1,N
           C=A(I,I)
           DO 70 J=1,N
70    A(I,J)=B(I,J)/C
      RETURN
      END
