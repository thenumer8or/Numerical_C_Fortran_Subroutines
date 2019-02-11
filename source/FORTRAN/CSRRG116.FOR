C******************************************************************************
C                                                                             *
C                CUBIC SPLINE RAYLEIGH-RITZ METHOD                            *
C                                                                             *
C******************************************************************************
C
C
C     TO APPROXIMATE THE SOLUTION TO THE BOUNDARY-VALUE PROBLEM
C
C         -D(P(X)Y')/DX + Q(X)Y = F(X), 0<=X<=1, Y(0)=Y(1)=0
C
C     WITH A SUM OF CUBIC SPLINES:
C
C     INPUT   INTEGER N
C
C     OUTPUT  COEFFICIENTS C(0),...,C(N+1) OF THE BASIS FUNCTIONS
C
C     (NOTE THAT C(I) USED HERE IS 6 TIMES AS LARGE AS THE C(I)
C     USED IN THE ALGORITHM, SINCE PHI(X) IS DEFINED TO BE 1/6
C     OF THE BASIS FUNCTION USED IN THE TEXT.)
C
C     TO CHANGE PROBLEMS DO THE FOLLOWING:
C         1. CHANGE FUNCTIONS P, Q AND F IN THE SUBPROGRAMS NAMED P,Q,F
C         2. CHANGE N
C         3. CHANGE DIMENSIONS IN MAIN PROGRAM TO VALUES GIVEN BY
C               DIMENSION A(N+2,N+3),X(N+2),C(N+2),CO(N+2,4,4),
C                         DCO(N+2,4,3),AF(N+1),BF(N+1),CF(N+1),DF(N+1),
C                         AP(N+1),BP(N+1),CP(N+1),DP(N+1),
C                         AQ(N+1),BQ(N+1),CQ(N+1),DQ(N+1)
C         4. CHANGE DIMENSIONS IN SUBROUTINE COEF TO
C               DIMENSION AA(N+2),BB(N+2),CC(N+2),DD(N+2),
C                         H(N+2),XA(N+2),XL(N+2),XU(N+2),XZ(N+2)
C
C
C     GENERAL OUTLINE
C
C         1. NODES LABELLED X(I)=(I-1)*H, 1 <= I <= N+2, WHERE
C            H=1/(N+1) SO THAT ZERO SUBSCRIPTS ARE AVOIDED
C         2. THE FUNCTIONS PHI(I) AND PHI'(I) ARE SHIFTED SO THAT
C            PHI(1) AND PHI'(1) ARE CENTERED AT X(1), PHI(2) AND PHI'(2)
C            ARE CENTERED AT X(2), . . . , PHI(N+2) AND
C            PHI'(N+2) ARE CENTERED AT X(N+2)---FOR EXAMPLE,
C                     PHI(3) = S((X-X(3))/H)
C                            = S(X/H + 2)
C        3. THE FUNCTIONS PHI(I) ARE REPRESENTED IN TERMS OF THEIR
C            COEFFICIENTS IN THE FOLLOWING WAY:
C            (PHI(2))(X) = CO(I,K,1) + CO(I,K,2)*X + CO(I,K,3)*X**2
C                       CO(I,K,4)*X**3
C            FOR X(J) <= X <= X(J+1) WHERE
C            K=1 IF J=I-2, K=2 IF J=I-1, K=3 IF J=I, K=4 IF J=I+1
C            SINCE PHI(I) IS NONZERO ONLY BETWEEN X(I-2) AND X(I+2)
C            UNLESS I = 1, 2, N+1 OR N+2
C            (SEE SUBROUTINE PHICO)
C         4. THE DERIVATIVE OF PHI(I) DENOTED PHI'(I) IS REPRESENTED
C            AS IN 3. BY ITS COEFFICIENTS DCO(I,K,L), L = 1, 2, 3
C            (SEE SUBROUTINE DPHICO).
C         5. THE FUNCTIONS P,Q AND F ARE REPRESENTED BY THEIR CUBIC
C            SPLINE INTERPOLANTS USING CLAMPLED BOUNDARY CONDITIONS
C            (SEE ALGORITHM 3.3).  THUS, FOR X(I) <= X <= X(I+1) WE
C            USE AF(I) + BF(I)*X + CF(I)*X**2 + DF(I)*X**3 TO
C            REPRESENT F(X).  SIMILARLY, AP,BP,CP,DP ARE USED FOR P AND
C            AQ,BQ,CQ,DQ ARE USED FOR Q.  (SEE SUBROUTINE COEF).
C         6. THE INTEGRANDS IS STEPS 6 AND 9 ARE REPLACED BY PRODUCTS
C            OF CUBIC POLYNOMIAL APPROXIMATIONS ON EACH SUBINTERVAL OF
C            LENGTH H AND THE INTEGRALS OF THE RESULTING POLYNOMIALS
C            ARE COMPUTED EXACTLY.  (SEE SUBROUTINE XINT).
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL P,Q,F
      DIMENSION A(11,12),X(11),C(11),CO(11,4,4),DCO(11,4,3),AF(10)
      DIMENSION DF(10),AP(10),BP(10),CP(10),DP(10),AQ(10),BQ(10),DQ(10)
      DIMENSION BF(10),CF(10),CQ(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      COMMON N
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Cubic Spline Rayleigh-Ritz Method.'
      WRITE(6,*) 'Have the functions P, Q and F been created? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input a positive integer n where'
            WRITE(6,*) 'x(0) = 0, ..., x(n+1) = 1.'
            WRITE(6,*) ' '
            READ(5,*) N
            IF ( N .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       WRITE(6,*) 'Input the derivative of F at 0 and at 1.'
         WRITE(6,*) ' '
         READ(5,*) FPL, FPR
         WRITE(6,*) 'Input the derivative of P at 0 and at 1.'
         WRITE(6,*) ' '
         READ(5,*) PPL, PPR
         WRITE(6,*) 'Input the derivative of Q at 0 and at 1.'
         WRITE(6,*) ' '
         READ(5,*) QPL, QPR
      ELSE
         WRITE(6,*) 'The program will end so that P, Q and F '
         WRITE(6,*) 'can be created '
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
      WRITE(OUP,*) 'CUBIC SPLINE RAYLEIGH-RITZ METHOD'
C     STEP 1
      H=1.0/(N+1)
      N1=N+1
      N2=N+2
      N3=N+3
C     INITIALIZE MATRIX A AT ZERO, NOTE THAT A(I,N+3) = B(I)
      DO 20 I=1,N2
           DO 20 J=1,N3
20              A(I,J)=0
C
C     STEP 2
C
C     X(1)=0,...,X(I)=(I-1)*H,...,X(N+2)=1
      DO 30 I=1,N2
30         X(I)=(I-1)*H
C
C     STEPS 3 AND 4 ARE IMPLEMENTED IN WHAT FOLLOWS
C
C     INITIALIZE COEFFICIENTS CO(I,J,K) AND DCO(I,J,K)
      DO 40 I=1,N2
           DO 50 J=1,4
                DO 60 K=1,4
                     CO(I,J,K)=0
                     IF(K.NE.4) THEN
                          DCO(I,J,K)=0
                     END IF
60              CONTINUE
C     JJ CORRESPONDS THE COEFFICIENTS OF PHI AND PHI' TO THE PROPER
C     INTERVAL INVOLVING J
                JJ=I+J-3
                CALL PHICO(I,JJ,CO(I,J,1),CO(I,J,2),CO(I,J,3),CO(I,J,4))
                CALL DPHICO(I,JJ,DCO(I,J,1),DCO(I,J,2),DCO(I,J,3))
50         CONTINUE
40    CONTINUE
C     OUTPUT THE BASIS FUNCTIONS
      WRITE(OUP,1)
      WRITE(OUP,2)
      DO 70 I=1,N2
           WRITE(OUP,5) I,I
           DO 70 J=1,4
                IF(I.NE.1.OR.(J.NE.1.AND.J.NE.2)) THEN
                     IF(I.NE.2.OR.J.NE.2) THEN
                          IF(I.NE.N1.OR.J.NE.4) THEN
                              IF(I.NE.N2.OR.(J.NE.3.AND.J.NE.4)) THEN
                                   JJ1=I-3+J
                                   JJ2=I-2+J
       WRITE(OUP,3) JJ1,JJ2,(CO(I,J,K),K=1,4),(DCO(I,J,K),K=1,3)
                               END IF
                          END IF
                     END IF
                END IF
70    CONTINUE
C
C     OBTAIN COEFFICIENTS FOR F, P AND Q--NOTE THAT THE 2ND AND
C     3RD ARGUMENTS ARE THE DERIVATIVES AT 0 AND 1 RESP.
C
      CALL COEF(F,FPL,FPR,X,AF,BF,CF,DF,N2,N1)
      CALL COEF(P,PPL,PPR,X,AP,BP,CP,DP,N2,N1)
      CALL COEF(Q,QPL,QPR,X,AQ,BQ,CQ,DQ,N2,N1)
C
C     STEPS 5, 6, 7, 8, 9 ARE IMPLEMENTED IN WHAT FOLLOWS
C
      WRITE(OUP,6)
      DO 80 I=1,N2
C     INDICES OF LIMITS OF INTEGRATION FOR A(I,I) AND B(I)
           J1=MIN0(I+2,N+2)
           JO=MAX0(I-2,1)
           J2=J1-1
C     INTEGRATE OVER EACH SUBINTERVAL WHERE PHI(I) NONZERO
           DO 90 JJ=JO,J2
C     LIMITS OF INTEGRATION FOR EACH CALL
                XU=X(JJ+1)
                XL=X(JJ)
C     COEFFICIENTS OF BASES
                K=INTY(I,JJ)
                A1=DCO(I,K,1)
                B1=DCO(I,K,2)
                C1=DCO(I,K,3)
                D1=0
                A2=CO(I,K,1)
                B2=CO(I,K,2)
                C2=CO(I,K,3)
                D2=CO(I,K,4)
C     CALL SUBPROGRAM FOR INTEGRATION
      A(I,I)=A(I,I)+XINT(XU,XL,AP(JJ),BP(JJ),CP(JJ),DP(JJ),A1,B1,C1,D1,
     1A1,B1,C1,D1)+XINT(XU,XL,AQ(JJ),BQ(JJ),CQ(JJ),DQ(JJ),A2,B2,C2,D2,
     1A2,B2,C2,D2)
90    A(I,N+3)=A(I,N+3)+XINT(XU,XL,AF(JJ),BF(JJ),CF(JJ),DF(JJ),A2,B2,C2,
     1D2,1.0D+00,0.0D+00,0.0D+00,0.0D+00)
C     COMPUTE A(I,J) FOR J = I+1, . . . , MIN(I+3,N+2)
      K3=I+1
      IF(K3.LE.N2) THEN
           K2=MIN0(I+3,N+2)
           DO 100 J=K3,K2
                JO=MAX0(J-2,1)
                DO 110 JJ=JO,J2
                     XU=X(JJ+1)
                     XL=X(JJ)
                     K=INTY(I,JJ)
                     A1=DCO(I,K,1)
                     B1=DCO(I,K,2)
                     C1=DCO(I,K,3)
                     D1=0
                     A2=CO(I,K,1)
                     B2=CO(I,K,2)
                     C2=CO(I,K,3)
                     D2=CO(I,K,4)
                     K=INTY(J,JJ)
                     A3=DCO(J,K,1)
                     B3=DCO(J,K,2)
                     C3=DCO(J,K,3)
                     D3=0
                     A4=CO(J,K,1)
                     B4=CO(J,K,2)
                     C4=CO(J,K,3)
                     D4=CO(J,K,4)
110   A(I,J)=A(I,J)+XINT(XU,XL,AP(JJ),BP(JJ),CP(JJ),DP(JJ),A1,B1,C1,D1,
     1A3,B3,C3,D3)+XINT(XU,XL,AQ(JJ),BQ(JJ),CQ(JJ),DQ(JJ),A2,B2,C2,D2,A4
     1,B4,C4,D4)
100   A(J,I)=A(I,J)
C     OUTPUT A(I,J) FOR J=I,I+1,...,MIN(I+3,N+2) AND J=N+3
      WRITE(OUP,7) (I,J,A(I,J),J=I,K2),I,N3,A(I,N3)
      END IF
80    CONTINUE
C
C     STEP 10
C
      DO 120 I=1,N1
           II=I+1
           DO 130 J=II,N2
                CC=A(J,I)/A(I,I)
                DO 140 K=II,N3
140             A(J,K)=A(J,K)-CC*A(I,K)
130        A(J,I)=0
120   CONTINUE
      C(N2)=A(N2,N3)/A(N2,N2)
      DO 150 I=1,N1
           J=N1-I+1
           C(J)=A(J,N3)
           JJ=J+1
           DO 160 KK=JJ,N2
160        C(J)=C(J)-A(J,KK)*C(KK)
150   C(J)=C(J)/A(J,J)
C
C     STEP 11
C
C     OUTPUT THE COEFFICIENTS C(I)
      WRITE(OUP,8)
      WRITE(OUP,9) (I,C(I),I=1,N2)
C     OBTAIN VALUES OF THE APPROX. AT THE NODES
      DO 170 I=1,N2
           S=0
           DO 180 J=1,N2
                JO=MAX0(J-2,1)
                J1=MIN0(J+2,N+2)
                SS=0
                IF(I.LT.JO.OR.I.GE.J1) THEN
                     S=S+C(J)*SS
                ELSE
                     K=INTY(J,I)
                     SS=((CO(J,K,4)*X(I)+CO(J,K,3))*X(I)+CO(J,K,2))*
     *                   X(I)+CO(J,K,1)
                     S=S+C(J)*SS
               END IF
180   CONTINUE
      WRITE(OUP,10) I,X(I),S
170   CONTINUE
C
C     STEP 12
C     PROCEDURE IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'BASIS FUNCTION: A + B*X + C*X**2 + D*X**3',21X,
     1 'DER. OF BASIS FUNCTION: A + B*X + C*X**2 ')
2     FORMAT(1X,'A',14X,'B',15X,'C',15X,'D',18X,'A',15X,'B',14X,
     1 'C')
3     FORMAT(1X,'ON  (X(',I2,'),X(',I2,'))',5X,4(D13.6,2X),3X,3(D13.6,
     1 2X))
5     FORMAT(1X,'PHI(',I2,')',77X,'PHI''(',I2,')')
6     FORMAT(1X,'NONZERO ENTRIES IN MATRIX A')
7     FORMAT(1X,5('A(',I2,',',I2,') = ',D15.8))
8     FORMAT(1X,'SOLUTION IS C(1)*PHI(1) + ... + C(N+2)*PHI(N+2)
     1WHERE ')
9     FORMAT(5(1X,'C(',I2,') = ', D15.8))
10    FORMAT(1X,'VALUE OF APPROX. AT X(',I2,') = ',D15.8,3X,'IS ',D15.8)
      END
C
      SUBROUTINE COEF(F,FPO,FPN,X,A,B,C,D,N,M)
C
C     THIS IMPLEMENTS ALGORITHM 3.3--CLAMPED BOUNDARY CUBIC SPLINE
C     F = FUNCTION TO BE APPROXIMATED, FPO = F'(0), FPN = F'(1),
C     X IS THE SET OF NODES, A,B,C,D ARE THE COEFFICIENTS TO BE
C     RETURNED, N = NUMBER OF NODES, M = NUMBER OF SUBINTERVALS
C     THE APPROX. WILL BE
C          A(I) + B(I)*X + C(I)*X**2 + D(I)*X**3 ON (X(I),X(I+1))
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AA(11),BB(11),CC(11),DD(11)
      DIMENSION X(N),A(M),B(M),C(M),D(M),H(11),XA(11),XL(11),XU(11)
      DIMENSION XZ(11)
      DO 10 I=1,M
           AA(I)=F(X(I))
10         H(I+1)=X(I+1)-X(I)
      AA(N)=F(X(N))
      XA(1)=3*(AA(2)-AA(1))/H(2)-3*FPO
      XA(N)=3*FPN-3*(AA(N)-AA(N-1))/H(N)
      XL(1)=2*H(2)
      XU(1)=.5D+00
      XZ(1)=XA(1)/XL(1)
      DO 20 I=2,M
           XA(I)=3*(AA(I+1)*H(I)-AA(I)*(X(I+1)-X(I-1))+AA(I-1)*H(I+1))
     *     /(H(I+1)*H(I))
           XL(I)=2*(X(I+1)-X(I-1))-H(I)*XU(I-1)
           XU(I)=H(I+1)/XL(I)
20         XZ(I)=(XA(I)-H(I)*XZ(I-1))/XL(I)
      XL(N)=H(N)*(2-XU(N-1))
      XZ(N)=(XA(N)-H(N)*XZ(N-1))/XL(N)
      CC(N)=XZ(N)
      DO 30 I=1,M
           J=N-I
           CC(J)=XZ(J)-XU(J)*CC(J+1)
           BB(J)=(AA(J+1)-AA(J))/H(J+1)-H(J+1)*(CC(J+1)+2*CC(J))/3
30         DD(J)=(CC(J+1)-CC(J))/(3*H(J+1))
      DO 40 I=1,M
           A(I)=(((-DD(I)*X(I))+CC(I))*X(I)-BB(I))*X(I)+AA(I)
           B(I)=(3*DD(I)*X(I)-2*CC(I))*X(I)+BB(I)
           C(I)=CC(I)-3*DD(I)*X(I)
40         D(I)=DD(I)
      RETURN
      END
C
      FUNCTION XINT (XU,XL,A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3)
C
C     FORMS THE PRODUCT (A1+B1*X+C1*X**2+D1*X**3)*(A2+B2*X+C2*X**2+
C     D2*X**3)*(A3+B3*X+C3*X**2+D3*X**3) TO OBTAIN
C     10*C(1)*X**9 + 9*C(2)*X**8 + . . . + 2*C(9)*X + C(10)
C     WHICH IS INTEGRATED FROM XL TO XU
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(10)
      AA=A1*A2
      BB=A1*B2+A2*B1
      CC=A1*C2+B1*B2+C1*A2
      DD=A1*D2+B1*C2+C1*B2+D1*A2
      EE=B1*D2+C1*C2+D1*B2
      FF=C1*D2+D1*C2
      GG=D1*D2
      C(10)=AA*A3
      C(9)=(AA*B3+BB*A3)/2
      C(8)=(AA*C3+BB*B3+CC*A3)/3
      C(7)=(AA*D3+BB*C3+CC*B3+DD*A3)/4
      C(6)=(BB*D3+CC*C3+DD*B3+EE*A3)/5
      C(5)=(CC*D3+DD*C3+EE*B3+FF*A3)/6
      C(4)=(DD*D3+EE*C3+FF*B3+GG*A3)/7
      C(3)=(EE*D3+FF*C3+GG*B3)/8
      C(2)=(FF*D3+GG*C3)/9
      C(1)=(GG*D3)/10
      XHIGH=0
      XLOW=0
      DO 10 I=1,10
           XHIGH=(XHIGH+C(I))*XU
10         XLOW=(XLOW+C(I))*XL
      XINT=XHIGH-XLOW
      RETURN
      END
C
      SUBROUTINE PHICO(I,J,A,B,C,D)
C
C     COMPUTES PHI(I) AS A+B*X+C*X**2+D*X**3 FOR X IN (X(J),X(J+1)) BY
C     MULTIPLYING AND SIMPLIFYING THE EQUATIONS IN STEP 2
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON N
      H=1.0D+00/(N+1)
      A=0
      B=0
      C=0
      D=0
      E=I-1
      IF(J.LT.I-2.OR.J.GE.I+2) RETURN
      IF(I.EQ.1.AND.J.LT.I) RETURN
      IF(I.EQ.2.AND.J.LT.I-1) RETURN
      IF(I.EQ.N+1.AND.J.GT.N+1) RETURN
      IF(I.EQ.N+2.AND.J.GE.N+2) RETURN
      IF(J.LE.I-2) THEN
           A= (((-E+6)*E-12)*E+8)/24
           B=((E-4)*E+4)/(8*H)
           C=(-E+2)/(8*H*H)
           D=1/(24*H**3)
           RETURN
      ELSE
           IF(J.GT.I-1.AND.J.GT.I) THEN
                A=(((E+6)*E+12)*E+8)/24
                B=((-E-4)*E-4)/(8*H)
                C=(E+2)/(8*H*H)
                D=-1/(24*H**3)
                RETURN
           ELSE
                IF(J.GT.I-1) THEN
                     A=((-3*E-6)*E*E+4)/24
                     B=(3*E+4)*E/(8*H)
                     C=(-3*E-2)/(8*H*H)
                     D=1/(8*H**3)
                     IF(I.NE.1.AND.I.NE.N+1) RETURN
                ELSE
                     A=((3*E-6)*E*E+4)/24
                     B=(-3*E+4)*E/(8*H)
                     C=(3*E-2)/(8*H*H)
                     D=-1/(8*H**3)
                     IF(I.NE.2.AND.I.NE.N+2) RETURN
                END IF
           END IF
      END IF
      IF(I.LE.2) THEN
           AA=1.0D+00/24
           BB=-1/(8*H)
           CC=1/(8*H*H)
           DD=-1/(24*H**3)
           IF(I.EQ.2) THEN
                A=A-AA
                B=B-BB
                C=C-CC
                D=D-DD
                RETURN
           ELSE
                A=A-4*AA
                B=B-4*BB
                C=C-4*CC
                D=D-4*DD
                RETURN
           END IF
      ELSE
           EE=N+2
           AA=(((-EE+6)*EE-12)*EE+8)/24
           BB=((EE-4)*EE+4)/(8*H)
           CC=(-EE+2)/(8*H*H)
           DD=1/(24*H**3)
           IF(I.EQ.N+1) THEN
                A=A-AA
                B=B-BB
                C=C-CC
                D=D-DD
                RETURN
           ELSE
                A=A-4*AA
                B=B-4*BB
                C=C-4*CC
                D=D-4*DD
           END IF
      END IF
      RETURN
      END
C
      SUBROUTINE DPHICO(I,J,A,B,C)
C
C     SAME AS PHICO EXCEPT THE COEFFICIENTS ARE FOR PHI'(I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON N
      H=1.0D+00/(N+1)
      A=0
      B=0
      C=0
      E=I-1
      IF(J.LT.I-2 .OR. J.GE.I+2) RETURN
      IF(I.EQ.1 .AND. J.LT.I) RETURN
      IF(I.EQ.2 .AND. J.LT.I-1) RETURN
      IF(I.EQ.N+1 .AND. J.GT.N+1) RETURN
      IF(I.EQ.N+2 .AND. J.GE.N+2) RETURN
      IF(J.LE.I-2) THEN
           A=((E-4)*E+4)/(8*H)
           B=(-E+2)/(4*H*H)
           C=1/(8*H**3)
           RETURN
      ELSE
           IF(J.GT.I-1.AND.J.GT.I) THEN
                A=((-E-4)*E-4)/(8*H)
                B=(E+2)/(4*H*H)
                C=-1/(8*H**3)
                RETURN
           ELSE
                IF(J.GT.I-1) THEN
                     A=(3*E+4)*E/(8*H)
                     B=(-3*E-2)/(4*H*H)
                     C=3/(8*H**3)
                     IF(I.NE.1 .AND.I.NE.N+1) RETURN
                ELSE
                     A=(-3*E+4)*E/(8*H)
                     B=(3*E-2)/(4*H*H)
                     C=-3/(8*H**3)
                     IF(I.NE.2 .AND.I.NE.N+2) RETURN
                END IF
           END IF
      END IF
      IF(I.LE.2) THEN
           AA=-1/(8*H)
           BB=1/(4*H*H)
           CC=-1/(8*H**3)
           IF(I.EQ.2) THEN
                A=A-AA
                B=B-BB
                C=C-CC
                RETURN
           ELSE
                A=A-4*AA
                B=B-4*BB
                C=C-4*CC
                RETURN
           END IF
      ELSE
           EE=N+2
           AA=((EE-4)*EE+4)/(8*H)
           BB=(-EE+2)/(4*H*H)
           CC=1/(8*H**3)
           IF(I.EQ.N+1) THEN
                A=A-AA
                B=B-BB
                C=C-CC
                RETURN
           ELSE
                A=A-4*AA
                B=B-4*BB
                C=C-4*CC
           END IF
      END IF
      RETURN
      END
C
      FUNCTION INTY(J,JJ)
C
C     CORRESPONDS THE INTERVAL (X(JJ),X(JJ+1)) TO PROPER INDEX
C     K IN CO(J,K,L),L=1,...,4 AND IN DCO(J,K,L),L=1,2,3
      IMPLICIT REAL*8(A-H,O-Z)
      IF(JJ.LE.J-3 .OR. JJ.GE.J+2) THEN
           WRITE(6,2)
2          FORMAT(1X,'ERROR IN INTY')
           STOP
      ELSE
           IF(JJ.EQ.J-2) INTY=1
           IF(JJ.EQ.J-1) INTY=2
           IF(JJ .EQ. J) INTY=3
           IF(JJ.EQ.J+1) INTY=4
      END IF
      RETURN
      END
C
      FUNCTION P(X)
      IMPLICIT REAL*8(A-H,O-Z)
      P=1.0D+00
      RETURN
      END
C
      FUNCTION Q(X)
      IMPLICIT REAL*8(A-H,O-Z)
      Q=3.141593*3.141593
      RETURN
      END
C
      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      F=2*(3.141593)**2*DSIN(3.141593*X)
      RETURN
      END
