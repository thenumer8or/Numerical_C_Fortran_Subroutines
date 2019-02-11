C******************************************************************************
C                                                                             *
C                        FINITE ELEMENT METHOD                                *
C                                                                             *
C******************************************************************************
C
C     TO APPROXIMATE THE SOLUTION TO AN ELLIPTIC PARTIAL-DIFFERENTIAL
C     EQUATION SUBJECT TO DIRICHLET, MIXED, OR NEUMANN BOUNDARY
C     CONDITIONS:
C
C     INPUT   SEE STEP 0
C
C     OUTPUT  DESCRIPTION OF TRIANGLES, NODES, LINE INTEGRALS, BASIS
C             FUNCTIONS, LINEAR SYSTEM TO BE SOLVED, AND THE
C             COEFFICIENTS OF THE BASIS FUNCTIONS
C
C
C     STEP 0
C     GENERAL OUTLINE
C
C        1. TRIANGLES NUMBERED: 1 TO K FOR TRIANGLES WITH NO EDGES ON
C           SCRIPT-S-1 OR SCRIPT-S-2, K+1 TO N FOR TRIANGLES WITH
C           EDGES ON SCRIPT-S-2, N+1 TO M FOR REMAINING TRIANGLES.
C           NOTE: K=0 IMPLIES THAT NO TRIANGLE IS INTERIOR TO D.
C           NOTE: M=N IMPLIES THAT ALL TRIANGLES HAVE EDGES ON
C           SCRIPT-S-2.
C
C        2. NODES NUMBERED: 1 TO LN FOR INTERIOR NODES AND NODES ON
C           SCRIPT-S-2, LN+1 TO LM FOR NODES ON SCRIPT-S-1.
C           NOTE: LM AND LN REPRESENT LOWER CASE M AND N RESP.
C           NOTE: LN=LM IMPLIES THAT SCRIPT-S-1 CONTAINS NO NODES.
C           NOTE: IF A NODE IS ON BOTH SCRIPT-S-1 AND SCRIPT-S-2, THEN
C           IT SHOULD BE TREATED AS IF IT WERE ONLY ON SCRIPT-S-1.
C
C        3. NL=NUMBER OF LINE SEGMENTS ON SCRIPT-S-2
C           LINE(I,J) IS AN NL BY 2 ARRAY WHERE
C           LINE(I,1)= FIRST NODE ON LINE I AND
C           LINE(I,2)= SECOND NODE ON LINE I TAKEN
C           IN POSITIVE DIRECTION ALONG LINE I
C
C        4. FOR THE NODE LABELLED KK,KK=1,...,LM WE HAVE:
C           A) COORDINATES XX(KK),YY(KK)
C           B) NUMBER OF TRIANGLES IN WHICH KK IS A VERTEX= LL(KK)
C           C) II(KK,J) LABELS THE TRIANGLES KK IS IN AND
C              NV(KK,J) LABELS WHICH VERTEX NODE KK IS FOR
C              EACH J=1,...,LL(KK)
C
C        5. NTR IS AN M BY 3 ARRAY WHERE
C           NTR(I,1)=NODE NUMBER OF VERTEX 1 IN TRIANGLE I
C           NTR(I,2)=NODE NUMBER OF VERTEX 2 IN TRIANGLE I
C           NTR(I,3)=NODE NUMBER OF VERTEX 3 IN TRIANGLE I
C
C        6. FUNCTION SUBPROGRAMS:
C           A) P,Q,R,F,G,G1,G2 ARE THE FUNCTIONS SPECIFIED BY
C              THE PARTICULAR DIFFERENTIAL EQUATION
C           B) RR IS THE INTEGRAND
C              R*N(J)*N(K) ON TRIANGLE I IN STEP 4
C           C) FFF IS THE INTEGRAND F*N(J) ON TRIANGLE I IN STEP 4
C           D) GG1=G1*N(J)*N(K)
C              GG2=G2*N(J)
C              GG3=G2*N(K)
C              GG4=G1*N(J)*N(J)
C              GG5=G1*N(K)*N(K)
C              INTEGRANDS IN STEP 5
C           E) QQ(FF) COMPUTES THE DOUBLE INTEGRAL OF ANY
C              INTEGRAND FF OVER A TRIANGLE WITH VERTICES GIVEN BY
C              NODES J1,J2,J3 - THE METHOD IS AN O(H**2) APPROXIMATION
C              FOR TRIANGLES
C           F) SQ(PP) COMPUTES THE LINE INTEGRAL OF ANY INTEGRAND PP
C              ALONG THE LINE FROM (XX(J1),YY(J1)) TO (XX(J2),YY(J2))
C              BY USING A PARAMETERIZATION GIVEN BY:
C                X=XX(J1)+(XX(J2)-XX(J1))*T
C                Y=YY(J1)+(YY(J2)-YY(J1))*T
C              FOR 0 <= T <= 1
C              AND APPLYING SIMPSON'S COMPOSITE METHOD WITH H=.01
C
C        7. ARRAYS:
C           A) A,B,C ARE M BY 3 ARRAYS WHERE THE BASIS FUNCTION N
C              FOR THE ITH TRIANGLE, JTH VERTEX IS
C              N(X,Y)=A(I,J)+B(I,J)*X+C(I,J)*Y
C              FOR J=1,2,3 AND I=1,2,...,M
C           B) XX,YY ARE LM BY 1 ARRAYS TO HOLD COORDINATES OF NODES
C           C) LINE,LL,II,NV,NTR HAVE BEEN EXPLAINED ABOVE
C           D) GAMMA AND ALPHA ARE CLEAR
C
C        8. NOTE THAT A,B,C,XX,YY,I,I1,I2,J1,J2,J3,DELTA ARE LABELLED
C           COMMON STORAGE SO THAT IN ANY SUBPROGRAM WE KNOW THAT
C           TRIANGLE I HAS VERTICES (XX(J1),YY(J1)),(XX(J2),YY(J2)),
C           (XX(J3),YY(J3)). THAT IS, VERTEX 1 IS NODE J1, VERTEX 2 IS
C           NODE J2, VERTEX 3 IS NODE J3 UNLESS A LINE INTEGRAL IS
C           INVOLVED IN WHICH CASE THE LINE INTEGRAL GOES FROM NODE J1
C           TO NODE J2 IN TRIANGLE I OR UNLESS VERTEX I1 IS NODE J1
C           AND VERTEX I2 IS NODE J2 - THE BASIS FUNCTIONS INVOLVE
C           A(I,I1)+B(I,I1)*X+C(I,I1)**Y FOR VERTEX I1 IN TRIANGLE I
C           AND A(I,I2)+B(I,I2)*X+C(I,I2)*Y FOR VERTEX I2 IN TRIANGLE I
C           DELTA IS 1/2 THE AREA OF TRIANGLE I
C
C     TO CHANGE PROBLEMS:
C        1) CHANGE FUNCTION SUBPROGRAMS P,Q,R,F,G,G1,G2
C        2) CHANGE DATA INPUT FOR K,N,M,LN,LM,NL.
C        3) CHANGE DATA INPUT FOR XX,YY,LLL,II,NV.
C        4) CHANGE DATA INPUT FOR LINE.
C        5) CHANGE DIMENSION STATEMENTS TO READ :
C           DIMENSION A(M,3),B(M,3),C(M,3),XX(LM),YY(LM)
C                   THERE ARE 10 PLACES TO CHANGE.
C           DIMENSION LINE(NL,2),LL(LM),II(LM,MAX LL(LM)),
C                NV(LM,MAX LL(LM))--1 PLACE TO CHANGE
C           DIMENSION NTR(M,3),GAMMA(LM),ALPHA(LN,LN+1)
C                   THERE IS 1 PLACE TO CHANGE
C
      EXTERNAL P,Q,RR,FFF,GG1,GG2,GG3,GG4,GG5
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      DIMENSION LINE(5,2),LL(11),II(11,5),NV(11,5)
      DIMENSION NTR(10,3),GAMMA(11),ALPHA(5,6)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Finite Element Method.'
      OK = .FALSE.
      WRITE(6,*) 'The input will be from a text file in the following'
      WRITE(6,*) ' form:'
      WRITE(6,*) 'K     N     M     n     m     nl'
      WRITE(6,*) 'Each of the above is an integer -'
      WRITE(6,*) 'separate with at least one blank.'
      WRITE(6,*) 'Follow with the input for each node in the form:'
      WRITE(6,*) 'x-coor., y-coord., number of triangles in which the'
      WRITE(6,*) ' node is a vertex.'
      WRITE(6,*) 'Continue with the triangle number and vertex number'
      WRITE(6,*) 'for each triangle in which the node is a vertex.'
      WRITE(6,*) 'Separate each entry with at least one blank.  After'
      WRITE(6,*) 'all nodes have been entered follow with information'
      WRITE(6,*) 'on the lines over which line integrals must be'
      WRITE(6,*) 'computed. The format'
      WRITE(6,*) 'of this data will be the node number of the starting'
      WRITE(6,*) 'node, followed by the node number of the ending node'
      WRITE(6,*) 'for each line, taken in the positive direction.'
      WRITE(6,*) 'There should be 2 * nl such entries, each an integer'
      WRITE(6,*) 'separated by a blank.  Have all'
      WRITE(6,*) 'the functions P,Q,R,F,G,G1,G2 been created and'
      WRITE(6,*) 'has the input file been created.  Answer Y or N.'
      WRITE(6,*) ' '
      READ(5,*) AA
      IF((AA .EQ. 'y') .OR. (AA .EQ. 'Y')) THEN
         OK = .TRUE.
         WRITE(6,*) 'Input the file name in the form - drive: name.ext'
         WRITE(6,*) 'for example:  A:DATA.DTA'
         WRITE(6,*) ' '
         READ(5,*) NAME
         INP = 4
         OPEN(UNIT=INP,FILE=NAME,ACCESS='SEQUENTIAL')
         READ(4,*) K,N,M,LN,LM,NL
C        INPUT COORDINATES OF EACH NODE, NUMBER OF TRIANGLES NODE
C        IS IN (LLL) , WHICH TRIANGLES THE NODE IS IN, AND WHICH
C        VERTEX OF EACH TRIANGLE THE NODE IS IN
         DO 2 KK=1,LM
            READ(4,*) XX(KK),YY(KK),LLL,(II(KK,J),NV(KK,J),J=1,LLL)
3           FORMAT(2E15.8,1X,I2,1X,5(I2,1X,I2))
            LL(KK)=LLL
C           COMPUTE ENTRIES FOR NTR
            DO 4 J=1,LLL
                N1=II(KK,J)
                N2=NV(KK,J)
4               NTR(N1,N2)=KK
2        CONTINUE
C        INPUT LINE INFORMATION
         IF(NL.GT.0) THEN
           READ(4,*) ((LINE(I,J),J=1,2),I=1,NL)
         END IF
         CLOSE(UNIT=4)
      ELSE
         WRITE(6,*) 'The program will end so that the input file'
         WRITE(6,*) 'can be created or so the functions'
         WRITE(6,*) 'can be created.'
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
      WRITE(OUP,*) 'FINITE ELEMENT METHOD'
      K1=K+1
      N1=LN+1
C     OUTPUT INFORMATION ON TRIANGLES, VERTICES AND NODES
      WRITE(OUP,100)
100   FORMAT(1X,'TRIANGLES ARE AS FOLLOWS')
      WRITE(OUP,101)
101   FORMAT(1X,'TRIANGLE',20X,'VERTEX 1',20X,'VERTEX 2',20X,
     1'VERTEX 3',/)
      WRITE(OUP,102) (I,(NTR(I,J),J=1,3),I=1,M)
102   FORMAT(1X,I3,21X,'NODE ',I3,20X,'NODE ',I3,20X,'NODE ',I3)
      WRITE(OUP,103)
103   FORMAT(1X,'NODE',25X,'X-COORD.',25X,'Y-COORD.',/)
      WRITE(OUP,104) (KK,XX(KK),YY(KK),KK=1,LM)
104   FORMAT(1X,I3,22X,E15.8,18X,E15.8)
C     OUTPUT INFO ON LINES
      WRITE(OUP,109)
109   FORMAT(1X,'DESCRIPTION OF LINES FOR LINE INTEGRALS',//)
      WRITE(OUP,110) (L1,LINE(L1,1),LINE(L1,2),L1=1,NL)
110   FORMAT(1X,'LINE',I3,'  FROM NODE',I3,' TO NODE',I3)
C     STEP 1
      IF(LM.NE.LN) THEN
           DO 5 L=N1,LM
5          GAMMA(L)=G(XX(L),YY(L))
      END IF
C     STEP 2 - INITIALIZATION OF ALPHA AND NOTE THAT
C     ALPHA(I,LN+1)=BETA(I)
      DO 6 I=1,LN
            DO 6 J=1,N1
6           ALPHA(I,J)=0.0
C     STEPS 3,4 AND 6 THROUGH 12 ARE WITHIN THE NEXT LOOP
C     FOR EACH TRIANGLE I LET NODE J1 BE VERTEX 1,NODE J2 BE VERTEX 2,
C     AND NODE J3 BE VERTEX 3
C     STEP 3
      DO 7 I=1,M
           J1=NTR(I,1)
           J2=NTR(I,2)
           J3=NTR(I,3)
           DELTA=XX(J2)*YY(J3)-XX(J3)*YY(J2)-XX(J1)*(YY(J3)-YY(J2))
     *     +YY(J1)*(XX(J3)-XX(J2))
           A(I,1)=(XX(J2)*YY(J3)-YY(J2)*XX(J3))/DELTA
           B(I,1)=(YY(J2)-YY(J3))/DELTA
           C(I,1)=(XX(J3)-XX(J2))/DELTA
           A(I,2)=(XX(J3)*YY(J1)-YY(J3)*XX(J1))/DELTA
           B(I,2)=(YY(J3)-YY(J1))/DELTA
           C(I,2)=(XX(J1)-XX(J3))/DELTA
           A(I,3)=(XX(J1)*YY(J2)-YY(J1)*XX(J2))/DELTA
           B(I,3)=(YY(J1)-YY(J2))/DELTA
           C(I,3)=(XX(J2)-XX(J1))/DELTA
C          STEP 4
C          I1=J FOR STEP 4 AND I1=K FOR STEP 7
           DO 8 I1=1,3
C          STEP 8
                JJ1=NTR(I,I1)
C               I2=K FOR STEP 4 AND I2=J FOR STEP 9
                DO 9 I2=1,I1
C                    STEP 10 AND STEP 4
                     JJ2=NTR(I,I2)
C                    ZZ=Z(I1,I2)**I
                     ZZ=B(I,I1)*B(I,I2)*QQ(P)+C(I,I1)*C(I,I2)*QQ(Q)
     *                    -QQ(RR)
C                        STEP 11 AND 12
                     IF(JJ1.LE.LN) THEN
                         IF(JJ2.LE.LN) THEN
                              ALPHA(JJ1,JJ2)=ALPHA(JJ1,JJ2) + ZZ
                              IF(JJ1.NE.JJ2) ALPHA(JJ2,JJ1) =
     *                         ALPHA(JJ2,JJ1)+ZZ
                         ELSE
                              ALPHA(JJ1,N1)=ALPHA(JJ1,N1)-ZZ*GAMMA(JJ2)
                         END IF
                     ELSE
                         IF(JJ2.LE.LN) ALPHA(JJ2,N1)=ALPHA(JJ2,N1)
     *                       -ZZ*GAMMA(JJ1)
                     END IF
9               CONTINUE
           HH=-QQ(FFF)
           IF(JJ1.LE.LN) ALPHA(JJ1,N1)=ALPHA(JJ1,N1)+HH
8          CONTINUE
7     CONTINUE
C     OUTPUT BASIS FUNCTIONS
      WRITE(OUP,105)
105   FORMAT(1X,'BASIS FUNCTIONS ON EACH TRIANGLE')
      WRITE(OUP,106)
106   FORMAT(1X,'TRIANGLE',6X,'VERTEX',9X,'NODE',33X,'FUNCTION',/)
      DO 108 I=1,M
      DO 108 J=1,3
108        WRITE(OUP,107) I,J,NTR(I,J),A(I,J),B(I,J),C(I,J)
107        FORMAT(1X,I3,10X,I3,10X,I3,10X,E15.8,' + ',E15.8,' * X + '
     *     ,E15.8,' * Y')
C     STEP 5
C     FOR EACH LINE SEGMENT JI=1,...,NL AND FOR EACH TRIANGLE I,
C     I=K1,...,N WHICH MAY CONTAIN LINE JI. SEARCH ALL 3 VERTICES
C     FOR POSSIBLE CORRESPONDENCE
C     STEP 5 AND STEPS 13 THROUGH 19
      IF(NL.NE.0.AND.N.NE.K) THEN
           DO 11 JI=1,NL
                DO 12 I=K1,N
                     DO 13 I1=1,3
C                    I1=J IN STEP 5 AND I1=K IN STEP 14
C                    STEP 15
                     J1=NTR(I,I1)
                     IF(LINE(JI,1).EQ.J1) THEN
                         DO 14 I2=1,3
C                        I2=K IN STEP 5 AND I2=J IN STEP 16
C                        STEP 17
                         J2=NTR(I,I2)
                         IF(LINE(JI,2).EQ.J2) THEN
C                        WE HAVE CORRESPONDENCE OF VERTEX I1
C                        IN TRIANGLE I WITH NODE J1 AS START OF LINE JI
C                        AND VERTEX I2 WITH NODE J2 AS END OF LINE JI
C                            STEP 5
                             XJ=SQ(GG1)
                             XJ1=SQ(GG4)
                             XJ2=SQ(GG5)
                             XI1=SQ(GG2)
                             XI2=SQ(GG3)
C                            STEP 8 AND 19
                             IF(J1.LE.LN) THEN
                                    IF(J2.LE.LN) THEN
                                         ALPHA(J1,J1)=ALPHA(J1,J1)+XJ1
                                         ALPHA(J1,J2)=ALPHA(J1,J2)+XJ
                                         ALPHA(J2,J2)=ALPHA(J2,J2)+XJ2
                                         ALPHA(J2,J1)=ALPHA(J2,J1)+XJ
                                         ALPHA(J1,N1)=ALPHA(J1,N1)+XI1
                                         ALPHA(J2,N1)=ALPHA(J2,N1)+XI2
                                    ELSE
                                         ALPHA(J1,N1)=ALPHA(J1,N1)-
     *                                      XJ*GAMMA(J2)+XI1
                                         ALPHA(J1,J1)=ALPHA(J1,J1)+XJ1
                                    END IF
                             ELSE
                                    IF(J2.LE.LN) THEN
                                         ALPHA(J2,N1)=ALPHA(J2,N1)-
     *                                       XJ*GAMMA(J1)+XI2
                                         ALPHA(J2,J2)=ALPHA(J2,J2)+XJ2
                                    END IF
                             END IF
                         END IF
14                       CONTINUE
                     END IF
13                   CONTINUE
12              CONTINUE
11         CONTINUE
      END IF
C     OUTPUT ALPHA
      WRITE(OUP,111)
111   FORMAT(1X,'ENTRIES OF ALPHA',/)
      WRITE(OUP,112) ((I,J,ALPHA(I,J),J=1,N1),I=1,LN)
112   FORMAT(3(1X,'ALPHA(',I3,',',I3,') = ',E15.8))
C     STEP 20
      IF(LN.GT.1) THEN
      INN=LN-1
      DO 30 I=1,INN
           I1=I+1
           DO 31 J=I1,LN
                CCC=ALPHA(J,I)/ALPHA(I,I)
                DO 32 J1=I1,N1
32              ALPHA(J,J1)=ALPHA(J,J1)-CCC*ALPHA(I,J1)
31         ALPHA(J,I)=0.0
30    CONTINUE
      END IF
      GAMMA(LN)=ALPHA(LN,N1)/ALPHA(LN,LN)
      IF(LN.GT.1) THEN
      DO 33 I=1,INN
           J=LN-I
           CCC=ALPHA(J,N1)
           J1=J+1
           DO 34 KK=J1,LN
34              CCC=CCC-ALPHA(J,KK)*GAMMA(KK)
33         GAMMA(J)=CCC/ALPHA(J,J)
      END IF
C     STEP 21
C     OUTPUT GAMMA
      WRITE(OUP,113)
113   FORMAT(1X,'COEFFICIENTS GAMMA',/)
      DO 115 I=1,LM
      LLL=LL(I)
115   WRITE(OUP,114) I,GAMMA(I),I,(II(I,J),J=1,LLL)
114   FORMAT(1X,'GAMMA(',I3,') = ',E15.8,'  USED WITH NODE',I3,' IN TRI
     1ANGLES ',5(I3,2X))
C     STEP 23
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
      END
      FUNCTION P(X,Y)
      P=1
      RETURN
      END
      FUNCTION Q(X,Y)
      Q=1
      RETURN
      END
      FUNCTION R(X,Y)
      R=0
      RETURN
      END
      FUNCTION F(X,Y)
      F=0
      RETURN
      END
      FUNCTION G(X,Y)
      G=4
      RETURN
      END
      FUNCTION G1(X,Y)
      G1=0.0
      RETURN
      END
      FUNCTION G2(X,Y)
      G2=0
      T=1.0E-05
      IF((.2-T.LE.X.AND.X.LE.4+T).AND.ABS(Y-.2).LE.T) G2=X
      IF((.5-T.LE.X.AND.X.LE..6+T).AND.ABS(Y-.1).LE.T) G2=X
      IF((-T.LE.Y.AND.Y.LE..1+T).AND.ABS(X-.6).LE.T) G2=Y
      IF((-T.LE.X.AND.X.LE..2+T).AND.ABS(Y+X-.4).LE.T)G2=(X+Y)/SQRT(
     * 2.)
      IF((.4-T.LE.X.AND.X.LE..5+T).AND.ABS(Y+X-.6).LE.T)G2=(X+Y)/SQRT(
     * 2.)
      RETURN
      END
      FUNCTION RR(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      RR=R(X,Y)*(A(I,I1)+B(I,I1)*X+C(I,I1)*Y)*(A(I,I2)+B(I,I2)*X+C(I,I2)
     1*Y)
      RETURN
      END
      FUNCTION FFF(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      FFF=F(X,Y)*(A(I,I1)+B(I,I1)*X+C(I,I1)*Y)
      RETURN
      END
      FUNCTION GG1(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      GG1=G1(X,Y)*(A(I,I1)+B(I,I1)*X+C(I,I1)*Y)*(A(I,I2)+B(I,I2)*X+
     1C(I,I2)*Y)
      RETURN
      END
      FUNCTION GG2(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      GG2=G2(X,Y)*(A(I,I1)+B(I,I1)*X+C(I,I1)*Y)
      RETURN
      END
      FUNCTION GG3(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      GG3=G2(X,Y)*(A(I,I2)+B(I,I2)*X+C(I,I2)*Y)
      RETURN
      END
      FUNCTION GG4(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      GG4=G1(X,Y)*(A(I,I1)+B(I,I1)*X+C(I,I1)*Y)**2
      RETURN
      END
      FUNCTION GG5(X,Y)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      GG5=G1(X,Y)*(A(I,I2)+B(I,I2)*X+C(I,I2)*Y)**2
      RETURN
      END
      FUNCTION QQ(FF)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      A1=(XX(J1)+XX(J2)+XX(J3))/3
      A2=(YY(J1)+YY(J2)+YY(J3))/3
      QQ=3*(FF(XX(J1),YY(J1))+FF(XX(J2),YY(J2))+FF(XX(J3),YY(J3)))/60
     1+8*(FF((XX(J1)+XX(J2))/2,(YY(J1)+YY(J2))/2)+FF((XX(J1)+XX(J3))/2,
     1(YY(J1)+YY(J3))/2)+FF((XX(J2)+XX(J3))/2,(YY(J2)+YY(J3))/2))/60
     1+27*FF(A1,A2)/60
      QQ=QQ*ABS(DELTA)/2
      RETURN
      END
      FUNCTION SQ(PP)
      COMMON A(10,3),B(10,3),C(10,3),XX(11),YY(11),I,I1,I2,J1,J2,J3
      COMMON DELTA
      DIMENSION X(101)
      F1(T)=PP(T1*T+X1,T2*T+Y1)*SQRT(T1*T1+T2*T2)
      X1=XX(J1)
      Y1=YY(J1)
      X2=XX(J2)
      Y2=YY(J2)
      T1=X2-X1
      T2=Y2-Y1
      H=.01
      DO 2 II=1,101
2     X(II)=(II-1)*H
      SQ=(F1(X(1))+F1(X(101)))*H/3
      Q1=0.0
      Q2=0.0
      DO 3 II=1,49
      Q1=Q1+2*H*F1(X(2*II+1))/3
3     Q2=Q2+4*H*F1(X(2*II))/3
      Q2=Q2+4*H*F1(X(100))/3
      SQ=SQ+Q1+Q2
      RETURN
      END
