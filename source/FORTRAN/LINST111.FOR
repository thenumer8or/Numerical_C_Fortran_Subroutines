C***********************************************************************
C                                                                      *
C                    LINEAR SHOOTING METHOD                            *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION OF THE BOUNDARY-VALUE PROBLEM
C
C          -Y'' + P(X)Y' + Q(X)Y + R(X) = 0, A <= X <= B,
C                  Y(A) = ALPHA, Y(B) = BETA
C
C     NOTE: EQUATIONS (11.5),(11.6) ARE WRITTEN AS FIRST ORDER
C           SYSTEMS AND SOLVED.
C
C     INPUT:   ENDPOINTS A,B; BOUNDARY CONDITIONS ALPHA, BETA; NUMBER OF
C              SUBINTERVALS N.
C
C     OUTPUT:  APPROXIMATIONS W(1,I) TO Y(X(I)); W(2,I) TO Y'(X(I))
C              FOR EACH I=0,1,...N.
C
C     VECTORS WILL NOT BE USED FOR W1 AND W2
C     INITIALIZATION
      DIMENSION U(2,10),V(2,10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      P(X)=-2/X
      Q(X)=2/X**2
      R(X)=SIN(ALOG(X))/X**2
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Linear Shooting Method.'
      WRITE(6,*) 'Have the functions P, Q and R been created? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
110      IF (OK) GOTO 11
            WRITE(6,*) 'Input left and right endpoints separated by'
            WRITE(6,*) 'blank - include decimal point'
            WRITE(6,*) ' '
            READ(5,*) A, B
            IF (A.GE.B) THEN
               WRITE(6,*) 'Left endpoint must be less'
               WRITE(6,*) 'than right endpoint'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 110
11       OK = .FALSE.
         WRITE(6,*) 'Input Y(',A,')'
         WRITE(6,*) ' '
         READ(5,*) ALPHA
         WRITE(6,*) 'Input Y(',B,')'
         WRITE(6,*) ' '
         READ(5,*) BETA
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input a positive integer for the number'
            WRITE(6,*) 'of subinvervals '
            WRITE(6,*) ' '
            READ(5,*) N
            IF ( N .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'P, Q and R can be created '
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
      WRITE(OUP,*) 'LINEAR SHOOTING METHOD'
C     STEP 1
      H=(B-A)/N
C     U1=U(1,0), U2=U(2,0), V1=V(1,0), V2=V(2,0)
      U1=ALPHA
      U2=0.0
      V1=0.0
      V2=1.0
C     STEP 2
C     RUNGE-KUTTA METHOD FOR SYSTEMS IS USED IN STEPS 3, 4
C     THE INDEX I IS SHIFTED BY 1
      DO 10 I=1,N
C          STEP 3
           X=A+(I-1)*H
C          STEP 4
           T=X+H/2
           XK11=H*U2
           XK12=H*(P(X)*U2+Q(X)*U1+R(X))
           XK21=H*(U2+XK12/2)
           XK22=H*(P(T)*(U2+XK12/2)+Q(T)*(U1+XK11/2)+R(T))
           XK31=H*(U2+XK22/2)
           XK32=H*(P(T)*(U2+XK22/2)+Q(T)*(U1+XK21/2)+R(T))
           XK41=H*(U2+XK32)
           XK42=H*(P(X+H)*(U2+XK32)+Q(X+H)*(U1+XK31)+R(X+H))
           U1=U1+(XK11+2*XK21+2*XK31+XK41)/6
           U2=U2+(XK12+2*XK22+2*XK32+XK42)/6
           XK11=H*V2
           XK12=H*(P(X)*V2+Q(X)*V1)
           XK21=H*(V2+XK12/2)
           XK22=H*(P(T)*(V2+XK12/2)+Q(T)*(V1+XK11/2))
           XK31=H*(V2+XK22/2)
           XK32=H*(P(T)*(V2+XK22/2)+Q(T)*(V1+XK21/2))
           XK41=H*(V2+XK32)
           XK42=H*(P(X+H)*(V2+XK32)+Q(X+H)*(V1+XK31))
           V1=V1+(XK11+2*XK21+2*XK31+XK41)/6
           V2=V2+(XK12+2*XK22+2*XK32+XK42)/6
           U(1,I)=U1
           U(2,I)=U2
           V(1,I)=V1
           V(2,I)=V2
10    CONTINUE
      WRITE(OUP,1)
C     STEP 5
      W1= ALPHA
      Z=(BETA-U(1,N))/V(1,N)
      X=A
      I=0
      WRITE(OUP,2) I,X,W1,Z
C     STEP 6
      DO 20 I=1,N
           X=A+I*H
           W1=U(1,I)+Z*V(1,I)
           W2=U(2,I)+Z*V(2,I)
C     OUTPUT IS X(I), W(1,I), W(2,I)
20    WRITE(OUP,2) I,X,W1,W2
C     STEP 7
C     PROCESS IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'OUTPUT:I,X(I),W(1,I),W(2,I)',/)
2     FORMAT(1X,I2,3(3X,E15.8))
      END
