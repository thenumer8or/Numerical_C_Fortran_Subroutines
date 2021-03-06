      SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
      DIMENSION C(4,4),Y(4),Y1(4),Y2(4),Y12(4),CL(16),X(16),WT(16,16)
      DATA WT/1.,0.,-3.,2.,4*0.,-3.,0.,9.,-6.,2.,0.,-6.,
     *  4.,8*0.,3.,0.,-9.,6.,-2.,0.,6.,-4.,10*0.,9.,-6.,
     *  2*0.,-6.,4.,2*0.,3.,-2.,6*0.,-9.,6.,2*0.,6.,-4.,
     *  4*0.,1.,0.,-3.,2.,-2.,0.,6.,-4.,1.,0.,-3.,2.,8*0.,
     *  -1.,0.,3.,-2.,1.,0.,-3.,2.,10*0.,-3.,2.,2*0.,3.,
     *  -2.,6*0.,3.,-2.,2*0.,-6.,4.,2*0.,3.,-2.,0.,1.,-2.,
     *  1.,5*0.,-3.,6.,-3.,0.,2.,-4.,2.,9*0.,3.,-6.,3.,0.,
     *  -2.,4.,-2.,10*0.,-3.,3.,2*0.,2.,-2.,2*0.,-1.,1.,
     *  6*0.,3.,-3.,2*0.,-2.,2.,5*0.,1.,-2.,1.,0.,-2.,4.,
     *  -2.,0.,1.,-2.,1.,9*0.,-1.,2.,-1.,0.,1.,-2.,1.,10*0.,
     *  1.,-1.,2*0.,-1.,1.,6*0.,-1.,1.,2*0.,2.,-2.,2*0.,-1.,1./
      D1D2=D1*D2
      DO 11 I=1,4
        X(I)=Y(I)
        X(I+4)=Y1(I)*D1
        X(I+8)=Y2(I)*D2
        X(I+12)=Y12(I)*D1D2
11    CONTINUE
      DO 13 I=1,16
        XX=0.
        DO 12 K=1,16
          XX=XX+WT(I,K)*X(K)
12      CONTINUE
        CL(I)=XX
13    CONTINUE
      L=0
      DO 15 I=1,4
        DO 14 J=1,4
          L=L+1
          C(I,J)=CL(L)
14      CONTINUE
15    CONTINUE
      RETURN
      END
