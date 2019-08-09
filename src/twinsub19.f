      SUBROUTINE RA(M,MM,N,NN,NDAT,EIG,X,XX,X3,X4,X5,Y,
     1RTOT,CTOT,IIROW,IADDR,ROWWGT,COLWGT,IDAT)
C DOES RECIPROCAL AVERAGING, DERIVING EIGENVECTOR X AND EIGENVALUE
C EIG.  ROWWGT, COLWGT ARE ROW AND COLUMN WEIGHTS.  XX,X3,X4,X5,Y ARE
C WORKSPACE.
      REAL X(M),XX(M),Y(N),RTOT(M),CTOT(N),ROWWGT(M),COLWGT(N)
      REAL X3(M),X4(M),X5(M)
      INTEGER IIROW(M),IADDR(MM),IDAT(NDAT)
C---FOLLOWING COMMON BLOCK ADDED BY P.R.MINCHIN SEPT 1986
c      COMMON /LUNITS/ IUIN1,IUIN2,IUOUT1,IUOUT2,IUOUT3
C---Variable PRECIS added P.Minchin  June 1997
      COMMON/LIMS/RARE,FEEBLE,FRQLIM,TOL,RATLIM,REPLIM,PRECIS
      COMMON/SWITCH/IFAIL,IDIAGR,ISTOP,IREWT
      IFAIL=0
      DO 10 I=1,M
   10 XX(I)=1.0
      DO 20 J=1,N
   20 Y(J)=1.0
      CALL XYMULT(M,MM,N,NN,NDAT,XX,CTOT,IIROW,IADDR,ROWWGT,
     1COLWGT,IDAT)
      CALL YXMULT(M,MM,N,NN,NDAT,Y,RTOT,IIROW,IADDR,ROWWGT,
     1COLWGT,IDAT)
      TOT=0.0
      DO 30 J=1,N
   30 TOT=TOT+CTOT(J)
      DO 40 I=1,M
   40 X(I)=FLOAT(I)
      X(1)=1.1
C WE NOW HAVE TRIAL VECTOR X, DELIBERATELY MADE IRREGULAR TO AVOID
C NASTY COINCIDENCES
      TTOL=1.0E-5
      ICOUNT=0
      IITIM=0
   50 A=0.0
C FIRST CENTER TO ZERO MEAN AND UNIT LENGTH
      DO 60 I=1,M
   60 A=A+X(I)*RTOT(I)
      A=A/TOT
      EX=0.0
      DO 70 I=1,M
      AX=X(I)-A
      EX=EX+AX*AX*RTOT(I)
   70 X(I)=AX
      EX=SQRT(EX)
      DO 80 I=1,M
   80 X(I)=X(I)/EX
      A11=0.0
      A12=0.0
      A22=0.0
      A23=0.0
      A33=0.0
      A34=0.0
      A44=0.0
      B13=0.0
      B14=0.0
      B24=0.0
      CALL XYMULT(M,MM,N,NN,NDAT,X,Y,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      DO 81 J=1,N
   81 Y(J)=Y(J)/CTOT(J)
      CALL YXMULT(M,MM,N,NN,NDAT,Y,XX,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      A=0.0
      DO 82 I=1,M
      AX=XX(I)
      XX(I)=AX/RTOT(I)
      A=A+AX
   82 A11=A11+AX*X(I)
      IF(A11.GT.TTOL) GOTO 83
      IFAIL=1
      RETURN
   83 A=A/TOT
      DO 84 I=1,M
      AX=XX(I)-A-A11*X(I)
      A12=A12+AX*AX*RTOT(I)
   84 XX(I)=AX
      A12=SQRT(A12)
      DO 86 I=1,M
   86 XX(I)=XX(I)/A12
      IF(A12.LT.TOL) GOTO 200
C---Maximum iteration limit increased by P.Minchin Jan 1997
C      IF(ICOUNT.GT.5) GOTO 200
      IF(ICOUNT.GT.999) GOTO 200
      ICOUNT=ICOUNT+1
      CALL XYMULT(M,MM,N,NN,NDAT,XX,Y,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      DO 91 J=1,N
   91 Y(J)=Y(J)/CTOT(J)
      CALL YXMULT(M,MM,N,NN,NDAT,Y,X3,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      A=0.0
      DO 92 I=1,M
      AX=X3(I)
      X3(I)=AX/RTOT(I)
      A=A+AX
      A22=A22+AX*XX(I)
   92 B13=B13+AX*X(I)
      IF(A22.LT.TTOL) GOTO 120
      A=A/TOT
      DO 94 I=1,M
      AX=X3(I)-A-A22*XX(I)-B13*X(I)
      A23=A23+AX*AX*RTOT(I)
   94 X3(I)=AX
      A23=SQRT(A23)
      DO 96 I=1,M
   96 X3(I)=X3(I)/A23
      CALL XYMULT(M,MM,N,NN,NDAT,X3,Y,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      DO 101 J=1,N
  101 Y(J)=Y(J)/CTOT(J)
      CALL YXMULT(M,MM,N,NN,NDAT,Y,X4,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      A=0.0
      DO 102 I=1,M
      AX=X4(I)
      X4(I)=AX/RTOT(I)
      A=A+AX
      A33=A33+AX*X3(I)
      B14=B14+AX*X(I)
  102 B24=B24+AX*XX(I)
      IF(A33.LT.TTOL) GOTO 120
      A=A/TOT
      DO 104 I=1,M
      AX=X4(I)-(A+A33*X3(I)+B14*X(I)+B24*XX(I))
      A34=A34+AX*AX*RTOT(I)
  104 X4(I)=AX
      A34=SQRT(A34)
      DO 106 I=1,M
  106 X4(I)=X4(I)/A34
      CALL XYMULT(M,MM,N,NN,NDAT,X4,Y,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      DO 111 J=1,N
  111 Y(J)=Y(J)/CTOT(J)
      CALL YXMULT(M,MM,N,NN,NDAT,Y,X5,IIROW,
     1IADDR,ROWWGT,COLWGT,IDAT)
      DO 112 I=1,M
  112 A44=A44+X4(I)*X5(I)
  120 AX1=1.0
      AX2=0.1
      AX3=0.01
      AX4=0.001
      DO 130 ITIMES=1,100
      AXX1=A11*AX1+A12*AX2
      AXX2=A12*AX1+A22*AX2+A23*AX3
      AXX3=A23*AX2+A33*AX3+A34*AX4
      AXX4=A34*AX3+A44*AX4
      AX1=A11*AXX1+A12*AXX2
      AX2=A12*AXX1+A22*AXX2+A23*AXX3
      AX3=A23*AXX2+A33*AXX3+A34*AXX4
      AX4=A34*AXX3+A44*AXX4
      EX=SQRT(AX1**2+AX2**2+AX3**2+AX4**2)
      AX1=AX1/EX
      AX2=AX2/EX
      AX3=AX3/EX
      AX4=AX4/EX
      IF(ITIMES.NE.(ITIMES/5)*5) GOTO 130
      EXX=SQRT(EX)
      RESI=SQRT((AX1-AXX1/EXX)**2+(AX2-AXX2/EXX)**2+
     1(AX3-AXX3/EXX)**2+(AX4-AXX4/EXX)**2)
      IF(RESI.LT.TOL*0.05) GOTO 135
  130 CONTINUE
  135 IITIM=IITIM+ITIMES
      DO 140 I=1,M
  140 X(I)=AX1*X(I)+AX2*XX(I)+AX3*X3(I)+AX4*X4(I)
      GOTO 50
c  200 WRITE (IUOUT2,1000) A11,ICOUNT
c 1000 FORMAT('EIGENVALUE',F6.3,'  AT ITERATION',I4)
c      IF(A12.GT.TOL) WRITE (IUOUT2,1010) ICOUNT,A12,TOL
c 1010 FORMAT('RA TROUBLE',I4,' ITERATIONS, AND RESIDUAL IS STILL ',
c     1F6.3,'  INSTEAD OF ',F8.6,' (THE TOLERANCE)')
      EIG=A11
      RETURN
      END
C
