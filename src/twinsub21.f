      SUBROUTINE POLISH(M,MM,N,NN,NDAT,X,XX,Y,YY,
     1IIROW,JJCOL,RTOT,ROWWGT,CCWT,COLWGT,IADDR,IDAT)
C SUBROUTINE POLARIZES ORDINATION X SO THAT ORDINATION APPROXIMATES
C BETTER TO A CLASSIFICATION.  TO DO THIS, DERIVES TWO NEW
C ORDINATIONS, TAKES THEIR SUM, AND PUTS BACK IN X.
C ORDINATION 1 IS ADDITIVE, DERIVED BY ADDITNG PREFERENCE SCORES IN
C A SUMMATION IN WHICH ONLY GOOD (I.E. REASONABLY FREQUENT AND
C REASONABLY PREFERENTIAL) SPECIES ARE GIVEN APPRECIABLE WEIGHT.
C ORDINATION 2 TAKES THE MEAN PREFERENCE SCORE OF THE SPECIES IN
C     EACH STAND, UNWEIGHTED EXCEPT FOR BASIC WEIGHTS CCWT.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IIROW(M),JJCOL(N),IADDR(MM),IDAT(NDAT)
      REAL(8) X(M),XX(M),Y(N),YY(N),RTOT(M),ROWWGT(M),CCWT(NN),COLWGT(N)
C---Variable PRECIS added P.Minchin  June 1997
      COMMON/LIMS/RARE,FEEBLE,FRQLIM,TOL,RATLIM,REPLIM,PRECIS
      COMMON/IARBS/ICWEXP,IEND,MMIN,IPREXP,LEVMAX
      COMMON/ARBS/CWTMIN,CRLONG,CRCUT
      CALL XMAXMI(X,AXMAX,AXMIN,M)
      CUTMID=(AXMAX+AXMIN)*0.5
      CUTHLF=(AXMAX-AXMIN)*0.5*CRCUT
      CUT1=CUTMID-CUTHLF
      CUT2=CUTMID+CUTHLF
      CALL INDSCO(M,MM,N,NDAT,X,AXPOS,AXNEG,CUT1,CUT2,
     1Y,YY,IIROW,ROWWGT,IADDR,IDAT)
C PUT PREFERENCE SCORE IN YY AND CALCULATE APPROPRIATE COLUMN WEIGHTS
      PRLIM=(RATLIM-1.0)/(RATLIM+1.0)
      DO 10 J=1,N
      JJ=JJCOL(J)
      AY=Y(J)/AXNEG
      AYY=YY(J)/AXPOS
      PREF=(AYY-AY)/(AYY+AY)
      FREQ=AYY+AY
      IF(FREQ.GT.FRQLIM) FREQ=FRQLIM
      IF(PREF.GT.PRLIM)PREF=PRLIM
      IF(PREF.LT.-PRLIM)PREF=-PRLIM
      IF(ABS(PREF).LT.0.001)PREF=0.001
      COLWGT(J)=CCWT(JJ)*(FREQ/FRQLIM)**ICWEXP*(ABS(PREF)/PRLIM)**IPREXP
      Y(J)=PREF/PRLIM
   10 CONTINUE
      CALL YXMULT(M,MM,N,NDAT,Y,X,IIROW,IADDR,ROWWGT,COLWGT,IDAT)
      DO 20 I=1,M
   20 X(I)=X(I)/ROWWGT(I)
      CALL XMAXMI(X,AXMAX,AXMIN,M)
      IF(ABS(AXMIN).GT.AXMAX)AXMAX=ABS(AXMIN)
      DO 30 I=1,M
   30 X(I)=X(I)/AXMAX
C WE NOW ABOLISH COLUMN WEIGHTS
      DO 40 J=1,N
      YY(J)=1.0
      JJ=JJCOL(J)
   40 COLWGT(J)=CCWT(JJ)
      CALL YXMULT(M,MM,N,NDAT,YY,RTOT,IIROW,IADDR,
     1ROWWGT,COLWGT,IDAT)
      CALL YXMULT(M,MM,N,NDAT,Y,XX,IIROW,IADDR,ROWWGT,COLWGT,IDAT)
      DO 50 I=1,M
   50 X(I)=X(I)+XX(I)/RTOT(I)
      RETURN
      END
C
