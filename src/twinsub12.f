      SUBROUTINE CLASS(MM,NN,NDAT,MIND,MMZ,MMS,IX,ICLASS,
     1IIROW,IADDR,INDPOT,INDORD,IZONE,IY,JJCOL,IDAT,
     2INDSIG,IPICT,X,XX,RTOT,RRWT,ROWWGT,Y,YY,CTOT,CCWT,
     3COLWGT,INAME1,INAME2,JNAME1,JNAME2,JNAM,X3,X4,X5,LIND,
     4INFLAG,INLEVMAX,INMMIN, EIGEN, INDICS, LIMPOS, ISEC)
C THIS DOES THE CLASSIFICATION.  IF ISEC.EQ.2 IT DOES A SPECIES
C CLASSIFICATION, PRINTING OUT LESS INFORMATION THAN IN THE CASE WHERE
C ISEC.EQ.1. Passes to caller EIGEN (eigen values of division), INDICS
C (signed indicator species for divisions) and LIMPOS (minimum indicator
C score for positive group).
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IX(MM),ICLASS(MM),IIROW(MM),IADDR(MM),INFLAG(MM),
     1INDPOT(NN),INDORD(NN),IZONE(MM),IY(NN),
     2JJCOL(NN),IDAT(NDAT),INDSIG(MMS),
     3IPICT(MMZ,MMS)
C---Changed ITEM, PLUS, MINUS to CHARACTER: P.Minchin 1997
      CHARACTER PLUS, MINUS
c     We do not want to have characters since we do not want to print
c     anything: change ITEM at the moment, and see later how to handle this
c     CHARACTER ITEM(200)*4
      INTEGER ITEM(200)
      INTEGER LIND(10)
C---Changed INAME, JNAME to CHARACTER: P.Minchin 1997
      INTEGER JNAM(NN)
      INTEGER INAME1(MM),INAME2(MM),JNAME1(NN),JNAME2(NN)
C---Declare function I2CHAR to be CHARACTER*4: P.Minchin June 1997
      CHARACTER*4 I2CHAR
      REAL(8) X(MM),XX(MM),RTOT(MM),RRWT(MM),ROWWGT(MM),Y(NN),YY(NN),
     1CTOT(NN),CCWT(NN),COLWGT(NN)
      REAL(8) X3(MM),X4(MM),X5(MM)
c     eigen passes the eigenvalues of each division, indics the indices
c     of indicator species so that negative indicators have negative
c     sign
      REAL(8) EIGEN(2**INLEVMAX-1)
      INTEGER INDICS(MIND, 2**INLEVMAX-1), LIMPOS(2**INLEVMAX-1)
      COMMON/WORK/ITEM
      COMMON/PICT/MZ,MS,IZD,IIZD,ISD,ISHIFT,MZIND,
     1MZCRIT,MZOUT,NIND
C---Variable PRECIS added P.Minchin  June 1997
      COMMON/LIMS/RARE,FEEBLE,FRQLIM,TOL,RATLIM,REPLIM,PRECIS
C      COMMON/SECOND/ISEC
      COMMON/SWITCH/IFAIL,IDIAGR,ISTOP,IREWT
      COMMON/IARBS/ICWEXP,IEND,MMIN,IPREXP,LEVMAX
      COMMON/ARBS/CWTMIN,CRLONG,CRCUT
C     Values for COMMON block variables
C     LIMS
      RARE = 0.2
      FEEBLE = 0.1
      FRQLIM = 0.2
      TOL = 5E-6
      RATLIM = 3.0
      REPLIM = 2.0
      PRECIS = 1E-7
C     IARBS: MMIN and LEVMAX were added as input arguments
      ICWEXP = 1
      MMIN = INMMIN
      IPREXP = 4
      LEVMAX = INLEVMAX
      IEND = 2**(LEVMAX+2)
C     ARBS
      CWTMIN = 0.01
      CRLONG = 0.2
      CRCUT = 0.2
C     SWITCH
C     idiagr is an argument to turn on drawing diagrams of division: never do
C     ifail and istop will be re-set later in this subroutine
      IFAIL = 1
      IDIAGR = 0
      ISTOP = 0
      IREWT = 2
C     PICT
C     MZ calculated similarly as in calling R function
C     Presumably calculated later, but token values given here:
C     izd, iizd, isd, ishift
      MZCRIT = 8
      MZOUT = 4
      MZ = MZCRIT + 2 * MZOUT
      MS = 1
      IZD = 0
      IIZD = 0
      ISD = 0
      ISHIFT = 1
      MZIND = 4
C---FOLLOWING COMMON BLOCK ADDED BY P.R.MINCHIN SEPT 1986
c      COMMON /LUNITS/ IUIN1,IUIN2,IUOUT1,IUOUT2,IUOUT3
      DATA PLUS/'+'/,MINUS/'-'/
      JLEV=0
      IC=1
      LEVEL=1
      ISTOP=0
      IEND=2**(LEVMAX+2)
      ID=1
      DO 5 II=1,MM
      IADDR(II)=ID
    4 JJ=IDAT(ID)
      ID=ID+1
      IF(JJ.EQ.-1) GOTO 5
      GOTO 4
    5 CONTINUE
      DO 10 II=1,MM
   10 ICLASS(II)=IC
      DO 15 II=1,MM
      ID=IADDR(II)
      IF(IDAT(ID).EQ.-1) ICLASS(II)=IEND
   15 CONTINUE
      GOTO 30
   20 CALL DECODE(M,MM,N,NN,NDAT,IIROW,IADDR,IDAT,JJCOL)
      CALL UPDATE(M,MM,IC,LEVEL,ICLASS,IIROW,IX)
      IF(IPUNCH.EQ.0) GO TO 21
C      IF(ISEC.EQ.2.AND.IPUNCH.EQ.1) GO TO 21
      KLEV=LEVEL - 1
      IF(KLEV.EQ.JLEV) GO TO 21
      JLEV=KLEV
c      IF(ISEC.EQ.1) WRITE (IUOUT3,1021) JLEV
c 1021 FORMAT('TWINSPAN QUADRAT CLASSIFICATION LEVEL ',I3)
c      IF(ISEC.EQ.2) WRITE (IUOUT3,1022) JLEV
c 1022 FORMAT('TWINSPAN SPECIES CLASSIFICATION LEVEL ',I3)
c      WRITE (IUOUT3,1023) (ICLASS(II),II=1,MM)
c 1023 FORMAT(26I3)
   21 CONTINUE
      IF(ISTOP.EQ.1) GOTO 500
   30 CALL RECODE(M,MM,N,NN,NDAT,IC,ICLASS,IIROW,JJCOL,IADDR,IDAT)
      CALL BIN(ITEM,IC)
c      WRITE (IUOUT2,1010)
c 1010 FORMAT(30('****'))
c      WRITE (IUOUT2,1000) IC,M,(ITEM(IT),IT=1,30)
c      WRITE (IUOUT1,1000) IC,M,(ITEM(IT),IT=1,30)
c 1000 FORMAT(/'DIVISION',I5,'  (N=',I4,')',
c     19X,'I.E. GROUP ',30A1)
      IF(M.GE.MMIN) GOTO 40
c      WRITE (IUOUT2,1001)
c 1001 FORMAT('DIVISION FAILS - THERE ARE TOO FEW ITEMS')
      IFAIL=1
      GOTO 20
   40 CALL WEIGHT(M,MM,N,NN,NDAT,X,Y,IIROW,IADDR,
     1RRWT,ROWWGT,JJCOL,CCWT,COLWGT,IDAT)
      CALL RA(M,MM,N,NN,NDAT,EIG,X,XX,X3,X4,X5,Y,RTOT,
     1     CTOT,IIROW,IADDR,ROWWGT,COLWGT,IDAT)
      EIGEN(IC) = EIG
      IF (IFAIL.EQ.1) GOTO 20
      CALL XMAXMI(X,AXMAX,AXMIN,M)
      IF (AXMAX.GT.-AXMIN) GOTO 45
      DO 44 I=1,M
   44 X(I)=-X(I)
      AX=AXMAX
      AXMAX=-AXMIN
      AXMIN=-AX
   45 CRMIN=AXMIN*CRLONG
      CRMAX=AXMAX*CRLONG
      CRMID=0.5*(CRMIN+CRMAX)
      CRHALF=0.5*(CRMAX-CRMIN)
      IF(IREWT.EQ.0) GOTO 48
      DO 47 IRE=1,IREWT
      DO 46 I=1,M
      AX=(X(I)-CRMID)/CRHALF
      IF(AX.GT.1.0) AX=1.0
      IF(AX.LT.-1.0) AX=-1.0
   46 XX(I)=AX
      CALL POLISH(M,MM,N,NN,NDAT,X,XX,Y,YY,IIROW,JJCOL,RTOT,
     1ROWWGT,CCWT,COLWGT,IADDR,IDAT)
      CALL XMAXMI(X,AXMAX,AXMIN,M)
      CRMID=0.5*(AXMIN+AXMAX)
      CRHALF=0.5*CRLONG*(AXMAX-AXMIN)
      CRMIN=CRMID-CRHALF
   47 CRMAX=CRMID+CRHALF
   48 IF(IC.EQ.1) GOTO 63
C WE NOW FIND OUT WHICH OF THE DERIVED GROUPS MORE CLOSELY RESEMBLES
C THE BROTHER OF THE GROUP THAT HAS JUST BEEN DIVIDED
C TO DO THIS WE REQUIRE SPECIES TOTALS (IGNORING PSEUDOSPECIES WEIGHTS)
      SMALL=1.0E-7
      DO 49 JJ=1,NN
      Y(JJ)=SMALL
   49 YY(JJ)=SMALL
      TOT0=0.0
      TOT1=0.0
      DO 52 I=1,M
      II=IIROW(I)
      ID=IADDR(II)
      IF(X(I).GE.CRMID) TOT1=TOT1+RRWT(II)
      IF(X(I).LE.CRMID) TOT0=TOT0+RRWT(II)
   51 J=IDAT(ID)
      ID=ID+1
      IF(J.EQ.-1) GOTO 52
      JJ=JJCOL(J)
      JJJ=IABS(INDPOT(JJ))
      IF(X(I).GE.CRMID) YY(JJJ)=YY(JJJ)+RRWT(II)*CCWT(JJ)
      IF(X(I).LE.CRMID) Y(JJJ)=Y(JJJ)+RRWT(II)*CCWT(JJ)
      GOTO 51
   52 CONTINUE
      IIC=IC+1
      IF(IIC.EQ.(IIC/2)*2)  IIC=IIC-2
      CALL CLOSER(YSCOR1,IIC,MM,N,NN,NDAT,IEND,TOT0,TOT1,
     1IDAT,ICLASS,IADDR,INDPOT,JJCOL,RRWT,CCWT,COLWGT,Y,YY)
      IIIIC=IC/2+1
      IF(IIIIC.EQ.(IIIIC/2)*2) IIIIC=IIIIC-2
      YSCOR2=0.0
      IF(IC.LE.3) GOTO 54
      CALL CLOSER(YSCOR2,IIIIC,MM,N,NN,NDAT,IEND,TOT0,TOT1,
     1IDAT,ICLASS,IADDR,INDPOT,JJCOL,RRWT,CCWT,COLWGT,Y,YY)
   54 IW=IC-(IC/4)*4
      W=0.5
      IF(IW.EQ.1.OR.IW.EQ.2) W=-W
      YSCORE=YSCOR1+W*YSCOR2
      IF(IIC.NE.(IIC/2)*2) GOTO 559
C IIC IS EVEN, SO WE WANT THE NEGATIVE GROUP TO BE CLOSER TO IT
      IF (YSCORE.LE.0) GOTO 63
      GOTO 60
C IN THIS CASE IIC IS ODD
  559 IF (YSCORE.GE.0) GOTO 63
   60 DO 61 I=1,M
   61 X(I)=-X(I)
      AX=AXMAX
      AXMAX=-AXMIN
      AXMIN=-AX
      AX=CRMAX
      CRMAX=-CRMIN
      CRMIN=-AX
      CRMID=-CRMID
   63 CALL ZONEUP(X,AXMIN,CRMIN,CRMAX,AXMAX,M,MZ,MZOUT,MZCRIT,IZONE)
      IF (MIND.EQ.0) GOTO 100
      CUT1=CRMID-0.5*(CRMAX-CRMIN)*FLOAT(MZIND)/(FLOAT(MZCRIT)+0.001)
      CUT2=2.0*CRMID-CUT1
      CALL INDSCO(M,MM,N,NN,NDAT,X,AXPOS,AXNEG,
     1CUT1,CUT2,Y,YY,IIROW,ROWWGT,IADDR,IDAT)
      DO 65 J=1,N
   65 Y(J)=YY(J)/AXPOS-Y(J)/AXNEG
C---Argument PRECIS added by P.Minchin June 1997
      CALL TOPIND(N,NN,MIND,FEEBLE,Y,INDORD,INDSIG,
     1JJCOL,INDPOT,JNAM,NIND,LIND,PRECIS)
      IF(NIND.EQ.0) GOTO 100
      MIS=10000
      CALL CODESC(M,MM,N,NN,NDAT,NIND,INDORD,IX,IY,IIROW,IADDR,IDAT)
C NOW FIND BEST NUMBER OF INDICATORS AND BEST DIVISION POINT
      DO 70 IIND=1,NIND
      CALL TABLE(M,MMZ,MMS,MIND,IX,IZONE,INDSIG,IPICT,0,IIND)
      CALL FIND(MMZ,MMS,IPICT,IZD1,ISD1,MISCL)
      IF(MISCL.GE.MIS)GOTO 70
      MIS=MISCL
      ISD=ISD1
      IZD=IZD1
      IIZD=IZD-MZIND
      IND=IIND
   70 CONTINUE
      CALL TABLE(M,MMZ,MMS,MIND,IX,IZONE,INDSIG,IPICT,1,IND)
c      WRITE (IUOUT2,1003)
c 1003 FORMAT('INDICATORS, TOGETHER WITH THEIR SIGN')
      IT=0
      DO 80 IIND=1,IND
      J=INDORD(IIND)
      JJ=JJCOL(J)
c     save indices of indicators (with sign)
      IF (INDSIG(IIND) .GT. 0) THEN
         INDICS(IIND, IC) = JJ
      ELSE
         INDICS(IIND, IC) = -JJ
      ENDIF
      IT=IT+1
c     ITEM(IT)=I2CHAR(JNAME1(JJ))
      ITEM(IT) = JNAME1(JJ)
      IT=IT+1
c     ITEM(IT)=I2CHAR(JNAME2(JJ))
      ITEM(IT) = JNAME2(JJ)
      IT=IT+1
C---Convert values of JNAM to CHARACTER: P.Minchin June 1997
c     ITEM(IT)=I2CHAR(JNAM(JJ))
      ITEM(IT) = JNAM(JJ)
c      ITEM(IT)=ITEM(IT)(4:4)
      IT=IT+1
      IF(INDSIG(IIND).GT.0) GOTO 75
c     Original function holds printable character string in ITEM and
c     adds either + or - to show the side of indicator species. We
c     currently work with integer indices, and make negative indicators
c     negative
      ITEM(IT)= -1 * ITEM(IT)
      GOTO 80
c     75 ITEM(IT)=PLUS
 75   CONTINUE
   80 CONTINUE
c      WRITE (IUOUT2,1004) (ITEM(IIT),IIT=1,IT)
C---Altered format: P.Minchin June 1997
C 1004 FORMAT((7(A4,1X,A4,I1,'(',A1,')',3X)))
c 1004 FORMAT((7(A4,1X,A4,A1,'(',A1,')',3X)))
      GOTO 120
  100 DO 110 I=1,M
  110 IX(I)=0
      MS=1
      IZD=(1+MZ)/2
      IIZD=IZD
      ISD=0
      IF(MIND.EQ.0) NIND=0
c      IF(MIND.NE.0) WRITE (IUOUT2,1006)
c 1006 FORMAT('(NO INDICATORS)')
  120 CALL REPORT(M,MM,N,NN,NDAT,MMZ,MMS,IC,X,XX,ROWWGT,
     1Y,YY,COLWGT,IX,IZONE,IIROW,IY,JJCOL,IADDR,IDAT,IPICT,
     2INAME1,INAME2,JNAME1,JNAME2,JNAM, ISD2, ISEC)
      LIMPOS(IC) = ISD2
      DO 130 I=1,M
      AX=X(I)
      IF(AX.GT.3.5) IX(I)=1
      IF(AX.LT.3.5) IX(I)=0
  130 CONTINUE
      GOTO 20
 500  CONTINUE
c  500 WRITE (IUOUT2,1010)
      IF(IPUNCH.EQ.0) GO TO 560
c     NUMBER=0
C     DO 510 II=1,MM
C     IF(ICLASS(II).NE.0) NUMBER=NUMBER+1
C 510 CONTINUE
C     WRITE(7,1020) NUMBER
C1020 FORMAT('THERE ARE',I6,'  ITEMS')
C     DO 550 II=1,MM
C     IC=ICLASS(II)
C     IF(IC.EQ.0) GOTO 550
C     CALL BIN(ITEM,IC)
C     IF(ISEC.EQ.2) GOTO 530
C     WRITE(7,1030) INFLAG(II),INAME1(II),INAME2(II),IC,
C    1(ITEM(IT),IT=1,30)
C     GOTO 550
C 530 WRITE(7,1040) INFLAG(II),INAME1(II),INAME2(II),IC,
C    1(ITEM(IT),IT=1,30)
C1030 FORMAT(I4,1X,2A4,I7,3X,30A1)
C1040 FORMAT(I4,1X,A4,1X,A4,I6,3X,30A1)
C 550 CONTINUE
  560 RETURN
      END
C
