      SUBROUTINE REPORT(M,MM,N,NN,NDAT,MMZ,MMS,IC,X,XX,ROWWGT,
     1Y,YY,COLWGT,IX,IZONE,IIROW,IY,JJCOL,IADDR,IDAT,IPICT,
     2INAME1,INAME2,JNAME1,JNAME2,JNAM)
C REPORTS ON DICHOTOMY.  TRUE INDICATOR SCORES ARE NOW IN IX.
C MATRIX IPICT CONTAINS DIAGRAM OF DICHOTOMY, RELATING ZONES
C TO INDICATOR SCORES.
      REAL X(M),XX(M),ROWWGT(M),Y(N),YY(N),COLWGT(N)
C---Changed ITEM to CHARACTER: P.Minchin June 1997
      CHARACTER*4 ITEM(200)
      INTEGER JNAM(NN)
C---Changed INAME, JNAME to CHARACTER: P.Minchin June 1997
      CHARACTER*4 INAME1(MM),INAME2(MM),JNAME1(NN),JNAME2(NN)
      INTEGER IX(M),IZONE(M),IIROW(M),IY(N),JJCOL(N),
     1IADDR(MM),IDAT(NDAT),IPICT(MMZ,MMS),NCAT(6)
C---FOLLOWING COMMON BLOCK ADDED BY P.R.MINCHIN SEPT 1986
      COMMON /LUNITS/ IUIN1,IUIN2,IUOUT1,IUOUT2,IUOUT3
      COMMON/PICT/MZ,MS,IZD,IIZD,ISD,ISHIFT,MZIND,MZCRIT,MZOUT,NIND
C---Variable PRECIS added P.Minchin  June 1997
      COMMON/LIMS/RARE,FEEBLE,FRQLIM,TOL,RATLIM,REPLIM,PRECIS
      COMMON/SWITCH/IFAIL,IDIAGR,ISTOP,IREWT
      COMMON/IARBS/ICWEXP,IEND,MMIN,IPREXP,LEVMAX
      COMMON/WORK/ITEM
      COMMON/SECOND/ISEC
      DO 10 ICAT=1,6
   10 NCAT(ICAT)=0
      DO 20 I=1,M
      IZ=IZONE(I)
      IS=IX(I)
      IF(IZ.GT.IIZD) GOTO 13
      IF(IS.GT.ISD) GOTO 12
      ICAT=1
      GOTO 18
   12 ICAT=3
      GOTO 18
   13 IF(IZ.GT.IZD) GOTO 16
      IF(IS.GT.ISD) GOTO 15
      ICAT=2
      GOTO 18
   15 ICAT=5
      GOTO 18
   16 IF(IS.GT.ISD) GOTO 17
      ICAT=6
      GOTO 18
   17 ICAT=4
   18 X(I)=FLOAT(ICAT)
      NCAT(ICAT)=NCAT(ICAT)+1
   20 CONTINUE
      NNEG=NCAT(1)+NCAT(2)+NCAT(3)
      NPOS=NCAT(4)+NCAT(5)+NCAT(6)
      CALL INDSCO(M,MM,N,NN,NDAT,X,AXPOS,AXNEG,3.5,3.5,
     1Y,YY,IIROW,ROWWGT,IADDR,IDAT)
      PRLIM=(REPLIM-1.0)/(REPLIM+1.0)-0.0001
      DO 30 J=1,N
      AY=Y(J)/AXNEG
      AYY=YY(J)/AXPOS
      IIY=0
C--Comparisons changed to use PRECIS: P.Minchin June 1997
C      IF(AY.GE.RARE) IIY=1
C      IF(AYY.GE.RARE) IIY=1
      IF((AY-RARE).GE.PRECIS) IIY=1
      IF((AYY-RARE).GE.PRECIS) IIY=1
      PREF=(AYY-AY)/(AYY+AY)
C--Comparisons changed to use PRECIS: P.Minchin June 1997
C      IF(PREF.GT.PRLIM) GOTO 25
C      IF(PREF.LT.-PRLIM) GOTO 23
      IF((PREF-PRLIM).GT.PRECIS) GOTO 25
      IF((-PREF-PRLIM).GT.PRECIS) GOTO 23
      IY(J)=IIY
      GOTO 30
   23 IY(J)=2*IIY
      GOTO 30
   25 IY(J)=3*IIY
   30 CONTINUE
C WE NOW HAVE A VARIABLE IY THAT IS 0 FOR RARITIES, 3 FOR POSITIVE
C PREFERENTIALS, 2 FOR NEGATIVE PREFERENTIALSS, 1 FOR INDIFFERENTS.
      IF(NIND.EQ.0) GOTO 60
      ISD1=ISD-ISHIFT
      ISD2=ISD1+1
      WRITE (IUOUT2,1050) ISD1,ISD2
 1050 FORMAT('MAXIMUM INDICATOR SCORE FOR NEGATIVE GROUP ',I4,
     1'     MINIMUM INDICATOR SCORE FOR POSITIVE GROUP ',I4)
      IF(IDIAGR.EQ.1) CALL DIAGR(IPICT,MMZ,MMS)
   60 IIC=IC*2
      CALL BIN(ITEM,IIC)
      WRITE (IUOUT2,1002) IIC,NNEG,(ITEM(IT),IT=1,30)
 1002 FORMAT(/'ITEMS IN NEGATIVE GROUP',I4,'  (N=',
     1I5,')',9X,'I.E. GROUP ',30A1)
      CALL XLIST(M,MM,X,IIROW,1,2,3,INAME1,INAME2)
      IF(NIND.EQ.0) GOTO 70
      NOUT=NCAT(2)
      IF(NOUT.EQ.0) GOTO 65
      WRITE (IUOUT2,1003) NOUT
 1003 FORMAT(/'BORDERLINE NEGATIVES    (N=',I5,')')
      CALL XLIST(M,MM,X,IIROW,2,2,2,INAME1,INAME2)
   65 NOUT=NCAT(3)
      IF(NOUT.EQ.0) GOTO 70
      WRITE (IUOUT2,1004) NOUT
 1004 FORMAT(/'MISCLASSIFIED NEGATIVES (N=',I5,')')
      CALL XLIST(M,MM,X,IIROW,3,3,3,INAME1,INAME2)
   70 IIC=IC*2+1
      CALL BIN(ITEM,IIC)
      WRITE (IUOUT2,1005) IIC,NPOS,(ITEM(IT),IT=1,30)
 1005 FORMAT(/'ITEMS IN POSITIVE GROUP',I4,'  (N=',
     1I5,')',9X,'I.E. GROUP ',30A1)
      CALL XLIST(M,MM,X,IIROW,4,5,6,INAME1,INAME2)
      IF(NIND.EQ.0) GOTO 80
      NOUT=NCAT(5)
      IF(NOUT.EQ.0) GOTO 75
      WRITE (IUOUT2,1006) NOUT
 1006 FORMAT(/'BORDERLINE POSITIVES    (N=',I5,')')
      CALL XLIST(M,MM,X,IIROW,5,5,5,INAME1,INAME2)
   75 NOUT=NCAT(6)
      IF(NOUT.EQ.0) GOTO 80
      WRITE (IUOUT2,1007) NOUT
 1007 FORMAT(/'MISCLASSIFIED POSITIVES (N=',I5,')')
      CALL XLIST(M,MM,X,IIROW,6,6,6,INAME1,INAME2)
C WE NOW LIST OUT THE COMMONER SPECIES ACCORDING AS THEY ARE
C PREFERENTIAL OR NEUTRAL
   80 IF(ISEC.EQ.2) GOTO 100
      CALL INDSCO(M,MM,N,NN,NDAT,X,AXPOS,AXNEG,3.5,3.5,
     1Y,YY,IIROW,ROWWGT,IADDR,IDAT)
      WRITE (IUOUT2,1010)
 1010 FORMAT(/'NEGATIVE PREFERENTIALS')
      CALL YLIST(N,NN,Y,YY,IY,JJCOL,2,JNAME1,JNAME2,JNAM)
      WRITE (IUOUT2,1011)
 1011 FORMAT(/'POSITIVE PREFERENTIALS')
      CALL YLIST(N,NN,Y,YY,IY,JJCOL,3,JNAME1,JNAME2,JNAM)
      WRITE (IUOUT2,1012)
 1012 FORMAT(/'NON-PREFERENTIALS')
      CALL YLIST(N,NN,Y,YY,IY,JJCOL,1,JNAME1,JNAME2,JNAM)
  100 RETURN
      END
C
