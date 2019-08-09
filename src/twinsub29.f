      SUBROUTINE YXMULT(M,MM,N,NN,NDAT,Y,X,IIROW,IADDR,
     1ROWWGT,COLWGT,IDAT)
C FORMS MATRIX PRODUCT X=AY, WEIGHTED AS SUBROUTINE XYMULT.
      REAL Y(N),X(M),ROWWGT(M),COLWGT(N)
      INTEGER IIROW(M),IADDR(MM),IDAT(NDAT)
C---FOLLOWING COMMON BLOCK ADDED BY P.R.MINCHIN SEPT 1986
      COMMON /LUNITS/ IUIN1,IUIN2,IUOUT1,IUOUT2,IUOUT3
      DO 10 J=1,N
   10 Y(J)=Y(J)*COLWGT(J)
      DO 20 I=1,M
      II=IIROW(I)
      AX=0.0
      ID=IADDR(II)
   15 J=IDAT(ID)
      IF(J.EQ.-1) GOTO 16
      AX=AX+Y(J)
      ID=ID+1
      GOTO 15
   16 X(I)=AX*ROWWGT(I)
   20 CONTINUE
      DO 30 J=1,N
   30 Y(J)=Y(J)/COLWGT(J)
      RETURN
      END
C
