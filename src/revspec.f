      SUBROUTINE revspec(NSPEC, MM, NDAT, IX, IY, X, Y,IDAT, IREV)
C WE HAVE NOW ESSENTIALLY GOT SAMPLES AND SPECIES ORDERED, BUT
C THE SPECIES ORDER MAY NEED TO BE REVERSED TO GET THINGS ON THE
C POSITIVE DIAGONAL.  ALSO WE ACCUMULATE SPECIES TOTALS IN INDPOT.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER IX(MM)
      INTEGER IY(NSPEC)
      INTEGER IDAT(NDAT)
      REAL(8) X(MM), Y(NSPEC)
      DO 500 II=1,MM
         I=IX(II)
 500  X(I)=DBLE(II)
      DO 510 JJ=1,NSPEC
         J=IY(JJ)
 510  Y(J)=DBLE(JJ)
      TOT1=0.0
      TOTX=0.0
      TOTY=0.0
      TOTXY=0.0
      ID=0
      DO 520 II=1,MM
         AX=X(II)
 512     ID=ID+1
         J=IDAT(ID)
         IF(J.EQ.-1) GOTO 520
         AY=Y(J)
         ID=ID+1
         JQ=IDAT(ID)
         Q=DBLE(JQ)
         TOT1=TOT1+Q
         TOTX=TOTX+Q*AX
         TOTY=TOTY+Q*AY
         TOTXY=TOTXY+Q*AX*AY
         GOTO 512
 520  CONTINUE
      TOTXY=TOTXY-TOTX*TOTY/TOT1
      IREV=1
      IF(TOTXY.GT.0.0) IREV=0
C IN THIS CASE TOTXY IS NEGATIVE, I.E. IX AND IY ARE ORDERS
C THAT ARE NEGATIVELY CORRELATED, AND IY MUST BE REVERSED
C AND THEN IREV=1
      RETURN
      END
