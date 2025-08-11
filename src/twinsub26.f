      SUBROUTINE FIND(MMZ,MMS,IPICT,IZD1,ISD1,MISCL)
C     INPUT IS MATRIX IPICT, RELATING INDICATOR SCORES TO ZONES.
C     SUBROUTINE FINDS CLASSIFICATION CENTER IZD1,ISD1 THAT MINIMIZES
C     NUMBER OF MISCLASSIFICATIONS MISCL.
      INTEGER IPICT(MMZ,MMS)
      COMMON/PICT/MZ,MS,IZD,IIZD,ISD,ISHIFT,MZIND,MZCRIT,MZOUT,NIND
      MINZ=MZOUT+MZIND
      MAXZ=MZOUT+MZCRIT
C     WE FORM THE TOTALS IN IPICT OF THE ELEMENTS THAT LIE BELOW AND TO
C     THE LEFT OF A GIVEN ONE. I.E. THE NEW VALUE OF IPICT(IZ,IS) IS THE
C     SUM OF ALL THE OLD VALUES IF IPICT(IIZ,IIS) FOR WHICH IIZ IS LESS
C     THAN OR EQUAL TO IZ AND IIS IS LESS THAN OR EQUAL TO IS.
      DO 10 IS=2,MS
         IPICT(1,IS)=IPICT(1,IS-1)+IPICT(1,IS)
 10   CONTINUE
      DO 20 IZ=2,MZ
         IPICT(IZ,1)=IPICT(IZ-1,1)+IPICT(IZ,1)
 20   CONTINUE
      DO 31 IS=2,MS
         DO 30 IZ=2,MZ
            IIS=IS-1
            IIZ=IZ-1
            IPICT(IZ,IS)=IPICT(IIZ,IS)+IPICT(IZ,IIS)-
     1           IPICT(IIZ,IIS)+IPICT(IZ,IS)
 30      CONTINUE
 31   CONTINUE
c     miscl=10000 guarantees that the first branch (35) will be taken on the
c     first pass in the arithmetic if(miss-miscl) and this sets iia, iib &
c     cc. However, compiler warns that these may be used uninitialized. Set
c     them any value here: they will be set later
      MISCL=10000
      IIA = 0
      IIB = 0
      CC = 0.0
      DO 41 IZ=MINZ,MAXZ
         DO 40 IS=1,MS
            IIZ=IZ-MZIND
            MISS=IPICT(IIZ,MS)-IPICT(IIZ,IS)+IPICT(MZ,IS)-IPICT(IZ,IS)
            A=FLOAT(IPICT(IIZ,MS))
            B=FLOAT(IPICT(MZ,MS)-IPICT(IZ,MS))
            C=ABS((A-B)/(A+B))
            IA=IABS(IS-ISHIFT)
            IB=IABS(1+MZ-IZ-IIZ)
c     See comment above about iia, iib & cc
            IF (MISS .LT. MISCL) THEN
               GOTO 35
            ELSEIF (MISS .EQ. MISCL) THEN
               GOTO 36
            ELSE
               GOTO 40
            ENDIF
 35         MISCL=MISS
            IZD1=IZ
            ISD1=IS
            IIA=IA
            IIB=IB
            CC=C
            GOTO 40
 36         CONTINUE
            IF(C .LT. CC) THEN
               GOTO 35
            ELSEIF (C .EQ. CC) THEN
               GOTO 38
            ELSE
               GOTO 40
            ENDIF
 38         IF(IA.LT.IIA) GOTO 35
            IF(IB.LT.IIB) GOTO 35
 40      CONTINUE
 41   CONTINUE
      RETURN
      END
C
