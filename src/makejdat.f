      subroutine makejdat(mm, nn, nspec, ndat, nmax, iaddr, idat,
     . indpot ,iclass, y, ccwt, rrwt, jdat)
c After classifying sampling units (quadrats), twinspan makes
c corresponding classification for species. This is not based on the
c original data, but data on species fidelity. This piece of code builds
c up the data vector JDAT that CLASS uses to build species
c classification.
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      integer mm, nn, nspec, ndat, nmax
      integer iaddr(mm), indpot(nn), iclass(mm)
      integer idat(ndat), jdat(ndat/2)
      real(8) y(nn), ccwt(nn), rrwt(mm), tot(511), totj(511)
c     Arbitrary constants
      SPE1 = 0.8
      SPE2 = 2.0
      SPE3 = 6.0
      WTHIGH = 2.0
c Initialize MMAX
      MMAX = MAX(mm, nn)
c ID is the first observation in the last SU, and the 415 loop
c advances to the end of the data
      ID=IADDR(MM)
  415 IF(IDAT(ID).EQ.-1) GOTO 420
      ID=ID+1
      GOTO 415
c     NDAT must be twice the size of real length: used as workspace
  420 NDAT=NDAT/2
c     IF (ID.GT.NDAT) GOTO 999
      if(id.gt.ndat) call rexit("1: id > ndat")
      IID=NDAT
      DO 422 J=1,NSPEC
 422     Y(J)=0.0
c     425 is a while loop with GOTO
  425 IDAT(IID)=IDAT(ID)
      IID=IID-1
      ID=ID-1
      IF(ID.EQ.0) GOTO 426
      GOTO 425
c     end 425 loop
  426 JJJ=0
      ICMAX=0
      DO 431 II=1,MM
         IC=ICLASS(II)
c TOT is hardcoded do length 511, but the following can give ICMAX that
c is higher, and then segfaults in loop 432. Here we go to upper level
c of hierarchy so long that all class values IC are below 512 which
c corresponds 8 to levels of division. Original 427 commented out.

c     427     IF(IC*3.LE.NMAX) GOTO 428
 427     if (ic .le. 511) goto 428
c ic/2 is a mother class at higher level
         IC=IC/2
         GOTO 427
 428     ICLASS(II)=IC
 431     IF(IC.GT.ICMAX) ICMAX=IC
c     IF (NSPEC.GT.MMAX) GOTO 999
         if(nspec.gt.mmax) call rexit("2 nspec > mmax")
      DO 432 IC=1,ICMAX
 432     TOT(IC)=0.0
      JD=0
      IBIG=10000
c 440 loop contains loop 433: you get out of this when you hit the
c quadrat endmarker -1 and GOTO 437
      DO 440 II=1,MM
         IC=ICLASS(II)
 433     IID=IID+1
         JJ=IDAT(IID)
         IF(JJ.EQ.-1) GOTO 437
         J=IABS(INDPOT(JJ))
         JD=JD+1
c     IF(JD.GT.NDAT) GOTO 999
         if(jd.gt.ndat) call rexit("3 jd > ndat")
         JDAT(JD)=J*IBIG+IC
         TOT(IC)=TOT(IC)+1.0
         IF(J.EQ.JJJ) GOTO 435
         JJJ=J
         Y(J)=Y(J)+CCWT(J)*RRWT(II)
         ID=ID+1
         IDAT(ID)=J
         ID=ID+1
         IDAT(ID)=1
c     IF(ID.GE.IID) GOTO 999
         if(id.ge.iid) call rexit("4 id > iid")
         GOTO 433
 435     IDAT(ID)=IDAT(ID)+1
         GOTO 433
 437     ID=ID+1
         IDAT(ID)=JJ
         JJJ=0
 440  CONTINUE
      DO 441 J=1,NSPEC
 441     RRWT(J)=Y(J)+0.00001
C WE NOW WORK ON SPECIES CLASSIFICATION
      IIBIG=IBIG-1
      DO 442 J=1,NSPEC
         JD=JD+1
c     IF(JD.GT.NDAT) GOTO 999
         if(jd.gt.ndat) call rexit("5 jd > ndat")
  442 JDAT(JD)=J*IBIG+IIBIG
      CALL ISORT(JDAT,JD)
C MOVE MATRIX UP TO BACK OF ARRAY
      JJD=NDAT
c 445 is a loop
  445 JDAT(JJD)=JDAT(JD)
        JJD=JJD-1
        JD=JD-1
        IF(JD.EQ.0) GOTO 450
      GOTO 445
C WE NOW RECONSTITUTE THE MATRIX AND PROCEED TO ACCUMULATE
C THE TOTALS FOR THE SAMPLE CATEGORIES
 450  IC=ICMAX
 455  ICC=IC/2
        TOT(ICC)=TOT(ICC)+TOT(IC)
        IC=IC-1
      IF(IC.GT.1) GOTO 455
      DO 470 JJJ=1,NSPEC
         DO 456 IC=1,ICMAX
 456        TOTJ(IC)=0.0
 458        JJD=JJD+1
              JJ=JDAT(JJD)
              J=JJ/IBIG
              IC=JJ-J*IBIG
              IF(IC.EQ.IIBIG) GOTO 460
              TOTJ(IC)=TOTJ(IC)+1.0
            GOTO 458
 460        IC=ICMAX
 462        IIC=IC/2
              TOTJ(IIC)=TOTJ(IIC)+TOTJ(IC)
              IC=IC-1
            IF(IC.GT.1) GOTO 462
            DO 464 IC=1,ICMAX
               RAT=TOTJ(IC)*(TOT(1)-TOT(IC))/(TOT(IC)*(TOTJ(1)-
     1              TOTJ(IC))+1.0E-7)
               IF (RAT.LT.SPE1) GOTO 464
               JD=JD+1
               JDAT(JD)=IC*3-2
               IF (RAT.LT.SPE2) GOTO 464
               JD=JD+1
               JDAT(JD)=IC*3-1
               IF (RAT.LT.SPE3) GOTO 464
               JD=JD+1
               JDAT(JD)=IC*3
 464        CONTINUE
            JD=JD+1
c     IF(JD.GE.JJD) GOTO 999
            if(jd.ge.jjd) call rexit("6 jd > jjd")
            JDAT(JD)=-1
 470     CONTINUE
C JDAT NOW CONTAINS DESIRED SPECIES INFORMATION.  IT IS NOW A
C MATTER OF FIXING UP WEIGHTS AND CLASSIFYING SPECIES
      DO 480 IC=1,ICMAX
         WEIGHT=1.0
         IIC=IC
 473     IIC=IIC*2
           IF(IIC.GT.1023) GOTO 476
              WEIGHT=WEIGHT*1.414
         GOTO 473
 476     WEIGHT=WEIGHT*TOT(IC)
         IIC=3*IC
         CCWT(IIC-2)=WEIGHT
         INDPOT(IIC-2)=IIC-2
         CCWT(IIC-1)=WEIGHT*WTHIGH
         INDPOT(IIC-1)=IIC-1
         CCWT(IIC)=WEIGHT*WTHIGH
         INDPOT(IIC)=IIC
 480  CONTINUE
      return
      end
