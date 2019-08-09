      SUBROUTINE ISORT(IX,N)
      INTEGER IX(N)
C EFFICIENT (HEAP-SORT) IN SITU SORTING OF VECTOR IX(N)
      DO 10 I=1,N
      J=I
      JX=IX(J)
    5 IF(J.EQ.1) GOTO 8
      JJ=J/2
      JJX=IX(JJ)
      IF(JJX.GE.JX) GOTO 8
      IX(J)=JJX
      J=JJ
      GOTO 5
    8 IX(J)=JX
   10 CONTINUE
      I=N
      GOTO 14
   12 IX(J)=JX
   14 IF(I.EQ.1) RETURN
      JX=IX(I)
      IX(I)=IX(1)
      I=I-1
      J=1
      JJ=2
   15 IF(I-JJ) 12,17,16
   16 JJX=IX(JJ)
      JJJ=JJ+1
      JJJX=IX(JJJ)
      IF(JJX.GE.JJJX) GOTO 18
      IF(JX.GE.JJJX) GOTO 12
      IX(J)=JJJX
      J=JJJ
      JJ=J*2
      GOTO 15
   17 JJX=IX(JJ)
   18 IF(JX.GE.JJX) GOTO 12
      IX(J)=JJX
      J=JJ
      JJ=J*2
      GOTO 15
      END
C
