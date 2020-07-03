*#RUN * :TIME=3 
                    PROGRAM EXA1
C**************************************************************
C*                                                            *
C*      ONE DIMENSIONAL UNSTEADY HEAT CONDUCTION PROBLEM      *
C*                                                            *
C*          BY M. KATSUKI              APRIL 20, 1989         *
C**************************************************************
C
C========================= SYMBOL LIST ========================
C--------   TEMPERATURE
C   T(I)      :TEMPERATURE  (K)
C   TO(I)     :TEMPERATURE AT PREVIOUS TIME STEP  (K)
C   TH        :TEMPERATURE OF HEATED SIDE  (K)
C   TL        :ATMOSPHERIC TEMPERATURE (K) 
C--------   GEOMERTY
C   DG(I)     :DISTANCE BETWEEN GRIDS  (M)
C   DX(I)     :WIDTH OF A CONTROL VOLUME  (M)
C   X(I)      :DISTANCE FROM THE HEATED SURFACE TO A GRID  (M)
C   NI        :MAXIMUM NUMBER OF GRIDS OR CONTROL VOLUMES
C   NIM1      :NI-1
C--------   PROPERTIES
C   ALMD      :THERMAL CONDUCTIVITY  (W/(M*K))
C   CP        :SPECIFIC HEAT  (J/(KG*K))
C   RHO       :DENSITY  (KG/M**3)
C--------   TIME
C   DELT      :TIME STEP  (S)
C   TIME      :TIME   (S)
C   TIMEND    :TIME LIMIT OF CALCULATION  (S)
C--------   TDMA
C   AP(I)     :COEFFICIENT FOR POLE GRID VALUE
C   AE(I)     :COEFFICIENT FOR EAST GRID VALUE
C   AW(I)     :COEFFICIENT FOR WEST GRID VALUE
C   AO(I)     :COEFFICIENT FOR POLE GRID VALUE AT PREVIOUS TIME STEP
C   B(I)      :UNSTEADY TERM
C   P(I),Q(I) :COEFFICIENTS OF TDMA SUBSTITUTION
C--------   PRINT
C   IE        :I-INDEX FOR PRINT END LOCATION
C   IVAL      :NUMBER OF DATA PRINTED IN A LINE
C   ISTART    :START I-INDEX IN A LINE
C   IEND      :END I-INDEX IN A LINE
C==============================================================
C
C--------------------------------------------------------------
CHAPTER   0   0   0   0   0   PRELIMINARIES   0   0   0   0   0
C
      DIMENSION DG(20),DX(20),T(20),TO(20),AP(20),AE(20),
     &  AW(20),B(20),AO(20),P(20),Q(20),X(20)
C
C*********    CONTROL VOLUME & GEOMETRY
C
      DATA NI/11/
      DATA (DG(I),I=1,10)/10*0.020/
      DATA (DX(I),I=1,11)/0.0100, 9*0.020, 0.0100/
C
C*********    PROPERTIES
C
      DATA ALAMD,RHO,CP/ 237.0, 2688.0, 905.0/
C
C*********    PARAMETERS
C
      DATA TH,TL/500.0, 300.0/
      DATA DELT,TIMEND/1.00, 20.0/
C
C---------------------------------------------------------------
CHAPTER   1   1   1   1   1   INITIALIZATION   1   1   1   1   1
      NIM1=NI-1
      X(1)=0.0
      DO 10 I=2,NI
        X(I)=X(I-1)+DG(I-1)
   10 CONTINUE
C
      DO 20 I=1,NI
        TO(I)=TL
   20 CONTINUE
C
      TIME=0.0
C
C---------------------------------------------------------------
CHAPTER   2   2   2   ITERATION LOOP FOR A TIME STEP   2   2   2
C
  999 CONTINUE
C
      TIME=TIME+DELT
C
C---------------------------------------------------------------
CHAPTER   3   3   3   3   ASSEMBLE  COEFFICIENTS   3   3   3   3
C
      DO 30 I=2,NIM1
        AW(I)=ALAMD/DG(I-1)
        AE(I)=ALAMD/DG(I)
        AO(I)=(DX(I)*RHO*CP)/DELT
        B(I) =AO(I)*TO(I)
        AP(I)=AE(I)+AW(I)+AO(I)
   30 CONTINUE
C
C---------------------------------------------------------------
CHAPTER   4   4   4   4    BOUNDARY  CONDITIONS    4   4   4   4
C
C*********     I=1
      AW(1)=0.0
      AE(1)=0.0
      AO(1)=1.0
      TO(1)=TH
      B(1)=AO(1)*TO(1)
      AP(1)=1.0
C
C*********     I=NI      
      AW(NI)=0.0
      AE(NI)=0.0
      AO(NI)=1.0
      B(NI) =AO(NI)*TL
      AP(NI)=1.0
C
C--------------------------------------------------------------
CHAPTER   5   5   5    SOLVE  EQUATIONS  BY  TDMA     5   5   5
C
      P(1)=AE(1)/AP(1)
      Q(1)=B(1)/AP(1)
C
      DO 40 I=2,NI
        P(I)=AE(I)/(AP(I)-AW(I)*P(I-1))
        Q(I)=(B(I)+AW(I)*Q(I-1))/(AP(I)-AW(I)*P(I-1))
   40 CONTINUE
C
      T(NI)=Q(NI)
C
      DO 50 I=1,NIM1
        J=NI-I
        T(J)=P(J)*T(J+1)+Q(J)
   50 CONTINUE
C
C--------------------------------------------------------------
CHAPTER   6   6    PREPARATION FOR THE NEXT TIME STEP     6   6
C
      DO 60 I=1,NI
        TO(I)=T(I)
   60 CONTINUE
C
      IF(TIME.LT.TIMEND) GO TO 999
C
C********     END OF ITERATION LOOP FOR A TIME STEP
C--------------------------------------------------------------
CHAPTER   7   7   7   PRINT  OUT  AND  TERMINATION   7   7   7 
C
      WRITE(6,900) 
      WRITE(6,910) TH
      WRITE(6,920) TL
      WRITE(6,930) DELT
      WRITE(6,940) TIMEND
C
      ISTART=1
      IE=NI
      IVAL=11
C
   70 IEND=ISTART+IVAL-1
      IF(IEND.GT.IE) IEND=IE
      WRITE(6,950) (I,I=ISTART,IEND)
      WRITE(6,960) (X(I),I=ISTART,IEND)
      WRITE(6,970) (T(I),I=ISTART,IEND)
      IF(IEND.EQ.IE) GO TO 80
      ISTART=ISTART+IVAL
      GO TO 70
C
   80 CONTINUE
C
C********     FORMAT STATEMENTS
C
  900 FORMAT(//5X,'ONE DIMENSIONAL UNSTEADY HEAT CONDUCTIVITY')
  910 FORMAT(//15X,'HIGH TEMPERATURE SIDE    ',1F6.1,2X,'K')
  920 FORMAT(/15X,'LOW TEMPERATURE SIDE     ',1F6.1,2X,'K')
  930 FORMAT(//15X,'TIME STEP                ',1F6.3,2X,'SEC')
  940 FORMAT(//5X,'TEMPERATURE DISTRIBUTION AT TIME = ',
     &  1F6.1,2X,'SEC') 
  950 FORMAT(/7X,'I =',I6,10I9)
  960 FORMAT(7X,'R (M)=',11F9.3)
  970 FORMAT(/7X,'T (K)=',1P11E9.2)
C
      STOP
C
      END
