C                     PROGRAM EXA2
C******************************************************************
C*                                                                *
C*        ONE-DIMENSIONAL STEADY HEAT CONDUCTION PROBLEM          *
C*            WITH HEAT TRANSFER AND HEAT GENERATION              *
C*                                                                *
C*                                                                *
C*            BY M. KATSUKI            APRIL 20, 1989             *
C******************************************************************
C
C============================ SYMBOL LIST =========================
C------------ SYMBOLS IN COMMON
C----- TEMP (TEMPERATURE)
C   T(I)       :TEMPERATURE  (K)
C   ALAMD(I)   :THERMAL CONDUCTIVITY  (W/(M*K))
C   TENV       :TEMPERATURE OF ENVIRONMENT  (K)
C   QS         :HEAT GENERATION RATE OF INNER MATERIAL  (W/M**3)
C   RESORT     :SUM OF RESIDUAL SOURCES FOR T-VARIABLES
C   TSF        :SURFACE TEMPERATURE  (K)
C----- GEOM (GEOMETRY)
C   R(I)       :RADIUS FROM SYMMETRY AXIS TO A GRID  (M)
C   RW(I)      :RADIUS TO THE WEST-WALL OF A CONTROL VOLUME  (M)
C   SEW(I)     :WIDTH OF A CONTROL VOLUME  (M)
C   SEP(I)     :XW(I+1)-X(I)  (M)
C   SPW(I)     :X(I)-XW(I)  (M)
C   NI         :MAXIMUM NUMBER OF GRIDS OR CONTROL VOLUMES 
C   NI1        :MAXIMUM I FOR INNER MATERIAL
C   DR         :WIDTH OF A CONTROL VOLUME  (M)
C   R1         :OUTER RADIUS OF THE INNER MATERIAL  (M)
C   R2         :OUTER RADIUS OF THE OUTER MATERIAL  (M)
C----- PRP (PROPERTIES)
C   C1IN,C2IN  :CONSTANTS FOR ALAMD(I) OF INNER MATERIAL
C   C1OUT,C2OUT:CONSTANTS FOR ALAMD(I) OF OUTER MATERIAL
C----- BND (BOUNDARY) 
C   ALMD       :THERMAL CONDUCTIVITY OF AIR  (W/(M*K))
C   U          :VELOCITY OF ENVIRONMENTAL AIR STREAM  (M/S)
C   AK1,AK2    :CONSTANTS FOR NUSSELT-REYNOLDS RELATION
C   AL1,AL2    :CONSTANTS FOR ALMD OF AIR
C   AN1,AN2    :CONSTANTS FOR ANYU OF AIR
C----- DSIZE (DIMENSION SIZE)
C   IT         :I-INDEX MAXIMUM DIMENSION OF DEPENDENT VARIABLES
C----- TDMA (T.D.M.A. COEFFICIENTS)
C   AP(I)      :COEFFICIENT FOR POLE GRID VALUE
C   AE(I)      :COEFFICIENT FOR EAST GRID VALUE
C   AW(I)      :COEFFICIENT FOR WEST GRID VALUE
C   SU(I),SP(I):COEFFICIENTS FOR THE LINEARIZED SOURCE TERM
C-------------- GENERAL SYMBOLS IN MAIN- & SUB-ROUTINES
C   I          :INDEX FOR FOR DEPENDENT VARIABLES AND COORDINATES
C   VOL        :VOLUME OF A CONTROL VOLUME  (M**3)
C   AREAE      :AREA OF THE EAST-WALL OF A CONTROL VOLUME  (M**2)
C   AREAW      :AREA OF THE WEST-WALL OF A CONTROL VOLUME  (M**2)
C   RESOR      :RESIDUAL SOURCE FOR THE INDIVIDUAL CONTROL VOLUME
C
C------------ SYMBOLS IN MAIN ROUTINE
C   ITERT      :ITERATION COUNTER
C   MAXIT      :MAXIMUM ITERATION FOR A RUN
C   NITERT     :TOTAL ITERATION NUMBER INCLUDING PREVIOUS RUNS
C   IPRINT     :PRINT INTERVAL OF ITERT FOR MONITORING VALUES
C   JPRINT     :PRINT INTERVAL OF ITERT FOR FIELD VALUES
C   IPRT       :ITERATION COUNTER TO PRINT MONITORING VALUES
C   JPRT       :ITERATION COUNTER TO PRINT FIELD VALUES
C   INITAL     :LOGICAL PARAMETER TO PRINT INITIAL FIELD VALUES
C   IREAD      :LOGICAL PARAMETER TO READ INITIAL FIELD VALUES 
C               FROM THE FILE
C   IWRITE     :LOGICAL PARAMETER TO WRITE FINAL FIELD VALUES 
C               TO THE FILE
C   SORMAX     :MAXIMUM ACCEPTABLE RESIDUAL SORCE FOR CONVERGED
C               SOLUTION
C   IMON       :I-INDEX AT THE MONITORING LOCATION
C   TMONIN     :REFERENCE VALUE FOR RESORT   (=QS*X1*X1*0.5)
C
C------------ SYMBOLS IN SUBROUTINE SOLVT
C   ALDE       :OVERALL HEAT CONDUCTANCE BETWEEN P AND E  (W/(M**2*K))
C   ALDW       :OVERALL HEAT CONDUCTANCE BETWEEN W AND P  (W/(M**2*K))
C
C------------ SYMBOLS IN SUBROUTINE BOUND
C   NIM1       :NI-1
C   TFILM      :FILM TEMPERATURE  (K)
C   REY        :REYNOLDS NUMBER
C   ANYU       :KINEMATIC VISCOSITY OF AIR  (M**2/s)
C   ANUS       :NUSSELT NUMBER
C   ALPHA      :HEAT TRANSFER COEFFICIENT ON SURFACE  (W/(M**2*K))
C   DM1,DM2    :DUMMY VALUES FOR OVERALL HEAT TRANSFER COEFFICIENT
C   ALDEN      :OVERALL HEAT TRANSFER COEFFICIENT THROUGH BOUNDARY
C   TMULT      :COEFFICIENT OF SOURCE TERM IN NIM1 CONTROL VOLUME
C
C------------ SYMBOLS IN SUBROUTINE TDMA
C   PHI(I)     :GENERAL VARIABLE SOLVED BY T.D.M.A.
C   ISTART     :T.D.M.A START I-INDEX
C   IEND       :T.D.M.A END I-INDEX
C   A(I),B(I),C(I),D(I),E(I)
C              :ELEMENT ARRAYS OF TRI-DIAGONAL MATRIX (=AP,AE,AW,SU,SP)
C   P(I),Q(I)  :COEFFICIENTS FOR T.D.M.A. SUBSTITUTION
C   DUMMY      :DUMMY VARIABLE
C
C------------ SYMBOLS IN SUBROUTINE PRINTS
C   IP         :PARAMETER TO SELECT THE VARIABLE TO BE PRINTED
C   AL(I)      :ARRAY TO PRINT
C   RA(I)      :POSITION ARRAY
C   IS         :I-INDEX OF PRINT START LOCATION
C   IE         :I-INDEX OF PRINT END LOCATION
C   IVAL       :NUMBER OF DATA PRINTED IN A LINE
C   ISTART     :START I-INDEX IN A LINE
C   IEND       :END I-INDEX IN A LINE
C==================================================================
C
CHAPTER   0   0   0   PRELIMINARIES   0   0   0   0   0
C
      COMMON
     &  /TEMP/T(20),ALAMD(20),TENV,QS,RESORT,TSF
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
     &  /PRP/C1IN,C2IN,C1OUT,C2OUT
     &  /BND/U,ALMD,AK1,AK2,AL1,AL2,AN1,AN2
     &  /DSIZE/IT
C
C--------     SET PARAMETERS
C--------     MAIN ROUTINE CONTROLLER
C
      DATA IPRINT,JPRINT,MAXIT/2,10,50/
      DATA INITAL/0/
      DATA IREAD,IWRITE/0,0/
      DATA SORMAX/1.0E-4/
      DATA IMON/5/
C
C--------     GEOMETRY
C
      NI=20
      NI1=8
      DR=0.01
      IT=NI
C
C--------     PROPERTIES
C
      TENV=500.0
      QS=1.0E6
      C1IN=2.84
      C2IN=0.0105
      C1OUT=102.8
      C2OUT=-7.07E-2
C
C--------     BOUNDARY PROPERTIES
C
      U=10.0
      AK1=0.174
      AK2=0.618
      AL1=7.24E-3
      AL2=6.36E-5
      AN1=-2.692E-6
      AN2=1.368E-8
C
C---------    SET INITIAL PARAMETERS
C
      ITERT=0
      NITERT=ITERT+MAXIT
      IPRT=ITERT+IPRINT
      JPRT=ITERT+JPRINT
C
C-------------------------------------------------------------
CHAPTER   1   1   1   INITIALIZATION   1   1   1   1   1
C
      CALL CONVOL
C
      CALL INITIA
C
      TMONIN=QS*R1*R1*0.5
      RESORT=1.0
      TSF=T(NI-1)
C
C--------     READ INITIAL FIELD VALUES
C
      IF(IREAD.EQ.0) GO TO 100
      READ(10,900) ITERT
      READ(10,901) T
      NITERT=ITERT+MAXIT
      IPRT=ITERT+IPRINT
      JPRT=ITERT+JPRINT
  100 CONTINUE
C
      CALL PROPS
C
C--------     PRINT HEADINGS
C
      WRITE(6,905)
      WRITE(6,912) R1,R2
      WRITE(6,913) TENV,U
      WRITE(6,910) IMON
C
C--------     PRINT INITIAL FIELD VALUES
C
      IF(INITAL.EQ.0) GO TO 110
      CALL PRINTS(1,T,R,1,NI)
      CALL PRINTS(2,ALAMD,R,1,NI)
  110 CONTINUE
C
C-------------------------------------------------------------
CHAPTER   2   2   2   ITERATION LOOP   2   2   2   2   2
C
  200 ITERT=ITERT+1
      SORCE=RESORT
C
C--------     UPDATE MAIN DEPENDENT VARIABLES
C
      CALL SOLVT
C
C--------     UPDATE PROPERTIES
C
      CALL PROPS
C
C--------     INTERMEDIATE OUTPUT
C
      RESORT=RESORT/TMONIN
C
C--------     MONITORING VALUES & RESIDUAL SOURCE
C
      IF(ITERT.NE.IPRT) GO TO 260
      IPRT=IPRT+IPRINT
      WRITE(6,911) ITERT,T(IMON),RESORT
  260 CONTINUE
C
C--------     ALL FIELD VALUES
C
      IF(ITERT.NE.JPRT) GO TO 261
      JPRT=JPRT+JPRINT
      CALL PRINTS(1,T,R,1,NI)
  261 CONTINUE
C
C--------     TERMINATION TESTS
C
      IF((ITERT.GT.NITERT).OR.(RESORT.GT.1.0E3)) GO TO 270
      IF(RESORT.LT.SORMAX) GO TO 271
      GO TO 200
C
C------------------------------------------------------------
CHAPTER   3   3   3   TERMINATING OPERATION   3   3   3   3 
C
C--------     NOT CONVERGING
C
  270 WRITE(6,920)
      GO TO 300
C
C--------     CONVERGING
C
  271 WRITE(6,921)
C
C--------     FINAL OUTPUT
C
  300 WRITE(6,930) ITERT
C
      CALL PRINTS(1,T,R,1,NI)
      CALL PRINTS(2,ALAMD,R,1,NI)
C
C--------     WITHOUT FORMAT
C
      IF(IWRITE.EQ.0) GO TO 310
      WRITE(11,900) ITERT
      WRITE(11,901) T
  310 CONTINUE
C
      STOP
C
C--------     FORMAT STATEMENTS
C
  900 FORMAT(I5)
  901 FORMAT(1P10E13.6)
  905 FORMAT(//5X,'ONE DIMENSIONAL HEAT CONDUCTION PROBLEM',
     &       /6X,'WITH HEAT TRANSFER AND HEAT GENERATION'/)
  910 FORMAT(/10X,' I-----------  ITERATION START !! -------------I'
     &       /10X,' I-------- MONITORING LOCATION =',I4,' --------I'
     &       /17X,'NO',7X,'TMON',8X,'RESORT')
  911 FORMAT(/14X,I5,2E13.4)
  912 FORMAT(/11X,'INNER RADIUS = ',F5.3,2X,'MM'
     &       /11X,'OUTER RADIUS = ',F5.3,2X,'MM')
  913 FORMAT(/11X,'ENVIRONMENTAL TEMPERATURE = ',F6.1,2X,'K'
     &       /11X,'AIR VELOCITY = ',F6.1,2X,'M/S'/)
  920 FORMAT(/' <<<<<<<<< NOT CONVERGING !! >>>>>>>>>>>')
  921 FORMAT(/' >>>>>>>>> CONVERGING !! <<<<<<<<<<')
  930 FORMAT(//' NUMBER OF ITERATIONS = ',I5)
C
      END
C
C*****************************************************************
C
      SUBROUTINE CONVOL
C
C*****************************************************************
C
      COMMON
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
C
C-----------------------------------------------------------------
      R(1)=-DR*0.5
C
      DO 10 I=2,NI
        R(I)=R(I-1)+DR
   10 CONTINUE
C
      R1=(R(NI1)+R(NI1+1))*0.5
      R2=(R(NI-1)+R(NI))*0.5
C
      RW(1)=-DR
C
      DO 20 I=2,NI
        RW(I)=0.5*(R(I)+R(I-1))
   20 CONTINUE
C
      SEW(NI)=DR
C
      DO 30 I=1,NI-1
        SEW(I)=RW(I+1)-RW(I)
   30 CONTINUE
C
      SEP(NI)=DR*0.5
      SPW(NI)=DR*0.5
C
      DO 40 I=1,NI-1
        SEP(I)=RW(I+1)-R(I)
        SPW(I)=R(I)-RW(I)
   40 CONTINUE
C
      RETURN
      END
C
C*****************************************************************
C
      SUBROUTINE INITIA
C
C*****************************************************************
C
      COMMON
     &  /TEMP/T(20),ALAMD(20),TENV,QS,RESORT,TSF
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
     &  /TDMAB/AP(20),AE(20),AW(20),SU(20),SP(20)
C
C--------     SET VARIABLES TO ZERO     --------------------------
      DO 10 I=1,NI
        ALAMD(I)=0.0
        T(I)=TENV
        AP(I)=0.0
        AE(I)=0.0
        AW(I)=0.0
        SU(I)=0.0
        SP(I)=0.0
   10 CONTINUE
C
      RETURN
      END
C
C*****************************************************************
C
      SUBROUTINE PROPS
C
C*****************************************************************
C
      COMMON
     &  /TEMP/T(20),ALAMD(20),TENV,QS,RESORT,TSF
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
     &  /PRP/C1IN,C2IN,C1OUT,C2OUT
C
C--------     THERMAL CONDUCTIVITY     ---------------------------
C--------     INNER MATERIAL
      DO 10 I=1,NI1
        IF(T(I).LT.0.0) T(I)=0.0
        ALAMD(I)=C1IN+C2IN*T(I)
   10 CONTINUE
C--------     OUTER MATERIAL
      DO 20 I=NI1+1,NI
        ALAMD(I)=C1OUT+C2OUT*T(I)
   20 CONTINUE
C
      RETURN
      END
C
C*****************************************************************
C
      SUBROUTINE SOLVT
C
C*****************************************************************
C
      COMMON
     &  /TEMP/T(20),ALAMD(20),TENV,QS,RESORT,TSF
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
     &  /TDMAB/AP(20),AE(20),AW(20),SU(20),SP(20)
     &  /DSIZE/IT
C
C--------     ASSEMBLY OF COEFFICIENTS     -----------------------
C
      DO 10 I=2,NI-1
        VOL=R(I)*SEW(I)
        AREAE=RW(I+1)
        AREAW=RW(I)
        ALDE=1.0/(SPW(I+1)/ALAMD(I+1)+SEP(I)/ALAMD(I))
        AE(I)=ALDE*AREAE
        ALDW=1.0/(SPW(I)/ALAMD(I)+SEP(I-1)/ALAMD(I-1))
        AW(I)=ALDW*AREAW
        IF(I.LE.NI1) THEN
          SOURCE=QS*VOL
        ELSE
          SOURCE=0.0
        ENDIF
        SU(I)=SOURCE
        SP(I)=0.0
   10 CONTINUE
C
C--------     BOUNDARY CONDITIONS
C
      CALL BOUND
C
C--------     FINAL COEFFICIENTS ASSEMBLY
C
      RESORT=0.0
C
      DO 20 I=2,NI-1
        AP(I)=AE(I)+AW(I)-SP(I)
        RESOR=AE(I)*T(I+1)+AW(I)*T(I-1)-AP(I)*T(I)+SU(I)
        RESORT=RESORT+ABS(RESOR)
   20 CONTINUE
C
C--------     SOLVE EQUATIONS BY T.D.M.A 
C
      CALL TDMA(2,NI-1,T)
C
      RETURN
      END
C
C*****************************************************************
C
      SUBROUTINE BOUND
C
C*****************************************************************
C
      COMMON
     &  /TEMP/T(20),ALAMD(20),TENV,QS,RESORT,TSF
     &  /GEOM/R(20),RW(20),SEW(20),SEP(20),SPW(20),NI,NI1,DR,R1,R2
     &  /TDMAB/AP(20),AE(20),AW(20),SU(20),SP(20)
     &  /BND/U,ALMD,AK1,AK2,AL1,AL2,AN1,AN2
C
C--------      BOUNDARY CONDITIONS     ---------------------------
C--------      AT SYMMETRY AXIS (I=2)
c
      T(1)=T(2)
      AW(2)=0.0
C
C--------      AT OUTER SURFACE (I=NI-1)
c
      NIM1=NI-1
      TFILM=(TENV+TSF)*0.5
c
C--------      HEAT TRANSFER COEFFICIENT
c
      ALMD=AL1+AL2*TFILM
      ANYU=AN1+AN2*TFILM
      REY=2.0*U*R2/ANYU
      ANUS=AK1*REY**AK2
      ALPHA=0.5*ANUS/R2*ALMD
C
      DM1=ALPHA
      DM2=ALAMD(NIM1)/SEP(NIM1)
      TSF=(DM1*TENV+DM2*T(NIM1))/(DM1+DM2)
      T(NI)=TSF
C
      ALDEN=1.0/(1.0/DM1+1.0/DM2)
      AREAE=RW(NI)
      TMULT=ALDEN*AREAE
      SU(NIM1)=SU(NIM1)+TENV*TMULT
      SP(NIM1)=SP(NIM1)-TMULT
      AE(NIM1)=0.0
C
      RETURN
      END
C
C****************************************************************
C
      SUBROUTINE TDMA(ISTART,IEND,PHI)
C
C****************************************************************
C
       COMMON
     &  /TDMAB/A(20),B(20),C(20),D(20),E(20)
     &  /DSIZE/IT
C
      DIMENSION PHI(IT),P(20),Q(20)
C
C----------------------------------------------------------------
C
      DO 10 I=ISTART,IEND
        DUMMY=1.0/(A(I)-C(I)*P(I-1))
        P(I)=B(I)*DUMMY
        Q(I)=(D(I)+C(I)*Q(I-1))*DUMMY
   10 CONTINUE
C
      DO 20 I=IEND,ISTART,-1
        PHI(I)=P(I)*PHI(I+1)+Q(I)
   20 CONTINUE
C
      RETURN
      END
C
C***************************************************************
C
      SUBROUTINE PRINTS(IP,AL,RA,IS,IE)
C
C***************************************************************
C
      COMMON
     &  /DSIZE/IT
C
      DIMENSION AL(IT),RA(IT)
C
C---------------------------------------------------------------
C
      IF(IP.EQ.1) WRITE(6,810)
      IF(IP.EQ.2) WRITE(6,820)
C
      IVAL=10
      ISTART=IS
C
   20 IEND=ISTART+IVAL-1
      IF(IEND.GT.IE) IEND=IE
      WRITE(6,700) (I,I=ISTART,IEND)
      WRITE(6,710) (RA(I),I=ISTART,IEND)
      WRITE(6,720) (AL(I),I=ISTART,IEND)
      IF(IEND.EQ.IE) GO TO 100
      ISTART=ISTART+IVAL
      GO TO 20
C
  100 CONTINUE
      RETURN
C
  700 FORMAT(/7X,'I    =',I6,10I9)
  710 FORMAT(7X,'R (M)=',10F10.3)
  720 FORMAT(/15X,1P10E10.3)
  810 FORMAT(//' =====>> T PROFILE  (K)')

  820 FORMAT(//' =====>> LAMDA PROFILE (W/MK)')
C
      END
