C
C********************************************************************
C
C     THE OTHER SUBPROGRAMS, THOMAS AND GAMWAL.
C
C--------------------------------------------------------------------
      SUBROUTINE THOMAS(N,AI,BI,CI,DI)
C
C-----THOMAS' ALGORITHM FOR TRI-DIAGONAL MATRIX.
C
      DIMENSION AI(N),BI(N),CI(N),DI(N)
      NM1=N-1
C-----USE BI FOR PI AND DI FOR QI.
      BI(1)=BI(1)/AI(1)
      DI(1)=DI(1)/AI(1)
      DO 10 I=2,N
      IM1=I-1
      DN=AI(I)-CI(I)*BI(IM1)
      IF(DN.EQ.0.) WRITE(6,100) DN
      BI(I)=BI(I)/DN
   10 DI(I)=(DI(I)+CI(I)*DI(IM1))/DN
C-----NOTE, RESULTS ARE STORED IN DI.
      DO 20 I=1,NM1
      NMI=N-I
   20 DI(NMI)=BI(NMI)*DI(NMI+1)+DI(NMI)
  100 FORMAT(1X,10(1H*),22H ZERO DETERMINANT, DN=,E10.3)
      RETURN
      END
C
C-------------------------------------------------------------------
      FUNCTION GAMWAL(UP,DY,RENO,PRLAM,IPHI)
C
C-----THIS FUNCTION EVALUATES THE EFFECTIVE EDDY DIFFUSIVITY AT
C-----THE WALL, GAMWAL, USING THE LOGARITHMIC LAW OF THE WALL.
C
      AK=0.41
      B=5.
C-----SET LAMINAR DIFFUSION COEFFICIENT.
      GAMWAL=1./RENO/PRLAM
      U=ABS(UP)
      YWALL=0.5*DY
      RE=U*YWALL*RENO
C-----USE LAMINAR DIFFUSIVITY WHEN THE VELOCITY IS TOO LOW.
      IF(RE.LT.117.)  GO TO 100
C-----SHOOT TWICE TO FIND UTAU FOR GIVEN U.
      C0=ALOG(RE)/AK+B
      C=C0*(1.-ALOG(C0)/(1.+AK*C0))
      UTAU=U/C/(1.-(AK*(C-C0)+ALOG(C))/(1.+AK*C))
      PRTUR=1.
C-----SET TURBULENT PRANDTL NUMBER.
      IF(IPHI.EQ.4) PRTUR=0.9
      GAMWAL=UTAU**2/U*YWALL/PRTUR
C-----USE TEMPERATURE LAW OF THE WALL.
      IF(IPHI.NE.4) GO TO 100
      GAMWAL=GAMWAL/(1.+UTAU/U*9.24*((PRLAM/PRTUR)**0.75-1.))
  100 RETURN
      END
C
C********************************************************************
C                                                                   C
C      USER'S WORKING AREA                                          C
C                                                                   C
C     THIS PART CONSISTS OF BLOCK DATA AND 2 ENTRIES, DATAIN AND    C
C     GAMAIN IN SUBROUTINE USERIN. A USER SHOULD CONCENTRATE ON     C
C     THIS PART, AND SHOULD NOT ALTER ANY OTHER PART OF THIS        C
C     PROGRAM UNLESS HE WISHES TO CHANGE THE SIZES OF ARRAYS        C
C     OR ADD OTHER TRANSPORT EQUATIONS TO IT.                       C
C                                                                   C
C********************************************************************
C
      BLOCK DATA
C-----
C     SET NTIMST=1 FOR STEADY FLOW SOLUTIONS.
C     SET RENO=SQRT(GRNO) FOR PURE FREE CONVECTION.
C     NI AND NJ MUST BE LESS THAN OR EQUAL TO 39.
C     THE REFERENCE POINT (IREF,JREF) MUST BE IN THE FLOW FIELD.
C     REMEMBER, IF ICNTDF=1, NUMERICAL STABILITY IS NOT GUARANTEED.
C     ********************** CAUTION ******************************
C
C     NTIMST = NO. OF TIME STEPS    ; NITERT = NO. OF ITERATIONS
C     DTIME = TIME INCREMENT        ; RENO = REYNOLDS NUMBER
C     GRNO = GRASHOF NUMBER         ; DIRCOS = 1 WHEN G POINTS WEST
C     PR = LAMINAR PRANDTL NUMBER   ; ITURB = 1 FOR TURBULENT FLOW
C     IRAD = 1 FOR CYLIND. COORD.   ; ICNTDF = 1 FOR CENTRAL DIFF.
C     NI= NO. OF VERTICAL C.S.LINES ; NJ= NO. OF HORIZONTAL C.S.LINES
C     ISCAN= I ADDRESS FOR SCANNING ; JSCAN= J ADDRESS FOR SCANNING
C     IREF = I ADDRESS FOR PREF     ; JREF = J ADDRESS FOR PREF
C     INTPRI= INTERVAL FOR PRINTING ; IREAD= 1 WHEN RESULTS EXIST
C     NLUMP = NO. OF LUMPS          ; LUMPW-N = ADDRESSES FOR LUMP
C     X = X-COORD. FOR C.S. LINES   ; Y = Y-COORD. FOR C.S. LINES
C     ISOLVE = 1 IF SOLUTION NEEDED ; IPRINT = 1 IF PRINTING NEEDED
C-----
      COMMON/BLOCK0/ NITERT,ITERT,INTPRI,NLUMP,ISCAN,JSCAN,DTIME,
     &               ISOLVE(10),IPRINT(10),
     &               LUMPE(10),LUMPW(10),LUMPN(10),LUMPS(10),
     &               IRMAX,JRMAX,RESMAX
      COMMON/BLOCK1/ NTIMST,ITURB,IRAD,IREAD,NPHI,IPHI,IU,IV,IPC,IT,
     &               IW,IP,NI,NJ,JDIM,IREF,JREF,RENO,GRNO,DIRCOS,PR,
     &               TIME,IFORS(40,40),ALPHA(10),NIM1,NJM1,NIP1,NJP1,
     &               IEND,JEND,ICNTDF
      COMMON/BLOCK2/ X(40),Y(40),XP(40),YP(40),DXP(40),DYP(40),
     &               DXU(40),DYV(40),DX(40),DY(40),DELX(40),DELY(40),
     &               FACX(40),FACY(40),RX(40),RY(40),BIGNO,ZERO
C
      DATA NTIMST,NITERT,DTIME,RENO,GRNO,DIRCOS,PR,ITURB,IRAD/
     &         10,    30,   .1, 100.,  0.,  0., 1.,    0,   0/
      DATA NI, NJ, ISCAN, JSCAN, IREF, JREF, INTPRI, ICNTDF, IREAD/
     &     11, 11,     6,    11,    2,    2,      2,      0,     0/
      DATA NLUMP/4/
      DATA LUMPW(1),LUMPE(1),LUMPS(1),LUMPN(1)/
     &            1,       1,       1,      12/
      DATA LUMPW(2),LUMPE(2),LUMPS(2),LUMPN(2)/
     &            1,      12,      12,      12/
      DATA LUMPW(3),LUMPE(3),LUMPS(3),LUMPN(3)/
     &           12,      12,       1,      12/
      DATA LUMPW(4),LUMPE(4),LUMPS(4),LUMPN(4)/
     &            1,      12,       1,       1/
      DATA (X(I),I=1,11)/
     &  0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1./
      DATA (Y(J),J=1,11)/
     &  0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1./
      DATA (ISOLVE(I),I=1,5)/
     &     1, 1, 1, 1, 0/
      DATA (IPRINT(I),I=1,6)/
     &     1, 1, 0, 1, 0, 1/
C-----DO NOT ALTER THE DATA BELOW --------------------------------
      DATA (ALPHA(I),I=1,6)/
     &     0.3, 0.3, 1.0, 0.5, 0.5, 0.2/
      DATA NPHI,IU,IV,IPC,IT,IW,IP,JDIM/
     &        5, 1, 2,  3, 4, 5, 6,  40/
      END
C
C--------------------------------------------------------------------
      SUBROUTINE USERIN
C
C-----THIS SUBROUTINE HAS 2 ENTRIES, DATAIN AND GAMAIN.
C
      COMMON         PHI(40,240),PHIOLD(40,240),GAM(40,40),SC(40,40),
     &               SP(40,40),AE(40,40),AW(40,40),AN(40,40),
     &               AS(40,40),SDX(40,40),SDY(40,40)
      DIMENSION U(40,40),V(40,40),PC(40,40),T(40,40),W(40,40),P(40,40)
      EQUIVALENCE    (PHI(1,1),U(1,1)),(PHI(1,41),V(1,1)),
     &               (PHI(1,81),PC(1,1)),(PHI(1,121),T(1,1)),
     &               (PHI(1,161),W(1,1)),(PHI(1,201),P(1,1))
      COMMON/BLOCK1/ NTIMST,ITURB,IRAD,IREAD,NPHI,IPHI,IU,IV,IPC,IT,
     &               IW,IP,NI,NJ,JDIM,IREF,JREF,RENO,GRNO,DIRCOS,PR,
     &               TIME,IFORS(40,40),ALPHA(10),NIM1,NJM1,NIP1,NJP1,
     &               IEND,JEND,ICNTDF
      COMMON/BLOCK2/ X(40),Y(40),XP(40),YP(40),DXP(40),DYP(40),
     &               DXU(40),DYV(40),DX(40),DY(40),DELX(40),DELY(40),
     &               FACX(40),FACY(40),RX(40),RY(40),BIGNO,ZERO
C
C-------------------------------------------------------------------
      ENTRY DATAIN
C
C-----ALL PHI-S ALONG THE VELOCITY-KNOWN BOUNDARIES SHALL BE
C     ASSUMED TO BE ZERO UNLESS DISTRIBUTIONS ARE SPECIFIED
C-----HERE FOR NON-ZERO PHI-S.
C
      DO 100 I=1,NI
  100 U(I,NJP1)=1.
      DO 200 I=1,NIP1
  200 T(I,NJP1)=1.
      DO 300 J=2,NJ
      T(1,J)=T(2,J)
  300 T(NIP1,J)=T(NI,J)
      RETURN
C
C----------------------------------------------------------------
      ENTRY GAMAIN
C
C-----FEED GAMMAS IN. DESCRIBE TEMPORAL OR SPATIAL VARIATION OF
C-----GAMMA, OR SUPPLY A HIGH REYNOLDS NUMBER TURBULENCE MODEL.
C
      RETURN
      END
