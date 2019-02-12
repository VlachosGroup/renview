C///////////////////////////////////////////////////////////////////
C
      PROGRAM TRANFT
C
C     WRITTEN BY:
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C
C/////////////////////////////////////////////////////////////////////
C
C     VERSION 1.9
C
C     CHANGES FROM PREVIOUS VERSION:
C     1.  Restructured data in OMEG12
C     2.  Changed REAL*8 to DOUBLE PRECISION
C     3.  Changed POLFIT and PCOEF to call single and double precision
C         SLATEC subroutines POLFIT,DPOLFT and PCOEF,DPCOEF
C     4.  Change vms open statements
C     CHANGES FROM VERSION 1.3
C     1.  Change THIGH to 3500
C     2.  Change name to TRANFT to conform to ANSI standard
C     3.  Add "unix" and "snla" change blocks
C     CHANGES FROM VERSION 1.4
C     1.  modify OPEN statements for unix
C     CHANGES FOR VERSION 1.6
C     1.  Find THIGH and TLOW from species thermodynamic data
C     CHANGES FOR VERSION 1.9
C     1.  Change linking file to include version, precision, error
C         status, and required length of work arrays
C     2.  Allow user input from a file.
C     3.  Implement non-formatted transport data input
C
C
C/////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C.....Ashish modification: KDIM, LENICK, LENRCK, and LENCCK increased
C..........................for large gas mechs.
C
      PARAMETER (LINC=25, LINKTP=35, LTRAN=31, LDATA=10, LOUT=6,
     1           KDIM=200, NFDIM=165, NO=4, NT=50, NCHECK=1000,
     2           LENICK=12000, LENRCK=12000, LENCCK=200, MAXFIT=7,
     3           NRANGE=2, MAXTMP=3)
C
      DIMENSION EPS(KDIM), SIG(KDIM), DIP(KDIM), POL(KDIM), ZROT(KDIM),
     1          NLIN(KDIM), WT(KDIM), CV(KDIM), FITWT(NT), FITRES(NT),
     2          ALOGT(NT), XLA(NT), XETA(NT), XD(NT), FITWRK(NFDIM),
     3          COFLAM(NO,KDIM), COFETA(NO,KDIM), COFD(NO,KDIM,KDIM),
     4          COFTD(NO,KDIM,10), ICKWRK(LENICK), RCKWRK(LENRCK),
     5          VALUE(6), KTDIF(KDIM), NTEMP(KDIM), TEMP(MAXTMP,KDIM),
     6          THERM(MAXFIT, NRANGE, KDIM)
C
      CHARACTER CCKWRK(LENCCK)*16, KSYM(KDIM)*16, LINE*80, VERS*16,
     1          PREC*16, FNAME*30
      LOGICAL IERR, KERR
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM,
     1                DIPMIN
C
      DATA NLIN/KDIM*NCHECK/, KOUNT/0/, EMAXL,EMAXE,EMAXD,EMAXTD/4*0.0/,
     1     KERR/.FALSE./, FNAME/' '/
C
C
      PI     = 3.1415926535
      EPSIL  = 0.0
      FDTCGS = 1.0E-18
      FATCM  = 1.0E+08
      DIPMIN = 1.0E-20
C
C*****OPEN statement > vms
CC
CC         OPEN THE LINK FILE
C      OPEN (LINC, STATUS='OLD', FORM='UNFORMATTED')
C      OPEN (LTRAN,  STATUS='OLD', FORM='FORMATTED', SHARED, READONLY)
C      OPEN (LINKTP, STATUS='NEW', FORM='UNFORMATTED')
C      INQUIRE (LDATA, NAME=FNAME)
C      IF (FNAME.NE.' ') OPEN (LDATA, FILE=FNAME, STATUS='OLD',
C     1    FORM='FORMATTED')
C*****END OPEN statement > vms
C*****OPEN statement > unix
      OPEN (LINC, FORM='UNFORMATTED', FILE='INP.d/cklink')
      OPEN (LTRAN,  FORM='FORMATTED'  , FILE='INP.d/trandat')
      OPEN (LINKTP, FORM='UNFORMATTED', FILE='INP.d/tplink')
      INQUIRE (LDATA, NAME=FNAME)
      IF (FNAME.NE.' ') OPEN (LDATA, FILE=FNAME, FORM='FORMATTED')
C*****END OPEN statement > unix
C
C          WRITE VERSION NUMBER
C
      VERS = '1.9'
      WRITE (LOUT, 15) VERS(:3)
   15 FORMAT(
     1/' TRANFT:  Transport property fitting,',
     2/'           CHEMKIN-II Version ',A,', November 1990',
C*****precision > double
     3/'           DOUBLE PRECISION')
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C     3/'           SINGLE PRECISION')
C      PREC = 'SINGLE'
C*****END precision > single
C
      WRITE (LOUT, 8310)
C
C         INITIALIZE CHEMKIN
C
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINC, LOUT, ICKWRK,
     1             RCKWRK, CCKWRK)
      CALL CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
	write(*,*)'in tranfit',kk
      IF (KK .GT. KDIM) THEN
         WRITE (LOUT, 8300) KDIM
         STOP
      ENDIF
C
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      KERR = KERR.OR.IERR
      CALL CKWT   (ICKWRK, RCKWRK, WT)
      CALL CKRP   (ICKWRK, RCKWRK, RU, RUC, PATM)
      P = PATM
      BOLTZ = RU/6.02217E23
C
      CALL CKATHM (MAXFIT, NRANGE, ICKWRK, RCKWRK, MAXTMP,
     1             NTEMP, TEMP, THERM)
      TLOW  = TEMP(1,1)
      THIGH = TEMP(NTEMP(1),1)
      DO 10 K = 2, KK
         THIGH = MIN (THIGH, TEMP (NTEMP(K), K))
         TLOW  = MAX (TLOW,  TEMP (1,K))
   10 CONTINUE
      DT = (THIGH-TLOW) / (NT-1)
C
      WRITE (LOUT, 20) TLOW, THIGH
   20 FORMAT (/,' DATA HAS BEEN FIT OVER THE RANGE:  TLOW=',F9.2,
     1        ', THIGH=',F9.2)
C
C       READ THE TRANSPORT DATA BASE
C
      WRITE (LOUT, 7020)
      LIN = LTRAN
C
   50 CONTINUE
      READ (LIN, '(A)', END=500) LINE
      ILEN = IPPLEN(LINE)
      IF (ILEN .GT. 0) THEN
         K1 = IFIRCH(LINE(:ILEN))
         K2 = K1 + INDEX(LINE(K1:ILEN),' ') - 1
         CALL CKCOMP (LINE(K1:K2-1), KSYM, KK, KFIND)
C
         IF (KFIND .GT. 0) THEN
            CALL CKXNUM (LINE(K2:ILEN), 6, LOUT, NVAL, VALUE, IERR)
            KERR = KERR.OR.IERR
            IF (IERR) WRITE (LOUT, 8000)
            WRITE (LOUT, 7010) LINE
C
            NLIN(KFIND) = INT(VALUE(1))
            EPS(KFIND)  = VALUE(2)
            SIG(KFIND)  = VALUE(3)
            DIP(KFIND)  = VALUE(4)
            POL(KFIND)  = VALUE(5)
            ZROT(KFIND) = VALUE(6)
         ENDIF
      ENDIF
      GO TO 50
C
  500 CONTINUE
      IF (FNAME.NE.' ' .AND. LIN.EQ.LTRAN) THEN
         LIN = LDATA
         GO TO 50
      ENDIF
C
      DO 600 K = 1, KK
         IF (NLIN(K) .EQ. NCHECK) THEN
            DO 750 J = K, KK
               IF (NLIN(J) .EQ. NCHECK) WRITE (LOUT, 8010) KSYM(J)
  750       CONTINUE
            KERR = .TRUE.
         ENDIF
  600 CONTINUE
C
      IF (.NOT. KERR) THEN
C
C        FIT THE CONDUCTIVITIES AND VISCOSITIES
C
         CALL LAMFIT (KK, NT, NO, NFDIM, LOUT, WT, SIG, EPS, DIP,
     1                ZROT, NLIN, P, TLOW, DT, ALOGT, FITRES, FITWT,
     2                FITWRK, XLA, XETA, CV, ICKWRK, RCKWRK, COFLAM,
     3                COFETA, EMAXL, EMAXE)
C
C        FIT THE DIFFUSION COEFFICIENTS
C
         CALL DIFFIT (KK, NT, NO, NFDIM, KDIM, LOUT, WT, SIG, EPS,
     1                DIP, POL, P, TLOW, DT, ALOGT, FITRES, FITWT,
     2                FITWRK, XD, COFD, EMAXD)
C
C        FIT THE THERMAL DIFFUSION RATIOS
C
         CALL THMFIT (KK, NT, NO, KDIM, LOUT, WT, SIG, EPS, DIP,
     1                POL, TLOW, DT, ALOGT, FITRES, FITWT, FITWRK,
     2                XD, NLITE, KTDIF, COFTD, EMAXTD)
C
C        PRINT THE FITS
C
         WRITE (LOUT, 7030)
         WRITE (LOUT, 7035) EMAXL
         WRITE (LOUT, 8200)
         WRITE (LOUT, 8100) (KSYM(K), (COFLAM(N,K),N=1,NO), K=1,KK)
C
         WRITE (LOUT, 7050)
         WRITE (LOUT, 7035) EMAXE
         WRITE (LOUT, 8200)
         WRITE (LOUT, 8100) (KSYM(K), (COFETA(N,K),N=1,NO), K=1,KK)
C
         WRITE (LOUT, 7060)
         WRITE (LOUT, 7035) EMAXD
         DO 2300 J = 1, KK
            WRITE (LOUT, 8200)
            WRITE (LOUT, 8110)
     1            (KSYM(J), KSYM(K), (COFD(N,J,K),N=1,NO), K=1,J)
 2300    CONTINUE
C
         WRITE (LOUT, 7070)
         WRITE (LOUT, 7035) EMAXTD
         DO 2400 M = 1, NLITE
            K = KTDIF(M)
            WRITE (LOUT, 8200)
            WRITE (LOUT, 8110)
     1            (KSYM(K), KSYM(J), (COFTD(N,J,M),N=1,NO), J=1,KK)
 2400    CONTINUE
      ELSE
         WRITE (LOUT, '(/A)')
     1   ' WARNING...THERE IS AN ERROR IN THE TRANSPORT LINKING FILE'
      ENDIF
C
C        WRITE THE LINKING FILE
C
      REWIND LINKTP
      WRITE (LINKTP) VERS, PREC, KERR
C
      LENIMC = 4*KK + NLITE
      LENRMC = (19 + 2*NO + NO*NLITE)*KK + (15+NO)*KK**2
C
      WRITE (LINKTP) LENIMC, LENRMC, NO, KK, NLITE
      WRITE (LINKTP) PATM, (WT(K), EPS(K), SIG(K), DIP(K),
     1               POL(K), ZROT(K), NLIN(K), K=1,KK),
     2               ((COFLAM(N,K),N=1,NO),K=1,KK),
     4               ((COFETA(N,K),N=1,NO),K=1,KK),
     5               (((COFD(N,J,K),N=1,NO),J=1,KK),K=1,KK),
     6               (KTDIF(N),N=1,NLITE),
     7               (((COFTD(N,J,L),N=1,NO),J=1,KK),L=1,NLITE)
      REWIND LINKTP
C
      STOP
C
C       FORMATS
C
 7000 FORMAT (A)
 7010 FORMAT (1X,A)
 7020 FORMAT (///,' TRANSPORT PARAMETERS FROM DATA BASE:',/)
 7030 FORMAT (///,' COEFFICIENTS FOR SPECIES CONDUCTIVITIES',/)
 7035 FORMAT ( '  MAXIMUM FITTING ERROR = ', 1PE12.3)
 7040 FORMAT (///,' COEFFICIENTS FOR MONOTOMIC PARTS OF CONDUCTIVITIES'
     1, /)
 7050 FORMAT (///,' COEFFICIENTS FOR SPECIES VISCOSITIES',/)
 7060 FORMAT (///,' COEFFICIENTS FOR SPECIES DIFFUSION COEFFICIENTS',
     1  /)
 7070 FORMAT (///,' COEFFICIENTS FOR THERMAL DIFFUSION RATIOS',/)
 8000 FORMAT (10X, 'ERROR IN CKXNUM READING FROM TRANSPORT DATA BASE')
 8010 FORMAT (10X, 'ERROR...TRANSPORT DATA NOT GIVEN FOR ',A10)
C
 8100 FORMAT (2X, A10, 4E12.3)
 8110 FORMAT (2X, A10, A10, 4E12.3)
 8200 FORMAT (1X, ' ')
 8300 FORMAT (10X, 'ERROR...THE TRANSPORT FITTING CODE IS DIMENSIONED
     1FOR ONLY', I3, ' SPECIES')
 8310 FORMAT(///,' OUTPUT FROM THE TRANSPORT PROPERTY FITTING PACKAGE'
     1,//)
      END
C
      SUBROUTINE LAMFIT (KK, NT, NO, NFDIM, NOUT, WT, SIG, EPS, DIP,
     1                   ZROT, NLIN, P, TLOW, DT, ALOGT, FITRES, FITWT,
     2                   FITWRK, XLA,  XETA, CV, ICKWRK, RCKWRK,
     3                   COFLAM, COFETA, EMAXL, EMAXE)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), ZROT(*), NLIN(*), CV(*),
     1          FITRES(*), FITWRK(*), ALOGT(*), FITWT(*), XLA(*),
     2          XETA(*), ICKWRK(*), RCKWRK(*), COFLAM(NO,*),
     3          COFETA(NO,*)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      FCON = 0.5 * FDTCGS**2 * FATCM**3 / BOLTZ
C
      DO 1000 K = 1, KK
         DST = FCON * DIP(K)**2 / (EPS(K) * SIG(K)**3)
         HELPE = 2.6693E-5 * SQRT(WT(K)) / SIG(K)**2
         HELPD = 2.6280E-3 * PATM / (P * SQRT(WT(K)) * SIG(K)**2)
C
         DO 300 N = 1, NT
            T = TLOW + (N-1)*DT
            TR = T / EPS(K)
            CALL CKCVML (T, ICKWRK, RCKWRK, CV)
            ALOGT(N) = LOG(T)
            DII = SQRT(T**3)*HELPD / OMEG12(1,TR,DST)
            XETA(N) = SQRT(T)*HELPE / OMEG12(2,TR,DST)
            RODET = DII * P * WT(K) / (RU * T * XETA(N))
            AA = 2.5 - RODET
            CALL PARKER (T, EPS(K), PARK)
            IF (NLIN(K) .EQ. 2) THEN
               BB = ZROT(K)*PARK + 2.0*(5.0/2.0 + RODET)/PI
               CROCT = 1.0
            ELSE
               BB = ZROT(K)*PARK + 2.0*(5.0/3.0 + RODET)/PI
               CROCT = 2.0/3.0
            ENDIF
            FTRA = 2.5 * (1.0 - 2.0 * CROCT * AA / (BB*PI))
            FROT = RODET * (1.0 + 2.0 * AA / (BB*PI))
            FVIB = RODET
            IF (NLIN(K) .EQ. 0) THEN
               FAKTOR = 5.0/2.0 * 1.5*RU
            ELSEIF (NLIN(K) .EQ. 1) THEN
               FAKTOR = (FTRA*1.5 + FROT)*RU + FVIB*(CV(K)-2.5*RU)
            ELSEIF (NLIN(K) .EQ. 2) THEN
               FAKTOR = (FTRA+FROT)*1.5*RU + FVIB*(CV(K)-3.0*RU)
            ENDIF
            XLA(N)  = LOG( XETA(N)/WT(K) * FAKTOR)
            XETA(N) = LOG(XETA(N))
300      CONTINUE
C
C      FIT CONDUCTIVITY
C
         FITWT(1) = -1.0
         EPSL = EPSIL
C
C*****precision > double
         CALL DPOLFT (NT, ALOGT, XLA, FITWT, NO-1, NORD, EPSL, FITRES,
     1                IERR, FITWRK)
C*****END precision > double
C*****precision > single
C         CALL POLFIT (NT, ALOGT, XLA, FITWT, NO-1, NORD, EPSL, FITRES,
C     1                IERR, FITWRK)
C*****END precision > single
C
         EMAXL = MAX (EMAXL, EPSL)
         IF (IERR .NE. 1) THEN
            WRITE (NOUT,7000)
            STOP
         ENDIF
C
         CCC = 0.0
C*****precision > double
         CALL DPCOEF (NORD, CCC, COFLAM(1,K), FITWRK)
C*****END precision > double
C*****precision > single
C         CALL PCOEF (NORD, CCC, COFLAM(1,K), FITWRK)
C*****END precision > single
C
C      FIT VISCOSITY
C
         FITWT(1) = -1.0
         EPSL = EPSIL
C*****precision > double
         CALL DPOLFT (NT, ALOGT, XETA, FITWT, NO-1, NORD, EPSL,
     1                FITRES, IERR, FITWRK)
C*****END precision > double
C*****precision > single
C         CALL POLFIT (NT, ALOGT, XETA, FITWT, NO-1, NORD, EPSL,
C     1                FITRES, IERR, FITWRK)
C*****END precision > single
C
         EMAXE = MAX (EMAXE, EPSL)
         IF (IERR .NE. 1) THEN
            WRITE (NOUT,7000)
            STOP
         ENDIF
C
         CCC = 0.0
C*****precision > double
         CALL DPCOEF (NORD, CCC, COFETA(1,K), FITWRK)
C*****END precision > double
C*****precision > single
C         CALL PCOEF (NORD, CCC, COFETA(1,K), FITWRK)
C*****END precision > single
C
1000  CONTINUE
C
7000  FORMAT(9X,'ERROR IN POLFIT WHILE FITTING VISCOSITY OR ',
     1'CONDUCTIVITY')
      RETURN
      END
C
      SUBROUTINE DIFFIT (KK, NT, NO, NFDIM, KDIM, NOUT, WT, SIG, EPS,
     1                   DIP, POL, P, TLOW, DT, ALOGT, FITRES, FITWT,
     2                   FITWRK, XD, COFD, EMAXD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), POL(*), FITRES(*),
     1          FITWT(*), FITWRK(*), ALOGT(*), XD(*), COFD(NO,KDIM,*)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      DO 50 N = 1, NT
         T = TLOW + (N-1)*DT
         ALOGT(N) = LOG(T)
   50 CONTINUE
C
      FCON = 0.5 * FDTCGS**2 * FATCM**3 / BOLTZ
C
      DO 1000 J = 1, KK
         DO 1000 K = 1, J
            SIGM = 0.5 * (SIG(J)+SIG(K))
            EPSM = SQRT(EPS(J)*EPS(K))
            IF (DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN) THEN
C
C        K IS POLAR, J IS NONPOLAR
C
               DST = 0.
               XI = 1.0 +
     1              POL(J) * FCON * DIP(K)**2 * SQRT(EPS(K)/EPS(J)) /
     2              ( 2.0 * SIG(J)**3 * EPS(K) * SIG(K)**3)
C
            ELSEIF (DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN) THEN
C
C        J IS POLAR, K IS NONPOLAR
C
               DST = 0.
               XI = 1.0 +
     1              POL(K) * FCON * DIP(J)**2 * SQRT(EPS(J)/EPS(K)) /
     2              (2.0 * SIG(K)**3 * EPS(J) * SIG(J)**3)
C
C
            ELSE
C
C         NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
C
               DST = FCON * DIP(J) * DIP(K) / (EPSM * SIGM**3)
               XI = 1.0
            ENDIF
C
            SIGM = SIGM * XI**(-1.0/6.0)
            EPSM = EPSM * XI**2
            HELP1 = (WT(J)+WT(K)) / (2.0*WT(J)*WT(K))
            DO 500 N = 1, NT
               T = TLOW + (N-1)*DT
               TR = T / EPSM
               HELP2 = 0.002628 * SQRT(HELP1*T**3)
               XD(N) = HELP2*PATM / (OMEG12(1,TR,DST) * SIGM**2 * P)
               XD(N) = XD(N) *
     1            D12(WT(J),WT(K),T,EPS(J),EPS(K),SIG(J),SIG(K),DST)
               XD(N) = LOG(XD(N))
500         CONTINUE
C
            FITWT(1) = -1.0
            EPSL = EPSIL
C
C*****precision > double
            CALL DPOLFT (NT, ALOGT, XD, FITWT, NO-1, NORD, EPSL,
     1                   FITRES, IERR, FITWRK)
C*****END precision > double
C*****precision > single
C            CALL POLFIT (NT, ALOGT, XD, FITWT, NO-1, NORD, EPSL,
C     1                   FITRES, IERR, FITWRK)
C*****END precision > single
C
            EMAXD = MAX (EMAXD, EPSL)
            IF (IERR .NE. 1) THEN
               WRITE (NOUT,7000) J,K
               STOP
            ENDIF
C
            CCC = 0.0
C
C*****precision > double
            CALL DPCOEF (NORD, CCC, COFD(1,K,J), FITWRK)
C*****END precision > double
C*****precision > single
C            CALL PCOEF (NORD, CCC, COFD(1,K,J), FITWRK)
C*****END precision > single
1000  CONTINUE
C
C          FILL OUT THE LOWER TRIANGLE
C
      DO 2000 K = 1, KK-1
         DO 2000 J = K+1, KK
            DO 2000 N = 1, NO
               COFD(N,J,K) = COFD(N,K,J)
 2000 CONTINUE
      RETURN
7000  FORMAT (10X, 'ERROR IN FITTING', I3, ',', I3,
     1                              'DIFFUSION COEFFICIENT')
      END
C
      SUBROUTINE THMFIT (KK, NT, NO, KDIM, NOUT, WT, SIG, EPS, DIP,
     1                   POL, TLOW, DT, AT, FITRES, FITWT, FITWRK, TD,
     2                   NLITE, KTDIF, COFTD, EMAXTD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), POL(*), KTDIF(*),
     1          COFTD(NO,KDIM,*), FITRES(*), FITWT(*), FITWRK(*),
     2          AT(*), TD(*), FITAST(7), FITBST(7), FITCST(7)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      DATA FITAST/ .1106910525E+01, -.7065517161E-02,-.1671975393E-01,
     1             .1188708609E-01,  .7569367323E-03,-.1313998345E-02,
     2             .1720853282E-03/
C
      DATA FITBST/ .1199673577E+01, -.1140928763E+00,-.2147636665E-02,
     1             .2512965407E-01, -.3030372973E-02,-.1445009039E-02,
     2             .2492954809E-03/
C
      DATA FITCST/ .8386993788E+00,  .4748325276E-01, .3250097527E-01,
     1            -.1625859588E-01, -.2260153363E-02, .1844922811E-02,
     2            -.2115417788E-03/
C
      NLITE = 0
      DO 50 N = 1, NT
         T = TLOW + (N-1)*DT
         AT(N) = T
   50 CONTINUE
C
      DO 1100 J = 1, KK
C
         IF (WT(J) .LE. 5.0) THEN
            NLITE = NLITE + 1
            KTDIF(NLITE) = J
            EPSJ = EPS(J) * BOLTZ
            SIGJ = SIG(J) * 1.0E-8
            POLJ = POL(J) * 1.0E-24
            POLJST = POLJ / SIGJ**3
C
            DO 1000 K = 1, KK
               EPSK = EPS(K) * BOLTZ
               SIGK = SIG(K) * 1.0E-8
               DIPK = DIP(K) * 1.0E-18
               DIPKST = DIPK / SQRT(EPSK*SIGK**3)
               EKOEJ = EPSK / EPSJ
               TSE = 1.0 + 0.25*POLJST*DIPKST**2*SQRT(EKOEJ)
               EOK = TSE**2 * SQRT(EPS(J)*EPS(K))
               WTKJ= (WT(K)-WT(J)) / (WT(K)+WT(J))
C
               DO 500 N = 1, NT
                  TSLOG = LOG(AT(N) / EOK)
                  CALL HORNER (7, FITAST, TSLOG, ASTAR)
                  CALL HORNER (7, FITBST, TSLOG, BSTAR)
                  CALL HORNER (7, FITCST, TSLOG, CSTAR)
                  TD(N) = 7.5 * WTKJ * (2.0*ASTAR + 5.0) *
     1                    (6.0*CSTAR - 5.0)/
     2               (ASTAR * (16.0*ASTAR - 12.0*BSTAR + 55.0))
  500          CONTINUE
C
               FITWT(1) = -1.0
               EPSL = EPSIL
C
C*****precision > double
               CALL DPOLFT (NT, AT, TD, FITWT, NO-1, NORD, EPSL,
     1                      FITRES, IERR, FITWRK)
C*****END precision > double
C*****precision > single
C               CALL POLFIT (NT, AT, TD, FITWT, NO-1, NORD, EPSL,
C     1                      FITRES, IERR, FITWRK)
C*****END precision > single
C
               EMAXTD = MAX (EMAXTD, EPSL)
               IF (IERR .NE. 1) THEN
                  WRITE (NOUT, 7000) J,K
                  STOP
               ENDIF
C
               CCC = 0.0
C
C*****precision > double
               CALL DPCOEF (NORD, CCC, COFTD(1,K,NLITE), FITWRK)
C*****END precision > double
C*****precision > single
C               CALL PCOEF (NORD, CCC, COFTD(1,K,NLITE), FITWRK)
C*****END precision > single
C
 1000       CONTINUE
         ENDIF
 1100 CONTINUE
 7000 FORMAT (10X, 'ERROR IN FITTING THE ', I3, ',', I3,
     1                ' THERMAL DIFFUSION RATIO')
      RETURN
      END
C
      SUBROUTINE HORNER (N, A, X, P)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION A(*)
      NM1 = N-1
      B = A(N)
      DO 10 I = 1, NM1
         B = A(N-I) + B*X
   10 CONTINUE
      P = B
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION D12
     1                (W1, W2, T, EPS1, EPS2, SIG1, SIG2, DST)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      REAL FUNCTION D12 (W1, W2, T, EPS1, EPS2, SIG1, SIG2, DST)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C        CORRECTION OF D(1,2).
C         REFERENCE: MARRERO AND MASON,J.PHYS CHEM REF DAT.1,3(1972)
C
      SUMW = W1+W2
      SIG12 = 0.5 * (SIG1 + SIG2)
      TR11 = T / EPS1
      TR22 = T / EPS2
      TR12 = T / SQRT(EPS1*EPS2)
C
      CALL OMEGXX (2, TR11, OM2F11, DST)
      CALL OMEGXX (2, TR22, OM2F22, DST)
      CALL OMEGXX (1, TR12, OM1F12, DST)
      CALL OMEGXX (2, TR12, OM2F12, DST)
      CALL OMEGXX (3, TR12, BST12,  DST)
      CALL OMEGXX (4, TR12, CST12,  DST)
      AST12 = OM2F12 / OM1F12
C
      H = (SIG1/SIG2)**2 * SQRT(2.0*W2/SUMW) * 2.0*W1**2/(SUMW*W2)
      P1 = H * OM2F11/OM1F12
C
      H = (SIG2/SIG1)**2 * SQRT(2.0*W1/SUMW) * 2.0*W2**2/(SUMW*W1)
      P2 = H * OM2F22/OM1F12
C
      P12 = 15.0 * ((W1-W2)/SUMW)**2 + 8.0*W1*W2*AST12/SUMW**2
C
      H = 2.0/(W2*SUMW) * SQRT(2.0*W2/SUMW) * OM2F11/OM1F12 *
     1    (SIG1/SIG12)**2
      Q1 = ((2.5-1.2*BST12)*W1**2 +
     1                           3.0*W2**2 + 1.6*W1*W2*AST12)*H
C
      H = 2.0/(W1*SUMW) * SQRT(2.0*W1/SUMW) * OM2F22/OM1F12 *
     1    (SIG2/SIG12)**2
      Q2 = ((2.5-1.2*BST12)*W2**2 +
     1                           3.0*W1**2 + 1.6*W1*W2*AST12)*H
C
      H = ((W1-W2)/SUMW)**2 * (2.5-1.2*BST12)*15.0 +
     1      4.0*W1*W2*AST12/SUMW**2 * (11.0 - 2.4*BST12)
      Q12 = H + 1.6*SUMW*OM2F11*OM2F22 / (SQRT(W1*W2)*OM1F12**2) *
     1      SIG1**2*SIG2**2/SIG12**4
      D12 = ((6.0*CST12-5.0)**2/10.0) * (P1+P2+P12) /
     1      (Q1+Q2+Q12) + 1.0
C
      RETURN
      END
C
      SUBROUTINE PARKER (T, EPS, P)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C        TEMPERATURE DEPENDENCE OF THE ROTATIONAL COLLISION NUMBERS
C
C         REF: BRAU, C.A., JONKMAN, R.M., "CLASSICAL THEORY
C            OF ROTATIONAL RELAXATION IN DIATOMIC GASES",
C            JCP,VOL52,P477(1970)
C
      DATA PI32O2/2.7842/, P2O4P2/4.4674/, PI32/5.5683/
      D  = EPS / T
      DR = EPS / 298.0
      P = (1.0 + PI32O2*SQRT(DR) + P2O4P2*DR + PI32*DR**1.5) /
     1    (1.0 + PI32O2*SQRT(D)  + P2O4P2*D  + PI32*D**1.5)
      RETURN
      END
C
      SUBROUTINE INTERP (ARG, X, Y, VAL)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// QUADRATIC INTERPOLATION //////////////////////////////////////////
C
C
      DIMENSION X(3),Y(3)
      VAL1 = Y(1) + (ARG-X(1))*(Y(2)-Y(1)) / (X(2)-X(1))
      VAL2 = Y(2) + (ARG-X(2))*(Y(3)-Y(2)) / (X(3)-X(2))
      FAC1 = (ARG-X(1)) / (X(2)-X(1)) / 2.0
      FAC2 = (X(3)-ARG) / (X(3)-X(2)) / 2.0
      IF (ARG .GE. X(2)) THEN
         VAL = (VAL1*FAC2+VAL2) / (1.0+FAC2)
      ELSE
         VAL = (VAL1+VAL2*FAC1) / (1.0+FAC1)
      ENDIF
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION OMEG12 (N, TR, DR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      REAL FUNCTION OMEG12 (N, TR, DR)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// CALC. OF COLLISION INTEGRALS FOR A KRIEGER-12-6-3-POTENTIAL //////
C
      DIMENSION VERT(3), ARG(3), VAL(3), TSTERN(37), DELTA(8), O(37,8),
     1          P(37,8)
      DATA TSTERN/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,2.5,
     1            3.,3.5,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,25.,
     2            30.,35.,40.,50.,75.,100./
      DATA DELTA/0.,.25,.5,.75,1.,1.5,2.,2.5/
      DATA ((O(I,L),L=1,8),I=1,10) /
     1   4.008 , 4.002 , 4.655 , 5.52  , 6.454 , 8.214 , 9.824 ,11.31  ,
     2   3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160 ,
     3   2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483 ,
     4   2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539 ,
     5   2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918 ,
     6   1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747 ,
     7   1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137 ,
     8   1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872 ,
     9   1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657 ,
     *   1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780/
      DATA ((O(I,L),L=1,8),I=11,20) /
     1   1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199 ,
     2   1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99  ,
     3   1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827 ,
     4   1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7   ,
     5   1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592 ,
     6   1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397 ,
     7    .95  ,  .9538,  .9656,  .9852, 1.012 , 1.082 , 1.168 , 1.265 ,
     8    .9131,  .9162,  .9256,  .9413,  .9626, 1.019 , 1.09  , 1.17  ,
     9    .8845,  .8871,  .8948,  .9076,  .9252,  .972 , 1.03  , 1.098 ,
     *    .8428,  .8446,  .850 ,  .859 ,  .8716,  .9053,  .9483,  .9984/
      DATA ((O(I,L),L=1,8),I=21,30) /
     1    .813 ,  .8142,  .8183,  .825 ,  .8344,  .8598,  .8927,  .9316,
     2    .7898,  .791 ,  .794 ,  .7993,  .8066,  .8265,  .8526,  .8836,
     3    .7711,  .772 ,  .7745,  .7788,  .7846,  .8007,  .822 ,  .8474,
     4    .7555,  .7562,  .7584,  .7619,  .7667,  .78  ,  .7976,  .8189,
     5    .7422,  .743 ,  .7446,  .7475,  .7515,  .7627,  .7776,  .796 ,
     6    .72022, .7206,  .722 ,  .7241,  .7271,  .7354,  .7464,  .76  ,
     7    .7025,  .703 ,  .704 ,  .7055,  .7078,  .7142,  .7228,  .7334,
     8    .68776, .688,   .6888,  .6901,  .6919,  .697 ,  .704 ,  .7125,
     9    .6751,  .6753,  .676 ,  .677 ,  .6785,  .6827,  .6884,  .6955,
     *    .664 ,  .6642,  .6648,  .6657,  .6669,  .6704,  .6752,  .681 /
      DATA ((O(I,L),L=1,8),I=31,37) /
     1    .6414,  .6415,  .6418,  .6425,  .6433,  .6457,  .649 ,  .653 ,
     2    .6235,  .6236,  .6239,  .6243,  .6249,  .6267,  .629 ,  .632 ,
     3    .60882, .6089,  .6091,  .6094,  .61  ,  .6112,  .613 ,  .6154,
     4    .5964,  .5964,  .5966,  .597 ,  .5972,  .5983,  .600 ,  .6017,
     5    .5763,  .5763,  .5764,  .5766,  .5768,  .5775,  .5785,  .58  ,
     6    .5415,  .5415,  .5416,  .5416,  .5418,  .542 ,  .5424,  .543 ,
     7    .518 ,  .518 ,  .5182,  .5184,  .5184,  .5185,  .5186,  .5187/
C
      DATA ((P(I,L),L=1,8),I=1,10) /
     1   4.1   , 4.266 , 4.833 , 5.742 , 6.729 , 8.624 ,10.34  ,11.890 ,
     2   3.263 , 3.305 , 3.516 , 3.914 , 4.433 , 5.57  , 6.637 , 7.618 ,
     3   2.84  , 2.836 , 2.936 , 3.168 , 3.511 , 4.329 , 5.126 , 5.874 ,
     4   2.531 , 2.522 , 2.586 , 2.749 , 3.004 , 3.64  , 4.282 , 4.895 ,
     5   2.284 , 2.277 , 2.329 , 2.46  , 2.665 , 3.187 , 3.727 , 4.249 ,
     6   2.084 , 2.081 , 2.13  , 2.243 , 2.417 , 2.862 , 3.329 , 3.786 ,
     7   1.922 , 1.924 , 1.97  , 2.072 , 2.225 , 2.641 , 3.028 , 3.435 ,
     8   1.7902, 1.795 , 1.84  , 1.934 , 2.07  , 2.417 , 2.788 , 3.156 ,
     9   1.682 , 1.689 , 1.733 , 1.82  , 1.944 , 2.258 , 2.596 , 2.933 ,
     *   1.593 , 1.60  , 1.644 , 1.725 , 1.84  , 2.124 , 2.435 , 2.746 /
      DATA ((P(I,L),L=1,8),I=11,20) /
     1   1.455 , 1.465 , 1.504 , 1.574 , 1.67  , 1.913 , 2.181 , 2.45  ,
     2   1.355 , 1.365 , 1.4   , 1.461 , 1.544 , 1.754 , 1.989 , 2.228 ,
     3   1.28  , 1.289 , 1.321 , 1.374 , 1.447 , 1.63  , 1.838 , 2.053 ,
     4   1.222 , 1.231 , 1.26  , 1.306 , 1.37  , 1.532 , 1.718 , 1.912 ,
     5   1.176 , 1.184 , 1.209 , 1.25  , 1.307 , 1.45  , 1.618 , 1.795 ,
     6   1.0933, 1.1   , 1.119 , 1.15  , 1.193 , 1.304 , 1.435 , 1.578 ,
     7   1.039 , 1.044 , 1.06  , 1.083 , 1.117 , 1.204 , 1.31  , 1.428 ,
     8    .9996, 1.004 , 1.016 , 1.035 , 1.062 , 1.133 , 1.22  , 1.32  ,
     9    .9699,  .9732,  .983 ,  .9991, 1.021 , 1.08  , 1.153 , 1.236 ,
     *    .9268,  .9291,  .936 ,  .9473,  .9628, 1.005 , 1.058 , 1.12  /
      DATA ((P(I,L),L=1,8),I=21,30) /
     1    .8962,  .8979,  .903 ,  .9114,  .923 ,  .9545,  .9955, 1.044 ,
     2    .8727,  .8741,  .878 ,  .8845,  .8935,  .918 ,  .9505,  .9893,
     3    .8538,  .8549,  .858 ,  .8632,  .8703,  .890 ,  .9164,  .9482,
     4    .8379,  .8388,  .8414,  .8456,  .8515,  .868 ,  .8895,  .916 ,
     5    .8243,  .8251,  .8273,  .8308,  .8356,  .8493,  .8676,  .89  ,
     6    .8018,  .8024,  .8039,  .8065,  .810 ,  .820 ,  .8337,  .8504,
     7    .7836,  .784 ,  .7852,  .7872,  .7899,  .7976,  .808 ,  .8212,
     8    .7683,  .7687,  .7696,  .771 ,  .7733,  .7794,  .788 ,  .7983,
     9    .7552,  .7554,  .7562,  .7575,  .7592,  .764 ,  .771 ,  .7797,
     *    .7436,  .7438,  .7445,  .7455,  .747 ,  .7512,  .757 ,  .7642/
      DATA ((P(I,L),L=1,8),I=31,37) /
     1    .71982, .72  ,  .7204,  .7211,  .7221,  .725 ,  .7289,  .7339,
     2    .701 ,  .7011,  .7014,  .702 ,  .7026,  .7047,  .7076,  .7112,
     3    .68545, .6855,  .686 ,  .686 ,  .6867,  .6883,  .6905,  .693 ,
     4    .6723,  .6724,  .6726,  .673 ,  .6733,  .6745,  .676 ,  .6784,
     5    .651 ,  .651 ,  .6512,  .6513,  .6516,  .6524,  .6534,  .6546,
     6    .614 ,  .614 ,  .6143,  .6145,  .6147,  .6148,  .6148,  .6147,
     7    .5887,  .5889,  .5894,  .59  ,  .5903,  .5901,  .5895,  .5885/
C
      IF (DR.LT.-.00001 .OR. DR.GT.2.5 .OR. TR.LT..09 .OR. TR.GT.500.
     1 .OR. (ABS(DR).GT.1.0E-5 .AND. TR.GT.75.0)) THEN
         WRITE (6,'(A)') 'COLLISION INTEGRAL UNDEFINED'
         RETURN
      ENDIF
C
      IF (TR .GT. 75.0)  THEN
         IF (N .EQ. 1) THEN
            OMEG12 = 0.623 - 0.136E-2*TR + 0.346E-5*TR*TR
     1                         -0.343E-8*TR*TR*TR
         ELSEIF (N .EQ. 2) THEN
            OMEG12 = 0.703 - 0.146E-2*TR + 0.357E-5*TR*TR
     1                         -0.343E-8*TR*TR*TR
         ENDIF
         RETURN
      ENDIF
C
      IF (TR .LE. 0.2) THEN
         II = 2
      ELSE
         II = 37
         DO 30 I = 2, 36
            IF (TSTERN(I-1).LT.TR .AND. TSTERN(I).GE.TR)  II=I
   30    CONTINUE
      ENDIF
C
      IF (ABS(DR) .GE. 1.0E-5)  THEN
         IF (DR .LE. 0.25) THEN
            KK = 2
         ELSE
            KK = 7
            DO 10 K = 2, 7
               IF (DELTA(K-1).LT.DR .AND. DELTA(K) .GE. DR)  KK=K
   10       CONTINUE
         ENDIF
         DO 50 I = 1, 3
            DO 40 K = 1, 3
               ARG(K) = DELTA(KK-2+K)
               IF (N .EQ. 1)  THEN
                  VAL(K) = O(II-2+I, KK-2+K)
               ELSEIF (N .EQ. 2)  THEN
                  VAL(K) = P(II-2+I, KK-2+K)
               ENDIF
   40       CONTINUE
            CALL INTERP (DR, ARG, VAL, VERT(I))
   50    CONTINUE
         DO 60 I = 1, 3
            ARG(I) = TSTERN(II-2+I)
   60    CONTINUE
         CALL INTERP (TR, ARG, VERT, OMEG12)
C
      ELSE
         DO 80 I = 1, 3
            ARG(I) = TSTERN(II-2+I)
            IF (N .EQ. 1) THEN
               VAL(I) = O(II-2+I, 1)
            ELSEIF (N .EQ. 2) THEN
               VAL(I) = P(II-2+I, 1)
            ENDIF
   80    CONTINUE
         CALL INTERP (TR, ARG, VAL, OMEG12)
      ENDIF
      RETURN
      END
C
      SUBROUTINE OMEGXX (N, TR, OM, DST)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// CALCULATION OF THE REDUCED COLLISION INTEGRALS ///////////////////
C
      IF (N.EQ.1 .OR. N.EQ.2) THEN
         OM = OMEG12 (N, TR, DST)
      ELSEIF (N .EQ. 3) THEN
         IF (TR .LE. 5.0) THEN
            OM = 1.36 - 0.223*TR + 0.0613*TR**2 -0.00554*TR**3
         ELSE
            OM = 1.095
         ENDIF
      ELSEIF (N .EQ. 4) THEN
         IF (TR .LE. 5.0) THEN
            OM =0.823 + 0.0128*TR + 0.0112*TR**2 -0.00193*TR**3
         ELSE
            OM = .9483
         ENDIF
      ELSE
         WRITE (6,'(10X,A)') 'OMEGXX IS CALLED WITH WRONG PARAMETERS'
      ENDIF
      RETURN
      END

