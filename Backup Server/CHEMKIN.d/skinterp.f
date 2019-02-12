C     processor (interpreter) for the CHEMKIN collection of codes.
C@(#)
C@(#)
C@(#)                   FILE =  skinterp.f
C@(#)
C@(#)  ---------------  VERSION = 4.5
C@(#)  |  SCCS  FILE |
C@(#)  |   SUMMARY   |  CURRENT CHECKOUT DATE = 08/10/94
C@(#)  ---------------                           at 16:40:04
C@(#)                   DATE OF NEWEST DELTA = 08/10/94
C@(#)                                            at 16:40:03
C@(#)  SCCS file name = /users/chemkin/SCCS/s.skinterp.f
C@(#)===================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      PROGRAM SKINTP
C
C     SKINTERP is the surface phase, symbolic, chemical mechanism pre-
C     processor (interpreter) for the CHEMKIN collection of codes.
C     (The acronym CHEMKIN is a registered Trademark).
C
C----------------------------------------------------------------------C
C
C///////////////////////////////////////////////////////////////////
C
C            SKINTERP: SURFACE KINETICS INTERPRETER
C                      VERSION 4.5
C
C     WRITTEN BY:
C         FRAN M. RUPLEY
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3657
C
C       Copyright 1990, Sandia Corporation.
C       The U.S. Goverment retains a limited license in this
C       software.
C
C       The U.S. Government retains, in this software, a paid-up,
C       nonexclusive, irrevocable worldwide license to reproduce,
C       prepare derivative works, perform publicly and display
C       publicly by or for the Government, including the right to
C       distribute to other Government contractors.
C
C       Neither the United States, the U.S. Dept. of Energy, nor
C       any of their employees, makes any warranty, express or
C       implied, or assumes any legal liability or responsibility
C       for the accuracy, completeness, or usefulness of any
C       information, apparatus, product, or process disclosed, or
C       represents that its use would not infringe privately owned
C       rights.
C
C/////////////////////////////////////////////////////////////////////
C
C     VERSION 4.3
C       CHANGES SINCE VERSION 1.0
C         1. Replace "REAL*8" with "DOUBLE PRECISION"
C       CHANGES SINCE VERSION 1.1
C         1. Include thermodynamic properties for all species
C         2. Allow reversible/irreversible surface reactions
C         3. Allow reaction species to end with '=' or '-'
C         4. Allow real values of elemental composition in THERMO cards
C         5. Allow upper/upper case input
C       CHANGES SINCE VERSION 1.3
C         1. Change name to SINTRP to conform to ANSI standard
C         2. Change MAXSITE to MXSITE to conform to ANSI standard
C       CHANGES SINCE VERSION 1.4
C         1. Add "unix" change blocks
C       CHANGES SINCE VERSION 1.5
C         1. Replace DEPOSIT species with MIXSOL and PURE species
C       CHANGES TO VERSION 3.2
C         1. Eliminate MIXSOL and PURE in favor of BULK.
C         2. Implement slash-delimited options.
C       CHANGES TO VERSION 3.3
C         1. Add auxiliary reaction keyword "STICK" for sticking
C            coefficients, in which case there must be one gas-phase
C            reactant.
C       CHANGES TO VERSION 3.4
C         1. 6 reactants, 6 products
C         2. 'NONCON'servation of sites option
C         3. default phase names
C         4. continuation character '&' for reaction lines
C         5. phase names must be unique
C         6. site-phase species names may not duplicate gas-phase
C            species names, but MAY duplicate other site-phase species
C            names not in same site phase.
C         7. bulk-phase species names may not duplicate either gas-
C            or site-phase species names, but MAY duplicate other
C            bulk-phase species names not in same bulk phase.
C         8. site densities are now moles/cm**2 instead of number
C            densities
C         CHANGES FOR VERSION 3.5
C         1. Add IKCOV array of species for reactions with coverage
C            parameters
C         2. Add integer NIICON, the number of reactions which do
C            not conserve sites, to first record of binary file
C         3. Add first record to binary file, which contains a
C            character*16 string VERS to name the version number of
C            the binary file, and PRES to name its precision
C            (DOUBLE or SINGLE)
C         4. Correction to SKBULK and SKSURF to initialize KNUM2=0.
C         5. Check that pre-exponential factor is positive.
C         CHANGES TO VERSION 3.6
C         1. Correct unit conversions
C         2. Bring up to date with latest version of manual
C         CHANGES TO VERSION 3.61
C         1. Error if site phase or bulk phas has no species declared
C         2. Change SKUNIT to parse LINE instead of SUB(*) to
C            correct misinterpretation of unit strings containing
C            a slash
C         3. Increase length of real work space by NPHASE in order
C            to have scratch space for phases
C         4. Modify SKTHRM such that if a species occurs in more than
C            one site or bulk phase, thermodynamic data will be
C            stored for each occurrence.
C         5. With sticking coefficients, the one gas-phase reactant
C            must have a stoichiometric coefficient of one.
C         6. A site must have a density delcaration.
C         CHANGES TO VERSION 3.62
C         1. Need NU in argument list for SKAUXL
C         CHANGES TO VERSION 3.63
C         1. Set KERR=.TRUE. if coefficient greater than 1 for the
C            gas phase species in a reaction with sticking coeff.
C         CHANGES TO VERSION 3.64
C         1. Modify reaction interpretation to allow continuation
C            lines as an array (thus avoiding CHARACTER*160).
C         CHANGES TO VERSION 3.7
C         1. INCF(*,*) previously calculated for sites only, is now
C            calculated for all phases (SUBROUTINE SKBAL).
C         CHANGES TO VERSION 3.71
C         1. Error in UPCASE was causing auxiliary keywords to be
C            ignored.
C         CHANGES TO VERSION 3.72
C         1. Indexing in scaling of 3rd Reverse Arrhenius parameter
C            storage was wrong, should be RPAR(3, NREV)*EFAC,
C            not RPAR(3, II).
C         CHANGES TO VERSION 3.73
C         1. Initialize logical variable LPDEN
C         CHANGES TO VERSION 3.74
C         1. Correct error similar to V3.72 for RPAR(1, NREV)*AFAC,etc.
C         CHANGES TO VERSION 3.75
C         1. Additional thermodynamic checks:
C            a) if TLO,TMID,THI not given on THERMO CARDS, use values
C               from THERMO.DAT
C            b) Check TLO < THI, TLO <= TMID, TMID <= THI
C         CHANGES TO VERSION 3.76
C         1. Need to get TLO,THI,TMID from Thermodynamic database
C            BEFORE reading user's THERMO data
C         CHANGES TO VERSION 3.77
C         1. Change "LINE" to "ILINE" for reading Thermodynamic
C            database.
C         CHANGES TO VERSION 3.78
C         1. Increase length LENRSK by NIISUR+NIIREV for perturbation
C            factor.
C         CHANGES TO MAKE VERSION 4.0
C         1. Change units on surface-coverage modification of the the
C            rate of progress (concentration units --> site fractions)
C         2. Add an extra NIISUR to length of real work array required
C            to account for new multiplicative factor (EQFAC) in SKLIB
C         CHANGES TO MAKE VERSION 4.01
C         (7/11/91 F. Rupley per M. Coltrin)
C         1. Additional error checking:
C            no species on a site or bulk
C            site occupancy numbers > 0
C            site density required and > 0
C            bulk density > 0 if given
C            lack of Arhennius coefficients
C            duplicate phase names
C         CHANGES TO MAKE VERSION 4.02
C         (7/17/91 F. Rupley per M. Coltrin)
C         1. SKPRNT was checking sites, bulks, and species and
C            printing error messages for an empty surface mechanism;
C            corrected to skip these checks if there are no sites
C            or bulks.
C         CHANGES TO MAKE VERSION 4.03
C         1. TMID must be passed to SKSPEC
C         CHANGES FOR VERSION 4.04 (4/13/92 F. Rupley per M. Coltrin)
C         1. Correct logic for SKDUP, requires additional argument.
C         CHANGES FOR VERSION 4.05 (1/26/94 F. Rupley per R. Kee)
C         1. Allow real stoichometric coefficients; NIIRNU, IRNU(*),
C            RNU(MAXSPR,*)
C        CHANGES FOR VERSION 4.06 (3/15/94 F. Rupley)
C        1.  DOS/PC compatibility effort includes adding filenames to
C            OPEN statements, removing unused variables in CALL lists,
C            unusued but possibly initialized variables.
C        CHANGES FOR VERSION 4.07 (4/19/94 F. Rupley)
C        1.  correct indexing in SKBAL, SKRBAL
C        CHANGES FOR VERSION 4.08 (4/29/94 F. Rupley)
C        1.  Cannot change RORD for an irreversible reaction.
C        CHANGES FOR VERSION 4.08b (5/20/94 F. Rupley per E. Meeks)
C        1.  Incorporate plasma options.
C        CHANGES FOR VERSION 4.08c (6/3/94 F. Rupley per H. Moffat)
C        1.  Add error checks/messages upon opening therm.dat,
C            chem.bin, surf.inp
C        2.  Allow comments (!) in thermo data
C        3.  Correct phase name logic
C        CHANGES FOR VERSION 4.09 (6/28/94 E. Meeks)
C        1.  Add REACTION line keyword MWOFF to turn off Motz-Wise
C            correction to sticking coefficient formulation.
C            Addition integer flag MOTZWS passed through linking file.
C        CHANGES FOR VERSION 4.10 (7/14/94 F. Rupley)
C        1.  Implement MATERIAL keyword to enable multiple materials
C        2.  Subroutine SKSET to initialize values
C        CHANGES FOR VERSION 4.2
C        1.  Correct LENISK (NIIBOHM > NIIBHM)
C        CHANGES FOR VERSION 4.3
C        1.  Lengthen ISKWRK by 3
C        CHANGES FOR VERSION 4.31 (8/10/94 H. Moffat)
C        1.  Changed physical constants to conform to 1986 CODATA
C            recommendations
C        CHANGES FOR VERSION 4.5 (10/3/94 F. Rupley)
C        1.  Add NIIORD to argument list for SKSET and initialize to 0.
C/////////////////////////////////////////////////////////////////////
C
C     The surface chemistry interpreter reads a CHEMKIN binary file
C     LINC, containing information about the gas-phase chemistry for
C     a problem, then reads and interprets a surface mechanism data
C     file and creates a binary file LINKSK for the surface chemistry.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     INCLUDE 'sk.cmn'
C.....Ashish modification: LRWK and LIWK increased for large mechs.
C.....Ashish modification: KDIM increased for large mechs.
C.....JE Sutton modification: KDIM and IDIM increased again
C                             KDIM = Max # surface species
C                             IDIM = Max # surface reactions
C.....JE Sutton modification: LRWK and LIWK increased again
C
      PARAMETER (LRWK =50000, LIWK =50000, LCWK =1000, SKMIN=1.0E-3)
C
      PARAMETER (MDIM=10, KDIM=2000, IDIM=10000, NSPAR=3, NPC=5,
     1           MAXFIT=NPC+2, MAXTP=3, NTR=MAXTP-1, MAXSPR=12,
     2           NSCOV=3, MXPHSE=10, MAXORD=10, LIN=35, LOUT=36,
     3           LTHRM=17, LINC=25, LINKSK=26)
C
      DIMENSION WORK(LRWK), IWORK(LIWK), NT(KDIM), TMP(MAXTP,KDIM),
     *          A(MAXFIT, NTR, KDIM), KCHRG(KDIM), KPHSE(KDIM),
     1          AWT(MDIM), WTM(KDIM), DEN(KDIM), KNCF(MDIM,KDIM),
     2          KCOV(KDIM), KFIRST(MXPHSE), KLAST(MXPHSE),
     3          KKPHAS(MXPHSE), PDEN(MXPHSE), NR(IDIM),
     4          NU(MAXSPR,IDIM), NUNK(MAXSPR,IDIM),
     5          RNCF(MXPHSE,IDIM), NREAC(IDIM), NUSUMK(IDIM),
     6          SPAR(NSPAR,IDIM), IDUP(IDIM), IICOV(IDIM),
     7          IKCOV(IDIM), CPAR(NSCOV,IDIM), IIREV(IDIM),
     8          RSPAR(NSPAR,IDIM), IISTK(IDIM), IRNU(IDIM),
     9          RNU(MAXSPR,IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     *          RORD(MAXORD,IDIM), IBOHM(IDIM), IBK(IDIM),
     1          IBT(IDIM)
C
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM), CWORK(LCWK),
     1             PNAME(MXPHSE), VERS, PREC, MATNAM
C     CHARACTER infile*50
      LOGICAL KERR, ITHRM(KDIM), NONCON, LPDEN(MXPHSE),
     1        LKDEN(KDIM), LMTZWS
C
C----------------------------------------------------------------------C
C
C ### J.M: define names for input-output files
C     infile = 'allin'
C     CALL Filnam(infile) 
      OPEN (LIN, FORM='FORMATTED', STATUS='OLD', FILE='INP.d/surf.inp',
     1      ERR=22222)
      READ (LIN,'(A)',END=22222)
      REWIND (LIN)
C
      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN', 
     1   FILE='OUT.d/surf.out')
      OPEN (LINKSK, FORM='UNFORMATTED', STATUS='UNKNOWN',
     1      FILE='INP.d/sklink')
C
      VERS = '4.5'
      WRITE (LOUT, 15) VERS(:5)
   15 FORMAT(/' SURFACE INTERPRETER OUTPUT: ',
     *       /' Copyright 1990, Sandia Corporation.',
     1       /' The U.S. Government retains a limited license',
     2        ' in this software.',
     3       /' CHEMKIN-II Version ',A,' January, 1995'
C*****precision > double
     4       /' DOUBLE PRECISION')
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C     4       /' SINGLE PRECISION')
C      PREC = 'SINGLE'
C*****END precision > single
C----------------------------------------------------------------------C
C
      NMAT = 0
  111 CONTINUE
C
C     Initialize variables
C
      NMAT = NMAT + 1
      CALL SKSET (NMAT, LRWK, WORK, LIWK, IWORK, KDIM, NT, MAXTP, TMP,
     1                  MAXFIT, NTR, A, KCHRG, KPHSE, AWT, WTM, DEN,
     2                  MDIM, KNCF, KCOV, MXPHSE, KFIRST, KLAST, KKPHAS,
     3                  PDEN, IDIM, NR, MAXSPR, NU, NUNK, RNCF, NREAC,
     4                  NUSUMK, NSPAR, SPAR, IDUP, IICOV, IKCOV, NSCOV,
     5                  CPAR, IIREV, RSPAR, IISTK, IRNU, RNU, IORD,
     6                  MAXORD, KORD, RORD, IBOHM, IBK, IBT, KERR,
     7                  ITHRM, NONCON, LPDEN, LKDEN, LMTZWS, KNAME,
     8                  ENAME, LCWK, CWORK, PNAME, NPHASE, NNSUR,
     9                  NNBLK, NKKGAS, NKKSUR, NKKBLK, NKKTOT, NIISUR,
     *                  NIICOV, NIIREV, NIISTK, NIICON, NIIRNU, NIIORD,
     1                  NIIBHM)
C
      IF (NMAT .EQ. 1) THEN
C
C        Initialize gas-phase chemistry
C
         OPEN (LINC, FORM='UNFORMATTED', STATUS='OLD',
     1                 FILE='INP.d/cklink', ERR=11111)
         READ (LINC, END=11111)
c	 write(*,*)'i read cklink'
         REWIND (LINC)
         CALL SKSTRT (LIWK, LRWK, LCWK, LINC, LOUT, MDIM, KDIM,
     1                MAXTP, MAXFIT, NTR, IWORK, WORK, CWORK, NELEM,
     2                AWT, ENAME, NKKGAS, KNAME, WTM, ITHRM, KCHRG,
     3                KPHSE, NT, TMP, KNCF, NFIT, A, KERR)
         CLOSE (LINC)
      ENDIF
C
C     Interpret surface-phase mechanism
C
      CALL SKKEY  (LIN, LTHRM, LOUT, MDIM, NELEM, ENAME, AWT, KDIM,
     1             KNAME, ITHRM, WTM, KNCF, NONCON, NIICON, KCOV,
     2             KPHSE, KCHRG, MAXTP, NT, NTR, TMP, MAXFIT, A,
     3             MXPHSE, NKKGAS, IDIM, NSPAR, MAXSPR, NSCOV, NPHASE,
     4             PNAME, PDEN, LPDEN, KFIRST, KLAST, KKPHAS, NFSUR,
     5             NLSUR, NNSUR, NKKSUR, NFBLK, NLBLK, NNBLK, NKKBLK,
     6             NKKTOT, DEN, LKDEN, NIISUR, SPAR, NR, NU, NUNK,
     7             RNCF, NREAC, NUSUMK, IDUP, NIIREV, IIREV, RSPAR,
     8             NIICOV, IICOV, IKCOV, CPAR, NIISTK, IISTK, KERR,
     9             NIIRNU, IRNU, RNU, SKMIN, NIIORD, MAXORD, IORD, KORD,
     *             RORD, NIIBHM, IBOHM, IBK, IBT, LMTZWS, NMAT,
     1             MATNAM, MORE)
C
C------version and precision of surface binary file
C
      WRITE (LINKSK) VERS, PREC, KERR, MORE
C
      IF (KERR) THEN
         WRITE (LOUT, '(/A)')
     1   ' WARNING...THERE IS AN ERROR IN THE SURFACE LINKING FILE'
         CLOSE (LINKSK)
         CLOSE (LOUT)
         STOP
      ENDIF
C
      LENISK = (4 + NELEM)*NKKTOT + 3*NPHASE + (3 + 2*MAXSPR)*NIISUR
     1       + 2*NIICOV + NIIREV + NIISTK + NIIRNU
     2       + NIIORD*(1 + MAXORD) + 3*NIIBHM + 75
      LENRSK = 4 + NELEM + NKKTOT*(4 + MAXTP + NFIT*(MAXTP-1))
     1           + 2*NPHASE + NSPAR*(NIISUR + NIIREV + NIISTK)
     2           + (7+NPHASE)*NIISUR + NSCOV*NIICOV + NIIREV
     3           + NIIRNU*MAXSPR + NIIORD*MAXORD
      LENCSK = NELEM + NKKTOT + NPHASE + 1
C
      NCP   = NFIT - 2
      MOTZWS = 0
      IF (LMTZWS) MOTZWS = 1
C
C------integer constants
C
      WRITE (LINKSK)  LENISK, LENRSK, LENCSK, MAXSPR, MAXTP, NCP,
     1                NELEM,  NKKGAS, NKKSUR, NKKBLK, NKKTOT,
     2                NPHASE, NFSUR,  NLSUR,  NNSUR,  NFBLK, NLBLK,
     3                NNBLK,  NIISUR, NSPAR,  NSCOV,  NIICOV, NIIREV,
     4                NIISTK, NIICON, NIIBHM, NIIRNU, NIIORD,
     5                MAXORD, SKMIN,  MOTZWS
C
C------element data
      WRITE (LINKSK) (ENAME(M), AWT(M), M=1,NELEM )
C------species data
      WRITE (LINKSK) (KNAME(K), WTM(K), KPHSE(K), KCHRG(K), NT(K),
     1               KCOV(K), DEN(K), (KNCF(M,K),M=1,NELEM),
     2               (TMP(L,K),L=1,MAXTP),
     3               ( (A(M,L,K), M=1,NFIT), L=1,NTR), K=1,NKKTOT)
C-------phase data
      WRITE (LINKSK) MATNAM, (PNAME(N), KFIRST(N), KLAST(N), KKPHAS(N),
     1               PDEN(N), (RNCF(N,I),I=1,NIISUR), N=1, NPHASE)
C-------reaction data
      IF (NIISUR .GT. 0) THEN
         WRITE (LINKSK)
     1   (NR(I), NREAC(I), NUSUMK(I), (SPAR(N,I),N=1,NSPAR),
     2   (NUNK(N,I), NU(N,I), N=1,MAXSPR), I=1,NIISUR)
         IF (NIICOV .GT. 0) WRITE (LINKSK)
     1   (IICOV(I), IKCOV(I), (CPAR(N,I),N=1,NSCOV), I=1,NIICOV)
         IF (NIIREV .GT. 0) WRITE (LINKSK)
     1   (IIREV(I), (RSPAR(N,I),N=1,NSPAR), I=1,NIIREV)
         IF (NIISTK .GT. 0) WRITE (LINKSK) (IISTK(I), I=1,NIISTK)
         IF (NIIBHM .GT. 0) WRITE (LINKSK)
     1   (IBOHM(I), IBK(I), IBT(I), I = 1, NIIBHM)
         IF (NIIRNU .GT. 0) WRITE (LINKSK)
     1   (IRNU(I), (RNU(N,I), N=1,MAXSPR), I=1,NIIRNU)
         IF (NIIORD .GT. 0) WRITE (LINKSK)
     1   (IORD(I), (KORD(N,I), RORD(N,I), N=1,MAXORD), I=1,NIIORD)
C
      ELSE
         WRITE (LOUT, '(/A)')
     1       ' WARNING...NO SURFACE REACTIONS FOUND; ',
     2       ' SURFACE LINKING FILE HAS NO REACTION INFORMATION ON IT.'
      ENDIF
C
      WRITE (LOUT, '(/A)')
     1   ' NO ERRORS FOUND ON INPUT...SURFACE LINKING FILE WRITTEN.'
      WRITE (LOUT, '(/A,3(/A,I6))')
     1   ' WORKING SPACE REQUIREMENTS ARE',
     2   '    INTEGER:   ',LENISK,
     3   '    REAL:      ',LENRSK,
     4   '    CHARACTER: ',LENCSK
C
      IF (MORE .GT. 0) GO TO 111
C
      CLOSE (LIN)
      CLOSE (LOUT)
      CLOSE (LINKSK)
      STOP
C
11111 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read cklink...'
      CLOSE (LINC)
      STOP 2
C
22222 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read surf.inp...'
      CLOSE (LIN)
      STOP 2
C
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKSET (NMAT, LRWK, WORK, LIWK, IWORK, KDIM, NT, MAXTP,
     1                  TMP,MAXFIT, NTR, A, KCHRG, KPHSE, AWT, WTM, DEN,
     2                  MDIM, KNCF, KCOV, MXPHSE, KFIRST, KLAST, KKPHAS,
     3                  PDEN, IDIM, NR, MAXSPR, NU, NUNK, RNCF, NREAC,
     4                  NUSUMK, NSPAR, SPAR, IDUP, IICOV, IKCOV, NSCOV,
     5                  CPAR, IIREV, RSPAR, IISTK, IRNU, RNU, IORD,
     6                  MAXORD, KORD, RORD, IBOHM, IBK, IBT, KERR,
     7                  ITHRM, NONCON, LPDEN, LKDEN, LMTZWS, KNAME,
     8                  ENAME, LCWK, CWORK, PNAME, NPHASE, NNSURF,
     9                  NNBULK, KKGAS, KKSURF, KKBULK, KKTOT, NIISUR,
     *                  NIICOV, NIIREV, NIISTK, NIICON, NIIRNU, NIIORD,
     1                  NIIBHM)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C     SKSET set the initial values of the arrays and matrices
C
      DIMENSION WORK(*), IWORK(*), NT(*), TMP(MAXTP,*),
     *          A(MAXFIT, NTR, *), KCHRG(*), KPHSE(*), AWT(*), WTM(*),
     1          DEN(*), KNCF(MDIM,*), KCOV(*), KFIRST(*), KLAST(*),
     2          KKPHAS(*), PDEN(*), NR(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     3          RNCF(MXPHSE,*), NREAC(*), NUSUMK(*), SPAR(NSPAR,*),
     4          IDUP(*), IICOV(*), IKCOV(*), CPAR(NSCOV,*), IIREV(*),
     5          RSPAR(NSPAR,*), IISTK(*), IRNU(*), RNU(MAXSPR,*),
     6          IORD(*), KORD(MAXORD,*), RORD(MAXORD,*), IBOHM(*),
     7          IBK(*), IBT(*)
C
      LOGICAL KERR, ITHRM(*), NONCON, LPDEN(*), LKDEN(*), LMTZWS
      CHARACTER*16 KNAME(*), ENAME(*), CWORK(*), PNAME(*)
C
      IF (NMAT .EQ. 1) THEN
         KKGAS = 0
         DO 1 L = 1, LRWK
            WORK(L) = 0.0
    1    CONTINUE
         DO 2 L = 1, LIWK
            IWORK(L) = 0
    2    CONTINUE
         DO 3 L = 1, LCWK
            CWORK(L) = ' '
    3    CONTINUE
         DO 5 M = 1, MDIM
            AWT(M) = 0.0
            ENAME(M) = ' '
    5    CONTINUE
      ENDIF
C
      DO 100 K = KKGAS+1, KDIM
         KNAME(K) = ' '
         NT(K) = MAXTP
         DO 10 M = 1, MAXTP
            TMP(M,K) = -1.0
   10    CONTINUE
         DO 15 M = 1, MAXFIT
            DO 15 N = 1, NTR
               A(M, N, K) = 0.0
   15    CONTINUE
         KCHRG(K) = 0
         KPHSE(K) = 0
         WTM(K) = 0.0
         DEN(K) = -1.0
         DO 20 M = 1, MDIM
            KNCF(M,K) = 0
   20    CONTINUE
         KCOV(K) = 1
         ITHRM(K) = .FALSE.
         LKDEN(K) = .FALSE.
  100 CONTINUE
C
      DO 200 M = 1, MXPHSE
         PNAME(M) = ' '
         KFIRST(M) = 0
         KLAST(M) = 0
         KKPHAS(M) = 0
         PDEN(M) = -1.0
         LPDEN(M) = .FALSE.
  200 CONTINUE
C
      DO 300 I = 1, IDIM
         NR(I) = 0
         DO 205 M = 1, MAXSPR
            NU(M,I) = 0
            NUNK(M,I) = 0
  205    CONTINUE
         DO 210 M = 1, MXPHSE
            RNCF(M,I) = 0.0
  210    CONTINUE
         NREAC(I) = 0
         NUSUMK(I) = 0
         DO 215 N = 1, NSPAR
            SPAR(N,I) = 0.0
  215    CONTINUE
         IDUP(I) = 0
         IICOV(I) = 0
         IKCOV(I) = 0
         DO 220 N = 1, NSCOV
            CPAR(N,I) = 0.0
  220    CONTINUE
         IIREV(I) = 0
         DO 225 N = 1, NSPAR
            RSPAR(N,I) = 0.0
  225    CONTINUE
         IISTK(I) = 0
         IRNU(I) = 0
         DO 230 M = 1, MAXSPR
            RNU(M,I) = 0.0
  230    CONTINUE
         IORD(I) = 0
         DO 235 M = 1, MAXORD
            KORD(M,I) = 0.0
            RORD(M,I) = 0.0
  235    CONTINUE
         IBOHM(I) = 0
         IBK(I) = 0
         IBT(I) = 0
  300 CONTINUE
C
      KERR = .FALSE.
      LMTZWS = .TRUE.
      NONCON = .FALSE.
      NPHASE = 0
      NNSURF = 0
      NNBULK = 0
      KKSURF = 0
      KKBULK = 0
      KKTOT = 0
      NIISUR = 0
      NIICOV = 0
      NIIREV = 0
      NIISTK = 0
      NIICON = 0
      NIIRNU = 0
      NIIORD = 0
      NIIBHM = 0
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKSTRT (LIWK, LRWK, LCWK, LINC, LOUT, MDIM,
     1                   KDIM, MAXTP, MAXFIT, NTR, IWORK, WORK, CWORK,
     2                   NELEM, AWT, ENAME, KKGAS, KNAME, WTM, ITHRM,
     3                   KCHRG, KPHSE, NT, T, KNCF, NFIT, A, KERR)
C
C     SKSTRT initializes a CHEMKIN binary file to get the
C     gas-phase chemistry
C
C     INPUT:
C        LIWK     - length provided for the Chemkin integer work array
C        LRWK     - length provided for the Chemkin real work array
C        LCWK     - length provided for the Chemkin character work array
C        LINC   - unit number for the Chemkin binary file
C        LOUT     - unit number for output messages
C        MDIM     - maximum number of elements allowed
C        KDIM     - maximum number of species allowed
C        MAXTP    - maximum number of temperatures allowed for the
C                   fits of thermodynamic data for the species
C        MAXFIT   - maximum number of coefficients allowed for the
C                   fits of thermodynamic data for the species
C        NTR      - number of temperature ranges provided over which
C                   the fits of thermodynamic data are valid
C
C     OUTPUT:
C        IWORK(*) - Chemkin integer work array
C        RWORK(*) - Chemkin real work array
C        CWORK(*) - Chemkin character work array
C        NELEM    - actual number of elements declared
C        AWT(*)   - real array of element atomic weights
C        ENAME(*) - character array of element names
C        KKGAS    - total number of gas-phas species
C        KNAME(*) - character array of species names
C        WTM(*)   - real array of species molecular weights
C        ITHRM(*) - logical array of flags to indicate presence
C                   of thermodynamic data for the species
C        KCHRG(*) - integer array of electronic charges of the species
C        KPHSE(*) - integer array of phases of the species
C        NT(*)    - integer array of the number of temperatures used
C                   in the fits of thermodynamic data for the species
C        T(*,*)   - matrix of temperatures used in the fits of
C                   thermodynamic data for the species
C        KNCF(*,*)- matrix of elemental composition of the species
C        NFIT     - actual number of coefficients in the fits of
C                   thermodynamic data for the species
C        A(*,*,*) - three-dimensional array of fit coefficients
C                   over the temperature ranges for the species
C        KERR     - logical error flag
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IWORK(*), WORK(*), KNCF(MDIM,*), AWT(*), WTM(*),
     1          KCHRG(*), NT(*), T(MAXTP,*), A(MAXFIT, NTR, *),
     2          KPHSE(*)
      CHARACTER*16 KNAME(*), ENAME(*), CWORK(*)
      LOGICAL ITHRM(*), IERR, KERR
C
c	write(*,*)'i am in skstrt'
      CALL CKINIT (LIWK, LRWK, LCWK, LINC, LOUT, IWORK, WORK, CWORK)
c     1             IFLAG)
c	write(*,*)'i passed ckinit',LIWK
c      IF (IFLAG .GT. 0) STOP
      CALL CKINDX (IWORK, WORK, MM, KK, II, NFIT)
      NELEM = MM
      KKGAS = KK
C
      IF (MDIM .LT. MM) THEN
         WRITE (LOUT, *) ' Element dimension too small,',
     1                   ' MDIM must be at least...',MM
         KERR = .TRUE.
      ELSE
         CALL CKSYME (CWORK, LOUT, ENAME, IERR)
         KERR = KERR.OR.IERR
         CALL CKAWT  (IWORK, WORK, AWT)
      ENDIF
      IF (KDIM .LT. KK) THEN
         WRITE (LOUT, *) ' Species dimension too small,',
     1                   ' KDIM must be at least...',KK
         KERR = .TRUE.
      ELSE
         CALL CKSYMS (CWORK, LOUT, KNAME, IERR)
         KERR = KERR.OR.IERR
         CALL CKWT   (IWORK, WORK, WTM)
         CALL CKCHRG (IWORK, WORK, KCHRG)
         CALL CKPHAZ (IWORK, WORK, KPHSE)
         IF (MDIM .GE. MM) CALL CKNCF  (MDIM, IWORK, WORK, KNCF)
         CALL CKATHM (MAXFIT, NTR, IWORK, WORK, MAXTP, NT, T, A)
         DO 10 K = 1, KKGAS
            ITHRM(K) = .TRUE.
   10    CONTINUE
      ENDIF
      IF (KERR) STOP
C
      RETURN
      END
C----------------------------------------------------------------------C
C
      SUBROUTINE SKKEY (LIN, LTHRM, LOUT, MDIM, NELEM, ENAME, AWT, KDIM,
     1                  KNAME, ITHRM, WTM, KNCF, NONCON, NIICON, KCOV,
     2                  KPHSE, KCHRG, MAXTP, NT, NTR, T, MAXFIT, A,
     3                  MXPHSE, KKGAS, IDIM, NSPAR, MAXSPR, NSCOV,
     4                  NPHASE, PNAME, PDEN, LPDEN, KFIRST, KLAST,
     5                  KKPHAS, NFSURF, NLSURF, NNSURF, KKSURF, NFBULK,
     6                  NLBULK, NNBULK, KKBULK, KKTOT, DEN, LKDEN,
     7                  NIISUR, SPAR, NR, NU, NUNK, RNCF, NREAC, NUSUMK,
     8                  IDUP, NIIREV, IIREV, RSPAR, NIICOV, IICOV,
     9                  IKCOV, CPAR, NIISTK, IISTK, KERR, NIIRNU, IRNU,
     *                  RNU, SKMIN, NIIORD, MAXORD, IORD, KORD, RORD,
     1                  NIIBHM, IBOHM, IBK, IBT, LMTZWS, NMAT, MATNAM,
     2                  MORE)
C
C     SKKEY interprets and stores the surface mechanism data
C
C     INPUT:
C        LIN      - unit number assigned to the ASCII surface data
C        LTHRM    - unit number assigned to thermodynamic database
C        LOUT     - unit number assigned for ASCII output
C        MDIM     - maximum number of elements allowed
C        NELEM    - actual number of elements
C        ENAME(*) - character element names
C        AWT(*)   - real atomic weights of elements
C        KDIM     - maximum number of species allowed
C        KNAME(*) - character species names
C        ITHRM(*) - logical to indicate that thermodynamic
C                   properties of species has been completed
C        WTM(*)   - real molecular weights of species
C        KNCF(*,*)- matrix of elemental composition of species
C        NONCON   - logical error flag to permit non-conservation
C                   of sites in surface reactions
C        NIICON   - total number of reactions which do not
C                   conserve sites
C        KCOV(*)  - species site coverage
C        KPHSE(*) - integer phases of species
C        KCHRG(*) - integer ionic charges of species
C        MAXTP    - maximum number of temperatures used in fits
C                   of thermodynamic properties coefficients
C        NT(*)    - integer number of temperatures used for fits
C        NTR      - integer number of temperature ranges (NT-1)
C        T(*,*)   - matrix of real temperatures used for the fits
C        MAXFIT   - maximum number of fit coefficients for each
C                   temperature range
C        A(*,*,*) - three dimensional array of fit coefficients
C                   for the species and the temperature ranges
C        MXPHSE   - maximum number of phases allowed
C        KKGAS    - total number of gas-phase species
C        IDIM     - maximum number of surface reactions allowed
C        NSPAR    - number of surface reaction parameters required
C        MAXSPR   - maximum number of reactants+products allowed
C                   in any surface reaction
C        NSCOV    - total number of coverage parameters required
C                   for a coverage reaction
C
C     OUTPUT:
C        Information about the phases
C        ----------------------------
C        NPHASE   - total number of phases (gas,site,bulk)
C        PNAME(*) - character array of names for the NPHASE phases
C        PDEN(*)  - real densities for the NPHASE phases
C                   (only used for site phases)
C        KFIRST(*)- integer starting species index for a phase
C        KLAST(*) - integer final species index for a phase
C        KKPHAS(*)- integer number of species per phase
C        NFSURF   - the first phase which is a site
C        NLSURF   - the last  phase which is a site
C        NNSURF   - the total number of site phases
C        NFBULK   - the first phase which is a bulk phase
C        NLBULK   - the last phase which is a bulk phase
C        NNBULK   - the total number of bulk phases
C
C        Information about the species
C        -----------------------------
C        KKSURF   - the total number of surface species
C        KKBULK   - the total number of bulk species
C        KKTOT    - the total number of species
C                   (KKGAS + KKSURF + KKBULK)
C        DEN(*)   - real densities for the KKTOT species
C                   (only used with bulk species)
C
C        Information about the surface reactions
C        ---------------------------------------
C        NIISUR   - total number of surface reactions considered
C        SPAR(*,*)- reaction parameters for the NIISUR surface
C                   reactions
C        NR(*)    - total number of reactants + products in the
C                   NIISUR surface reactions;
C                   NR(I) < 0 = reaction I is irreversible,
C                         > 0 = reaction I is reversible
C        NU(*,*)  - the coefficients of the reactants
C                   and products in the NIISUR surface reactions
C        NUNK(*,*)- the species indices of the reactants
C                   and products in the NIISUR surface reactions
C        NREAC(*) - total number of reactants only in the NIISUR
C                   surface reactions
C        NUSUMK(*)- sums of the coefficients of the
C                   reactants and products in the surface reactions
C        IDUP(*)  - integer flag to identify duplicate reactions
C        NIIREV   - total number of reactions with reverse parameters
C        IIREV(*) - the NIIREV reaction indices
C        RSPAR(*,*)-real reverse parameters for the NIIREV reactions
C        NIICOV   - total number of surface reactions which have
C                   coverage parameters
C        IICOV(*) - the NIICOV reaction indices
C        IKCOV(*) - surface species index numbers for the NIICOV
C                   reactions
C        CPAR(*,*)- real coverage parameters for the NIICOV reactions
C        NIISTK   - total number of surface reactions which have
C                   sticking parameters
C        IISTK(*) - the NIISTK reaction indices
C        LMTZWS   - logical flag to allow modification of sticking
C                   coefficients according to Motz-Wise formula
C        KERR     - logical error flag
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (MXCONT=10)
      DIMENSION AWT(*), WTM(*), KNCF(MDIM,*), KCOV(*), KPHSE(*),
     1          KCHRG(*), NT(*), T(MAXTP,*), A(MAXFIT,NTR,*),
     2          PDEN(*), KKPHAS(*), KFIRST(*), KLAST(*), DEN(*),
     3          SPAR(NSPAR,*), NR(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     4          RNCF(MXPHSE,*), NREAC(*), NUSUMK(*), IDUP(*),
     5          IIREV(*), RSPAR(NSPAR,*), IICOV(*), IKCOV(*),
     6          CPAR(NSCOV,*), IISTK(*), VALUE(5), IRNU(*),
     7          RNU(MAXSPR,*), IORD(*), KORD(MAXORD,*),
     8          RORD(MAXORD,*), IBOHM(*), IBK(*), IBT(*)
C
      CHARACTER*16 KNAME(*), ENAME(*), PNAME(*), MATNAM
      CHARACTER ISTR*80, SUB(80)*80, IUNITS*80, LINE(MXCONT)*80,
     1          ILINE*80, UPCASE*4, SKEY*4, IST(20)*2, EUNITS*4,
     2          AUNITS*4, IAUXL*3
      LOGICAL KERR, IERR, THERMO, ITHRM(*), NONCON, LPDEN(*),
     1              LKDEN(*), LCONT, LTHERM, LMTZWS, LMATE
      DATA IST/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10',
     1         '11','12','13','14','15','16','17','18','19','20'/
C
      THERMO = .TRUE.
      LTHERM = .FALSE.
      IF (NMAT .EQ. 1) THEN
         WRITE (LOUT, 1500)
         WRITE (LOUT, 1600) (ENAME(M)(:2),M=1,NELEM)
         WRITE (LOUT, 1500)
C
         WRITE (LOUT, '(/1X,A)') 'Gas phase species:'
         DO 75 K = 1, KKGAS
            WRITE (LOUT, 1650)
     1      K, KNAME(K), WTM(K), (KNCF(M,K),M=1,NELEM)
   75    CONTINUE
      ENDIF
C
      NPHASE = 1
      KFIRST(NPHASE) = 1
      KLAST (NPHASE) = KKGAS
      KKPHAS(NPHASE) = KKGAS
      PNAME (NPHASE) = 'GAS'
      KKTOT = KKGAS
      ITASK = 0
      MORE = 0
      LMATE = .FALSE.
      MATNAM = 'MATERIAL'//IST(NMAT)
   15 CONTINUE
      NLINES = 0
      LCONT  = .FALSE.
   10 CONTINUE
C
      READ (LIN, '(A)', END=9999) ILINE
      ILEN = IPPLEN(ILINE)
      IF (ILEN .LE. 0) GO TO 10
C
      IS = INDEX (ILINE(:ILEN), '&')
      IF (IS .NE. 0) THEN
         ILEN = IS - 1
         IF (ILEN.LE.0 .OR. ILINE(:ILEN).EQ.' ') GO TO 10
         LCONT = .TRUE.
      ELSE
         LCONT = .FALSE.
      ENDIF
C
      CALL SKISUB (ILINE(:ILEN), SUB, NSUB)
C
C-----look for 'keyword' input
C
      SKEY = UPCASE(SUB(1), 4)
      IF (SKEY.EQ.'MATE') THEN
         IF (LMATE) THEN
            BACKSPACE (LIN)
            MORE = 1
            GO TO 9999
         ELSE
            LMATE = .TRUE.
            IF (NSUB .GT. 1) MATNAM = SUB(2)
            WRITE (LOUT, *)
            WRITE (LOUT, 1500)
            WRITE (LOUT, '(1X,A,A)') 'MATERIAL:',MATNAM
            WRITE (LOUT, 1600) (ENAME(M)(:2),M=1,NELEM)
            WRITE (LOUT,1500)
            GO TO 10
         ENDIF
      ENDIF
C
      IF (SKEY.EQ.'SITE' .OR. SKEY.EQ.'BULK') THEN
C
C--------new phase
C
         IF (NPHASE+1 .GT. MXPHSE) THEN
            WRITE (LOUT, '(A)')' Error...too many phases...'
            KERR = .TRUE.
C
         ELSE
C
            NPHASE = NPHASE + 1
            I1 = INDEX(SUB(1),'/')
            IF (I1 .GT. 0) THEN
               I2 = INDEX(SUB(1)(I1+1:),'/') + I1
               IF (I2 .GT. I1) PNAME(NPHASE) = SUB(1)(I1+1:I2-1)
            ENDIF
C
            IERR = .FALSE.
            KOLD = KKTOT
            IF (KFIRST(NPHASE) .EQ. 0) KFIRST(NPHASE) = KOLD+1
C
            IF (SKEY .EQ. 'SITE') THEN
C
               ITASK = 1
               NNSURF = NNSURF + 1
C
               IF (PNAME(NPHASE) .EQ. ' ')
     1             PNAME(NPHASE)='SITE'//IST(NNSURF)
C
               CALL SKCOMP (PNAME(NPHASE), PNAME, NPHASE-1, I, NF)
               IF (I .GT. 0) THEN
                  ILS = ILASCH(PNAME(NPHASE))
                  ISTR = ' '
                  ISTR = PNAME(NPHASE)(:ILS)
                  WRITE (LOUT, '(A)')
     1            ' Warning...duplicate phase name'//ISTR(:ILS)
                  KERR = .TRUE.
               ENDIF
C
               IF (NNSURF .EQ. 1) NFSURF = NPHASE
               NLSURF = NPHASE
C
               IF (NNBULK .GT. 0) THEN
                  WRITE (LOUT, '(A)')
     1 ' Error...sites must precede bulk phases...',PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (NIISUR .GT. 0) THEN
                  WRITE (LOUT, '(A)')
     1 ' Error...site must precede reactions...',PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (LTHERM) THEN
                  WRITE (LOUT, '(A)')
     1 ' Error...site must precede THERMO data...',PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               CALL SKSURF (SUB(2), NSUB-1, KDIM, KNAME, KCOV, KKTOT,
     1                      KKGAS, KFIRST(NPHASE), LPDEN(NPHASE),
     2                      PDEN(NPHASE), IERR, LOUT)

            ELSE
               ITASK = 2
               NNBULK = NNBULK + 1
C
               IF (PNAME(NPHASE) .EQ. ' ')
     1             PNAME(NPHASE)='BULK'//IST(NNBULK)
C
               CALL SKCOMP (PNAME(NPHASE), PNAME, NPHASE-1, I, NF)
               IF (I .GT. 0) THEN
                  ILS = ILASCH(PNAME(NPHASE))
                  ISTR = ' '
                  ISTR = PNAME(NPHASE)(:ILS)
                  WRITE (LOUT, '(A)')
     1            ' Warning...duplicate phase name',ISTR(:ILS)
                  KERR = .TRUE.
               ENDIF
C
               IF (NNBULK .EQ. 1) NFBULK = NPHASE
               NLBULK = NPHASE
C
               IF (NIISUR .GT. 0) THEN
                  WRITE (LOUT, '(A)')
     1' Error...bulk phase must precede reactions...',PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
               IF (LTHERM) THEN
                  WRITE (LOUT, '(A)')
     1' Error...bulk phase must precede THERMO data...',PNAME(NPHASE)
                  KERR = .TRUE.
               ENDIF
C
C--------------Get species for this bulk phase
C
               CALL SKBULK (SUB(2), NSUB-1, KDIM, KNAME, DEN, LKDEN,
     1                      KKTOT, KKGAS+KKSURF, KFIRST(NPHASE),
     2                      IERR, LOUT)
            ENDIF
C
            KNEW   = KKTOT - KOLD
            KKPHAS(NPHASE) = KKPHAS(NPHASE) + KNEW
            IF (ITASK .EQ. 1) THEN
               KKSURF = KKSURF + KNEW
            ELSE
               KKBULK = KKBULK + KNEW
            ENDIF
            KLAST(NPHASE) = KKTOT
            KERR = KERR.OR.IERR
         ENDIF
C
      ELSEIF (SKEY .EQ. 'REAC') THEN
C
C--------reactions
C
         ITASK = 3
C
C        UNITS or 'NONCON' OPTIONS CAN BE SPECIFIED ON 'REACTION' LINE
C        or 'MWOFF' OPTION (MOTZ-WISE CORRECTION OFF)
C
         CALL SKUNIT (ILINE(:ILEN), AUNITS, EUNITS, IUNITS, NONCON,
     1                LMTZWS)
C
         IF (THERMO) THEN
            OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',
     1            FILE='INP.d/thermdat', ERR=33333)
  311       CONTINUE
            READ (LTHRM,'(A)',END=33333) ILINE
            IF (IPPLEN(ILINE).LE.0 .OR. INDEX(ILINE,'THERM').GT.0
     1          .OR. INDEX(ILINE,'therm').GT.0) GO TO 311
C
            TLO = -1.0
            TMID= -1.0
            THI = -1.0
            CALL CKXNUM (ILINE, 3, LOUT, NVAL, VALUE, IERR)
            IF (IERR) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID= VALUE(2)
               THI = VALUE(3)
            ENDIF
            CALL SKTHRM (LTHRM, MDIM, ENAME, NELEM, AWT, KNAME,
     1                   KKGAS, KKTOT, KNCF, KPHSE, KCHRG, WTM,
     2                   MAXTP, NT, NTR, TLO, TMID, THI, T,
     3                   MAXFIT, A, ITHRM, KERR, LOUT, ILINE)
            CALL SKPRNT (LOUT, MAXTP, KKTOT, KNAME, ITHRM, T,
     1                   NT, TMID, WTM, PDEN, LPDEN, DEN, LKDEN,
     2                   NELEM, MDIM, KNCF, KCOV, NFSURF,
     3                   NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     4                   KFIRST, KLAST, KKPHAS, PNAME, KERR)
C
            WRITE (LOUT, *)
            WRITE (LOUT, 1500)
            THERMO = .FALSE.
            CLOSE (LTHRM)
         ENDIF
         WRITE (LOUT, 1800)
C
      ELSEIF (SKEY .EQ. 'THER') THEN
C
C--------thermodynamic data
C
         LTHERM = .TRUE.
         ITASK = 0
         TLO = -1.0
         TMID= -1.0
         THI = -1.0
         IF (NSUB .GT. 1) THEN
            IF (UPCASE(SUB(2), 3) .EQ. 'ALL') THEN
               THERMO = .FALSE.
  309          CONTINUE
               READ (LIN,'(A)') ILINE
               IF (IPPLEN(ILINE) .LE. 0) GO TO 309
C
               CALL CKXNUM (ILINE, 3, LOUT, NVAL, VALUE, IERR)
               IF (IERR) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 333)
               ELSE
                  TLO = VALUE(1)
                  TMID= VALUE(2)
                  THI = VALUE(3)
               ENDIF
            ENDIF
         ELSE
C
C           NEED TLO,TMID,THI FROM THERMODYNAMIC DATABASE
            OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',
     1                   FILE='INP.d/thermdat', ERR=33333)
312         CONTINUE
            READ (LTHRM,'(A)',END=33333) ILINE
            IF (IPPLEN(ILINE).LE.0 .OR. INDEX(ILINE,'THERM').GT.0
     1          .OR. INDEX(ILINE,'therm').GT.0) GO TO 312
C
            CALL IPPARR (ILINE, -1, 3, VALUE, NVAL, IER, LOUT)
            IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID = VALUE(2)
               THI = VALUE(3)
            ENDIF
            CLOSE (LTHRM)
         ENDIF
C
         CALL SKTHRM (LIN, MDIM, ENAME, NELEM, AWT, KNAME,
     1                KKGAS, KKTOT, KNCF, KPHSE, KCHRG, WTM,
     2                MAXTP, NT, NTR, TLO, TMID, THI, T,
     3                MAXFIT, A, ITHRM, KERR, LOUT, ILINE)
C
         IF (.NOT. THERMO) THEN
            CALL SKPRNT (LOUT, MAXTP, KKTOT, KNAME, ITHRM, T,
     1                   NT, TMID, WTM, PDEN, LPDEN, DEN,
     2                   LKDEN, NELEM, MDIM, KNCF, KCOV, NFSURF,
     3                   NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     4                   KFIRST, KLAST, KKPHAS, PNAME, KERR)
C
            WRITE (LOUT, *)
            WRITE (LOUT, 1500)
         ENDIF
         I1 = IFIRCH(ILINE)
         IF (UPCASE(ILINE(I1:), 4) .EQ. 'REAC') GO TO 15
C
      ELSEIF (INDEX(SKEY, 'END') .GT. 0) THEN
         ITASK = 0
C
      ELSE
C
         IF (ITASK .EQ. 3) THEN
C
C--------------reaction data
C
            IND = 0
            DO 400 N = IFIRCH(ILINE), ILEN-2
               IAUXL = UPCASE(ILINE(N:), 3)
               IF (IAUXL.EQ.'COV' .OR. IAUXL.EQ.'DUP' .OR.
     1             IAUXL.EQ.'REV' .OR. IAUXL.EQ.'STI' .OR.
     2             IAUXL.EQ.'ORD' .OR. IAUXL.EQ.'BOH')
     3             IND=MAX(IND,N)
  400       CONTINUE
C
            IF (IND .GT. 0) THEN
C
C              AUXILIARY REACTION DATA
C
               CALL SKAUXL (ILINE(:ILEN), SUB, NSUB, NIISUR, LOUT,
     1                      NR(NIISUR), NREAC(NIISUR), NU(1,NIISUR),
     2                      NUNK(1,NIISUR), KKGAS, KKTOT, KNAME,
     3                      NPHASE, PNAME, NFSURF, NLSURF,
     4                      KFIRST, KLAST, KKPHAS, NSPAR, NSCOV, IDUP,
     5                      NIIREV, IIREV, RSPAR, NIICOV, IICOV,
     6                      IKCOV, CPAR, NIISTK, IISTK, NIIRNU, IRNU,
     7                      MAXSPR, RNU, NIIORD, IORD,
     7                      MAXORD, KORD, RORD, KERR,
     9                      NIIBHM, IBOHM, IBK, IBT)
C
            ELSE
C
C              THIS IS A REACTION STRING
C
               IF (NIISUR .LT. IDIM) THEN
C
C                 NEW REACTION
C
C                  IF (LCONT) THEN
                     IF (NLINES+1 .GT. MXCONT) THEN
                        KERR = .TRUE.
                        WRITE (LOUT, *)
     1                  ' Error...more than 10 continuation lines...'
                     ELSE
                        NLINES = NLINES + 1
                        LINE(NLINES) = ' '
                        LINE(NLINES) = ILINE(:ILEN)
                     ENDIF
                     IF (LCONT) GO TO 10
C                  ENDIF
C
C              CHECK PREVIOUS REACTION FOR BALANCE AND DUPLICATION
C
                  IF (NIISUR .GT. 0)
     1               CALL SPREAC (NIISUR, MAXSPR, KKGAS, KKSURF, NR,
     2                            NREAC, NSPAR, SPAR, RSPAR, AUNITS,
     3                            EUNITS, NUNK, NU, RNCF(1,NIISUR),
     4                            KCHRG, NIISTK, IISTK, NIICOV, IICOV,
     5                            NSCOV, CPAR, KCOV, NIIREV, IIREV,
     6                            MDIM, NELEM, KNCF, NONCON, NIICON,
     7                            IDUP, NNSURF, NPHASE, KFIRST, KLAST,
     8                            LOUT, KERR, NIIRNU, IRNU, RNU,SKMIN)
C
                  NIISUR = NIISUR + 1
C
                  CALL SKREAC (NIISUR, LINE, NLINES, KKTOT,
     1                         KNAME, NPHASE, PNAME, KKPHAS,
     2                         NSPAR, LOUT, MAXSPR, NR(NIISUR),
     3                         NREAC(NIISUR), NUNK(1,NIISUR),
     4                         NU(1,NIISUR), SPAR(1,NIISUR), KERR,
     5                         NIIRNU, IRNU, RNU)
C
                  DO 500 N = 1, MAXSPR
                     IF (NUNK(N,NIISUR) .LE. KKGAS)
     1               NUSUMK(NIISUR) = NUSUMK(NIISUR) + NU(N,NIISUR)
  500             CONTINUE
               ENDIF
            ENDIF
C
         ELSEIF (ITASK .GT. 0) THEN
C
C-----------species data
C
            IERR = .FALSE.
            KOLD = KKTOT
            IF (KFIRST(NPHASE) .EQ. 0) KFIRST(NPHASE) = KOLD + 1
            IF (ITASK .EQ. 1) THEN
               CALL SKSURF (SUB, NSUB, KDIM, KNAME, KCOV, KKTOT,
     1                      KKGAS, KFIRST(NPHASE), LPDEN(NPHASE),
     2                      PDEN(NPHASE), IERR, LOUT)
            ELSE
               CALL SKBULK (SUB, NSUB, KDIM, KNAME, DEN, LKDEN,KKTOT,
     1                      KKGAS+KKBULK, KFIRST(NPHASE), IERR, LOUT)
            ENDIF
            KNEW   = KKTOT - KOLD
            KKPHAS(NPHASE) = KKPHAS(NPHASE) + KNEW
            IF (ITASK .EQ. 1) THEN
               KKSURF = KKSURF + KNEW
            ELSE
               KKBULK = KKBULK + KNEW
            ENDIF
            KLAST(NPHASE) = KKTOT
            KERR = KERR.OR.IERR
         ENDIF
      ENDIF
      GO TO 15
 9999 CONTINUE
C
C     END OF INPUT
C
      IF (NIISUR .GT. 0) THEN
C
C        CHECK FINAL REACTION
C
         CALL SPREAC (NIISUR, MAXSPR, KKGAS, KKSURF, NR, NREAC, NSPAR,
     1                SPAR, RSPAR, AUNITS, EUNITS, NUNK, NU,
     2                RNCF(1,NIISUR), KCHRG, NIISTK, IISTK, NIICOV,
     3                IICOV, NSCOV, CPAR, KCOV, NIIREV, IIREV, MDIM,
     4                NELEM, KNCF, NONCON, NIICON, IDUP, NNSURF,
     5                NPHASE, KFIRST, KLAST, LOUT, KERR, NIIRNU,
     6                IRNU, RNU, SKMIN)
C
C        CHECK REACTIONS DECLARED AS DUPLICATES
C
         DO 600 I = 1, NIISUR
            IF (IDUP(I) .LT. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1095) I
            ENDIF
  600    CONTINUE
C
         WRITE (LOUT, '(/1X,A)') ' NOTE:  '//IUNITS(:ILASCH(IUNITS))
         IF (NONCON) WRITE (LOUT, '(10X,A)')
     1               'Site non-conservation specified'
         IF (.NOT.LMTZWS) WRITE (LOUT, '(10X,A)')
     1     'Motz-Wise correction to sticking coefficients is turned off'
C
      ELSEIF (THERMO) THEN
C
C        THERE WAS NO REACTION DATA, MAKE SURE SPECIES DATA IS COMPLETE
C
         OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',
     1                FILE='INP.d/thermdat', ERR=33333)
C
         IF (KKTOT .GT. KKGAS) THEN
            TLO = -1.0
            TMID= -1.0
            THI = -1.0
  313       CONTINUE
            READ (LTHRM,'(A)',END=33333) ILINE
            IF (IPPLEN(ILINE).LE.0 .OR. INDEX(ILINE,'THERM').GT.0
     1          .OR. INDEX(ILINE,'therm').GT.0) GO TO 313
C
            CALL CKXNUM (ILINE, 3, LOUT, NVAL, VALUE, IERR)
            IF (IERR) THEN
               WRITE (LOUT, 333)
               KERR = .TRUE.
            ELSE
               TLO = VALUE(1)
               TMID= VALUE(2)
               THI = VALUE(3)
            ENDIF
            CALL SKTHRM (LTHRM, MDIM, ENAME, NELEM, AWT, KNAME,
     1                   KKGAS, KKTOT, KNCF, KPHSE, KCHRG, WTM,
     2                   MAXTP, NT, NTR, TLO, TMID, THI, T,
     3                   MAXFIT, A, ITHRM, KERR, LOUT, ILINE)
         ENDIF
C
         IF (NNSURF.GT.0 .OR. NNBULK.GT.0) THEN
            CALL SKPRNT (LOUT, MAXTP, KKTOT, KNAME, ITHRM, T,
     1                   NT, TMID, WTM, PDEN, LPDEN, DEN,
     2                   LKDEN, NELEM, MDIM, KNCF, KCOV, NFSURF,
     3                   NLSURF, NNSURF, NFBULK, NLBULK, NNBULK,
     4                   KFIRST, KLAST, KKPHAS, PNAME, KERR)
         ENDIF
         WRITE (LOUT, *)
         WRITE (LOUT, 1500)
         CLOSE (LTHRM)
      ENDIF
C
  333 FORMAT (/6X,'Error...no TLO,TMID,THI given for THERMO...'/)
 1020 FORMAT (6X, 'Error in site declaration...',A)
 1500 FORMAT (1X, 79('-'))
 1600 FORMAT (1X, 'SPECIES',16X,'MOLECULAR',24X,'ELEMENT COUNT',/
     1        1X, 'CONSIDERED',13X,'WEIGHT',7X,'Density',4X,'Nsites',
     2        3X, 15A3)
 1650 FORMAT (I4,'. ', A16, F11.5, 23X, 10I3)
 1800 FORMAT (54X, '(k = A T**b exp(-E/RT))',/,
     1         6X, 'SURFACE REACTIONS CONSIDERED',
     2        22X, 'A',8X,'b',8X,'E',/)
 1095 FORMAT (6X,'Error...no duplicate declared for reaction no.',I3)
      RETURN
C
33333 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read therm.dat...'
      CLOSE (LTHRM)
      STOP 2
C
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKREAC (II, LINE, NLINES, KKTOT, KNAME, NPHASE,
     1                   PNAME, KKPHAS, NPAR, LOUT, MAXSPR,
     2                   NSPEC, NREAC, KSPEC, KCOEF, PAR, KERR, NIIRNU,
     3                   IRNU, RNU)
C
C     Accepts a reaction character string, finds reactants, products,
C     and their coefficients, and the Arrhenius parameters
C
C     INPUT:
C        II       - the index number of the reaction
C        LINE(*)  - array of character strings describing the reaction
C        NLINES   - total integer number of lines
C        KKTOT    - total integer number of species
C        KNAME(*) - character species names
C        NPHASE   - total integer number of phases
C        PNAME(*) - character phase names
C        KKPHAS(*)- integer array of number of species of the phases
C        NPAR     - number of parameters required
C        LOUT     - output unit for error messages
C        MAXSPR   - maximum number of species allowed in reaction
C
C     OUTPUT:
C        NSPEC    - total number of reactants+products in reaction
C        NREAC    - total number of reactants
C        KSPEC(*) - the NSPEC species indices
C        KCOEF(*) - the NSPEC stoichiometric coefficients
C        PAR(*)   - real parameters for this reaction
C        KERR     - logical error flag
C
C                                      F. Rupley, Div. 8245, 5/13/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C----------------------------------------------------------------------C
      DIMENSION PAR(*), KSPEC(*), KCOEF(*), KKPHAS(*),
     1          IPLUS(20), RNU(MAXSPR,*), IRNU(*)
      CHARACTER*16 KNAME(*), PNAME(*)
      CHARACTER*(*) LINE(*) 
      CHARACTER*1 CNUM(11)
      CHARACTER*80 ISTR, IREAC, IPROD, ISPEC, INAME, ITEMP
      LOGICAL KERR, IERR, LRSTO
      DATA CNUM/'.','0','1','2','3','4','5','6','7','8','9'/
C
      NSPEC = 0
      NREAC = 0
C
C----------Find NPAR real parameters------------------------
C
      NFOUND = NPAR
      NEND   = NLINES
C
C----------Loop from last line to first line
C
      DO 50 N = NEND, 1, -1
C
C----------Loop from last character on a line to first character
C
   10    CONTINUE
         I1 = IFIRCH(LINE(N))
         I2 = ILASCH(LINE(N))
C
         DO 30 L = I2, I1, -1
C
C--------------Right-most non-blank character found
C
            IF (LINE(N)(L:L) .NE. ' ') THEN
               JEND = L
               DO 20 J = L, I1, -1
C
C-----------------Left-most non-blank character found
C
                  IF (LINE(N)(J:J).EQ.' ' .OR. J.EQ.I1) THEN
                     IF (LINE(N)(J:J) .EQ. ' ') THEN
                        JSTART =  J + 1
                     ELSE
                        JSTART = I1
                     ENDIF
C
C--------------------Convert to real value
C
                     ISTR = ' '
                     ISTR = LINE(N)(JSTART:JEND)
C
C--------------------If Arhennius coefficients are missing,
C                    the string found is some portion of
C                    the reaction -
C
C                    a) portion of reaction string with delimiter
                     KNUM = INDEX(ISTR,'=')
C
C                    b) species
                     IF (KNUM .LE. 0)
     1               CALL SKCOMP (ISTR, KNAME, KKTOT, KNUM, NF)
C
C                    c) contains '/'

                     IF (KNUM .LE. 0) KNUM = INDEX(ISTR,'/')
C     1               CALL SKPCMP (ISTR, KNAME, KKTOT, PNAME, NPHASE,
C     2                            KKPHAS, KNUM, NF)
C
C                    d) portion of reaction string with '+'
                     IF (KNUM .LE. 0) THEN
                        ILAST = INDEX(ISTR,' ')
                        IND = -1
                        DO 99 M = ILAST-1, 1
                           IF (ISTR(M:M).EQ.'+' .AND. IND.LT.0) IND=M
   99                   CONTINUE
                        IF (IND .GT. 0)
     1                  CALL SKPCMP (ISTR(IND+1:), KNAME, KKTOT, PNAME,
     2                  NPHASE, KKPHAS, KNUM, NF)
                     ENDIF
C
                     IF (KNUM .GT. 0) THEN
                        KERR = .TRUE.
                        GO TO 52
                     ENDIF
C
                     LINE(N)(JSTART:) = ' '
                     CALL IPPARR (ISTR, -1, 1, VAL, NVAL, IERR, LOUT)
                     KERR = KERR.OR.IERR
                     IF (NVAL .EQ. 1) THEN
                        PAR(NFOUND) = VAL
C
C-----------------------Have all parameters been found?
C
                        NFOUND = NFOUND -1
                        IF (NFOUND .EQ. 0) THEN
                           IF (J.EQ.I1 .OR. LINE(N).EQ.' ')
     1                         NLINES = N - 1
                           GO TO 52
                        ELSE
                           IF (J .EQ. I1) GO TO 50
                           GO TO 10
                        ENDIF
                     ENDIF
                  ENDIF
   20          CONTINUE
            ENDIF
   30    CONTINUE
   50 CONTINUE
C
   52 CONTINUE
C
C-----Remove blanks from reaction strings
C
      DO 60 N = 1, NLINES
         INAME = ' '
         L = 0
         DO 55 I = IFIRCH(LINE(N)), ILASCH(LINE(N))
            IF (LINE(N)(I:I) .NE. ' ') THEN
               L = L + 1
               INAME(L:L) = LINE(N)(I:I)
            ENDIF
   55    CONTINUE
         LINE(N) = ' '
         LINE(N) = INAME
   60 CONTINUE
C
C-------Find delimiter '<=>", '=', or '=>'
C
      IR   = 1
      DO 75 N = 1, NLINES
C
         IF (INDEX(LINE(N), '<=>') .GT. 0) THEN
            IND = INDEX(LINE(N), '<=>')
            NNEXT = IND + 3
            LINNO = N
         ELSEIF (INDEX(LINE(N), '=>') .GT. 0) THEN
            IND = INDEX(LINE(N), '=>')
            IR  = -1
            NNEXT = IND + 2
            LINNO = N
         ELSEIF (INDEX(LINE(N), '=') .GT. 0) THEN
            IND = INDEX(LINE(N), '=')
            IF ((IND.EQ.1 .AND. LINE(N)(IND+1:IND+1).NE.'=') .OR.
     1         (IND.GT.1 .AND. LINE(N)(IND-1:IND-1).NE.'=')) THEN
               NNEXT = IND + 1
               LINNO = N
            ENDIF
         ENDIF
   75 CONTINUE
C
      IF (IND .LE. 0) THEN
C
C        DELIMETER WAS NOT FOUND
C
         WRITE (LOUT, 800)
         KERR = .TRUE.
      ELSE
C---------------------store reactant string, product string
C
         IREAC = ' '
         IPROD = ' '
         NCHAR = IND - 1
         DO 66 J = 1, LINNO-1
            IF (IREAC .EQ. ' ') THEN
               IREAC = LINE(J)
            ELSE
               IREAC(ILASCH(IREAC)+1:) = LINE(J)
            ENDIF
   66    CONTINUE
         IF (NCHAR .GT. 0) THEN
            IF (IREAC .EQ. ' ') THEN
               IREAC = LINE(LINNO)(:NCHAR)
            ELSE
               IREAC(ILASCH(IREAC)+1:) = LINE(LINNO)(:NCHAR)
            ENDIF
         ENDIF
C
         IPROD = LINE(LINNO)(NNEXT:)
         DO 67 J = LINNO+1, NLINES
            IF (IPROD .EQ. ' ') THEN
               IPROD = LINE(J)
            ELSE
               IPROD(ILASCH(IPROD)+1:) = LINE(J)
            ENDIF
   67    CONTINUE
C
         INAME = ' '
         IN = ILASCH(IREAC)
         IP = ILASCH(IPROD)
         IF (IN+IP .GE. 38) THEN
            WRITE (LOUT, 1900) II, IREAC(:IN), (PAR(N),N=1,NPAR)
            IF (IR .GT. 0) THEN
               WRITE (LOUT, 1920) '<=>'//IPROD(:IP)
            ELSE
               WRITE (LOUT, 1920) '=>'//IPROD(:IP)
            ENDIF
         ELSE
            IF (IR .GT. 0) THEN
               WRITE (LOUT, 1900) II, IREAC(:IN)//'<=>'//IPROD(:IP),
     1                            (PAR(N), N = 1, NPAR)
            ELSE
               WRITE (LOUT, 1900) II, IREAC(:IN)//'=>'//IPROD(:IP),
     1                            (PAR(N), N = 1, NPAR)
            ENDIF
         ENDIF
C
         IF (NFOUND .NE. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...NPAR Arrhenius coefficients not found...'
         ENDIF
         IF (PAR(1) .LT. 0.0) THEN
            WRITE (LOUT, 4010)
            KERR = .TRUE.
         ENDIF

      ENDIF
C
C
C-----Find reactants, products
C
      LRSTO = ((INDEX(IREAC,'.').GT.0) .OR.INDEX(IPROD,'.').GT.0)
      IF (LRSTO) THEN
         NIIRNU = NIIRNU + 1
         IRNU(NIIRNU) = II
      ENDIF
C
      DO 600 J = 1, 2
         ISTR = ' '
         IF (J .EQ. 1) THEN
            ISTR = IREAC
            NS = 0
         ELSE
            ISTR = IPROD
            NS = MAXSPR/2
         ENDIF
C
C-----------store pointers to '+'-signs
C
         NPLUS = 1
         IPLUS(NPLUS) = 0
         DO 500 L = 2, ILASCH(ISTR)-1
            IF (ISTR(L:L) .EQ. '+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = L
            ENDIF
  500    CONTINUE
         NPLUS = NPLUS + 1
         IPLUS(NPLUS) = ILASCH(ISTR)+1
C
         NSTART = 1
  505    CONTINUE
         N1 = NSTART
         DO 510 N = NPLUS, N1, -1
            ISPEC = ' '
            ISPEC = ISTR ( IPLUS(N1)+1 : IPLUS(N)-1)
C
C-----------does this string start with a number?
C
            IND = 0
            DO 334 L = 1, LEN(ISPEC)
               NTEST = 0
               DO 333 M = 1, 11
                  IF (ISPEC(L:L) .EQ. CNUM(M)) THEN
                     NTEST=M
                     IND = L
                  ENDIF
  333          CONTINUE
               IF (NTEST .EQ. 0) GO TO 335
  334       CONTINUE
  335       CONTINUE
C
            RVAL = 1.0
            IVAL = 1
            IF (IND .GT. 0) THEN
               IF (LRSTO) THEN
                  CALL IPPARR (ISPEC(:IND), 1, 1, RVAL, NVAL,
     1                         IER, LOUT)
               ELSE
                  CALL IPPARI (ISPEC(:IND), 1, 1, IVAL, NVAL,
     1                         IER, LOUT)
               ENDIF
               ITEMP = ' '
               ITEMP = ISPEC(IND+1:)
               ISPEC = ' '
               ISPEC = ITEMP
            ENDIF
C
            CALL SKPCMP (ISPEC, KNAME, KKTOT, PNAME, NPHASE,
     1                   KKPHAS, KNUM, NF)
            IF (NF .GT. 1) THEN
               WRITE (LOUT, 675) ISPEC(:ILASCH(ISPEC))
               KERR = .TRUE.
               GO TO 600
            ENDIF
C
            IF (KNUM .EQ. 0) THEN
               IF ((N-N1) .GT. 1) GO TO 510
               WRITE (LOUT, 680) ISPEC(:ILASCH(ISPEC))
               KERR = .TRUE.
               GO TO 600
            ELSE
C
C--------------a species has been found
C
               IF (J .EQ. 1) THEN
                  IVAL = -IVAL
                  RVAL = -RVAL
               ENDIF
C
               NNUM = 0
               IF (LRSTO) THEN
                  DO 110 K = 1, NS
                     IF (KNUM.EQ.KSPEC(K) .AND.
     1                  RNU(K,NIIRNU)/RVAL.GT.0) THEN
                        NNUM = K
                        RNU(K,NIIRNU)=RNU(NNUM,NIIRNU) + RVAL
                     ENDIF
  110             CONTINUE
               ELSE
                  DO 111 K = 1, NS
                     IF (KNUM.EQ.KSPEC(K) .AND.
     1                   KCOEF(K)/IVAL.GT.0) THEN
C
C-----------------increment species coefficient count
C
                         NNUM = K
                         KCOEF(NNUM) = KCOEF(NNUM) + IVAL
                     ENDIF
  111             CONTINUE
               ENDIF
C
               IF (NNUM .LE. 0) THEN
C
C-----------------are there too many species?
C
                  IF (J.EQ.1 .AND. NS.EQ.MAXSPR/2) THEN
                     WRITE (LOUT, 690)
                     KERR = .TRUE.
                  ELSEIF (J.EQ.2 .AND. NS.EQ.MAXSPR) THEN
                     WRITE (LOUT, 700)
                     KERR = .TRUE.
                  ELSE
C
C--------------------increment species count
C
                     NS = NS + 1
                     NSPEC = NSPEC+1
                     IF (J .EQ. 1) NREAC = NS
                     KSPEC(NS) = KNUM
                     IF (LRSTO) THEN
                        RNU(NS,NIIRNU) = RVAL
                     ELSE
                        KCOEF(NS) = IVAL
                     ENDIF
                  ENDIF
               ENDIF
               IF (N .EQ. NPLUS) GO TO 600
               NSTART = N
               GO TO 505
            ENDIF
C
  510    CONTINUE
  600 CONTINUE
C
      NSPEC = IR*NSPEC
C
  650 FORMAT (6X,'Error...reaction string not found...')
  675 FORMAT (6X,'Error...ambiguous duplicate species...',A,
     1       /6X,'possible correction is "speciesname/phasename/"...')
  680 FORMAT (6X,'Error...undeclared species...',A)
  690 FORMAT (6X,'Error...more than 6 reactants...')
  700 FORMAT (6X,'Error...more than 6 products...')
  800 FORMAT (6X,'Error in reaction delimiter...')
 1900 FORMAT (I4,'. ', A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 1920 FORMAT (6X,A)
 4010 FORMAT (6X,'Error...pre-exponential factor must be positive...')
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKTHRM (LUNIT, MDIM, MCHAR, MM, AWT, KCHAR, KK, KKTOT,
     1                   KNCF, KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO,
     2                   TMID, THI, T, MAXFIT, A, ITHRM, KERR, LOUT,
     3                   ISTR)
C
C     Finds thermodynamic data and elemental composition for species
C     INPUT:
C        LUNIT    - unit number for input of thermo properties
C        MDIM     - maximum number of elements allowed
C        MCHAR(*) - character element names
C        MM       - total number of elements declared
C        AWT(*)   - atomic weights of elements
C        KCHAR(*) - character species names
C        KK       - total number of gas-phase species
C        KKTOT    - total number of species
C        LOUT     - output unit for messages
C        NT(*)    - number of temperature values used for fits
C        NTR      - number of temperature ranges
C
C     OUTPUT:
C        KNCF(*,*)- elemental composition of species
C        KPHSE(*) - species phases
C        KCHRG(*) - species charges, -1)*number of electrons
C        WTM(*)   - molecular weights of species
C        A(*,*,*) - thermodynamic coefficients for fits
C        T(*)     - temperatures used for fits
C        ITHRM(*) - logical; ITHRM(K)=.TRUE. if
C                         thermodynamic properties found
C        KERR     - logical error flag
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WTM(*), NT(*), T(MAXTP,*), KPHSE(*), KNCF(MDIM,*),
     1          KCHRG(*), A(MAXFIT,NTR,*), AWT(*)
      CHARACTER*16 MCHAR(*), KCHAR(*), ELEM
      CHARACTER*80 LINE(4), ISTR, SUB(80)
      CHARACTER UPCASE*4
      LOGICAL KERR, IERR, ITHRM(*)
C
      IF (MM.LE.0 .OR. (KKTOT-KK).LE.0) WRITE (LOUT, 80)
C
      GO TO 20
   10 CONTINUE
      ISTR = ' '
      READ (LUNIT, '(A)', END=70) ISTR
   20 CONTINUE
      ILEN = IPPLEN(ISTR)
      IF (ILEN .LE. 0) GO TO 10
C
      CALL SKISUB (ISTR(:ILEN), SUB, NSUB)
      IF (UPCASE(SUB(1), 3).EQ.'END' .OR.
     1    UPCASE(SUB(1), 4).EQ.'REAC') RETURN
C
      IF (ILEN.LT.80 .OR. ISTR(80:80).NE.'1') GO TO 10
C
      CALL SKCOMP (SUB(1), KCHAR, KKTOT, KNUM, NK)
      IF (KNUM.LE.0 .OR. (KNUM.GT.0 .AND. ITHRM(KNUM))) GO TO 10
      ITHRM(KNUM) = .TRUE.
      LINE(1) = ' '
      LINE(1) = ISTR
      L = 2
  111 CONTINUE
      READ (LUNIT,'(A)',END=70) LINE(L)
      IF (IPPLEN(LINE(L)) .GE. 80) THEN
         IF (LINE(L)(80:80) .EQ. '4') THEN
            GO TO 25
         ELSEIF (LINE(L)(80:80).EQ.'2' .OR.
     1           LINE(L)(80:80).EQ.'3') THEN
            L = L + 1
         ENDIF
      ENDIF
      GO TO 111
C
   25 CONTINUE
C
      ICOL = 20
      DO 30 I = 1, 5
         ICOL = ICOL + 5
         IF (I .EQ. 5) ICOL = 74
         ELEM  = LINE(1)(ICOL:ICOL+1)
         IELEM = 0
C
         IF (LINE(1)(ICOL+2:ICOL+4) .NE. ' ') THEN
            CALL CKXNUM (LINE(1)(ICOL+2:ICOL+4), 1, LOUT, NVAL, RVAL,
     1                   IERR)
            KERR = KERR.OR.IERR
            IELEM = INT(RVAL)
         ENDIF
C
         IF (ELEM.NE.' ' .AND. IELEM.NE.0) THEN
            IF (UPCASE(ELEM, 1) .EQ. 'E')
     1         KCHRG(KNUM) = KCHRG(KNUM) + IELEM*(-1)
            CALL SKCOMP (ELEM, MCHAR, MM, M, NE)
            IF (M .GT. 0) THEN
               KNCF(M,KNUM) = IELEM
               WTM(KNUM) = WTM(KNUM) + AWT(M)*FLOAT(IELEM)
            ELSE
               WRITE (LOUT, 100) ELEM,KCHAR(KNUM)(:10)
               KERR = .TRUE.
            ENDIF
         ENDIF
   30 CONTINUE
C
      KPHSE(KNUM) = 0
      IF (UPCASE(LINE(1)(45:), 1) .EQ. 'L') KPHSE(KNUM)=1
      IF (UPCASE(LINE(1)(45:), 1) .EQ. 'S') KPHSE(KNUM)=-1
C
C-----Currently allows for three temperatures, two ranges;
C     in future, NT(K) may vary, NTR = NT(K)-1
C
      NT(KNUM)  = 3
      T(1,KNUM) = TLO
      IF (LINE(1)(46:55) .NE. ' ')
     1   CALL CKXNUM (LINE(1)(46:55), 1, LOUT, NVAL, T(1,KNUM), IERR)
      KERR = KERR.OR.IERR
C
      T(2,KNUM) = TMID
      IF (LINE(1)(66:73) .NE. ' ')
     1   CALL CKXNUM (LINE(1)(66:73), 1, LOUT, NVAL, T(2,KNUM), IERR)
      KERR = KERR.OR.IERR
C
      T(NT(KNUM),KNUM) = THI
      IF (LINE(1)(56:65) .NE. ' ')
     1   CALL CKXNUM (LINE(1)(56:65), 1, LOUT, NVAL, T(NT(KNUM),KNUM),
     2                IERR)
      KERR = KERR.OR.IERR
C
      READ (LINE(2)(:75),'(5E15.8)') (A(I,NTR,KNUM),I=1,5)
      READ (LINE(3)(:75),'(5E15.8)')
     1            (A(I,NTR,KNUM),I=6,7),(A(I,1,KNUM),I=1,3)
      READ (LINE(4)(:60),'(4E15.8)') (A(I,1,KNUM),I=4,7)
C
      IF (NK .GT. 1) THEN
C
C        SPECIES NAME OCCURS AGAIN (IN ANOTHER PHASE)
C
         DO 60 K = KNUM+1, KKTOT
            IF (KCHAR(K) .EQ. KCHAR(KNUM)) THEN
               ITHRM(K) = ITHRM(KNUM)
               KCHRG(K) = KCHRG(KNUM)
               KPHSE(K) = KPHSE(KNUM)
               WTM(K)   = WTM(KNUM)
               NT(K)    = NT(KNUM)
               DO 40 M = 1, MM
                  KNCF(M,K) = KNCF(M,KNUM)
   40          CONTINUE
               DO 50 N = 1, NT(K)
                  T(N,K) = T(N,KNUM)
   50          CONTINUE
               DO 55 N = 1, 7
                  A(N,1,  K) = A(N,1,  KNUM)
                  A(N,NTR,K) = A(N,NTR,KNUM)
   55          CONTINUE
            ENDIF
   60    CONTINUE
      ENDIF
      GO TO 10
C
   70 RETURN
   80 FORMAT (6X,'Warning...THERMO cards misplaced will be ignored...')
  100 FORMAT (6X,'Error...element...',A,'not declared for...',A)
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKBAL (MXSPEC, KSPEC, KCOEF, RNCF, MDIM, MM, KNCF,
     1                  KCOV, NONCON, NCON, KCHRG, NNSURF, NPHASE,
     1                  KFIRST, KLAST, KERR, SKMIN)
C
C     Checks balance of ionic charge, element, and sites for
C     a surface reaction
C
C     INPUT:
C        MXSPEC   - number of species allowed in a reaction
C        KSPEC(*) - the MXSPEC species indices
C        KCOEF(*) - the MXSPEC stoichiometric coefficients
C                   (reactants < 0, products > 0)
C        RNCF(*)  - net change of phases
C        MDIM     - maximum number of elements allowed
C        MM       - actual integer number of elements
C        KNCF(*,*)- integer elemental composition of species
C        KCOV(*)  - integer species site coverage
C        NONCON   - logical flag to indicate whether or not non-
C                   conservation of sites reactions are allowed
C        KCHRG(*) - electronic charge of species
C        NNSURF   - total number of sites
C        NPHASE   - total number of phases
C        KFIRST(*)- index number of the first species of phases
C        KLAST(*) - index number of the last species of phases
C
C     OUTPUT:
C        NCON     - total number of surface reactions which do not
C                   conserve sites
C        KERR     - logical, .TRUE. if reaction does not balance
C
C-----------------------------------------------------------------------
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KCOEF(*), KSPEC(*), KNCF(MDIM,*), KCOV(*), KCHRG(*),
     1          KFIRST(*), KLAST(*), RNCF(*)
      LOGICAL KERR,NONCON
C
      KBAL = 0
      MBAL = 0
      DO 60 N = 1, MXSPEC
C
         IF (KSPEC(N) .NE. 0) THEN
C
            KBAL = KBAL + KCOEF(N)*KCHRG(KSPEC(N))
C
            DO 40 M = 1, MM
               MBAL = MBAL + KCOEF(N)*KNCF(M,KSPEC(N))
   40       CONTINUE
C
C           rncf(j) is the net change of phase j for reaction I
C
            DO 50 J = 1, NPHASE
               IF (KSPEC(N).GE.KFIRST(J) .AND. KSPEC(N).LE.KLAST(J))
     1              RNCF(J) = RNCF(J) + KCOEF(N)*KCOV(KSPEC(N))
   50       CONTINUE
C
         ENDIF
C
   60 CONTINUE
      KERR = KERR .OR. (MBAL.NE.0) .OR. (KBAL.NE.0)
C
      ICON = 0
      DO 70 J = 1, NNSURF
         IF (ABS(RNCF(J+1)) .GT. SKMIN) ICON = ICON + 1
   70 CONTINUE
      IF (ICON .NE. 0) THEN
         IF (.NOT. NONCON) THEN
            KERR = .TRUE.
         ELSE
C           This is a legitimate non-conservation of sites reaction
            NCON = NCON + 1
         ENDIF
      ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKRBAL (MXSPEC, KSPEC, RCOEF, RNCF, MDIM, MM, KNCF,
     1                   KCOV, NONCON, NCON, KCHRG, NNSURF, NPHASE,
     2                   KFIRST, KLAST, KERR, SKMIN)
C
C     Checks balance of ionic charge, element, and sites for
C     a surface reaction
C
C     INPUT:
C        MXSPEC   - total number of species in this reaction
C        KSPEC(*) - the MXSPEC species indices
C        RCOEF(*) - the MXSPEC stoichiometric coefficients
C                   (reactants < 0, products > 0)
C        RNCF(*)  - net change of phases
C        MDIM     - maximum number of elements allowed
C        MM       - actual integer number of elements
C        KNCF(*,*)- integer elemental composition of species
C        KCOV(*)  - integer species site coverage
C        NONCON   - logical flag to indicate whether or not non-
C                   conservation of sites reactions are allowed
C        KCHRG(*) - electronic charge of species
C        NNSURF   - total number of sites
C        NPHASE   - total number of phases
C        KFIRST(*)- index number of the first species of phases
C        KLAST(*) - index number of the last species of phases
C
C     OUTPUT:
C        NCON     - total number of surface reactions which do not
C                   conserve sites
C        KERR     - logical, .TRUE. if reaction does not balance
C
C-----------------------------------------------------------------------
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RCOEF(*), KSPEC(*), KNCF(MDIM,*), KCOV(*), KCHRG(*),
     1          KFIRST(*), KLAST(*), RNCF(*)
      LOGICAL KERR,NONCON
C
      SKBAL = 0.0
      SMBAL = 0.0
      DO 60 N = 1, MXSPEC
C
         IF (KSPEC(N) .NE. 0) THEN
C
            SKBAL = SKBAL + RCOEF(N)*KCHRG(KSPEC(N))
C
            DO 40 M = 1, MM
               SMBAL = SMBAL + RCOEF(N)*KNCF(M,KSPEC(N))
   40       CONTINUE
C
C           rncf(j) is the net change of phase j for reaction I
C
            DO 50 J = 1, NPHASE
               IF (KSPEC(N).GE.KFIRST(J) .AND. KSPEC(N).LE.KLAST(J))
     1              RNCF(J) = RNCF(J) + RCOEF(N)*KCOV(KSPEC(N))
   50       CONTINUE
C
         ENDIF
C
   60 CONTINUE
      KERR = KERR .OR. (ABS(SMBAL).GT.SKMIN)
     1            .OR. (ABS(SKBAL).GT.SKMIN)
C
      ICON = 0
      DO 70 J = 1, NNSURF
         IF (ABS(RNCF(J+1)) .GT. SKMIN) ICON = ICON + 1
   70 CONTINUE
      IF (ICON .NE. 0) THEN
         IF (.NOT. NONCON) THEN
            KERR = .TRUE.
         ELSE
C           This is a legitimate non-conservation of sites reaction
            NCON = NCON + 1
         ENDIF
      ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKDUP (I, MAXSP, NS, NR, NU, NUNK, ISAME)
C
C     Checks reaction I against the (I-1) reactions for duplication
C
C     INPUT:
C        I         - index number of a reaction
C        MAXSP     - maximum number of species allowed in a reaction
C        NS(*)     - array of the number of species in the reactions
C        NR(*)     - array of the number of reactants in the reactions
C        NU(*,*)   - matrix of the stoichiometric coefficients of the
C                    species in the reactions
C        NUNK(*,*) - matrix of the species indices for the species in
C                    the reactions
C     OUTPUT:
C        ISAME     - integer number which is 0 if there is no duplicate
C                    reaction to I, else ISAME is the index number of
C                    the duplicte reaction
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NS(*), NR(*), NU(MAXSP,*), NUNK(MAXSP,*)
C
      ISAME = 0
      NRI = NR(I)
      NPI = ABS(NS(I)) - NR(I)
C
      DO 500 J = 1, I-1
C
         NRJ = NR(J)
         NPJ = ABS(NS(J)) - NR(J)
C
         IF (NRJ.EQ.NRI .AND. NPJ.EQ.NPI) THEN
C
            NSAME = 0
            DO 20 N = 1, MAXSP
C
               KI = NUNK(N,I)
               NI = NU(N,I)
C
               DO 15 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. NJ.EQ.NI)
     1            NSAME = NSAME + 1
   15          CONTINUE
   20       CONTINUE
C
            IF (NSAME .EQ. ABS(NS(J))) THEN
C
C           same products, reactants, coefficients
C
               ISAME = J
               RETURN
            ENDIF
         ENDIF
C
         IF (NPI.EQ.NRJ .AND. NPJ.EQ.NRI) THEN
C
            NSAME = 0
            DO 30 N = 1, MAXSP
               KI = NUNK(N,I)
               NI = NU(N,I)
C
               DO 25 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. -NJ.EQ.NI)
     1            NSAME = NSAME + 1
   25          CONTINUE
   30       CONTINUE
C
            IF (NSAME.EQ.ABS(NS(J)) .AND.
     1          (NS(J).GT.0 .OR. NS(I).GT.0)) THEN
C
C              same products as J reactants, and vice-versa
C
               ISAME = J
               RETURN
            ENDIF
         ENDIF
C
  500 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKPRNT (LOUT, MAXTP, KKTOT, KNAME, ITHRM, T, NT, TMID,
     1                   WTM, PDEN, LPDEN, DEN, LKDEN, NELEM, MDIM,
     2                   KNCF, KCOV, NFSURF, NLSURF, NNSURF, NFBULK,
     3                   NLBULK, NNBULK, KFIRST, KLAST, KKPHAS, PNAME,
     4                   KERR)
C
C     Print surface and solid species data
C
C     INPUT:
C        LOUT     - unit number for output messages
C        MAXTP    - maximum number of temperatures allowed for
C                   thermodynamic fits
C        KKTOT    - total number of species
C        KNAME(*) - character array of species names
C        ITHRM(*) - logical array to indicate presence of thermodynamic
C                   data of the species
C        T(*,*)   - real array of temperatures used in thermodynamic
C                   fits for the species
C        NT(*)    - integer array of the number of temperatures used
C        WTM(*)   - real array of species molecular weights
C        PDEN(*)  - real array of phase densities
C        DEN(*)   - real array of species densities
C        NELEM    - total number of elements declared
C        MDIM     - maximum number of elements allowed
C        KNCF(*,*)- matrix of elemental composition of the species
C        KCOV(*)  - integer array of site coverage of the species
C        NFSURF   - the number of the first phase which is a site
C        NLSURF   - the number of the last phase which is a site
C        NFBULK   - the number of the first phase which is a bulk
C        NLBULK   - the number of the last phase which is a bulk
C        KFIRST(*)- integer array of starting species of the phases
C        KLAST(*) - integer array of ending species of the phases
C        PNAME(*) - character array of phase names
C
C     OUTPUT:
C        KERR     - logical error flag
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WTM(*), PDEN(*), DEN(*), KNCF(MDIM,*), KCOV(*),
     1          T(MAXTP,*), NT(*), KFIRST(*), KLAST(*), KKPHAS(*)
      CHARACTER*(*) PNAME(*), KNAME(*)
      LOGICAL KERR, ITHRM(*), LPDEN(*), LKDEN(*)
C
      IF (NNSURF .GT. 0) THEN
         DO 10 J = NFSURF,NLSURF
C
            IF (LPDEN(J)) THEN
               WRITE (LOUT,'(/2(1X,A),10X,E10.3,A)')
     1         'SITE:', PNAME(J), PDEN(J), ' moles/cm**2'
               IF (PDEN(J) .LE. 0.0) THEN
                  WRITE (LOUT, 1700)
                  KERR = .TRUE.
               ENDIF
            ELSE
               WRITE (LOUT,'(/2(1X,A))') 'SITE:', PNAME(J)
               WRITE (LOUT, 1710)
               KERR = .TRUE.
            ENDIF
C
            IF (KKPHAS(J) .LE. 0) THEN
               WRITE (LOUT, 1720)
               KERR = .TRUE.
            ELSE
               DO 5 K = KFIRST(J), KLAST(J)
                  WRITE (LOUT, 1655) K, KNAME(K), WTM(K), KCOV(K),
     1                               (KNCF(M,K),M=1,NELEM)
                  IF (KCOV(K) .LE. 0) THEN
                     WRITE (LOUT, 1660)
                     KERR = .TRUE.
                  ENDIF
                  CALL SKSPEC (K, KNAME, KKTOT, ITHRM, MAXTP, T, NT,
     1                         TMID, LOUT, KERR)
    5          CONTINUE
            ENDIF
   10    CONTINUE
      ENDIF
C
      IF (NNBULK .GT. 0) THEN
         DO 20 J = NFBULK, NLBULK
            WRITE (LOUT, '(/2(1X,A))') 'BULK:', PNAME(J)
            IF (KKPHAS(J) .LE. 0) THEN
               WRITE (LOUT, 1730)
               KERR = .TRUE.
            ELSE
               DO 15 K = KFIRST(J), KLAST(J)
                  WRITE (LOUT, 1660) K, KNAME(K), WTM(K),
     1                               DEN(K), (KNCF(M,K), M=1,NELEM)
                  IF (LKDEN(K) .AND. DEN(K).LE.0.0) THEN
                     WRITE (LOUT, 1740)
                     KERR = .TRUE.
                  ENDIF
                  CALL SKSPEC (K, KNAME, KKTOT, ITHRM, MAXTP, T, NT,
     1                         TMID, LOUT, KERR)
   15          CONTINUE
            ENDIF
   20    CONTINUE
      ENDIF
C
      RETURN
C
 1655 FORMAT (I4,'. ', A16, F11.5, 17X, I3, 3X, 10I3)
 1660 FORMAT (I4,'. ', A16, F11.5, E10.3,' g/cm**3', 5X, 10I3)
 1700 FORMAT (6X,'Error...site phase density must be > 0')
 1710 FORMAT (6X,'Error...site phase density is required input')
 1720 FORMAT (6X,'Error...site phase has no species')
 1730 FORMAT (6X,'Error...bulk phase has no species')
 1740 FORMAT (6X,'Error...bulk species density must be > 0')
C
      END
C
      SUBROUTINE SKUNIT (LINE, AUNITS, EUNITS, IUNITS, NONCON,
     1                   LMTZWS)
C
C     Determines units of Arrhenius parameters, and establishes
C     whether non-conservation of sites reactions are allowed.
C
C     INPUT:
C        LINE      - character string
C
C     OUTPUT:
C        AUNITS    - character flag to indicate units of A,
C                    the pre-exponential factor
C        EUNITS    - character flag to indicate units of E,
C                    the activation energy
C        IUNITS    - character string describing the input units
C        NONCON    - logical flag to indicate whether or not non-
C                    conservation of sites reactions are allowed
C        LMTZWS    - logical flag to allow modification of sticking
C                    coefficients according to Motz-Wise formula
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) LINE, IUNITS, EUNITS, AUNITS
      CHARACTER*4 UPCASE
      LOGICAL NONCON, LMTZWS
C
      IUNITS = ' '
      EUNITS = ' '
      AUNITS = ' '
      DO 85 N = 1, ILASCH(LINE)-3
         IF (UPCASE(LINE(N:), 4) .EQ. 'NONC') THEN
            NONCON = .TRUE.
         ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'MWOF') THEN
            LMTZWS = .FALSE.
         ENDIF
         IND = ILASCH(IUNITS)
         IF (EUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4) .EQ. 'CAL/') THEN
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units cal/mole'
               ELSE
                  IUNITS(IND:) = ', E units cal/mole'
               ENDIF
               EUNITS = 'CAL/'
c JM: change line, use four first letters and not five
c            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KCAL/') THEN
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KCAL') THEN
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kcal/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kcal/mole'
               ENDIF
               EUNITS = 'KCAL'
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'JOUL') THEN
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Joules/mole'
               ELSE
                  IUNITS(IND:) = ', E units Joules/mole'
               ENDIF
               EUNITS = 'JOUL'
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KJOU') THEN
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kjoule/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kjoule/mole'
               ENDIF
               EUNITS = 'KJOU'
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KELV') THEN
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kelvins'
               ELSE
                  IUNITS(IND:) = ', E units Kelvins'
               ENDIF
               EUNITS = 'KELV'
            ENDIF
         ENDIF
         IF (AUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4) .EQ. 'MOLE') THEN
               IF (N+4.LE.ILASCH(LINE) .AND.
     1             UPCASE(LINE(N+4:), 1) .EQ. 'S') THEN
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units moles'
                  ELSE
                     IUNITS(IND:) = ', A units moles'
                  ENDIF
                  AUNITS = 'MOLE'
               ELSEIF (N+7.LE.ILASCH(LINE) .AND.
     1            UPCASE(LINE(N+4:), 4) .EQ. 'CULE') THEN
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units molecules'
                  ELSE
                     IUNITS(IND:) = ', A units molecules'
                  ENDIF
                  AUNITS = 'MOLC'
               ENDIF
            ENDIF
         ENDIF
   85 CONTINUE
C
      IF (AUNITS .EQ. ' ') THEN
         AUNITS = 'MOLE'
         IND = MAX(ILASCH(IUNITS), 1)
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', A units moles'
         ELSE
            IUNITS(IND:) = ' A units moles'
         ENDIF
      ENDIF
C
      IF (EUNITS .EQ. ' ') THEN
         EUNITS = 'CAL/'
         IND = MAX(ILASCH(IUNITS), 1)
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', E units cal/mole'
         ELSE
            IUNITS(IND:) = ' E units cal/mole'
         ENDIF
      ENDIF
C
      RETURN
      END
C
C----------------------------------------------------------------------C
      SUBROUTINE SKAUXL (LINE, SUB, NSUB, II, LOUT, NSPEC, NREAC, NU,
     1                   NUNK, KKGAS, KKTOT, KNAME, NPHASE, PNAME,
     2                   NFSURF, NLSURF, KFIRST, KLAST, KKPHAS,
     3                   NPAR, NSCOV, IDUP, NREV, IREV, RPAR, NCOV,
     4                   ICOV, IKCOV, CPAR, NSTICK, ISTICK, NIIRNU,
     5                   IRNU, MAXSP, RNU, NIIORD,
     6                   IORD, MAXORD, KORD, RORD, KERR,
     7                   NIIBHM, IBOHM, IBK, IBT)
C
C     SKAUXL parses the auxiliary CHAR*(*) lines representing
C     additional options for a surface reaction; data is stored
C     based on finding a 'keyword' followed by its required
C     parameters:
C
C     REV[ERSE]/val1 val2 val3/ :
C        if IREV(NREV)=I, this is an error, REV already declared;
C        if NSPEC(I)<0, this an error, as reaction I is irreversible;
C        else increment NREV, the number of reactions with reverse
C             parameters given,
C             IREV(NREV)=I, RPAR(N,NREV)=val(N),N=1,3;
C
C     DUP[LICATE]:
C        This reaction is allowed to be duplicated.
C
C     COV[ERAGE]/species val1 val2 val3/:
C
C     STICK:
C
C     INPUT:
C        LINE      - character string
C        SUB(*)    - character array of substrings
C        NSUB      - total number of substring in the line
C        II        - total number of reactions, and the index number
C                    of the current reaction
C        LOUT      - unit number for output messages
C        NSPEC     - total number of species in the current reaction
C        NREAC     - total number of reactants in the current reaction
C        NU(*)     - the NREAC species stoichiometric
C        NUNK(*)   - the NREAC species indices
C        KKGAS     - total number of gas-phase species in the problem
C        KKTOT     - total number of species in the problem
C        KNAME(*)  - character array of species names
C        NFSURF    - index number of the first phase which is a site
C        NLSURF    - index number of the last phase which is a site
C        KFIRST(*) - integer array of starting species of the phases
C        KLAST(*)  - integer array of ending species of the phases
C        KKPHAS(*) - integer array of total species in each phases
C        NPAR      - the number of Arrhenius parameters required for
C                    a reaction
C        NSCOV     - the number of coverage parameters required in
C                    a coverage declaration
C     OUTPUT:
C        IDUP(*)   - integer array of duplicate status of the reactions;
C                    = -1 indicates a legal duplicate
C        NREV      - total number of reactions which declared reverse
C                    Arrhenius parameters
C        IREV(*)   - the NREV reaction indices
C        RPAR(*,*) - real matrix of Arhennius parameters for the NREV
C                    reactions
C        NCOV      - total number of coverage declarations
C        ICOV(*)   - the NCOV reaction indices
C        IKCOV(*)  - integer array of species indices for the
C                    coverage declarations
C        CPAR(*,*) - real matrix of coverage parameters for the
C                    species in the NCOV coverage declarations
C        NSTICK    - total number of reactions with sticking coeff'nts
C        ISTICK(*) - the NSTICK reaction indices
C        KERR      - logical error flag
C
C                                        F. Rupley, Div. 8245, 5/27/87
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NU(*), NUNK(*), IDUP(*), IREV(*), RPAR(NPAR,*),
     1          ICOV(*), IKCOV(*), CPAR(NSCOV,*), ISTICK(*),
     2          KFIRST(*), KLAST(*), KKPHAS(*), IORD(*),
     3          KORD(MAXORD,*), RORD(MAXORD,*), IRNU(*),
     4          RNU(MAXSP,*), IBOHM(*), IBK(*), IBT(*)
      CHARACTER*(*) SUB(*), KNAME(*), PNAME(*), LINE
      CHARACTER RSTR*80, UPCASE*4
      LOGICAL KERR, IERR, LREV, LFORD, LRORD
C
      LREV = (NREV.GT.0 .AND. IREV(NREV).EQ.II)
C
      LSTART = IFIRCH(LINE)
      LEND   = ILASCH(LINE)
      DO 500 L = LSTART, LEND-2
C
         IF (UPCASE(LINE(L:), 3) .EQ. 'DUP') THEN
            IDUP(II) = -1
            WRITE (LOUT, 4000)
C
         ELSEIF (UPCASE(LINE(L:), 3) .EQ. 'REV') THEN
C
C        REVERSE ARRHENIUS PARAMETERS
C
            IF (LREV .OR. NSPEC.LT.0) THEN
               KERR = .TRUE.
               IF (LREV) WRITE (LOUT, 2050) LINE(L:LEND)
               IF (NSPEC .LT. 0) WRITE (LOUT, 2060) LINE(L:LEND)
            ELSE
               LREV = .TRUE.
               NREV = NREV + 1
               IREV(NREV) = II
               I1 = INDEX(LINE(L:), '/')
               I2 = INDEX(LINE(L+I1+1:), '/')
               IF (I1.LE.0 .OR. I2.LE.0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2085) LINE(L:LEND)
                  GO TO 500
               ENDIF
C
               RSTR = LINE(L+I1:L+I1+I2-1)
               CALL CKXNUM (RSTR, NPAR, LOUT, NVAL, RPAR(1,NREV), IERR)
               KERR = KERR.OR.IERR
               WRITE (LOUT, 1900) '   Reverse Arrhenius coefficients:',
     1                            (RPAR(J,NREV), J=1,3)
               IF (RPAR(1,NREV) .LT. 0) THEN
                  WRITE (LOUT, 4010)
                  KERR = .TRUE.
               ENDIF
C
            ENDIF
C
         ELSEIF (UPCASE(LINE(L:), 3) .EQ. 'COV') THEN
               NCOV = NCOV+1
               ICOV(NCOV) = II
               I1 = INDEX(LINE(L:), '/')
C
               IF (I1 .LE. 0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2090) LINE(L:LEND)
                  GO TO 500
               ENDIF
C
               CALL SKISUB (LINE(L+I1:), SUB, NSUB)
               IF (NSUB .LT. 4) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2090) LINE(L:LEND)
                  GO TO 500
               ENDIF
C
               CALL SKPCMP (SUB(1), KNAME, KKTOT, PNAME,
     1                      NPHASE, KKPHAS, KNUM, NF)
C
               IF (KNUM .LE. 0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 4001) LINE(L:LEND)
               ELSEIF (KNUM.LT.KFIRST(NFSURF) .OR.
     1                 KNUM.GT.KLAST(NLSURF)) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 4005) LINE(L:LEND)
               ENDIF
C
               IND = ILASCH(SUB(4))
               IF (SUB(4)(IND:IND) .EQ. '/') SUB(4)(IND:) = ' '
               RSTR = ' '
               RSTR = SUB(2)
               DO 25 J = 3, 4
                  IND = ILASCH(RSTR) + 2
                  RSTR(IND:) = SUB(J)
   25          CONTINUE
C
               CALL IPPARR (RSTR, 1, NSCOV, CPAR(1,NCOV), NVAL, IERR,
     1                      LOUT)
               IF (IERR) THEN
                  KERR = .TRUE.
               ELSE
                  IKCOV(NCOV) = KNUM
                  WRITE (LOUT,3035) KNAME(KNUM),(CPAR(J,NCOV),J=1,NSCOV)
               ENDIF
C            ENDIF
C
         ELSEIF (UPCASE(LINE(L:), 3) .EQ. 'STI') THEN
C
            WRITE (LOUT, 3050)
            NSTICK = NSTICK + 1
            ISTICK(NSTICK) = II
C
C           there must be only one gas-phase reactant
C
            NGAS = 0
            DO 400 J = 1, NREAC
               IF (NUNK(J) .LE. KKGAS) THEN
                  NGAS = NGAS + 1
                  IF (ABS(NU(J)) .GT. 1) THEN
                     WRITE (LOUT, 3042) KNAME(NUNK(J))
                     KERR = .TRUE.
                  ENDIF
                  IF (NIIRNU .GT. 0) THEN
                      IF (IRNU(NIIRNU).EQ.II .AND.
     1                RNU(J,NIIRNU).GT.1.0) THEN
                         WRITE (LOUT, 3042) KNAME(NUNK(J))
                         KERR = .TRUE.
                      ENDIF
                  ENDIF
               ENDIF
  400       CONTINUE
            IF (NGAS .NE. 1) THEN
               KERR = .TRUE.
               IF (NGAS .EQ. 0) THEN
                  WRITE (LOUT, 3045)
               ELSE
                  WRITE (LOUT, 3040)
               ENDIF
            ENDIF
C
         ELSEIF (UPCASE(LINE(L:), 3) .EQ. 'BOH') THEN
C
            WRITE (LOUT, 3052)
            NIIBHM = NIIBHM + 1
            IBOHM(NIIBHM) = II
            I1 = INDEX(LINE(L:), '/')
C
            IF (I1 .LE. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 2096) LINE(L:LEND)
               GO TO 500
            ENDIF
C
            CALL SKISUB (LINE(L+I1:), SUB, NSUB)
            IF (NSUB .LT. 2) THEN
               KERR = .TRUE.
               WRITE (LOUT, 2096) LINE(L:LEND)
               GO TO 500
            ENDIF
C
            CALL SKPCMP (SUB(1), KNAME, KKTOT, PNAME,
     1                      NPHASE, KKPHAS, KNUM, NF)
C
            IF (KNUM .LE. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 4015) LINE(L:LEND)
            ELSEIF (KNUM.LT.1 .OR.
     1                 KNUM.GT.KKGAS) THEN
               KERR = .TRUE.
               WRITE (LOUT, 4020) LINE(L:LEND)
            ENDIF
C
            IND = ILASCH(SUB(2))
            IF (SUB(2)(IND:IND) .EQ. '/') SUB(2)(IND:) = ' '
            CALL CKXNUM (SUB(2), 1, LOUT, NVAL, VALUE, IERR)
            IBT(NIIBHM) = INT(VALUE)
            IF (IERR .OR. IBT(NIIBHM) .GT. KKGAS) THEN
               KERR = .TRUE.
               WRITE (LOUT,3060) KKGAS
            ELSE
               IBK(NIIBHM) = KNUM
               WRITE (LOUT,3055) KNAME(KNUM), IBT(NIIBHM)
            ENDIF
C
C    the species specified to control the rate must be a reactant
C    with a coefficient of 1.
C
            LHS = 0
            DO 40 J = 1, NREAC
               IF (NUNK(J) .EQ. KNUM) THEN
                  LHS = 1
                  IF (ABS(NU(J)) .NE. 1) THEN
                     WRITE (LOUT, 3058) KNAME(KNUM)
                     KERR = .TRUE.
                  ENDIF
                  IF (NU(J) .GT.0) THEN
                     WRITE (LOUT, 3059) KNAME(KNUM)
                     KERR = .TRUE.
                  ENDIF
               ENDIF
 40         CONTINUE
            IF (LHS .NE. 1) THEN
               WRITE (LOUT, 3059) KNAME(KNUM)
               KERR = .TRUE.
            ENDIF
C
         ELSE
            LFORD = (UPCASE(LINE(L:),4) .EQ. 'FORD')
            LRORD = (UPCASE(LINE(L:),4) .EQ. 'RORD')
            IF (LFORD .OR. LRORD) THEN
               IF (LRORD .AND. NSPEC.LT.0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2065)
               ELSE
               IF (NIIORD.EQ.0 .OR.(NIIORD.GT.0.AND.IORD(NIIORD).NE.II))
     1         THEN
                  NIIORD = NIIORD + 1
                  IORD(NIIORD) = II
                  NKORD = 0
                  IF (NIIRNU.GT.0 .AND. IRNU(NIIRNU).EQ.II) THEN
                     DO 111 M = 1, 12
                        IF (NUNK(M).NE.0) THEN
                           NKORD = NKORD + 1
                           IF (RNU(M,NIIRNU) .LT. 0) THEN
                              KORD(NKORD,NIIORD) = -NUNK(M)
                              RORD(NKORD,NIIORD) = ABS(RNU(M,NIIRNU))
                           ELSE
                              KORD(NKORD,NIIORD) = NUNK(M)
                              RORD(NKORD,NIIORD) = RNU(M,NIIRNU)
                           ENDIF
                        ENDIF
  111                CONTINUE
                  ELSE
                     DO 113 M = 1, 12
                        IF (NUNK(M) .NE. 0) THEN
                           NKORD = NKORD + 1
                           IF (NU(M) .LT. 0) THEN
                              KORD(NKORD,NIIORD) = -NUNK(M)
                              RORD(NKORD,NIIORD) = IABS(NU(M))
                           ELSE
                              KORD(NKORD,NIIORD) = NUNK(M)
                              RORD(NKORD,NIIORD) = NU(M)
                           ENDIF
                        ENDIF
  113                CONTINUE
                  ENDIF
               ENDIF
               ENDIF
C
               I1 = L + INDEX(LINE(L:),'/')
               IF (I1 .LE. L) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2095) LINE(L:LEND)
                  GO TO 500
               ENDIF
               I2 = I1 + INDEX(LINE(I1+1:),'/') - 1
               CALL SKISUB (LINE(I1:I2), SUB, NSUB)
               IF (NSUB .LT. 2) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 2095) LINE(L:LEND)
                  GO TO 500
               ENDIF
C
C JM insert:  Change option to have change of raction order
C             only for gas phase species and extend it to 
C             all species, gaseous and surface ones.  
C              CALL CKCOMP (SUB(1), KNAME, KKGAS, K)
               CALL SKCOMP (SUB(1), KNAME, KKTOT, K, K_N)
               IF (K .LE. 0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 4002) LINE(L:LEND)
               ENDIF
C
               IF (LFORD) K = -K
               NK = 0
               DO 121 M = 1, MAXORD
                  IF (KORD(M,NIIORD).EQ.0) THEN
                     NK = M
                     GO TO 122
                  ELSEIF (KORD(M,NIIORD).EQ.K) THEN
                     IF (LFORD) THEN
                        WRITE (LOUT,*)
     1'      Warning...changing order for reactant...',
     2                  KNAME(-K)
                     ELSE
                        WRITE (LOUT,*)
     1'      Warning...changing order for product...',
     2                  KNAME(K)
                     ENDIF
                     NK = M
                     GO TO 122
                  ENDIF
  121          CONTINUE
  122          CONTINUE
C
               IND = ILASCH(SUB(2))
               IF (SUB(2)(IND:IND) .EQ. '/') SUB(2)(IND:)=' '
               CALL IPPARR (SUB(2), 1, 1, VAL, NVAL, IER, LOUT)
C
               KORD(NK,NIIORD) = K
               RORD(NK,NIIORD) = VAL
               IF (LFORD) THEN
                  WRITE (LOUT, 3015) KNAME(-K),VAL
               ELSE
                  WRITE (LOUT, 3016) KNAME(K),VAL
               ENDIF
            ENDIF
         ENDIF
  500 CONTINUE
C
C     FORMATS
C
 2050 FORMAT (6X,'Error...REV declared more than once...',A)
 2060 FORMAT (6X,'Error...REV declared for irreversible reaction...',A)
 2065 FORMAT (6X,'Error...RORD declared for irreversible reaction...')
 2085 FORMAT (6X,'Error in REV declaration...',A)
 2090 FORMAT (6X,'Error in COV declaration...',A)
 2095 FORMAT (6X,'Error in ORDER declaration...',A)
 2096 FORMAT (6X,'Error in BOHM declaration...',A)
 1900 FORMAT (6X,A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 3015 FORMAT (7X,A16,' Forward order ',1PE12.3)
 3016 FORMAT (7X,A16,' Reverse order ',1PE12.3)
 3035 FORMAT (9X,'Coverage parameters for species ',A,': ',/,10X,
     1            3(1PE10.3))
 3040 FORMAT (6X,'Error...only one gas-phase reactant allowed with',
     1           ' sticking option...')
 3042 FORMAT (6X,'Error...gas-phase reactant ',A,' must not have a ',
     1           'coefficient greater than 1...')
 3045 FORMAT (6X,'Error...no gas-phase reactant with ',
     1           'sticking option...')
 3050 FORMAT (9X,'Coefficients are sticking parameters...')
 3052 FORMAT (9X,'Coefficients are Bohm sticking parameters...')
 3055 FORMAT (9X,'Bohm temperature dependence for species ',A,':',I3)
 3058 FORMAT (9X,'Error...gas-phase reactant ',A,' must have a',
     1           'coefficient of 1...')
 3059 FORMAT (9X,'Error...the Bohm species ',A,' is not a reactant')
 3060 FORMAT (9x,'Error...the temperature dependence flag must be ',
     1           'less than ',I5)
 4000 FORMAT (6X,'Declared duplicate reaction...')
 4001 FORMAT (6X,'Error...unrecognized COV species...',A)
 4002 FORMAT (6X,'Error...unrecognized ORDER species...',A)
 4005 FORMAT (6X,'Error...COV species must be site-phase...',A)
 4010 FORMAT (6X,'Error...pre-exponential factor must be positive...')
 4015 FORMAT (6x,'Error...unrecognized BOH species...',A)
 4020 FORMAT (6x,'Error...BOH species must be gas-phase...',A)
C
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SPREAC (II, MAXSP, KKGAS, KKSUR, NSPEC, NREAC, NPAR,
     1                   PAR, RPAR, AUNITS, EUNITS, NUNK, NU, RNCF,
     2                   KCHRG, NSTK, ISTK, NCOV, ICOV, NSCOV, CPAR,
     3                   KCOV, NREV, IREV, MDIM, MM, KNCF, NONCON, NCON,
     4                   IDUP, NNSURF, NPHASE, KFIRST, KLAST, LOUT,
     5                   KERR, NIIRNU, IRNU, RNU, SKMIN)
C
C     Checks for reaction balance and duplication, conversion of units
C
C     INPUT:
C        II        - index of the current reaction
C        MAXSP     - maximum number of reactants+products allowed
C        KKGAS     - total number of gas-phase species
C        KKSUR     - total number of site-phase species
C        NSPEC(*)  - array of the number of species in the reactions
C        NREAC(*)  - array of the number of reactants in the reactions
C        NPAR      - number of Arrhenius parameters, reverse
C                    Arrhenius parameters, or sticking coefficients
C                    required for a reaction
C        PAR(*,*)  - matrix of Arrhenius parameters (or sticking
C                    coefficients) for the reactions
C        RPAR(*,*) - matrix of reverse parameters for the reactions
C                    which have declared them
C        AUNITS    - character string which describes the input
C                    units of A, the pre-exponential factor PAR(1,I)
C        EUNITS    - character string which describes the input
C                    units of E, the activation energy PAR(3,I)
C        NUNK(*,*) - matrix of the species indices of the reactants
C                    and products in the reactions
C        NU(*,*)   - matrix of the stoichiometric coefficients of
C                    the reactants and products in the reactions
C        RNCF(*)   - array of the net change in sites of reaction II
C                    for the phases
C        KCHRG(*)  - array of electronic charges for the species
C        NSTK      - total number of reactions with sticking coeff'nts
C        ISTK(*)   - the NSTK reaction indices
C        NCOV      - total number of coverage declarations
C        ICOV(*)   - the NCOV reaction indices
C        NSCOV     - number of coverage parameters required
C        CPAR(*,*) - matrix of the coverage parameters
C        KCOV(*)   - array of site coverage parameters for the species
C        NREV      - total number of reactions with reverse
C                    parameters defined
C        IREV(*)   - the NREV reaction indices
C        MDIM      - maximum number of elements allowed
C        MM        - actual number of elements declared
C        KNCF(*,*) - matrix of elemental composition of the species
C        NONCON    - logical flag to allow non-conservation of sites
C        NCON      - total number of reactions which do not conserve
C                    sites
C        IDUP(*)   - integer array of flags to indicate duplicate
C                    reactions
C        NNSURF    - total number of sites
C        NPHASE    - total number of phases
C        KFIRST(*) - array of the first species indices of the phases
C        KLAST(*)  - array of the final species indices of the phases
C        LOUT      - unit number for output messages
C     OUTPUT:
C        KERR      - logical error flag
C
C           (Value of Avrogadro's Constant and Gas Constant
C            from 1986 CODATA recommended values (1993 CRC)
C            J. Research National Bureau of Standards, 92, 95, 1987
C            6.0221367(39)E23 mol-1
C            8.314510 Joules / (Mole K)       )
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
       DOUBLE PRECISION    RU_JOUL,              AVAG
       PARAMETER          (RU_JOUL = 8.314510D0, AVAG = 6.0221367D23)
C*****END precision > double
C*****precision > single
C       IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C       REAL               RU_JOUL,              AVAG
C       PARAMETER         (RU_JOUL = 8.314510E0, AVAG = 6.0221367E23)
C*****END precision > single
C
      DIMENSION PAR(NPAR,*), RPAR(NPAR,*), NSPEC(*), NREAC(*),
     1          NUNK(MAXSP,*), NU(MAXSP,*), RNCF(*), KCHRG(*), ISTK(*),
     2          ICOV(*), CPAR(NSCOV,*), KCOV(*), IREV(*), KNCF(MDIM,*),
     3          IDUP(*), KFIRST(*), KLAST(*), IRNU(*), RNU(MAXSP,*)
      CHARACTER*(*) AUNITS, EUNITS
      LOGICAL IERR, KERR, NONCON
C
      IERR = .FALSE.
      IF (NIIRNU .GT. 0) THEN
         IF (IRNU(NIIRNU) .EQ. II)
     *   CALL SKRBAL (MAXSP, NUNK(1,II), RNU(1,NIIRNU), RNCF, MDIM, MM,
     1                KNCF, KCOV, NONCON, NCON, KCHRG, NNSURF, NPHASE,
     2                KFIRST, KLAST, IERR, SKMIN)
      ELSE
         CALL SKBAL (MAXSP, NUNK(1,II), NU(1,II), RNCF, MDIM, MM,
     1         KNCF, KCOV, NONCON, NCON, KCHRG, NNSURF, NPHASE,
     2         KFIRST, KLAST, IERR, SKMIN)
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1060)
      ENDIF
C
      CALL SKDUP (II, MAXSP, NSPEC, NREAC, NU, NUNK, ISAME)
C
      IF (ISAME .GT. 0) THEN
         IF (IDUP(ISAME).NE.0 .AND. IDUP(II).NE.0) THEN
            IDUP(ISAME) = ABS(IDUP(ISAME))
            IDUP(II)    = ABS(IDUP(II))
         ELSE
            KERR = .TRUE.
            WRITE (LOUT, 1050) ISAME
         ENDIF
      ENDIF
C
      IF (EUNITS .EQ. 'KELV') THEN
         EFAC = 1.0
      ELSEIF (EUNITS .EQ. 'CAL/') THEN
C        convert E from cal/mole to Kelvin
         EFAC = 4.184  / RU_JOUL
      ELSEIF (EUNITS .EQ. 'KCAL') THEN
C        convert E from kcal/mole to Kelvin
         EFAC = 4184.0 / RU_JOUL
      ELSEIF (EUNITS .EQ. 'JOUL') THEN
C        convert E from Joules/mole to Kelvin
         EFAC = 1.0    / RU_JOUL
      ELSEIF (EUNITS .EQ. 'KJOU') THEN
C        convert E from Kjoule/mole to Kelvin
         EFAC = 1000.0 / RU_JOUL
      ENDIF
C
C JM insert
      PAR(3,II) = PAR(3,II) * EFAC
c     write(6,*) 'Reaction, E(K)  ',II, PAR(3,II)
      IF (NREV .GT. 0) THEN
        IF (IREV(NREV) .EQ. II) RPAR(3,NREV) = RPAR(3,NREV) * EFAC
      ENDIF
C
      NSTOR = 0
      NSTOP = 0
      DO 50 N = 1, MAXSP
         IF (NUNK(N,II) .LE. KKGAS+KKSUR) THEN
            IF (NU(N,II) .LT. 0) THEN
C              sum of stoichiometric coefficients of reactants
               NSTOR = NSTOR + ABS(NU(N,II))
            ELSEIF (NU(N,II) .GT. 0) THEN
C              sum of stoichiometric coefficients of products
               NSTOP = NSTOP + NU(N,II)
            ENDIF
         ENDIF
   50 CONTINUE
C
      IF (AUNITS .EQ. 'MOLC') THEN
         IF (NSTK.GT.0 .AND. ISTK(NSTK).EQ.II) THEN
C           no unit conversion of A if sticking coefficient (no units)
         ELSE
C           convert A from molecules to moles
            IF (NSTOR .GT. 0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
         ENDIF

         IF (NREV.GT.0 .AND. IREV(NREV).EQ.II .AND. NSTOP.GT.0)
     1      RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)

      ENDIF
C
      DO 100 N = 1, NCOV
C
C        conversion of coverage parameter units
C
         IF (ICOV(N) .EQ. II) THEN
            CPAR(3,N) = CPAR(3,N) * EFAC
         ENDIF
  100 CONTINUE
C
      RETURN
C
 1050 FORMAT (6X,'Error...undeclared duplicate to reaction number ',I3)
 1060 FORMAT (6X,'Error...reaction does not balance...')
      END
C
C----------------------------------------------------------------------C
      SUBROUTINE SKISUB (LINE, SUB, NSUB)
C
C     Generates an array of CHAR*(*) substrings from a CHAR*(*) line,
C     using blanks as a delimiter
C
C     Input:  LINE  - a CHAR*(*) line
C     Output: SUB   - a CHAR*(*) array of substrings
C             NSUB  - number of substrings found
C     A '!' will comment out a line, or remainder of the line.
C                                      F. Rupley, Div. 8245, 5/15/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) SUB(*), LINE
      NSUB = 0
      IF (IPPLEN(LINE) .LE. 0) RETURN
C
      ILEN = ILASCH(LINE)
C
      NSTART = IFIRCH(LINE)
   10 CONTINUE
      ISTART = NSTART
      NSUB = NSUB + 1
      SUB(NSUB) = ' '
C
c      DO 100 I = ISTART, ILEN
         ILAST = INDEX(LINE(ISTART:),' ') - 1
         IF (ILAST .GT. 0) THEN
            ILAST = ISTART + ILAST - 1
         ELSE
            ILAST = ILEN
         ENDIF
         SUB(NSUB) = LINE(ISTART:ILAST)
         IF (ILAST .EQ. ILEN) RETURN
C
         NSTART = ILAST + IFIRCH(LINE(ILAST+1:))
C
C        Does SUB have any slashes?
C
         I1 = INDEX(SUB(NSUB),'/')
         IF (I1 .LE. 0) THEN
            IF (LINE(NSTART:NSTART) .NE. '/') GO TO 10
            NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
            IND = INDEX(SUB(NSUB),' ')
            SUB(NSUB)(IND:) = LINE(NSTART:NEND)
            IF (NEND .EQ. ILEN) RETURN
            NSTART = NEND + IFIRCH(LINE(NEND+1:))
            GO TO 10
         ENDIF
C
C        Does SUB have 2 slashes?
C
         I2 = INDEX(SUB(NSUB)(I1+1:),'/')
         IF (I2 .GT. 0) GO TO 10
C
         NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
         IND = INDEX(SUB(NSUB),' ') + 1
         SUB(NSUB)(IND:) = LINE(NSTART:NEND)
         IF (NEND .EQ. ILEN) RETURN
         NSTART = NEND + IFIRCH(LINE(NEND+1:))
         GO TO 10
c  100 CONTINUE
c      RETURN
      END
C----------------------------------------------------------------------C
      CHARACTER*(*) FUNCTION UPCASE(STR, ILEN)
      CHARACTER*(*) STR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      UPCASE = ' '
      UPCASE = STR(:ILEN)
      JJ = MIN (LEN(UPCASE), LEN(STR), ILEN)
      DO 10 J = 1, JJ
         DO 10 N = 1, 26
            IF (STR(J:J) .EQ. LCASE(N)) UPCASE(J:J) = UCASE(N)
   10 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKBULK (SUB, NSUB, NDIM, STRAY, RAY, LRAY, NN, KCHECK,
     1                   KSTART, KERR, LOUT)
C
C     Extracts names and real values from an array of CHAR*(*)
C     substrings; stores names in STRAY array, real values in RAY;
C     i.e. can be used to store element and atomic weight data,
C     species names, etc.
C
C     Input:   SUB(N),N=1,NSUB  - array of CHAR*(*) substrings
C              NSUB             - number of substrings
C              NDIM             - size of STRAY,RAY arrays
C              NN               - actual number of STRAY found
C              STRAY(N),N=1,NN  - CHAR*(*) array
C              RAY(N),N=1,NN    - Real array
C              LOUT             - output unit for error messages
C     Output:  NN               - incremented if more STRAY found
C              STRAY(N),N=1,NN  - incremented array of STRAY
C              RAY(N),N=1,NN    - incremented array of reals
C              KERR             - logical, .TRUE. = error in data
C
C                                       F. Rupley, Div. 8245, 2/5/88
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RAY(*), PAR(1)
      CHARACTER*(*) SUB(*), STRAY(*)
      CHARACTER ISTR*80, UPCASE*4
      LOGICAL KERR, LRAY(*)
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF (UPCASE(SUB(N), 3) .EQ. 'END') RETURN
C
         ISTR = ' '
         I1 = INDEX(SUB(N),'/')
         IF (I1 .EQ .1) THEN
            KERR = .TRUE.
            WRITE (LOUT, 130) SUB(N)(:ILASCH(SUB(N)))
         ELSE
            IF (I1 .LE. 0) THEN
               ISTR = SUB(N)
            ELSE
               ISTR = SUB(N)(:I1-1)
            ENDIF
            CALL SKCOMP (ISTR, STRAY, NN, KNUM1, NF)
C
            IF (NN .GT. KSTART) THEN
               CALL SKCOMP (ISTR, STRAY(KSTART), NN-KSTART, KNUM2, NF)
            ELSE
               KNUM2 = 0
            ENDIF
C
            IF (KNUM1 .GT. 0) THEN
               WRITE (LOUT, 100) SUB(N)(:ILASCH(SUB(N)))
               KERR = .TRUE.
            ELSEIF (KNUM2 .GT. 0) THEN
               WRITE (LOUT, 105) SUB(N)(:ILASCH(SUB(N)))
            ELSE
               IF (NN .LT. NDIM) THEN
                  IF (ISTR(ILEN+1:) .NE. ' ') THEN
                     WRITE (LOUT, 120) SUB(N)(:ILASCH(SUB(N)))
                     KERR = .TRUE.
                  ELSE
                     NN = NN + 1
                     STRAY(NN) = ' '
                     STRAY(NN) = ISTR(:ILEN)
                     IF (I1 .GT. 0) THEN
                        I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                        ISTR = ' '
                        ISTR = SUB(N)(I1+1:I2-1)
                        CALL IPPARR (ISTR, 1, 1, PAR, NVAL, IER, LOUT)
                        KERR = KERR .OR. (IER.NE.0)
                        RAY(NN) = PAR(1)
                        LRAY(NN) = .TRUE.
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (LOUT, 110) SUB(N)(:ILASCH(SUB(N)))
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
C
  100 FORMAT (6X,
     1'Error...bulk species name duplicates gas or site species...',A)
  105 FORMAT (6X,'Warning...duplicate bulk species name ignored...',A)
  110 FORMAT (6X,'Error...species array size too small for  ...',A)
  120 FORMAT (6X,'Error...bulk species name too long...',A)
  130 FORMAT (6X,'Error...misplaced value...',A)
      END
C----------------------------------------------------------------------C
      SUBROUTINE SKSURF (SUB, NSUB, NDIM, STRAY, IRAY, NN, KCHECK,
     1                   KSTART, LPDEN, PDEN, KERR, LOUT)
C
C     Extracts names and real values from an array of CHAR*(*)
C     substrings; stores names in STRAY array, real values in RAY;
C     i.e. can be used to store element and atomic weight data,
C     species names, etc.
C
C     Input:   SUB(N),N=1,NSUB  - array of CHAR*(*) substrings
C              NSUB             - number of substrings
C              NDIM             - size of STRAY,RAY arrays
C              NN               - actual number of STRAY found
C              STRAY(N),N=1,NN  - CHAR*(*) array
C              IRAY(N),N=1,NN   - Integer array
C              LOUT             - output unit for error messages
C     Output:  NN               - incremented if more STRAY found
C              STRAY(N),N=1,NN  - incremented array of STRAY
C              IRAY(N),N=1,NN   - incremented array of integers
C              KERR             - logical, .TRUE. = error in data
C
C                                       F. Rupley, Div. 8245, 2/5/88
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IRAY(*), IPAR(1)
      CHARACTER*(*) SUB(*), STRAY(*)
      CHARACTER ISTR*80, UPCASE*4
      LOGICAL KERR, IERR, LPDEN
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF (UPCASE(SUB(N), 3) .EQ. 'END') RETURN
         IF (UPCASE(SUB(N), 4) .EQ. 'SDEN') THEN
            I1 = INDEX(SUB(N),'/')
            IF (I1 .GT. 0) THEN
               I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
               IF (I2 .GT. I1) THEN
                  ISTR = ' '
                  ISTR = SUB(N)(I1+1:I2-1)
                  CALL CKXNUM (ISTR, 1, LOUT, NVAL, PDEN, IERR)
               ENDIF
            ENDIF
            IF (IERR) WRITE (LOUT, 1010) SUB(N)(:ILASCH(SUB(N)))
            KERR = KERR.OR.IERR
            IF (LPDEN) THEN
               WRITE (LOUT, 1020) SUB(N)(:ILASCH(SUB(N)))
            ELSE
               LPDEN = .TRUE.
            ENDIF
         ELSE
            ISTR = ' '
            I1 = INDEX(SUB(N),'/')
            IF (I1 .EQ .1) THEN
               KERR = .TRUE.
               WRITE (LOUT, 130) SUB(N)(:ILASCH(SUB(N)))
            ELSE
               IF (I1 .LE. 0) THEN
                  ISTR = SUB(N)
               ELSE
                  ISTR = SUB(N)(:I1-1)
               ENDIF
               CALL SKCOMP (ISTR, STRAY, KCHECK, KNUM1, NF)
               IF (NN .GT. KSTART) THEN
                  CALL SKCOMP (ISTR, STRAY(KSTART), NN-KSTART, KNUM2,
     1                         NF)
               ELSE
                  KNUM2 = 0
               ENDIF
C
               IF (KNUM1 .GT. 0) THEN
                  WRITE (LOUT, 100) SUB(N)(:ILASCH(SUB(N)))
                  KERR = .TRUE.
               ELSEIF (KNUM2 .GT. 0) THEN
                  WRITE (LOUT, 105) SUB(N)(:ILASCH(SUB(N)))
               ELSE
                  IF (NN .LT. NDIM) THEN
                     IF (ISTR(ILEN+1:) .NE. ' ') THEN
                        WRITE (LOUT, 120) SUB(N)(:ILASCH(SUB(N)))
                        KERR = .TRUE.
                     ELSE
                        NN = NN + 1
                        STRAY(NN) = ' '
                        STRAY(NN) = ISTR(:ILEN)
                        IF (I1 .GT. 0) THEN
                           I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                           ISTR = ' '
                           ISTR = SUB(N)(I1+1:I2-1)
                           CALL IPPARI (ISTR,1,1,IPAR,NVAL,IER,LOUT)
                           KERR = KERR .OR. (IER.NE.0)
                           IRAY(NN) = IPAR(1)
                        ENDIF
                     ENDIF
                  ELSE
                     WRITE (LOUT, 110) SUB(N)(:ILASCH(SUB(N)))
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
C
  100 FORMAT (6X,
     1'Error...site species name duplicates gas species name...',A)
  105 FORMAT (6X,
     1'Warning...duplicate species name on a site ignored...',A)
  110 FORMAT (6X,'Error...species array size too small for  ...',A)
  120 FORMAT (6X,'Error...site species name too long...',A)
  130 FORMAT (6X,'Error...misplaced value...',A)
 1010 FORMAT (6X,'Error in density declaration...',A)
 1020 FORMAT (6X,'Error...more than one density declaration...',A)
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NF)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NF)
C     Returns the first location of a character string if it also
C     occurs in a given array of character strings; also returns
C     the total number of times the character string appears in
C     the array.
C
C  INPUT
C     ISTR   - A reference character string that is to be found.
C     IRAY   - Array of character strings.
C                   Data type - character array
C                   Dimension IRAY(*) at least NN.
C     NN     - Length of IRAY(*).
C                   Data type - integer scalar
C
C  OUTPUT
C     IND    - Index of the position in IRAY(I) that matches with ISTR.
C              If no match is found, IND=0.
C                   Data type - integer scalar
C     NF     - Total number of times ISTR is found in IRAY.
C                   Data type - integer scalar
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) ISTR,IRAY(*)
      IND = 0
      NF = 0
      DO 10 I = NN, 1, -1
         IF (ISTR .EQ. IRAY(I)) THEN
            IND = I
            NF = NF + 1
         ENDIF
   10 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NF)
C
C  START PROLOGUE
C
C    SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NF)
C       Given an array of character strings IRAY with NN entries,
C       an array of character strings SETS with NSET members, and
C       a correlation between IRAY and SETS as follows:
C
C       If NSETS=2, with
C          SETS      = {"COLORS", "STONES"},
C          "COLORS"  = {"RED", "BLUE", "JADE", "RUBY"}, and
C          "STONES"  = {"TOPAZ", "JADE"},
C       then
C          IRAY = {"RED", "BLUE", "JADE", "RUBY", "TOPAZ", "JADE"},
C          and NN=6, the total number of entries.
C
C       ISET is an integer array of the number of entries of IRAY
C       in a SET;  ISET(1)=4 and ISET(2)=2.
C
C       If a given character string ISTR also occurs in IRAY, IND
C       is returned as the integer index in IRAY where it occurs;
C       if ISTR occurs more than once, IND is the integer index in
C       IRAY of its first occurrence;  NF is the total number of
C       times ISTR appears in IRAY.
C       For example, if ISTR="BLUE" then IND=2 and NF=1,
C                    if ISTR="PINK" then IND=0 and NF=0,
C                and if ISTR="JADE", IND=3 and NF=2.
C
C       However, if ISTR ends with a slash-delimited substring,
C       an additional condition for finding IND is that the
C       substring be a SET and the preceding portion of ISTR
C       be a member of that SET.
C       For example, if ISTR="BLUE/COLORS/" then IND=2 and NF=1,
C                    if ISTR="BLUE/STONES/" then IND=0 and NF=0,
C                and if ISTR="JADE/STONES/", IND=7 and NF=1.
C
C  INPUT
C     ISTR     - A character string which may or may not end with a
C                slash-delimited substring.
C     IRAY(*)  - A reference array of character strings.
C                    Data type - character array
C                    Dimension IRAY(*) at least NN
C     NN       - Number of entries in IRAY(*).
C                    Data type - integer scalar
C     SETS(*)  - A cross-reference array of character strings
C                with which a subset of IRAY(*) is associated.
C                    Data type - character array
C                    Dimension SETS(*) at least NSETS
C     NSETS   - Number of entries in NSETS(*)
C                    Data type - integer scalar
C     ISET(*)  - Integer total number of a subset of IRAY.
C                    Data type - integer array
C
C  OUTPUT
C     IND      - An integer index whereby either ISTR is identical
C                to IRAY(IND), or the slash-delimited substring of
C                ISTR is identical to an element NS of SETS and the
C                first part of the substring of ISTR is identical
C                to an IRAY(IND) associated with SETS(NS);
C                If none of these conditions are met, IND=0.
C                   Data type - integer scalar
C     NF       - Total number of times ISTR occurs in IRAY,
C                or total number of times ISTR occurs in a subset
C                of IRAY.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISET(*)
      CHARACTER*(*) ISTR, IRAY(*), SETS(*)
C
      IND = 0
      NF  = 0
C
      IF (ISTR .EQ. ' ') RETURN
C
      I2 = ILASCH(ISTR)
      IF (ISTR(I2:I2) .NE. '/') THEN
         CALL SKCOMP (ISTR, IRAY, NN, IND, NF)
         RETURN
      ENDIF
C
      I1 = 0
      DO 10 L = I2-1, 1, -1
         IF (ISTR(L:L).EQ.'/' .AND. I1.EQ.0) I1 = L
   10 CONTINUE
      IF (I1 .LE. 0) RETURN
C
      K = 0
      DO 50 N = 1, NSETS
         DO 50 J = 1, ISET(N)
            K = K + 1
            IF (ISTR(I1+1:I2-1) .EQ. SETS(N) .AND.
     1          ISTR(:I1-1) .EQ. IRAY(K) .AND. IND.EQ.0) THEN
                IND = K
                NF = NF + 1
            ENDIF
   50 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKSPEC (K, KNAME, KKTOT, ITHRM, MAXTP, T, NT,
     1                   TMID, LOUT, KERR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NT(*), T(MAXTP,*), IPLUS(10)
      CHARACTER*(*) KNAME(*)
      CHARACTER INUM(10)*1
      LOGICAL ITHRM(*), KERR
C
      DATA INUM/'0','1','2','3','4','5','6','7','8','9'/
C
C        each species must have thermodynamic data
C
      IF (.NOT.ITHRM(K)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1110) KNAME(K)
      ENDIF
C
C        thermodynamic temperatures must be in order
C
      IF (T(1,K) .LT. 0.0) THEN
         KERR = .TRUE.
         WRITE (LOUT, 575)
      ENDIF
      IF (T(NT(K),K) .LT. 0.0) THEN
         KERR = .TRUE.
         WRITE (LOUT, 580)
      ENDIF
      IF (T(2,K) .LT. 0.0) T(2,K) = TMID
      IF (T(1,K) .GE. T(NT(K),K)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 500)
      ENDIF
      IF (T(1,K) .GT. T(2,K)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 510)
      ENDIF
      IF (T(NT(K),K) .LT. T(2,K)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 520)
      ENDIF
C
C        a species cannot start with a number
C
      CALL SKCOMP (KNAME(K)(:1), INUM, 10, I, NF)
      IF (I .GT. 0) THEN
         KERR = .TRUE.
         WRITE (LOUT, 210)
      ENDIF
C
C        is there another species name after a +
C
      NPLUS = 0
      DO 30 N = 1, ILASCH(KNAME(K))
         IF (KNAME(K)(N:N) .EQ. '+') THEN
            NPLUS = NPLUS + 1
            IPLUS(NPLUS) = N
         ENDIF
   30 CONTINUE
      DO 40 N = 1, NPLUS
         I1 = IPLUS(N)
         IF (I1 .EQ. 1) THEN
            WRITE (LOUT, 220)
            KERR = .TRUE.
         ELSE
C
C        is there another species name after a +
C
            I1 = I1 + 1
            IF (N .LT. NPLUS) THEN
               DO 35 L = N+1, NPLUS
                  I2 = IPLUS(L)
                  IF (I2 .GT. I1) THEN
                     CALL SKCOMP (KNAME(K)(I1:I2-1),KNAME,
     1                            KKTOT, KNUM, NF)
                     IF (KNUM .GT. 0) THEN
                        WRITE (LOUT, 230)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
   35          CONTINUE
            ENDIF
            I2 = ILASCH(KNAME(K))
            IF (I2 .GE. I1) THEN
               CALL SKCOMP (KNAME(K)(I1:I2), KNAME, KKTOT, KNUM, NF)
               IF (KNUM .GT. 0) THEN
                  WRITE (LOUT, 230)
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
   40 CONTINUE
C
  210 FORMAT (6X,'Error...species starts with a number')
  220 FORMAT (6X,'Error...species starts with a plus')
  230 FORMAT (6X,'Error...illegal + in species name')
  500 FORMAT (6X,'Error...High temperature must be > Low temperature')
  510 FORMAT (6X,'Error...Low temperature must be <= Mid temperature')
  520 FORMAT (6X,'Error...High temperature must be => Mid temperature')
  575 FORMAT (6X,'Error...No Low thermodynamics temperature given')
  580 FORMAT (6X,'Error...No High thermodynamics temperature given')
 1110 FORMAT (6X,'Error...no thermodynamic properties for species ',A)
C
      RETURN
      END
