C@(#)===================================================================
C@(#)
C@(#)
C@(#)                   FILE =  sklib.f
C@(#)
C@(#)  ---------------  VERSION = 4.15
C@(#)  |  SCCS  FILE |
C@(#)  |   SUMMARY   |  CURRENT CHECKOUT DATE = 08/10/94
C@(#)  ---------------                           at 17:21:57
C@(#)                   DATE OF NEWEST DELTA = 08/10/94
C@(#)                                            at 17:21:54
C@(#)  SCCS file name = /users/chemkin/SCCS/s.sklib.f
C@(#)===================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C ---------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKABS
C     SKLIB is the library of surface kinetics subroutines for the
C     collection of codes known as CHEMKIN.  The acronym CHEMKIN is a
C     registered Trademark.
C///////////////////////////////////////////////////////////////////
C
C            SKLIB: SURFACE KINETICS SUBROUTINE LIBRARY
C                   VERSION 5.0
C
C     WRITTEN BY:
C         FRAN M. RUPLEY AND
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C       AND
C         MICHAEL E. COLTRIN
C         SANDIA NATIONAL LABORATORIES
C         SURFACE PROCESSING SCIENCES DIVISION
C         ALBUQUERQUE, NM 87185
C         (505) 844-7843
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
C     VERSION 4.5
C
C     CHANGES FROM VERSION 1.0:
C       1. Added subroutine SKRAEX
C     CHANGES FROM VERSION 1.1:
C       1. Changed from individual site densities to SUMDEN in
C          Subroutines SKRROP and SKRAT
C     CHANGES FROM VERSION 1.2:
C       1. Replaced "REAL*8" statements with "DOUBLE PRECISION"
C     CHANGES FROM VERSION 1.3:
C       1. Added subroutines SKDRDA, SKDRDC, and SKRDWC
C     CHANGES FROM VERSION 1.4:
C       1. Added thermodynamic properties subroutines SKTHM, SKHML,
C          SKHMS, SKHRT
C       2. Changed from SCOPY for integer arrays to DO loops
C       3. Changed order of call list for SKWT, SKKNAM
C     CHANGES FROM VERSION 1.5:
C       1. Correct call lists in some subroutines using NKK...'s as
C          an argument instead of defined in COMMON /SKSTRT/
C       2. Add double precision change blocks with DBLE(I) to
C          parallel REAL(I) used in calculations
C       3. Replace Subroutine SKRTK
C     CHANGES FROM VERSION 1.6:
C       1. Added thermodynamic properties subroutines -
C             SKTHM3  SKHML3  SKHMS3 SKHRT3
C             SKAML   SKAML3  SKAMS  SKAMS3
C             SKCPML  SKCPL3  SKCPMS SKCPS3  SKCPOR
C             SKGML   SKGML3  SKGMS  SKGMS3
C             SKSML   SKSML3  SKSMS  SKSMS3  SKSMH  SKSMH3 SKSOR
C             SKUML   SKUML3  SKUMS  SKUMS3
C     CHANGES FROM VERSION 1.7
C       1. Added PROLOGUE modifications
C       2. Added subroutine SKNU, SKSTOI and SKCOVI
C       3. Eliminate SKSEL
C       4. Switch old SKRAT  to SKRAT3 and add SKRAT,
C          switch old SKRATI to SKRTI3 and add SKRTI,
C          switch old SKRTI  to SKRATI
C       5. Added subroutines SKEQYP and SKYTCZ
C       6. Correct call list for CALL SKSTOI in subroutine SKNU
C       7. Added DO 20 loop to subroutine SKRROP
C       8. Pointer in 2nd call to DCOPY in SKYTCZ
C       9. Eliminated SKWT, since essentially a duplicate of SKKWT
C      10. Added missing DIMENSION ISKWRK(*) to subroutine SKSITP
C     CHANGES FROM VERSION 1.8
C       1. Changed DDOT/SDOT to DO loops
C       2. Allow reversible Arhennius coefficients, with changes
C          to /SKSTRT/ of NIIREV, the number of reactions with
C          reverse Arrhenius parameters defined, ISKWRK pointer IiIREV
C          for the array of reaction numbers, RSKWRK pointer IrRPAR
C          for the array of reverse parameters, and additional
C          reaction workspace pointer IrIT3
C       3. New SKRROP to allow reverse Arrhenius parameters.
C     CHANGES FROM VERSION 1.9
C       1. Character manipulation subroutines have additional
C          error handling and arguments
C       2. Change IckNAME to IckNAM in SKDNAM
C     CHANGES FROM VERSION 2.0
C       1. Get rid of non-standard comments
C     CHANGES FOR VERSION 2.5
C       1. Implement new mixture and pure solids.
C     CHANGES FOR VERSION 3.0
C       1. Implement phases, bulk species
C     CHANGES FOR VERSION 3.1
C       1. vax/cray change blocks for machine constants changed to
C          smallexp/bigexp blocks
C     CHANGES FOR VERSION 3.1
C       1. Mixture/pure species eliminated in favor of bulk species
C          only; binary file and common skstrt restructured
C     CHANGES FOR VERSION 3.2
C       1. Phase pointer for SDEN was incorrect in loop 210 of SKDEN
C          and loop 130 of SKATCZ
C     CHANGES FOR VERSION 3.3
C       1. CZ(K) for site species multiplied by their site coverage
C                parameter
C       2. Add sticking coefficient, and 'NONCON' options to allow
C          non-balancing site information for reactions;
C          RNCF(N,I), N=1,NPHASE is site information for the species
C          in reaction I
C     CHANGES FOR VERSION 3.4
C       1. Remove pointer and read statement for RSKWRK(IrSTK); real
C          sticking coefficients are now stored in the Arrhenius
C          parameter work space.
C       2. Site densities now have units of moles/cm**2 instead of
C          # densities.
C       3. SKRROP modified (DO 120 loop) for 'sticking' coefficients
C          of reaction I:
C          If equivalent first Arrhenius coefficient is
C             PAR1 = PAR(1,I) * 3637.601 / SQRT(WT1),
C JM insert:
C 3637.601 = SQRT(R/2*PI) i.e. this forms the 'preexponential factor '
C
C             where WT1 is the molecular wt. of the gas-phase reactant,
C          and equivalent second Arrhenius coefficient is
C             PAR2 = 0.5 + PAR(2,I),
C JM insert: this is the SQRT(T) term as PAR(2,I) is usually zero
C
C          then the forward rate expression
C             RKF(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
C          becomes
C             RKF(I) = PAR1 * EXP(PAR2*ALOGT - PAR(3,I)/T) / (SDTOT**M)
C JM insert: this is equation 46 in SURF. CHEMKIN manual
C
C          where SDTOT is the total of the site densities (moles/cm**2),
C          and M is the sum of the stoichiometric coefficients of the
C          site-phase (surface) reactants.
C          This affects all subroutines that CALL SKRROP, since SDTOT,
C          WT(*), NIISTK, IISTK(*), NKKGAS and NKKSUR must be known.
C      4.  Subroutine SKATCZ and all subroutines which call it now
C          require additional input argument SDEN, the site densities.
C      5.  Subroutine SKRATS is renamed SKRAT to replace previous SKRAT
C          and returns additional array SDOT.
C      6.  Correct SKATCZ: divide by coverage parameter rather than
C          multiply
C      CHANGES FOR VERSION 3.5
C      1.  Additional pointer IiKCOV and modified read for the species
C          numbers of the surface species associated with the coverage
C          parameters in a surface reaction; additional subroutine
C          SKICOV provides information for surface reaction I as to
C          it total number of coverage species, their species numbers,
C          and their coverage parameter numbers.
C      2.  Modified SKRAT to correct SDOT loop
C      3.  Modified SKDEN to correct units and calculation of site
C          species densities
C      4.  Modified SKRROP for coverage factors in calculations of
C          forward rates.
C      5.  Subroutine SKLEN reads first record of surface binary
C          file in order to provide lengths required for work arrays.
C      6.  Removed arrays KFIRST, KLAST, and KKPHAS from argument
C          list of Subroutine SKINDX.
C      7.  Change Subroutine SKPKK to return arrays KFIRST, KLAST,
C          and KKPHAS for all the phases, instead of just for
C          a specified phase.
C      8.  Additional input SDEN required for Subroutine SKDEN.
C      9.  Added integer constant NIICON to first record of binary
C          file, to COMMON/SKSTRT/, and to Subroutine ???,
C           for the total number of reactions which do not conserve
C           sites.
C     10.  No longer need argument ISKWRK for SKINDX.
C     11.  New Subroutine SKCONT
C     12.  First record of binary file contains character strings
C          VERS and PREC to indication binary file version and
C          precision, and the logical flag KERR to indicate an
C          error in the binary file; KERR previously was the
C          first item in the binary file
C     13.  Correct SKRAT DO 50, SKRATI DO 100, and SKRDWC DO 400 to
C          make sure NUNK(N,I) .ne. 0.
C     14.  Delete SKRPAR, since somewhat duplicates SKABE.
C     15.  Change SKKEL to SKNCF to parallel Chemkin.
C     16.  Change SKHRT to SKHORT     "        "    .
C     17.  CHANGE SKTHM TO SKATHM     "
C     18.  Delete SKRDWC and include it's calculation in SKDRDC,
C          since SKDRDC is the only subroutine that calls SKRDWC.
C     19.  COMMON /MACH/ renamed to COMMON /MACHN/ to avoid conflict
C          with CKLIB
C     20.  Additional factor in forward rate calculation
C     CHANGES TO VERSION 3.6
C      1.  Bring up to date with manual changes
C      2.  SKCOMP has additional argument NT to indicate how many
C          occurrences of a character string occurs in an array
C          of character strings.
C      3.  New Subroutine SKPCMP to find and species index number
C          for a SPECNAME/PHASENAME/ character string.
C      4.  Subroutine SKABE must convert sticking
C          coefficients to Arrhenius coefficients
C      5.  Subroutine SKSNUM
C      6.  Subroutine SKISTK
C      7.  User subroutine SKRTI renamed to SKRATI, and internal
C          subroutine SKRATI becomes SKRTI
C      8.  Eliminate SKRTI, since same logic is in SKRATI; the
C          only subroutine that called SKRTI is SKDRDA; SKDRDA
C          now calls SKRATI instead.
C      9.  Correct DO 20 loop in SKRAT to go from NFSUR to NLSUR
C          instead of 1 to NNSUR
C      CHANGES TO VERSION 3.61
C      1.  Correct SKABE conversion of sticking coefficients
C      2.  Change SKDEN to return densities of surface species
C          in gm/cm**2 instead of moles/cm**2
C      3.  Check that there site-phase species, bulk-phase species,
C          sites, or site phases before filling appropriate arrays.
C      4.  Correct SUMS in SKABE and SKRROP to include only
C          stoichiometric coefficients of Surface Species, and
C          check that SUMS is greater than 0.0
C      5.  SIGN function in SKRROP to preserve sign of rates
C      6.  Add VERS 3.61 to list of correct versions in SKLEN
C      7.  Correct DO 50 loop in SKICOV
C      8.  DO 115 in SKABE  only if NNSUR.GT.0
C          DO  60 in SKDRDC  "       "
C          DO  60 in SKEQ    "       "
C          DO  20 in SKRAT   "       "
C             (and combine two loops over NNSUR into one loop)
C          DO  60 in SKROP   "       "
C      9.  In SKSDEN, if NNSUR.LE.0 RETURN
C     10.  In SKDRDA, CALL SKABE instead of SKRAEX to get RA for
C          reaction IR.
C     11.  In SKDRDC, ROPNU is a sum if species present as reactant
C          and as product
C     12.  In SKRAT, replace the variable SDOT with SITDOT,
C          then replace WDOT with SDOT
C     13.  In SKRATI, replace the variable WDOTI with SDOTI, and
C          add the array SITDTI (phase change due to a reaction)
C     14.  Add pointer IrPT1 to COMMON/SKSTRT for real dummy work
C          space of length NPHASE, in order to accept the array
C          SITDTI(*) where SUBROUTINE SKDRDA has a CALL SKRATI
C     15.  In SKSYMR, include slash-delimited phase name in a
C          reaction string where a species name is not unique
C     16.  New subroutine SKPNT, SKSAVE to read, write binary
C          file information and work arrays.
C     17.  SKSYMR corrected to start the reaction string in the
C          first character space.
C     CHANGES FOR VERSION 3.62
C      1.  New versions of SKDRDA and SKDRDC by M. Coltrin.
C      2.  New SKDRDA reflects reversal of change #10 above.
C      3.  SKDRDA call list changed from (IR, ROP, ISKWRK, ...) to
C          (IR, P, T, ACT, SDEN, ISKWRK, ...)
C     CHANGES FOR VERSION 3.63
C      1.  Accept binary file versions 3.62 and 3.63
C      2.  Modify SKRROP according to new Equation 35, need
C          NFSUR, NLSUR, KFIRST, KLAST, SDEN, and RNCF in call list.
C     CHANGES FOR VERSION 3.64
C      1.  Accept Surface binary file 3.64
C     CHANGES FOR VERSION 3.7
C      1.  SKNCON requires ISKWRK input, and returns array of
C          length NPHASE
C      2.  SKNU returns array KSTOIC, stoichiometric coefficients
C          for the surface reactions, and array NSTOIC, the net change
C          in phases for the surface reactions.
C      3.  SKSDEN return array SDEN changed to SDEN0
C      4.  Correct indexing for ISKWRK(IiNCF) in SKINIT, SKRAT and
C          SKRATI
C      5.  Initialize NSIG=0 in Subroutine SKRROP after the start of
C          the DO 80 loop, not before it.
C      CHANGES FOR VERSION 3.71
C      1.  Accept Surface binary file 3.71.
C      CHANGES FOR VERSION 3.72
C      1.  In SKRROP, if SDTOT=0 then RKF(I)=0.
C      CHANGES FOR VERSION 3.73
C      1.  Accept Surface binary file 3.72 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.74
C      1.  Correction to SKRROP
C      CHANGES FOR VERSION 3.75
C      1.  Accept Surface binary file 3.73 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.76
C      1.  Accept Surface binary file 3.74 (correction made to
C          interpreter)
C      CHANGES FOR VERSION 3.77
C      1.  Accept Surface binary file 3.75 (correction)
C      CHANGES FOR VERSION 3.78 (F. Rupley, per M. Coltrin 1/3/91)
C      1.  Correct bug in Subroutine SKNU
C      CHANGES FOR VERSION 3.79 (F. Rupley, per M. Coltrin 1/10/91)
C      1.  Expand logic in SKCOMP to ignore leading/trailing spaces
C          in strings being compared.
C      CHANGES FOR VERSION 3.80 (F. Rupley, per M. Coltrin 1/11/91)
C      1.  Add COMMON MACHN and set values of BIG, SMALL, EXPARG in
C          Subroutine SKPNT.
C      CHANGES FOR VERSION 3.81 (F. Rupley)
C      1.  Accept surface binary file V.3.76 (correction)
C      CHANGES FOR VERSION 3.82 (F. Rupley)
C      1.  Accept surface binary file V.3.77 (correction)
C      CHANGES FOR VERSION 3.83 (2/18/91 F. Rupley, per R. Kee)
C      1.  Modify equation using sticking parameters for rates in
C          SKRROP, per Peter Glarborg.
C       2. Add a fourth parameter to the array of Arhennius
C          coefficients for the IISUR reactions;
C          increase the value of NSPAR in COMMON /SKSTRT/ by one.
C          (this also increases the length of the array of reverse
C          Arhennius parameters);
C          initialize the value of the fourth parameter to 1.0 in
C          SKINIT;
C          use this value as a "perturbation factor" for the forward
C          rates in SKRROP;
C          add SUBROUTINE SKRDEX to allow applications codes to change
C          the perturbation factor RD(I) in sensitivity calculations.
C       3. Accept binary file V.3.78 (LENRSK was increased by NIISUR
C          + NIIREV to reflect above changes in RSKWRK array.
C       CHANGES FOR VERSION 3.84 (4/2/91 F. Rupley)
C       1. Add Subroutine SKRHEX to get/put thermodynamic polynomial
C          coefficient a6, to enable applications codes to perturb
C          heat of formation for species.
C       2. Add Subroutine SKMAXTP to find the maximum number of
C          temperatures used in fitting thermodynamic data for the
C          species.
C       3. Change Subroutine SKABE to return either Arhennius or
C          sticking coefficients for the reactions, and an array
C          of integer flags to indicate the type of the coefficients.
C       CHANGES FOR VERSION 3.85 (4/29/91 F. Rupley per M. Coltrin)
C       1. Correct indexing in SKRDEX
C       CHANGES FOR VERSION 3.86 (5/9/91 F. Rupley)
C       1. In SKRROP, perform the "d" (PAR(4,I)) perturbation
C          before the checking for change of sign of the rates instead
C          of after.
C       CHANGES FOR VERSION 3.87 (5/23/91 H. Moffat)
C       1. In SKRROP, check to see if sticking coefficient is greater
C          than one. If it is, set it to one.
C       CHANGES FOR VERSION 4.0 (6/5/91 M. Coltrin)
C       1. Add new subroutines SKFLGS, SKIREV, SKNUF (H. Moffat)
C       2. Change units on surface-coverage modification of the the
C          rate of progress (concentration units --> site fractions)
C       3. Fix error in calculation of equilibrium constant. Added
C          array of length IISUR to hold correction factor for each
C          reaction. Added pointer IrEQ to common block SKSTRT.
C          Added EQFAC to arguments of SKRROP.
C       4. In SKATCZ and SKDEN we were dividing a real number by
C          ISKWRK(IiNSCV). Put in change blocks to first convert
C          the number to floating point (either REAL or DBLE).
C       CHANGES FOR VERSION 4.01 (7/11/91 F. Rupley)
C       1. Accept Interpreter binary file V.4.01.
C       CHANGES FOR VERSION 4.02 (7/17/91 F. Rupley)
C       1. Accept Interpreter binary file V.4.02.
C       CHANGES FOR VERSION 4.03 (8/1/91 F. Rupley per M. Coltrin)
C       1. New SKSTRT pointer IrEQ added to record read in SKPNT and
C          added to record written in SKSAVE.
C       CHANGES FOR VERSION 4.04 (8/28/91 F. Rupley per J. Grcar)
C       1. Instead of using NSPAR+1 in the calls to SKROP for the
C          dimension of SPAR, pass NSPAR, then dimension
C          SPAR(NSPAR+1,*) in SKROP, to avoid confusion as to the
C          value of NSPAR.
C       CHANGES FOR VERSION 4.06 (4/13/92 F. Rupley)
C       1. Accept Interpreter binary file V.4.04. (correction to
C          Subroutine SKDUP)
C       CHANGES FOR VERSION 4.07 (10/1/92 F. Rupley per M. Coltrin)
C       1. COMMON /SKCONS/ VERS, PREC, KERR, LENI, LENR, LENC
C          eliminates the need for LINKSK argument to SKSAVE
C       CHANGES FOR VERSION 4.08 (1/27/94 F. Rupley per R. Kee)
C       1. Real stoichometric coefficients added; NRNU,  IRNU,
C          RNU (additional pointers)
C       2. Integer phase (site) balance INCF becomes real RNCF;
C          require RCKWRK be added to calls to SKNCON, SKNU
C      CHANGES FOR VERSION 4.09 (3/15/94 F. Rupley)
C      1.  DOS/PC compatibility effort includes adding file names to
C          OPEN statements, removing unused variables in CALL lists,
C          unusued but possibly initialized variables.
C      CHANGES FOR VERSION 4.10 (4/14/94 F. Rupley, per E. Meeks)
C      1.  use INCLUDE 'skstrt.inc' instead of having COMMON /SKSTRT/ in
C          every subroutine
C      CHANGES FOR VERSION 4.11 (4/19/94 F. Rupley)
C      1.  accept linking file 4.07 (correction in SKBAL, SKRBAL)
C      CHANGES FOR VERSION 4.12 (4/28/94 F. Rupley, per M. Coltrin)
C      1.  New Subroutines SKINU, SKIRNU, SKIORD for real
C          stoichiometric coefficients and change of order.
C      CHANGES FOR VERSION 4.12b (5/20/94 F. Rupley per E. Meeks)
C      1.  Incorporate plasma options (linking file 4.08b)
C      CHANGES FOR VERSION 4.12c (6/3/94 F. Rupley)
C      1.  Accept linking file 4.08c (bugfixes per H. Moffat)
C      CHANGES FOR VERSION 4.13 (6/28/94 E. Meeks)
C      1.  Add subroutine SKIBHM and correct Bohm reaction rate formula
C      2.  Correct common block /SKSTRT/ to include IiIBHM
C      3.  Added option for correction to sticking parameter; Option
C          is turned off in the interpreter input using REACTION line
C          keyword MWOFF and passed through linking file to SKRROP.
C      CHANGES FOR VERSION 4.14 (7/14/94 F.Rupley)
C      1.  Changes to allow multiple materials include:
C          pointers are now stored in ISKWRK
C          most subroutines now require ISKWRK input
C          backspaces instead of rewind in SKLEN
C          etc.
C      CHANGES FOR VERSION 4.2 (7/27/94 F. Rupley per E. Meeks)
C      1.  Correction NIBOHM > NIIBHM; accept linking file 4.2.
C      2.  Correct DO 70 in SKRROP for BOHM rates.
C      CHANGES FOR VERSION 4.3 (8/2/94 F. Rupley)
C      1.  Bugfixes regarding using old constants vs. newer
C          ISKWRK pointers
C      2.  Add pointers to SKSTRT in order to store work array lengths
C      3.  Accept linking file V.4.3.
C      CHANGES FOR VERSION 4.4 (8/4/94 F. Rupley)
C      1.  Correct READ statement in SKINIT for 'IF (NIICOV'.
C      2.  Correct 'IND =' index in SKABE
C      3.  Use REWIND in SKINIT if last material in linking file
C      CHANGES FOR VERSION 4.5 (8/6/94 F. Rupley)
C      1.  Assign local variables before CALL SKRROP in order to
C          shorten call list.
C      CHANGES FOR VERSION 4.51 (8/10/94 H. Moffat)
C      1.  Changed gas constant to conform to 1986 CODATA
C          recommendations. RUC, obtained from RU by dividing
C          by 4.184, is now compatible with RU up to DP machine
C          precision.
C      CHANGES FOR VERSION 4.52 (8/26 F. Rupley)
C      1.  Correct above value for RUC (RSKWRK(ISKWRK(IrRUC)))
C          in SKINIT.
C      CHANGES FOR VERSION 4.7 (9/7/94 F. Rupley per M. Coltrin)
C      1.  Loop 10 in SKEQ contained a wrong pointer (IrEQ should have
C          been IRIT3).  This subroutine has not returned the
C          correct value since version 4.2.
C      CHANGES FOR VERSION 4.8 (9/30/94 F. Rupley per G. Evans)
C      1.  Correct SKRROP (DO 70); variable SDOT should be SDTOT
C      CHANGES FOR VERSION 4.9 (12/15/94 F. Rupley per E. Meeks)
C      1.  Correct value of TEMP in PSKHML.
C      CHANGES FOR VERSION 5.0 (1/19/95 F. Rupley per M. Coltrin)
C      1.  Add integer error flag to SKLEN, SKINIT call lists)
C/////////////////////////////////////////////////////////////////////
C
C  START PROLOGUE
C
C  SUBROUTINE SKABS
C
C  The work arrays contain all the pertinent information about the
C  species and the reaction mechanism.  They also contain some work
C  space needed by various routines for internal manipulations.
C  If a user wishes to modify an SKLIB subroutine or to write new
C  routines, he will probably want to use the work arrays directly.
C  The starting addresses for information stored in the work arrays
C  are found in the labeled common block, COMMON /SKSTRT/, and are
C  explained below.
C
C   COMMON /SKSTRT/ MAXSPR, NELEM,   NKKGAS, NKKSUR, NKKBLK, NKKTOT,
C  1                NPHASE, NFSUR,  NLSUR,  NNSUR,  NFBLK,  NLBLK,
C  2                NNBLK,  NIISUR, NSPAR,  NSCOV,  NIICOV, NIIREV,
C  3                NIISTK, NIICON, MAXTP,   NCP,    NCP1,   NCP2,
C  4                NCP2T,  IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG,
C  5                IiKCMP, IiNSCV, IiKNT,  IiNRPP, IiNREA, IiNUNK,
C  6                IiNU,   IiNCF,  IiNSUM, IiICOV,  IiKCOV,  IiIREV,
C  7                IiISTK,  IrPATM, IrRU,   IrRUC,  IrSDEN, IrKTMP,
C  8                IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR,  IrCOV,
C  9                IrRPAR,  IrEQ,   IrKT1,  IrKT2,  IrPT1,  IrIT1,
C  *                IrIT2,  IrIT3,  IcENAM, IcKNAM, IcPNAM
C
C  INDEX CONSTANTS.
C
C     MAXSPR - Maximum number of species in any surface reaction.
C              Unless changed in the interpreter MAXSPR=12.
C     NELEM   - Number of elements.
C     NKKGAS - Number of gas-phase species.
C     NKKSUR - Number of surface species.
C     NKKBLK - Number of bulk species.
C     NKKTOT - Total number of species.  NKKTOT=NKKGAS+NKKSUR+NKKBLK
C     NPHASE - Number of phases; gas + sites + bulk.
C     NFSUR  - Pointer to the first surface phase.
C     NLSUR  - Pointer to the last surface phase.
C     NNSUR  - Number of surface phases.
C     NFBLK  - Pointer to the first bulk phase.
C     NLBLK  - Pointer to the last bulk phase.
C     NNBLK  - Number of bulk phases.
C     NIISUR - Number of surface reactions.
C     NSPAR  - Number of parameters in the rate expression.
C              In the current formulation NSPAR=3.
C     NSCOV  - Number of parameters in a coverage expression.
C     NIICOV - Number of surface reactions with coverage parameters.
C     NIIREV - Number of surface reactions with reverse parameters.
C     NIISTK - Number of surface reactions with sticking coefficients.
C     NIICON - Number of surface reactions with non-conserved sites.
C     MAXTP   - Maximum number of temperatures used to fit thermodynamic
C              properties of species;
C              unless changed in the interpreter MAXTP=3.
C     NCP    - Number of polynomial coefficients to fits of CP/R.
C              Unless the interpreter and the thermodynamic data base
C              are changed NCP=5.
C     NCP1   - NCP + 1.
C     NCP2   - NCP + 2.
C     NCP2T  - (MAXTP-1) * NCP2.  Total number of thermodynamic fit
C              coefficients for (MAXTP-1) temperature ranges.
C              Unless changed NCP2T=14.
C
C  STARTING ADDRESSES FOR THE INTEGER WORK SPACE, ISKWRK.
C
C     IiPKST - Starting address of the starting species numbers of
C              the NPHASE phases.  ISKWRK(IiPKST+N-1) is the first
C              species in the Nth phase.
C     IiPKND - Starting address of the final species numbers of
C              the NPHASE phases.  ISKWRK(IiPKND+N-1) is the last
C              species in the Nth phase.
C     IiPTOT - Starting address of the total number of species in
C              the NPHASE phases.  ISKWRK(IiPTOT+N-1) is the total
C              number of species in the Nth phase.
C     IiKPHS - Starting address of the phases of the NKKTOT species.
C              ISKWRK(IiKPHS+K-1) = -1, the Kth species is a solid
C                                 =  0, the Kth species is a gas
C                                 = +1, the Kth species is a liquid
C     IiKCHG - Starting address of the electronic charges of the
C              NKKTOT species.
C              ISKWRK(IiKCHG+K-1) = -2, the Kth species has two excess
C                                       electrons.
C     IiKCMP - Starting address of the elemental content of the NELEM
C              elements in the NKKTOT species.
C              ISKWRK(IiKCMP + (K-1)*NELEM + M - 1) is the number of
C              atoms of the Mth element in the Kth species.
C     IiNSCV - Starting address of site coverage for the NKKTOT species.
C              ISKWRK(IiNSCV+K-1) is site coverage for the Kth species.
C     IiKTMP - Starting address of the number of temperatures used
C              to fit thermodynamic coefficients for the NKKTOT species.
C              ISKWRK(IiKTMP+K-1) = N, N temperatures were used in the
C                                      fit for the Kth species.
C     IiNRPP - Starting address of the total number of participant
C              species for the NIISUR reactions, and the reversibility
C              of the reactions.
C              ISKWRK(IiNRPP+I-1) = +N, the Ith reaction is reversible
C                                       and has N participant species
C                                       (reactants + products)
C                                 = -N, the Ith reaction is irreversible
C                                       and has N participant species
C                                       (reactants + products)
C     IiNREA - Starting address of the number of reactants only for
C              the NIISUR reactions.  ISKWRK(IiNREA+I-1) is the total
C              number of reactants (not including products) in the Ith
C              reaction.
C     IiNUNK - Starting address of a matrix of species index numbers for
C              the MAXSPR species in the NIISUR reactions.
C              ISKWRK(IiNUNK+(I-1)*MAXSPR+N-1) = K, the species number
C              of the Nth participant species in the Ith reaction.
C     IiNU   - Starting address of a matrix of stoichiometric
C              coefficients of the MAXSPR species in the NIISUR
C              reactions.  ISKWRK(IiNU+(I-1)*MAXSPR+N-1) is the
C              coefficient of the Nth participant species in the
C              Ith reaction.
C     IiNCF  - Starting address of the net change in sites for the
C              NPHASE phases.
C     IiNSUM - Starting address of the sum of the stoichiometric
C              coefficients of the gas-phase reactants and products
C              in the NIISUR surface reactions.  ISKWRK(IiNSUM+I-1)
C              is the sum of the stoichiometric coefficients of the
C              gas-phase species in surface reaction I.
C     IiICOV  - Starting address of the NIICOV reaction numbers for
C              surface reactions with coverage parameters.
C              ISKWRK(IiICOV+N-1) is the surface reaction index I
C              of the Nth coverage declaration.
C     IiKCOV  - Starting address of the NIICOV species index numbers
C              for surface reactions with coverage parameters.
C              ISKWRK(IiKCOV+N-1) is the species index K of the Nth
C              coverage declaration; the reaction number is
C              ISKWRK(IiICOV+N-1).
C     IiIREV  - Starting address of the NIIREV reaction numbers for
C              surface reactions with reverse parameters.
C              ISKWRK(IiIREV+N-1) is the surface reaction index I
C              of the Nth reaction with reverse parameters.
C     IiISTK  - Starting address of the NIISTK reaction numbers for
C              surface reactions with sticking coefficients.
C              ISKWRK(IiISTK+N-1) is the surface reaction index I
C              of the Nth reaction with sticking coefficients.
C
C  STARTING ADDRESSES FOR THE REAL WORK SPACE, RSKWRK.
C
C     IrPATM - RSKWRK(IrPATM) is the pressure of one standard atmosphere
C              (dynes/cm**2).
C     IrRU   - RSKWRK(IrRU) is the universal gas constant (ergs/mole-K).
C     IrRUC  - RSKWRK(IrRUC) is the universal gas constant (cal/mole-K).
C     IrSDEN - Starting address of densities for the NPHASE phases.
C              RSKWRK(IrSDEN+N-1) is the density of the Nth phase.
C     IrKTMP - Starting address of a matrix of the MAXTP temperatures
C              used in fits of thermodynamic properties for the NKKTOT
C              species.  RSKWRK(IrKTMP+(K-1)*MAXTP+N-1) is the Nth
C              temperature for the Kth species.
C     IrKTHM - Starting address of a three-dimensional array of
C              coefficients for the NCP2 fits to the thermodynamic
C              properties for the NKKTOT species, for (MAXTP-1)
C              temperature ranges.
C              RSKWRK(IrKTHM+(L-1)*NCP2+(K-1)*NCP2T+N-1) = A(N,L,K;
C              A(N,L,K),N=1,NCP2T = polynomial coefficients in the fits
C              for the Kth species and the Lth temperature range, where
C              the total number of temperature ranges for the
C              Kth species is ISKWRK(IiKNT+K-1) - 1.
C     IrKDEN - Starting address of the densities (gm/cm**3 for gas or
C              bulk species, gm/cm**2 for surface species) of the
C              NKKTOT species.  RSKWRK(IrKDEN+K-1) is the density of
C              the Kth species.
C     IrAWT  - Starting address of the atomic weights for the NELEM
C              elements.  RSKWRK(IrAWT+M-1) is the atomic weight
C              of the Mth element.
C     IrKWT  - Starting address of the molecular weights for the
C              NKKTOT species.  RSKWRK(IrKWT+K-1) is the molecular
C              weight of the Kth species.
C     IrPAR  - Starting address of a matrix of the NSPAR Arrhenius
C              parameters for the NIISUR surface reactions.
C              RSKWRK(IrPAR+(I-1)*NSPAR+N-1) is the Nth parameter
C              of the Ith surface reaction, where
C                 N=1 is the pre-exponential factor (mole-cm-sec-K),
C                 N=2 is the temperature exponent, and
C                 N=3 is the activation energy (kelvins).
C     IrCOV  - Starting address of a matrix of the NSCOV coverage
C              parameters for the NIICOV coverage declarations.
C              RSKWRK(IrCOV+(N-1)*NSCOV+L-1) is the Lth coverage
C              parameter for the Nth coverage declaration;
C              the reaction number is ISKWRK(IiICOV+N-1);
C              the species number is ISKWRK(IiKCOV+N-1).
C     IrRPAR  - Starting address of a matrix of the NSPAR reverse
C              parameters for the NIIREV surface reactions with
C              reverse parameters.  RSKWRK(IrRPAR+(N-1)*NSPAR+L-1)
C              is the Lth reverse parameter for the Nth reverse
C              reaction; the reaction number is ISKWRK(IiIREV+N-1).
C     IrEQ   - Starting address of multiplicative constant for the
C              calculation of NIISUR surface equilibrium constants.
C              RSKWRK(IrEQ+I-1) is the constant for the Ith reaction.
C     IrKT1  - Starting addresses of arrays of internal work space
C     IrKT2  -                of length NKKTOT.
C     IrIT1  - Starting addresses of arrays of internal work space
C     IrIT2  -                of length NIISUR.
C     IrIT3  -
C
C  STARTING ADDRESSES FOR THE CHARACTER WORK SPACE, CSKWRK.
C
C     IcENAM - Starting address of the NELEM element names.
C              CSKWRK(IcENAM+M-1) is the name of the Mth element.
C     IcKNAM - Starting address of the NKKTOT species names.
C              CSKWRK(IcKNAM+K-1) is the name of the Kth species.
C     IcPNAM - Starting address of the NPHASE phase names.
C              CSKWRK(IcPNAM+N-1) is the name of the Nth phase.
C
C The binary file consists of the following binary records:
C
C 1) Information about the binary file:  VERS, PREC, KERR
C    Where VERS   = character*16 string representing the version
C                   number of the interpreter which created the
C                   the binary file.
C          PREC   = character*16 string representing the machine
C                   precision of the binary file (SINGLE, DOUBLE).
C          KERR   = logical which indicates whether or not
C                   an error occurred in the interpreter input.
C 2) Index constants:
C    LENISK, LENRSK, LENCSK, MAXSPR, MAXTP, NCP, NELEM,
C    NKKGAS, NKKSUR, NKKBLK, NKKTOT, NPHASE, NFSUR, NLSUR,
C    NNSUR, NFBLK, NLBLK, NNBLK, NIISUR, NSPAR, NSCOV, NIICOV,
C    NIIREV, NIISTK, NIICOM
C    Where LENISK = required length of ISKWRK.
C          LENRSK = required length of RSKWRK.
C          LENCSK = required length of CSKWRK.
C
C 3) Element information:
C    (CSKWRK(IcENAM+M-1),                         !element names
C     RSKWRK(IrAWT+M-1),                          !atomic weights
C     M=1,NELEM)
C
C 4) Species information:
C    (CSKWRK(IcKNAM+K-1),                         !species names
C     RSKWRK(IrKWT+K-1),                          !molecular weight
C     ISKWRK(IiKPHS+K-1),                         !phase
C     ISKWRK(IiKCHG+K-1),                         !charge
C     ISKWRK(IiKNT+K-1),                          !# of fit temps
C     ISKWRK(IiNSCV+K-1),                         !site coverage
C     RSKWRK(IrKDEN+K-1),                         !species densities
C     (ISKWRK(IiKCMP+(K-1)*NELEM+M-1),M=1,NELEM),   !composition
C     (RSKWRK(IrKTMP+(K-1)*MAXTP+L-1), L=1,MAXTP),  !array of temps
C     ((RSKWRK(IrKTHM+(L-1)*NCP2+(K-1)*NCP2T+N-1),!fit coeff'nts
C              N=1,NCP2), L=1,(MAXTP-1)),
C     K=1,NKKTOT),
C
C 5) Phase Information:
C    (CSKWRK(IcPNAM+N-1),                         !phase names
C     ISKWRK(IiPKST+N-1),                         !start species
C     ISKWRK(IiPKND+N-1),                         !end species
C     ISKWRK(IiPTOT+N-1),                         !# of species
C     RSKWRK(IrSDEN+N-1),                         !phase densities
C     N=1,NPHASE),
C
C 6) Reaction information (if NIISUR > 0)
C    (ISKWRK(IiNRPP+I-1),                       !# of species
C     ISKWRK(IiNREA+I-1),                       !# of reactants
C     ISKWRK(IiNSUM+I-1),                       !sum gas-phase coeff.
C     (RSKWRK(IrPAR+(I-1)*NSPAR+N-1),N=1,NSPAR),!Arr. coefficients
C     (ISKWRK(IiNUNK+(I-1)*MAXSPR+N-1),         !stoic coef
C     ISKWRK(IiNU+(I-1)*MAXSPR+N-1),N=1,MAXSPR),!species numbers
C     I=1,NIISUR),
C    (ISKWRK(IiICOV+I-1), ISKWRK(IiKCOV+I-1),   !cov. reactions,species
C    (RSKWRK(IrCOV+(I-1)*NSCOV+L-1),L=1,NSCOV), !cov. parameters
C     I=1,NIICOV),
C    (ISKWRK(IiIREV+I-1),                          !rev. reactions
C    (RSKWRK(IrRPAR+(I-1)*NSPAR+L-1),L=1,NSPAR),   !rev. parameters
C     I=1,NIIREV),
C    (ISKWRK(IiISTK+I-1),I=1,NIISTK)               !sticking reactions
C
C  END PROLOGUE
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKABE  (ISKWRK, RSKWRK, RA, RB, RE, ISTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKABE  (ISKWRK, RSKWRK, RA, RB, RE, ISTFL)
C     Returns the Arrhenius coefficients or the sticking coefficients
C     of the surface reactions, and integer flags to indicate the type
C     of the coefficients.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real work space.
C                   Data type - real array
C
C  OUTPUT
C     RA     - Pre-exponential constants for the reactions.
C                   cgs units - mole-cm-sec-K
C                   Data type - real array
C                   Dimension RA(*) at least IISUR, the total
C                   number of surface reactions.
C     RB     - Temperature dependence exponents for the reactions.
C                   cgs units - none
C                   Data type - real array
C                   Dimension RB(*) at least IISUR, the total
C                   number of surface reactions.
C     RE     - Activation energies for the reactions.
C                   cgs units - Kelvins
C                   Data type - real array
C                   Dimension RE(*) at least IISUR, the total
C                   number of surface reactions.
C     ISTFL  - 0 if a reaction does not use sticking coefficients.
C              1 if a reaction does use sticking coefficients.
C                   Data type - integer array
C                   Dimension ISTFL(*) at least IISUR, the
C                   total number of surface reactions.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RA(*), RB(*), RE(*), ISKWRK(*), RSKWRK(*), ISTFL(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      DO 100 I = 1, ISKWRK(IiNIIS)
         IND   = ISKWRK(IrPAR) + (I-1)*(NSPAR+1)
         RA(I) = RSKWRK(IND)
         RB(I) = RSKWRK(IND+1)
         RE(I) = RSKWRK(IND+2)
C         WRITE(*,*) RE(I)
         ISTFL(I) = 0
  100 CONTINUE
C
      IF (ISKWRK(IiNSTK) .LE. 0) RETURN
C
      DO 150 N = 1, ISKWRK(IiNSTK)
         I     = ISKWRK(ISKWRK(IiISTK) + N - 1)
         ISTFL(I) = 1
  150 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKAML  (T, ISKWRK, RSKWRK, AML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKAML  (T, ISKWRK, RSKWRK, AML)
C     Returns an array of the standard state Helmholtz free energies
C     in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     AML    - Standard state Helmholtz free energies in molar units
C              for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C                   Dimension AML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), AML(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKSML (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT1)))
      CALL SKHML (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT2)))
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         AML(K) = RSKWRK(ISKWRK(IrKT2) + K - 1) - T *
     1           (RSKWRK(ISKWRK(IrRU)) + RSKWRK(ISKWRK(IrKT1) + K - 1))
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKAMS  (T, ISKWRK, RSKWRK, AMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKAMS  (T, ISKWRK, RSKWRK, AMS)
C     Returns an array of the standard state Helmholtz free energies
C     in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     AMS    - Standard state Helmholtz free energies in mass units
C              for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension AMS(*) at least KKTOT, the total number
C                   of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), AMS(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKSMS (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT1)))
      CALL SKHMS (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT2)))
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         AMS(K) = RSKWRK(ISKWRK(IrKT2) + K - 1) - T *
     1           (RSKWRK(ISKWRK(IrRU))/
     2            RSKWRK(ISKWRK(IrKWT) + K - 1) +
     3            RSKWRK(ISKWRK(IrKT1) + K - 1))
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C
C  START PROLOGUE
C
C  SUBROUTINE SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C     Returns the concentrations of the species, given the pressure,
C     temperature and activities.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     CZ     - Matrix of the concentrations of the gas-phase and
C              surface species in the problem, and the activities
C              of the bulk species.  The first KKGAS entries of CZ
C              are the gas-phase molar concentrations (moles/cm**3).
C              The next KKSURF entries are the surface species molar
C              concentrations (moles/cm**2).  The final KKBULK entries
C              are the activities of the bulk species.
C                   Data type - real array
C                   Dimension CZ(*) at least KKTOT, the total
C                   number of gas-phase + surface + bulk species.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), CZ(*)
      INCLUDE 'skstrt.inc'
C
C     COMPUTE THE MOLAR CONCENTRATIONS (C) FROM MOLE FRACTIONS
C     the multiplication by (t/298.15) is done by ashish to make sure
C     that gas phase temp (298.15K) is used and not Tsurface.
      PRUT = P/(RSKWRK(ISKWRK(IrRU))*T)
      DO 100 K = 1, NKKGAS
c         CZ(K) = ACT(K)*PRUT*T/298.15
         CZ(K) = ACT(K)*PRUT
  100 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
         KCOV = ISKWRK(IiNSCV)
         DO 210 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(ISKWRK(IiPKST) + N - 1)
            KLAST  = ISKWRK(ISKWRK(IiPKND) + N - 1)
            DO 210 K = KFIRST, KLAST
C*****precision > double
               CZ(K) = ACT(K) * SDEN(N) / DBLE(ISKWRK(KCOV + K - 1))
C*****END precision > double
C*****precision > single
C               CZ(K) = ACT(K) * SDEN(N) / REAL(ISKWRK(KCOV + K - 1))
C*****END precision > single
  210    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .GT. 0) THEN
         KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
         KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
         DO 220 K = KFIRST, KLAST
            CZ(K) = ACT(K)
  220    CONTINUE
      ENDIF
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKATHM (MDIM, NDIM1, NDIM2, ISKWRK, RSKWRK, NT, TMP,
     1                   A)
C
C  START PROLOGUE
C
C  SUBROUTINE SKATHM (MDIM, NDIM1, NDIM2, ISKWRK, RSKWRK, NT, TMP,
C                     A)
C     Returns the polynomial coefficients of the fits for
C     thermodynamic properties of all of the species.
C
C  INPUT
C
C     MDIM   - First dimension of an array of temperatures used
C              in the thermodynamic fits for the species; MDIM
C              must be at least MAXTP, the maximum number of
C              temperatures used to fit the thermodynamics.
C                 Data type - integer scalar
C     NDIM1  - First dimension of an array of thermodynamic fit
C              coefficients; NDIM1 must be at least NCP2, the total
C              number of coefficients for one temperature range.
C                   Data type - integer scalar
C     NDIM2  - Second dimension of the array of thermodynamic fit
C              coefficients; NDIM2 must be at least (MAXTP-1), the
C              number of temperature ranges.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C
C     Where NT(K) is the number of temperatures used in fitting the
C     thermodynamic properties of species K, TMP(N) is the Nth
C     temperature, NT(K)-1 is the number of temperature ranges for
C     which the polynomial coefficients are valid, then
C     A (L, N, K) is the Lth polynomial coefficient, for the Nth
C     temperature range, and the Kth species; i.e.,
C
C          | <   N = 1  >. <N=2> .                    .<  N = NT - 1>
C    P  E  |    .        .       .                    .        .
C    O  X  |    .        .       .                    .        .
C    L  P  |    .        .       .                    .        .
C    Y  R  |    .        .       .                    .        .
C    N  E  |    .        .       .                    .        .
C    O  S  |    .        .       .                    .        .
C    M  S  |    .        .       .                    .        .
C    I  I  |    .        .       .                    .        .
C    A  O  |    .        .       .                    .        .
C    L  N  |____.________._______.____________________.________.______
C             TMP(1)  TMP(2)  TMP(3) . .  .  .  . TMP(NT-1)  TMP(NT)
C
C
C
C     NT     - Number of temperatures used for fitting the
C              thermodynamic properties of the species.
C                   Data type - integer array
C                   Dimension NT(*) at least KKTOT, the total number of
C                   species.
C     TMP    - The temperatures which divide the temperature ranges
C              over which the polynomial coefficients are valid.
C                   cgs units - K
C                   Data type - real array
C                   Dimension TMP(*,*) exactly MAXTP (the maximum
C                   number of temperatures allowed) for the first
C                   dimension and at least KKTOT (the total number of
C                   species) for the second.
C     A      - Three-dimensional array of fit coefficients to the
C              thermodynamic data for the species.
C              The indicies in  A(N,L,K) mean-
C              N = 1,NN represent polynomial coefficients in CP/R
C                CP/R(K)=A(1,L,K) + A(2,L,K)*T + A(3,L,K)*T**2 + ...
C              N = NN+1 is for the formation enthalpies, i.e.,
C                HO/R = A(NN+1,L,K)
C              N = NN+2 is for the formation entropies, i.e.,
C                SO/R = A(NN+2,L,K)
C              L = 1 is for temperature <= TMP(2,K)
C              L = 2 is for TMP(2,K) < temperature <= TMP(3)
C                :
C              L = (NTMP-1) is for TMP(NTMP-1) <= temperature;
C              K  is  the  species index
C                   Data type - real array
C                   Dimension A(*,*,*) exactly NPCP2 (the total number
C                   of coefficients for each temperature range) for
C                   the first dimension, MAXTP-1 (the maximum number of
C                   temperature ranges) for the second dimension, and at
C                   least KKTOT (the total number of species) for the
C                   third.
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
      DIMENSION NT(*),TMP(MDIM,*),A(NDIM1,NDIM2,*),ISKWRK(*),
     1          RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         NT(K) = ISKWRK(ISKWRK(IiKNT) + K - 1)
  100 CONTINUE
C
      IF (MDIM .LT. MAXTP) RETURN
      DO 140 L = 1, MAXTP
         DO 140 K = 1, ISKWRK(IiKTOT)
            TMP(L,K) = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + L - 1)
  140 CONTINUE
C
      DO 150 K = 1, ISKWRK(IiKTOT)
         DO 150 L = 1, MAXTP-1
            NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
            DO 150 M = 1, NCP2
               A(M, L, K) = RSKWRK(NA1 + M - 1)
150   CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCHRG (ISKWRK, RSKWRK, KCHARG)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCHRG (ISKWRK, RSKWRK, KCHARG)
C     Returns an array containing electronic charges of the species.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     KCHARG - Electronic charges of the species.
C              KCHARG(K)=-2 indicates that the Kth species has two
C              excess electrons.
C                   Data type - integer array
C                   Dimension KCHARG(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), KCHARG(*)
      INCLUDE 'skstrt.inc'
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         KCHARG(K) = ISKWRK(ISKWRK(IiKCHG) + K - 1)
  100 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCOMP (ISTR, IRAY, NN, IND, NT)
C     Given a character string vector IRAY of vector length NN, this
C     subroutine will search for the first occurrence of the string
C     ISTR.  The position of ISTR in IRAY is returned in the integer
C     argument NT.  If the string ISTR is not found in IRAY, IND and
C     NT are assigned a zero value.
C
C     Consider the following example,
C        IRAY = {"BOOK","BLUE","BEAR","BOOK"}
C        NN=4.
C
C     If ISTR="BLUE" then IND=2 and NT=1;
C     if ISTR="RED"  then IND=0 and NT=0; and
C     if ISTR="BOOK",then IND=1 and NT=2.
C
C  INPUT
C     ISTR   - A reference character string.
C     IRAY   - Array of character strings.
C                   Data type - character array
C                   Dimension IRAY(*) at least NN.
C     NN     - Length of IRAY(*).
C                   Data type - integer scalar
C
C  OUTPUT
C     IND    - Index of the position in IRAY(*) containing ISTR.
C              If ISTR is not in IRAY(*), IND=0.
C                   Data type - integer scalar
C     NT     - Total number of times ISTR occurs in IRAY.
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
      NT  = 0
      DO 10 I = NN, 1, -1
         IS1 = IFIRCH(ISTR)
         IS2 = ILASCH(ISTR)
         IR1 = IFIRCH(IRAY(I))
         IR2 = ILASCH(IRAY(I))
         IF ( IS2.GE.IS1 .AND. IS2.GT.0 .AND.
     1        IR2.GE.IR1 .AND. IR2.GT.0 .AND.
     2        ISTR(IS1:IS2).EQ.IRAY(I)(IR1:IR2) ) THEN
            IND = I
            NT  = NT + 1
         ENDIF
   10 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCONT (KSPEC, ROP, ISKWRK, RSKWRK, CIK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCONT (KSPEC, ROP, ISKWRK, RSKWRK, CIK)
C     Returns the contributions of each of the surface reactions to
C     the molar production rate of species KSPEC.
C
C  INPUT
C     KSPEC  - Integer species number.
C                   Data type - integer scalar
C     ROP    - Rates of progress for the reactions.
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension ROP(*) at least IISUR, the total number
C                   of SURFACE reactions.
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C
C  OUTPUT
C     CIK    - Contributions of the surface reactions to the molar
C              production rate of species KSPEC.
C                   cgs units - mole/(cm**2*sec)
C                   Data type - real array
C                   Dimension CIK(*) at least IISUR, the total number
C                   of surface reactions.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), ROP(*), CIK(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      DO 100 I = 1, ISKWRK(IiNIIS)
         CIK(I) = 0.0
  100 CONTINUE
C
      DO 200 N = 1, MAXSPR
         DO 200 I = 1, ISKWRK(IiNIIS)
            NK = ISKWRK(ISKWRK(IiNUNK) + MAXSPR*(I-1) + N - 1)
            IF (NK .EQ. KSPEC) THEN
C*****precision > double
               RNC = DBLE (ISKWRK(ISKWRK(IiNU) + MAXSPR*(I-1) + N - 1))
C*****END precision > double
C*****precision > single
C               RNC = REAL (ISKWRK(ISKWRK(IiNU) + MAXSPR*(I-1) + N - 1))
C*****END precision > single
               CIK(I) = CIK(I) + RNC*ROP(I)
            ENDIF
200   CONTINUE
C
      IF (ISKWRK(IiNRNU) .LE. 0) RETURN
      DO 300 N = 1, ISKWRK(IiNRNU)
         I = ISKWRK(ISKWRK(IiIRNU) + N - 1)
         DO 300 L = 1, MAXSPR
            NK = ISKWRK(ISKWRK(IiNUNK) + MAXSPR*(I-1) + L - 1)
            IF (NK .EQ. KSPEC) THEN
               RNC = RSKWRK(ISKWRK(IrRNU) + MAXSPR*(N-1) + L - 1)
               CIK(I) = CIK(I) + RNC*ROP(I)
            ENDIF
  300 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCOV  (ISKWRK, KOCC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCOV  (ISKWRK, KOCC)
C     Returns an array of site occupancy numbers for the species.
C
C  INPUT
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C
C     KOCC   - Site occupancy numbers for the species.
C                   Data type - integer array
C                   Dimension KOCC(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), KOCC(*)
      INCLUDE 'skstrt.inc'
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         KOCC(K) = ISKWRK(ISKWRK(IiNSCV) + K - 1)
   50 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPML (T, ISKWRK, RSKWRK, CPML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPML (T, ISKWRK, RSKWRK, CPML)
C     Returns an array of the specific heats at constant pressure
C     in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     CPML   - Specific heats at constant pressure in molar units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CPML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), CPML(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         CPML(K) = RSKWRK(ISKWRK(IrRU))*SUM
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPMS (T, ISKWRK, RSKWRK, CPMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPMS (T, ISKWRK, RSKWRK, CPMS)
C     Returns an array of the specific heats at constant pressure
C     in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C   OUTPUT
C     CPMS   - Specific heats at constant pressure in mass units
C              for the species.
C                   cgs units - ergs/(gm*K)
C                   Data type - real array
C                   Dimension CPMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), CPMS(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         CPMS(K) = RSKWRK(ISKWRK(IrRU))*SUM /
     1             RSKWRK(ISKWRK(IrKWT) + K - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKCPOR (T, ISKWRK, RSKWRK, CPOR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKCPOR (T, ISKWRK, RSKWRK, CPOR)
C     Returns an array of the nondimensional specific heats at constant
C     pressure.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C   OUTPUT
C     CPOR   - Nondimensional specific heats at constant pressure
C              for the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension CPOR(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), TN(10), CPOR(*)
      INCLUDE 'skstrt.inc'
C
      TN(1) = 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         CPOR(K) = 0.0
         DO 250 N = 1, NCP
            CPOR(K) = CPOR(K) + TN(N)*RSKWRK(NA1 + N - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDEN  (P, T, ACT, SDEN, ISKWRK, RSKWRK, DEN)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDEN  (P, T, ACT, SDEN, ISKWRK, RSKWRK, DEN)
C     Returns a real array of species densities.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     DEN    - Densities for the species.
C              NOTE:  mass densities are not required to be input to
C                     the Interpreter for bulk-phase species.
C                     If they are input, they are returned by this
C                     subroutine.  If not, DEN = -1.0 for the bulk
C                     species
C                   cgs units - gm/(cm**3) for gas-phase phase species
C                               gm/(cm**2) for surface species
C                               gm/(cm**3) for bulk species
C
C                   Data type - real array
C                   Dimension DEN(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), DEN(*)
      INCLUDE 'skstrt.inc'
C     the mulitiplication by t/298.15 is done by ashish to make sure
c     that gas phase temp (298.15K) is used and not the Tsurface.
      PRUT = P / (T * RSKWRK(ISKWRK(IrRU)))
      DO 120 K = 1, NKKGAS
c         DEN(K) = ACT(K)*RSKWRK(ISKWRK(IrKWT)+K-1)*PRUT*T/298.15
         DEN(K) = ACT(K)*RSKWRK(ISKWRK(IrKWT)+K-1)*PRUT
  120 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
         DO 130 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(ISKWRK(IiPKST) + N - 1)
            KLAST  = ISKWRK(ISKWRK(IiPKND) + N - 1)
            DO 130 K = KFIRST, KLAST
C*****precision > double
               RCOV = DBLE (ISKWRK(ISKWRK(IiNSCV) + K - 1))
C*****END precision > double
C*****precision > single
C               RCOV = REAL (ISKWRK(ISKWRK(IiNSCV) + K - 1))
C*****END precision > single
               DEN(K) =  ACT(K)*SDEN(N)*RSKWRK(ISKWRK(IrKWT) + K - 1)
     1                   /RCOV
  130    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .GT. 0) THEN
         KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
         KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
         DO 140 K = KFIRST, KLAST
            DEN(K) = RSKWRK(ISKWRK(IrKDEN) + K - 1)
  140    CONTINUE
      ENDIF
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDRDA (IR, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDAI,
     1                   Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDRDA (IR, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDAI)
C     Returns the partial of the rates of production for each of the
C     species with respect to the pre-exponential constant of surface
C     reaction IR.
C
C  INPUT
C     IR     - Reaction index
C                   Data type - integer scalar
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C  OUTPUT
C     DKDAI  - Array of the partial of the production rates of the
C              species with respect to the pre-exponential
C              constant for reaction IR.
C                   cgs units - moles/(cm**2*sec) / (units of A)
C                   Data type - real array
C                   Dimension DKDAI(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), DKDAI(*)
      LOGICAL IFLAG
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDAI(K) = 0.0
   50 CONTINUE
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta)
C
C     PROCESS REACTION IR
      CALL SKRAEX  (IR, ISKWRK, RSKWRK, RA)
C
C     FIND OUT IF REVERSE ARRHENIUS PARAMETERS WERE SPECFIED FOR
C     THIS REACTION
      IFLAG = .TRUE.
      DO 70 N = 1, NIIREV
         IF (ISKWRK(ISKWRK(IiIREV)+N-1). EQ. IR) IFLAG = .FALSE.
   70 CONTINUE
C
      IRNU = 0
      DO 75 N = 1, NIIRNU
         IF (ISKWRK(ISKWRK(IiIRNU)+N-1) .EQ. IR)
     1      IRNU = ISKWRK(IrRNU)+MAXSPR*(N-1)
   75 CONTINUE
C
      IND = (IR-1) *MAXSPR
      RKF = RSKWRK(ISKWRK(IrIT1) + IR - 1)
      RKR = RSKWRK(ISKWRK(IrIT2) + IR - 1)
      DO 80 N = 1, MAXSPR
         NK = ISKWRK(ISKWRK(IiNUNK) + IND + N - 1)
         IF (NK .NE. 0) THEN
C*****precision > double
            RNU = DBLE (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > double
C*****precision > single
C           RNU = REAL (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > single
            IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
            IF (IFLAG) THEN
C              REVERSE PARAMETERS WERE NOT GIVEN FOR THIS REACTION
               DKDAI(NK) = DKDAI(NK) - RNU*RKR / RA
            ENDIF
            DKDAI(NK) = DKDAI(NK) + RNU*RKF / RA
         ENDIF
   80 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKDRDC (KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC,
     1                   Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKDRDC (KSPEC, P, T, ACT, SDEN, ISKWRK, RSKWRK, DKDC)
C     Returns the partial derivative of the production rates for
C     each of the species with respect to the concentration of species
C     KSPEC.
C
C  INPUT
C     KSPEC  - Species index
C                   Data type - integer scalar
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     DKDC   - Array of the partial of the production rates of the
C              species with respect to the concentration
C              of species KSPEC.
C                   cgs units - moles/(cm**2*sec) / (units of KSPEC)
C                   Data type - real array
C                   Dimension DKDC(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ACT(*), ISKWRK(*), RSKWRK(*), SDEN(*), DKDC(*)
      LOGICAL IFLAG
      INCLUDE 'skstrt.inc'
C
      DATA TEN/10.0/
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDC(K) = 0.0
   50 CONTINUE
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      CZK = RSKWRK(ISKWRK(IrKT1) + KSPEC - 1)
C
C     PROCESS EACH REACTION
C
      DO 1000 I = 1, NIISUR
C
C        FIND OUT IF REVERSE ARRHENIUS PARAMETERS WERE SPECFIED FOR
C        THIS REACTION
         IFLAG = .TRUE.
         DO 70 N = 1, NIIREV
            IF ( ISKWRK(ISKWRK(IiIREV)+N-1). EQ. I) IFLAG = .FALSE.
   70    CONTINUE
         IRNU = 0
         DO 75 N = 1, NIIRNU
            IF (ISKWRK(ISKWRK(IiIRNU) + N - 1) .EQ. I)
     1         IRNU = ISKWRK(IrRNU) + MAXSPR*(N-1)
   75    CONTINUE
C
         FFACT = 0.0
         RFACT = 0.0
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR
         INU = ISKWRK(IiNU)   + (I-1)*MAXSPR
         DO 350 N = 1, MAXSPR
            NK = ISKWRK(INK + N - 1)
            IF (NK .EQ. KSPEC) THEN
C*****precision > double
               RNU = DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C               RNU = REAL (ISKWRK(INU + N - 1))
C*****END precision > single
               IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
               IF (RNU .LT. 0) THEN
                   FFACT = -RNU / CZK
               ELSE
                  RFACT =  RNU / CZK
               ENDIF
            ENDIF
  350    CONTINUE
C
         DO 200 N = 1, NIICOV
            IF ( ISKWRK(ISKWRK(IiICOV)+N-1). EQ. I) THEN
C              COVERAGE PARAMETERS WERE SPECIFIED FOR THIS REACTION
               IF ( ISKWRK(ISKWRK(IiKCOV)+N-1). EQ. KSPEC) THEN
C                 SPECIES NUMBER KSPEC MODIFIES THE RATE EXPRESSION
C                 FOR REACTION I.
                  IC = ISKWRK(IrKCOV) + (N-1)*NSCOV - 1
                  CPAR1 = RSKWRK(IC + 1)
                  CPAR2 = RSKWRK(IC + 2)
                  CPAR3 = RSKWRK(IC + 3)
C                 ADD THE COV. MODIFICATION FOR FORWARD DIRECTION
                  CMOD = CPAR2/ACT(KSPEC) + LOG(TEN)*CPAR1 - CPAR3/T
C                 SCALE BY THE NUMBER OF SITES COVERED BY THIS SPECIES
C                 DIVIDED BY SITE DENSITY OF THIS SURFACE PHASE.
C                 THIS NUMBER TURNS OUT TO BE THE SITE FRACTION DIVIDED
C                 BY THE MOLAR CONCENTRATION FOR THE SURFACE SPECIES.
                  CMOD = CMOD * ACT(KSPEC) / CZK
                  FFACT = FFACT + CMOD
                  IF (IFLAG) THEN
C                    REVERSE ARRHENIUS PARAMETERS WERE NOT GIVEN FOR
C                    THIS REACTION, SO IT IS OK TO ADD-IN THE
C                    COVERAGE MODIFICATION TO THE REVERSE DIRECTION
                     RFACT = RFACT + CMOD
                  ENDIF
               ENDIF
            ENDIF
  200    CONTINUE
C
         RKF = RSKWRK(ISKWRK(IrIT1) + I - 1)
         RKR = RSKWRK(ISKWRK(IrIT2) + I - 1)
         DO 400 N = 1, MAXSPR
            NK = ISKWRK(INK + N - 1)
            IF (NK .NE. 0) THEN
C*****precision > double
               RNU = DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C               RNU = REAL (ISKWRK(INU + N - 1))
C*****END precision > single
               IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
               DKDC(NK) = DKDC(NK) + RNU*(RKF*FFACT - RKR*RFACT)
            ENDIF
  400    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKEQ   (P, T, ACT, SDEN, ISKWRK, RSKWRK, EQKC,
     1                   Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKEQ    (P, T, ACT, SDEN, ISKWRK, RSKWRK, EQKC)
C     Returns the equilibrium constants for the surface reactions
C     given pressure, temperature, species activities, and the site
C     densities.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     EQKC   - Equilibrium constants in concentration units
C              for the reactions.
C                   cgs units - (moles,cm), depends on reaction
C                   Data type - real array
C                   Dimension EQKC(*) at least IISUR, the total
C                   number of surface reactions.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), ACT(*), SDEN(*), EQKC(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      DO 10 I = 1, NIISUR
         EQKC(I) = RSKWRK(ISKWRK(IrIT3) + I - 1)
   10 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKFLGS (ISKWRK, NRPP, IREV, ISTFL, ICOV)
C
C  START PROLOGUE
C
C  SUBROUTINE SKFLGS (ISKWRK, NRPP, IREV, ISTFL, ICOV)
C     Returns several integer arrays describing the type of
C     surface reactions in the mechanism
C
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     NRPP   - Reversibility flags and number of species
C              (reactants plus products) for reactions.
C              NRPP(I)=1, reaction I is reversible
C              NRPP(I)=0, reaction I is irreversible
C                   Data type - integer array
C                   Dimension NRPP(*) at least IISUR, the total number
C                   of surface reactions.
C     IREV   - Flag to indicate that reaction has reverse Arrhenius
C              parameters defined
C              IREV(I)=1, Reaction I has reverse Arrhenius parameters
C                         defined
C              IREV(I)=0, Reaction I does not have reverse Arrhenius
C                         parameters defined.  It may or may not
C                         be a reversible reaction (check NRPP for this
C                         information.
C                   Data type - integer array
C                   Dimension IREV(*) at least IISUR, the total number
C                   of surface reactions.
C     ISTFL  - 0 if reaction does not use sticking coefficients
C              1 if reaction does use sticking coefficients
C                   Data type - integer array
C                   Dimension ISTFL(*) at least IISUR, the total number
C                   of reactions.
C     ICOV   - 0 if reaction does not have coverage dependent parameters
C              1 if reaction does have coverage dependent parameters
C                   Data type - integer array
C                   Dimension ICOV(*) at least IISUR, the total number
C                   of reactions.
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
      DIMENSION ISKWRK(*), NRPP(*), IREV(*), ISTFL(*), ICOV(*)
      INCLUDE 'skstrt.inc'
C
      DO 5 I = 1, ISKWRK(IiNIIS)
         ICOV(I)  = 0
         ISTFL(I) = 0
         IREV(I)  = 0
5     CONTINUE
C
      DO 10 N = 1, ISKWRK(IiNSTK)
         ISTFL( ISKWRK(ISKWRK(IiISTK) + N - 1) ) = 1
   10 CONTINUE
C
      DO 20 I = 1, ISKWRK(IiNIIS)
         IF (ISKWRK(ISKWRK(IiNRPP) + I - 1) .GT. 0) THEN
           NRPP(I) = 1
         ELSE
           NRPP(I) = 0
         ENDIF
   20 CONTINUE
C
      DO 30 N = 1, ISKWRK(IiNCOV)
         ICOV( ISKWRK(ISKWRK(IiICOV) + N - 1) ) = 1
   30 CONTINUE
C
      DO 40 N = 1, ISKWRK(IiNREV)
         IREV( ISKWRK(ISKWRK(IiIREV) + N - 1) ) = 1
   40 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKGML  (T, ISKWRK, RSKWRK, GML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKGML  (T, ISKWRK, RSKWRK, GML)
C     Returns an array of the standard state Gibbs free energies
C     in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     GML    - Standard state Gibbs free energies in molar units
C              for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C                   Dimension GML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), GML(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKHML (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT1)))
      CALL SKSML (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT2)))
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         GML(K) = RSKWRK(ISKWRK(IrKT1) + K - 1) - T*
     1            RSKWRK(ISKWRK(IrKT2) + K - 1)
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKGMS  (T, ISKWRK, RSKWRK, GMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKGMS  (T, ISKWRK, RSKWRK, GMS)
C     Returns an array of the standard state Gibbs free energies
C     in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     GMS    - Standard state Gibbs free energies in mass units
C              for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension GMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), GMS(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKHMS (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT1)))
      CALL SKSMS (T, ISKWRK, RSKWRK, RSKWRK(ISKWRK(IrKT2)))
C
      DO 100 K = 1, ISKWRK(IiKTOT)
         GMS(K) = RSKWRK(ISKWRK(IrKT1) + K - 1) - T*
     1            RSKWRK(ISKWRK(IrKT2) + K - 1)
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHML  (T, ISKWRK, RSKWRK, HML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHML  (T, ISKWRK, RSKWRK, HML)
C     Returns an array of the enthalpies in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     HML    - Enthalpies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C                   Dimension HML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), HML(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      RUT = T*RSKWRK(ISKWRK(IrRU))
      TN(1) = 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)/N
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         HML(K) = RUT*(SUM + RSKWRK(NA1 + NCP1 - 1)/T)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHMS  (T, ISKWRK, RSKWRK, HMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHMS  (T, ISKWRK, RSKWRK, HMS)
C     Returns an array of the enthalpies in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C  OUTPUT
C     HMS    - Enthalpies in mass units for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension HMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), HMS(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      RUT = T*RSKWRK(ISKWRK(IrRU))
      TN(1) = 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)/N
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         HMS(K) = RUT * (SUM + RSKWRK(NA1 + NCP1 - 1)/T)
     1                / RSKWRK(ISKWRK(IrKWT) + K - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKHORT (T, ISKWRK, RSKWRK, HORT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKHORT (T, ISKWRK, RSKWRK, HORT)
C     Returns an array of the nondimensional enthalpies.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C  OUTPUT
C     HORT   - Nondimensional enthalpies for the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension HORT(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), HORT(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      DO 100 N = 1, NCP
         TN(N) = T**(N-1)/N
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         HORT(K) = SUM + RSKWRK(NA1 + NCP1 - 1)/T
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C
      SUBROUTINE SKIBHM (IR, ISKWRK, IBMFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIBHM (IR, ISKWRK, IBMFL)
C     Returns an integer flag to indicate whether reaction IR uses
C     sticking coefficients.
C
C  INPUT
C     IR     - Surface reaction index number.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     IBMFL  - 0 if reaction IR does not use BOHM sticking coefficients
C              1 if reaction IR does use BOHM sticking coefficients
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
      DIMENSION ISKWRK(*)
      INCLUDE 'skstrt.inc'
C
      IBMFL = 0
      DO 10 N = 1, ISKWRK(IiNBHM)
         IF (IR .EQ. ISKWRK(ISKWRK(IiIBHM) + N - 1) ) IBMFL = 1
   10 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKICOV (IR, NDIM, ISKWRK, RSKWRK, NCOVI, KCOVI, CPARI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKICOV (IR, NDIM, ISKWRK, RSKWRK, NCOVI, KCOVI, CPARI)
C     Returns the coverage species index numbers and their coverage
C     parameters for reaction IR.
C
C  INPUT
C     IR     - Surface reaction index number.
C                   Data type - integer scalar
C     NDIM   - Actual first dimension of CPAR.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     NCOVI  - Total number of species that modify the rate of
C              reaction IR through coverage dependence.
C                   Data type - integer scalar
C     KCOVI  - Index numbers for the NCOVI species that modify the
C              rate of reaction IR through coverage dependence.
C                   Data type - integer array
C                   Dimension KCOVI(*) at least KKTOT, the
C                   total number of species.
C     CPARI - Coverage parameters for the coverage species
C              of reaction IR.
C                   Data type - real array
C                   Dimension CPARI(*,*) exactly NSCOV for the
C                   first dimension, the total number of coverage
C                   parameters required for each coverage species,
C                   and at least KKTOT for the second dimension,
C                   the total number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), KCOVI(*), CPARI(NDIM,*)
      INCLUDE 'skstrt.inc'
C
      NCOVI = 0
      DO 10 K = 1, ISKWRK(IiKTOT)
         KCOVI(K) = 0
         DO 10 N = 1, NSCOV
            CPARI(N, K) = 0.0
   10 CONTINUE
C
      DO 50 N = 1, ISKWRK(IiNCOV)
         IF (IR .EQ. ISKWRK(ISKWRK(IiICOV) + N - 1)) THEN
            NCOVI = NCOVI + 1
            KCOVI(NCOVI) = ISKWRK(ISKWRK(IiKCOV) + N - 1)
            DO 25 J = 1, NSCOV
               CPARI(J, NCOVI) =
     1            RSKWRK(ISKWRK(IrKCOV) + NSCOV*(N-1) + J - 1)
   25       CONTINUE
         ENDIF
   50 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINDX (ISKWRK, NELM, KKGAS, KKSUR, KKBULK, KKTOT,
     1                   NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                   NLBULK, IISUR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINDX (ISKWRK, NELM, KKGAS, KKSUR, KKBULK, KKTOT,
C                     NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
C                     NLBULK, IISUR)
C     Returns a group of indices defining the size of the surface
C     reaction mechanism.
C
C  INPUT
C     ISKWRK -
C
C  OUTPUT
C     NELM   - Number of elements.
C                   Data type - integer scalar
C     KKGAS  - Number of gas-phase species.
C                   Data type - integer scalar
C     KKSUR  - Number of surface species.
C                   Data type - integer scalar
C     KKBULK - Total number of bulk species.
C                   Data type - integer scalar
C     KKTOT  - Total number of species.  KKTOT=KKGAS+KKSUR+KKBULK
C                   Data type - integer scalar
C     NNPHAS - Number of phases; gas + sites + bulk.
C                   Data type - integer scalar
C     NNSURF - Number of surface phases.
C                   Data type - integer scalar
C     NFSURF - Pointer to the first surface phase.
C                   Data type - integer scalar
C     NLSURF - Pointer to the last surface phase.
C                   Data type - integer scalar
C     NNBULK - Number of bulk phases.
C                   Data type - integer scalar
C     NFBULK - Pointer to the first bulk phase.
C                   Data type - integer scalar
C     NLBULK - Pointer to the last bulk phase.
C                   Data type - integer scalar
C     IISUR  - Number of surface reactions.
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
      DIMENSION ISKWRK(*)
      INCLUDE 'skstrt.inc'
C
      NELM  = NELEM
      KKGAS = NKKGAS
      KKSUR = ISKWRK(IiKSUR)
      KKBULK= ISKWRK(IiKBLK)
      KKTOT = ISKWRK(IiKTOT)
      NNPHAS= ISKWRK(IiNPHA)
      NNSURF= ISKWRK(IiNSUR)
      NFSURF= ISKWRK(IiFSUR)
      NLSURF= ISKWRK(IiLSUR)
      NNBULK= ISKWRK(IiNBLK)
      NFBULK= ISKWRK(IiFBLK)
      NLBULK= ISKWRK(IiLBLK)
      IISUR = ISKWRK(IiNIIS)
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINIT (LSIWK, LSRWK, LSCWK, LINKSK, LOUT,
     1                   ISKWRK, RSKWRK, CSKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINIT (LSIWK, LSRWK, LSCWK, LINKSK, LOUT,
C                     ISKWRK, RSKWRK, CSKWRK, IFLAG)
C     Reads the surface binary file and creates the internal work
C     arrays ISKWRK, RSKWRK, and CSKWRK.  SKINIT must be called before
C     any other Surface Chemkin subroutine is called.  The work arrays
C     must then be made available as input to the other Surface Chemkin
C     subroutines.
C
C  INPUT
C     LSIWK  - Dimension of ISKWRK
C                   Data type - integer scalar
C     LSRWK  - Dimension of RSKWRK
C                   Data type - integer scalar
C     LSCWK  - Dimension of CSKWRK
C                   Data type - integer scalar
C     LINKSK  - Unit number assigned to binary file
C                   Data type - integer scalar
C     LOUT   - Unit number assigned for output
C                   Data type - integer scalar
C
C  OUTPUT
C     ISKWRK - Array of integer workspace containing integer data.
C                   Data type - integer array
C     RSKWRK - Array of real workspace containing real data.
C                   Data type - real array
C     CSKWRK - Array of character workspace containing character data.
C                   Data type - CHARACTER*16 array
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      CHARACTER*16 CSKWRK(*), VERS, PREC
C*****precision > double
      DOUBLE PRECISION RU,               PA
      PARAMETER       (RU = 8.314510D7,  PA = 1.01325D6)
C*****END precision > double
C*****precision > single
C      REAL             RU,               PA
C      PARAMETER       (RU = 8.314510E7,  PA = 1.01325E6)
C*****END precision > single
      LOGICAL IOK,ROK,COK,KERR
      COMMON /SKCONS/ VERS, PREC, KERR
      COMMON /MACHN/ SMALL, BIG, EXPARG
      INCLUDE 'skstrt.inc'
C
C----------------------------------------------------------------------C
C     Data about the machine dependent constants is carried in         C
C     COMMON/MACHN/SMALL,BIG,EXPARG
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      EXPARG = LOG(BIG)
C
C----------------------------------------------------------------------C
C
C          WRITE VERSION NUMBER
C
c     WRITE (LOUT,15)
   15 FORMAT(/1X,' SKLIB:  Surface kinetics library',
     *       /1X,'         Copyright 1990, Sandia Corporation.',
     1       /1X,'         The U.S. Government retains a limited ',
     2                     'license in this software.',
     3       /1X,'         CHEMKIN-II Version 5.0, January 1995',
C*****precision > double
     4       /1X,'         DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     4       /1X,'         SINGLE PRECISION')
C*****END precision > single
C
c     ashish insert
      open(linksk, file='INP.d/sklink', form='unformatted')
      CALL SKLEN (LINKSK, LOUT, LENI, LENR, LENC)
c      IF (IFLAG .GT. 0) RETURN
C
      IOK = (LSIWK .GE. LENI)
      ROK = (LSRWK .GE. LENR)
      COK = (LSCWK .GE. LENC)
      IF (.NOT.IOK) WRITE (LOUT,300) LENI
      IF (.NOT.ROK) WRITE (LOUT,350) LENR
      IF (.NOT.COK) WRITE (LOUT,375) LENC
      IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
c         IFLAG = 20
         RETURN
      ENDIF
C
      MORE = 0
      READ (LINKSK,err=999) VERS, PREC, KERR, MORE
      READ (LINKSK,err=999) LENISK, LENRSK, LENCSK, MAXSPR, MAXTP, NCP,
     1                     NELEM, NKKGAS, NKKSUR, NKKBLK, NKKTOT,
     2                     NPHASE, NFSUR, NLSUR, NNSUR, NFBLK, NLBLK,
     3                     NNBLK, NIISUR, NSPAR, NSCOV, NIICOV, NIIREV,
     4                     NIISTK, NIICON, NIIBHM, NIIRNU, NIIORD,
     5                     MAXORD, SKMIN, MOTZWS
C
      NCP1  = NCP+1
      NCP2  = NCP+2
      NCP2T = NCP2*(MAXTP-1)
C
C     ISKWRK pointers to integer variables, array locations
C
      IiLENI = 1
      IiLENR = IiLENI + 1
      IiLENC = IiLENR + 1
      IiKSUR = IiLENC + 1
      IiKBLK = IiKSUR + 1
      IiKTOT = IiKBLK + 1
      IiNPHA = IiKTOT + 1
      IiFSUR = IiNPHA + 1
      IiLSUR = IiFSUR + 1
      IiNSUR = IiLSUR + 1
      IiFBLK = IiNSUR + 1
      IiLBLK = IiFBLK + 1
      IiNBLK = IiLBLK + 1
      IiNIIS = IiNBLK + 1
      IiNCOV = IiNIIS + 1
      IiNREV = IiNCOV + 1
      IiNSTK = IiNREV + 1
      IiNCON = IiNSTK + 1
      IiNBHM = IiNCON + 1
      IiNRNU = IiNBHM + 1
      IiNORD = IiNRNU + 1
      IiMOTZ = IiNORD + 1
      IiPKST = IiMOTZ + 1
      IiPKND = IiPKST + 1
      IiPTOT = IiPKND + 1
      IiKPHS = IiPTOT + 1
      IiKCHG = IiKPHS + 1
      IiKCMP = IiKCHG + 1
      IiNSCV = IiKCMP + 1
      IiKNT  = IiNSCV + 1
      IiNRPP = IiKNT  + 1
      IiNREA = IiNRPP + 1
      IiNUNK = IiNREA + 1
      IiNU   = IiNUNK + 1
      IiNSUM = IiNU   + 1
      IiICOV = IiNSUM + 1
      IiKCOV = IiICOV + 1
      IiIREV = IiKCOV + 1
      IiISTK = IiIREV + 1
      IiIBHM = IiISTK + 1
      IiKBHM = IiIBHM + 1
      IiTBHM = IiKBHM + 1
      IiIRNU = IiTBHM + 1
      IiIORD = IiIRNU + 1
      IiKORD = IiIORD + 1
C
C     ISKWRK pointers for real variables, arrays
C
      IrSKMN = IiKORD + 1
      IrPATM = IrSKMN + 1
      IrRU   = IrPATM + 1
      IrRUC  = IrRU   + 1
      IrSDEN = IrRUC  + 1
      IrKTMP = IrSDEN + 1
      IrKTHM = IrKTMP + 1
      IrKDEN = IrKTHM + 1
      IrAWT  = IrKDEN + 1
      IrKWT  = IrAWT  + 1
      IrPAR  = IrKWT  + 1
      IrKCOV = IrPAR  + 1
      IrRPAR = IrKCOV + 1
      IrEQ   = IrRPAR + 1
      IrRNU  = IrEQ   + 1
      IrNCF  = IrRNU  + 1
      IrKORD = IrNCF  + 1
      IrKFT  = IrKORD + 1
      IrKRT  = IrKFT  + 1
      IrKT1  = IrKRT  + 1
      IrKT2  = IrKT1  + 1
      IrPT1  = IrKT2  + 1
      IrIT1  = IrPT1  + 1
      IrIT2  = IrIT1  + 1
      IrIT3  = IrIT2  + 1
C
C     ISKWRK pointers to character variables, arrays
C
      IcENAM = IrIT3  + 1
      IcKNAM = IcENAM + 1
      IcMNAM = IcKNAM + 1
      IcPNAM = IcMNAM + 1
C
      NPNTS = IcPNAM
C
      ISKWRK(IiLENI) = LENI
      ISKWRK(IiLENR) = LENR
      ISKWRK(IiLENC) = LENC
C! number of site species
      ISKWRK(IiKSUR) = NKKSUR
C! number of bulk species
      ISKWRK(IiKBLK) = NKKBLK
C! total number of species
      ISKWRK(IiKTOT) = NKKTOT
C! number of phases
      ISKWRK(IiNPHA) = NPHASE
C! first site phase
      ISKWRK(IiFSUR) = NFSUR
C! last site phase
      ISKWRK(IiLSUR) = NLSUR
C! number of site phases
      ISKWRK(IiNSUR) = NNSUR
C! first bulk phase
      ISKWRK(IiFBLK) = NFBLK
C! last bulk phase
      ISKWRK(IiLBLK) = NLBLK
C! number of bulk phases
      ISKWRK(IiNBLK) = NNBLK
C! number of surface reactions
      ISKWRK(IiNIIS) = NIISUR
C! number of reactions with coverage parameters
      ISKWRK(IiNCOV) = NIICOV
C! number of reactions with reverse parameters
      ISKWRK(IiNREV) = NIIREV
C! number of reactions with sticking coefficients
      ISKWRK(IiNSTK) = NIISTK
C! number of reactions with non-conservation of sites
      ISKWRK(IiNCON) = NIICON
C! number of reactions with Bohm parameters
      ISKWRK(IiNBHM) = NIIBHM
C! number of reactions with real stoichiometry
      ISKWRK(IiNRNU) = NIIRNU
C! number of reactions with change-of-order
      ISKWRK(IiNORD) = NIIORD
C! MOTZWS flag (0/1)
      ISKWRK(IiMOTZ) = MOTZWS
C! starting species location of phases
      ISKWRK(IiPKST) = NPNTS + 1
C! ending species location of phases
      ISKWRK(IiPKND) = ISKWRK(IiPKST) + NPHASE
C! number of species each phase
      ISKWRK(IiPTOT) = ISKWRK(IiPKND) + NPHASE
C! species phases
      ISKWRK(IiKPHS) = ISKWRK(IiPTOT) + NPHASE
C! species charges
      ISKWRK(IiKCHG) = ISKWRK(IiKPHS) + NKKTOT
C! species elemental composition
      ISKWRK(IiKCMP) = ISKWRK(IiKCHG) + NKKTOT
C! site coverage of species
      ISKWRK(IiNSCV) = ISKWRK(IiKCMP) + NKKTOT*NELEM
C! # of temperatures for fit
      ISKWRK(IiKNT)  = ISKWRK(IiNSCV) + NKKTOT
C! num of species in reactions
      ISKWRK(IiNRPP) = ISKWRK(IiKNT) + NKKTOT
C! num of reactants in reactions
      ISKWRK(IiNREA) = ISKWRK(IiNRPP) + NIISUR
C! spec. numbers in reactions
      ISKWRK(IiNUNK) = ISKWRK(IiNREA) + NIISUR
C! species coefficients
      ISKWRK(IiNU)   = ISKWRK(IiNUNK) + NIISUR*MAXSPR
C! phase (site) balances for reactions
C      IiNCF   = IiNU   + NIISUR*MAXSPR
C! total gas-phase coefficients
C      IiNSUM  = IiNCF  + NIISUR*NPHASE
      ISKWRK(IiNSUM) = ISKWRK(IiNU) + NIISUR*MAXSPR
C! reaction num. for coverage parameters
      ISKWRK(IiICOV)  = ISKWRK(IiNSUM) + NIISUR
C! species num. for coverage parameters
      ISKWRK(IiKCOV)  = ISKWRK(IiICOV)  + NIICOV
C!  " with reverse parameters
      ISKWRK(IiIREV)  = ISKWRK(IiKCOV)  + NIICOV
C!  " with sticking parameters
      ISKWRK(IiISTK)  = ISKWRK(IiIREV)  + NIIREV
C!  reactions with Bohm velocity conditions
      ISKWRK(IiIBHM)  = ISKWRK(IiISTK)  + NIISTK
C! species numbers used in Bohm velocity formulation
      ISKWRK(IiKBHM)   = ISKWRK(IiIBHM)  + NIIBHM
C! species temperature flags for Bohm velocity
      ISKWRK(IiTBHM)   = ISKWRK(IiKBHM)   + NIIBHM
C!  " with real stoichometric coeff's
      ISKWRK(IiIRNU)  = ISKWRK(IiTBHM)  + NIIBHM
C!  " with change of order
      ISKWRK(IiIORD)  = ISKWRK(IiIRNU) + NIIRNU
C!  change of order species
      ISKWRK(IiKORD)  = ISKWRK(IiIORD) + NIIORD
C! total size req. for ISKWRK
      IiTOT   = ISKWRK(IiKORD) + NIIORD*MAXORD
      ISKWRK(IiTOT) = MORE
C
C! minimum difference allowed for balancing stoichiometry,
C  sites, etc.
      ISKWRK(IrSKMN) = 1
      RSKWRK(ISKWRK(IrSKMN)) = SKMIN
C! pressure 1 atm
      ISKWRK(IrPATM) = ISKWRK(IrSKMN) + 1
      RSKWRK(ISKWRK(IrPATM)) = PA
C
      ISKWRK(IrRU)   = ISKWRK(IrPATM) + 1
      RSKWRK(ISKWRK(IrRU)) = RU
C
      ISKWRK(IrRUC)  = ISKWRK(IrRU) + 1
      RSKWRK(ISKWRK(IrRUC)) = RU / 4.184E7
C! densities for phases
      ISKWRK(IrSDEN) = ISKWRK(IrRUC) + 1
C! temperatures for fit
      ISKWRK(IrKTMP)  = ISKWRK(IrSDEN) + NPHASE
C! thermodynamic properties
      ISKWRK(IrKTHM)  = ISKWRK(IrKTMP) + MAXTP *NKKTOT
C! densities for the species
      ISKWRK(IrKDEN)  = ISKWRK(IrKTHM) + NCP2T*NKKTOT
C! atomic weights of elements
      ISKWRK(IrAWT)   = ISKWRK(IrKDEN) + NKKTOT
C! molecular weights of species
      ISKWRK(IrKWT)   = ISKWRK(IrAWT)  + NELEM
C! Arrhenius parameters
      ISKWRK(IrPAR)   = ISKWRK(IrKWT)  + NKKTOT
C! coverage parameters
      ISKWRK(IrKCOV)   = ISKWRK(IrPAR)  + NIISUR*(NSPAR+1)
C! reverse parameters
      ISKWRK(IrRPAR)   = ISKWRK(IrKCOV)  + NIICOV*NSCOV
C! multiplicative constant in surface equilibrium constants
      ISKWRK(IrEQ)    = ISKWRK(IrRPAR) + NIIREV*(NSPAR+1)
C! real stoichometric coefficients
      ISKWRK(IrRNU)   = ISKWRK(IrEQ) + NIISUR
C! phase (site) balances
      ISKWRK(IrNCF)   = ISKWRK(IrRNU) + NIIRNU*MAXSPR
C! species order coefficients
      ISKWRK(IrKORD)   = ISKWRK(IrNCF) + NIISUR*NPHASE
C!
      ISKWRK(IrKFT)   = ISKWRK(IrKORD) + NIIORD*MAXORD
C!
      ISKWRK(IrKRT)   = ISKWRK(IrKFT) + NIISUR
C! dummy species work space
      ISKWRK(IrKT1)   = ISKWRK(IrKRT) + NIISUR
C! dummy species work space
      ISKWRK(IrKT2)   = ISKWRK(IrKT1)  + NKKTOT
C! dummy phase work space
      ISKWRK(IrPT1)   = ISKWRK(IrKT2)  + NKKTOT
C! dummy reaction work space
      ISKWRK(IrIT1)   = ISKWRK(IrPT1)  + NPHASE
C! dummy reaction work space
      ISKWRK(IrIT2)   = ISKWRK(IrIT1)  + NIISUR
C! dummy reaction work space
      ISKWRK(IrIT3)   = ISKWRK(IrIT2)  + NIISUR
C! total size req. for RSKWRK
      IrTOT   = ISKWRK(IrIT3)  + NIISUR-1
C
C! element names
      ISKWRK(IcENAM) = 1
C! species names
      ISKWRK(IcKNAM) = ISKWRK(IcENAM) + NELEM
C! material name
      ISKWRK(IcMNAM) = ISKWRK(IcKNAM) + NKKTOT
C! phase names
      ISKWRK(IcPNAM) = ISKWRK(IcMNAM) + 1
C! total size req. for CSKWRK
      IcTOT  = ISKWRK(IcPNAM) + NPHASE - 1
C
C     RSKWRK(ISKWRK(IrSKMN)) = SKMIN
C! dynes per sq. cm.
      RSKWRK(ISKWRK(IrPATM)) = PA
C! calories/ mole K
      RSKWRK(ISKWRK(IrRUC))  = RU / 4.184E7
C! ergs/ mole K
      RSKWRK(ISKWRK(IrRU))   = RU
C
      READ (LINKSK,err=444)
     1 ( CSKWRK(ISKWRK(IcENAM)+M-1),
     2   RSKWRK(ISKWRK(IrAWT) +M-1), M=1,NELEM)
C
      READ (LINKSK,err=555)
     1 ( CSKWRK(ISKWRK(IcKNAM)+K-1),
     2   RSKWRK(ISKWRK(IrKWT) +K-1),
     3   ISKWRK(ISKWRK(IiKPHS)+K-1),
     2   ISKWRK(ISKWRK(IiKCHG)+K-1),
     5   ISKWRK(ISKWRK(IiKNT) +K-1),
     6   ISKWRK(ISKWRK(IiNSCV)+K-1),
     3   RSKWRK(ISKWRK(IrKDEN)+K-1),
     4   (ISKWRK(ISKWRK(IiKCMP)+(K-1)*NELEM+M-1),M=1,NELEM),
     5   (RSKWRK(ISKWRK(IrKTMP)+(K-1)*MAXTP + L - 1), L=1,MAXTP),
     6   ( (RSKWRK(ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T + N-1),
     7       N=1,NCP2), L=1,(MAXTP-1)) , K=1,NKKTOT)
C
      READ (LINKSK,err=666) CSKWRK(ISKWRK(IcMNAM)),
     1 ( CSKWRK(ISKWRK(IcPNAM)+N-1),
     2   ISKWRK(ISKWRK(IiPKST)+N-1),
     3   ISKWRK(ISKWRK(IiPKND)+N-1),
     2   ISKWRK(ISKWRK(IiPTOT)+N-1),
     5   RSKWRK(ISKWRK(IrSDEN)+N-1),
     3   (RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR+I-1), I=1,NIISUR),
     7   N=1,NPHASE)
C
      IF (NIISUR .LE. 0) THEN
         IF (MORE .LE. 0) REWIND (LINKSK)
         RETURN
      ENDIF
C
      READ (LINKSK,err=888)
     1 ( ISKWRK(ISKWRK(IiNRPP)+I-1), ISKWRK(ISKWRK(IiNREA)+I-1),
     2   ISKWRK(ISKWRK(IiNSUM)+I-1),
     3   (RSKWRK(ISKWRK(IrPAR)+(I-1)*(NSPAR+1)+N-1),N=1,NSPAR),
     4   (ISKWRK(ISKWRK(IiNUNK)+(I-1)*MAXSPR+N-1),
     5   ISKWRK(ISKWRK(IiNU)+(I-1)*MAXSPR+N-1), N=1,MAXSPR),
     6   I=1,NIISUR)
      IF (NIICOV .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiICOV)+I-1), ISKWRK(ISKWRK(IiKCOV)+I-1),
     2    (RSKWRK(ISKWRK(IrKCOV)+(I-1)*NSCOV+L-1),L=1,NSCOV),
     3   I=1,NIICOV)
      IF (NIIREV .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiIREV)+I-1),
     2    (RSKWRK(ISKWRK(IrRPAR)+(I-1)*(NSPAR+1)+L-1),L=1,NSPAR),
     3   I=1,NIIREV)
      IF (NIISTK .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiISTK)+I-1), I=1,NIISTK)
      IF (NIIBHM .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiIBHM)+I-1), ISKWRK(ISKWRK(IiKBHM)+I-1),
     2   ISKWRK(ISKWRK(IiTBHM)+I-1), I=1,NIIBHM)
      IF (NIIRNU .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiIRNU)+I-1),
     2   (RSKWRK(ISKWRK(IrRNU) + (I-1)*MAXSPR + N-1),N=1,MAXSPR),
     3   I=1,NIIRNU)
      IF (NIIORD .GT. 0) READ (LINKSK,err=888)
     1   (ISKWRK(ISKWRK(IiIORD)+I-1),
     2   (ISKWRK(ISKWRK(IiKORD) + (I-1)*MAXORD + N-1),
     3    RSKWRK(ISKWRK(IrKORD) + (I-1)*MAXORD + N-1), N=1,MAXORD),
     4    I=1,NIIORD)
C
      DO 10 I = 1, NIISUR
C
C        PERTURBATION FATOR
         IPAR = ISKWRK(IrPAR) + (I-1)*(NSPAR+1) + NSPAR
         RSKWRK(IPAR) = 1.0
C
C        MULTIPLICATIVE FACTOR IN SURFACE EQLIBRIUM CONSTANT
         IEQ = ISKWRK(IrEQ) + I - 1
         RSKWRK(IEQ) = 1.0
   10 CONTINUE
C
      IF (NNSUR .LE. 0) THEN
         IF (MORE .LE. 0) REWIND (LINKSK)
         RETURN
      ENDIF
C
      DO 50 I = 1, NIISUR
         INU = ISKWRK(IiNU) + (I-1)*MAXSPR - 1
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR - 1
         IEQ = ISKWRK(IrEQ) + I - 1
         DO 40 N = NFSUR, NLSUR
            NUSUM = 0
            KST = ISKWRK (ISKWRK(IiPKST) + N - 1)
            KND = ISKWRK (ISKWRK(IiPKND) + N - 1)
            SDEN = RSKWRK(ISKWRK(IrSDEN) + N - 1)
            DO 30 K = KST, KND
               DO 20 M = 1, MAXSPR
                  NSTOIC = ISKWRK(INU + M)
                  KSPEC =   ISKWRK(INK + M)
                  IF (NSTOIC.NE.0 .AND. KSPEC.EQ.K) THEN
                     KCOV = ISKWRK (ISKWRK(IiNSCV) + K - 1)
C*****precision > double
                     RSKWRK(IEQ) = RSKWRK(IEQ) * DBLE(KCOV)**(-NSTOIC)
C*****END precision > double
C*****precision > single
C                     RSKWRK(IEQ) = RSKWRK(IEQ) * REAL(KCOV)**(-NSTOIC)
C*****END precision > single
                     NUSUM = NUSUM + NSTOIC
                  ENDIF
   20          CONTINUE
   30       CONTINUE
            RSKWRK(IEQ) = RSKWRK(IEQ) * SDEN**NUSUM
   40    CONTINUE
   50 CONTINUE
C
      IF (NIIRNU .LE. 0) THEN
         IF (MORE .LE. 0) REWIND (LINKSK)
         RETURN
      ENDIF
C
      DO 150 L = 1, NIIRNU
         I = ISKWRK (ISKWRK(IiIRNU) + L - 1)
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR - 1
         INU = ISKWRK(IrRNU)  + (L-1)*MAXSPR - 1
         IEQ = ISKWRK(IrEQ) + I - 1
C
         DO 140 N = NFSUR, NLSUR
            KST = ISKWRK (ISKWRK(IiPKST) + N - 1)
            KND = ISKWRK (ISKWRK(IiPKND) + N - 1)
            SDEN = RSKWRK (ISKWRK(IrSDEN) + N - 1)
            RNUSUM = 0.0
            DO 130 K = KST, KND
               DO 120 M = 1, MAXSPR
                  STOICH = RSKWRK(INU + M)
                  KSPEC  = ISKWRK(INK + M)
                  IF (ABS(STOICH).GT.SKMIN .AND. KSPEC.EQ.K) THEN
                     KCOV = ISKWRK (ISKWRK(IiNSCV) + K - 1)
C*****precision > double
                     RSKWRK(IEQ) = RSKWRK(IEQ) * DBLE(KCOV)**(-STOICH)
C*****END precision > double
C*****precision > single
C                     RSKWRK(IEQ) = RSKWRK(IEQ) * REAL(KCOV)**(-STOICH)
C*****END precision > single
                     RNUSUM = RNUSUM + STOICH
                  ENDIF
  120          CONTINUE
  130       CONTINUE
            RSKWRK(IEQ) = RSKWRK(IEQ) * SDEN**RNUSUM
  140    CONTINUE
  150 CONTINUE
C
      IF (MORE .LE. 0) REWIND (LINKSK)
      RETURN
C
  444 WRITE (LOUT,*) ' Error reading element data...'
c      IFLAG = 3
      RETURN
  555 WRITE (LOUT,*) ' Error reading species data...'
c      IFLAG = 4
      RETURN
  666 WRITE (LOUT,*) ' Error reading phase data...'
c      IFLAG = 5
      RETURN
  888 WRITE (LOUT,*) ' Error reading reaction data...'
c      IFLAG = 6
      RETURN
  999 WRITE (LOUT,*) ' Error reading Surface binary file...'
c      IFLAG = 1
      RETURN
C
  300 FORMAT (10X,'ISKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  350 FORMAT (10X,'RSKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  375 FORMAT (10X,'CSKWRK MUST BE DIMENSIONED AT LEAST ',I5)
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKINU   (IS, NDIM, ISKWRK, RSKWRK, NSPEC, KI, NU)
C
C  START PROLOGUE
C
C  SUBROUTINE SKINU   (IS, NDIM, ISKWRK, RSKWRK, NSPEC, KI, NU)
C     Returns the number of species in a surface reaction, and the
C     species indices and stoichiometric coefficients.
C
C  INPUT
C     IS     - Index number of a surface reaction;  IS must be greater
C              than 0, and less than or equal to NIISUR, the number of
C              reactions in the surface mechanism.
C                   Data type - integer scalar
C     NDIM   - Dimension of the arrays KI and NU;  NDIM must be at
C              MAXSPR, the number of species allowed in a surface
C              reaction.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C
C  OUTPUT
C     NSPEC  - Number of species in surface reaction IS.
C                   Data type - integer scalar
C     KI     - Array of species indices for the species in surface
C              reaction IS;  KI(N) is the index of the Nth species
C              in reaction IS.
C                   Data type - integer array
C                   Dimension KI(*) at least MAXSPR, the number of
C                   species allowed in a surface reaction.
C     NU     - Array of stoichiometric coefficients for the species
C              in surface reaction IS.  NU(N) is the stoichiometric
C              coefficient of the Nth species in surface reaction IS;
C              NU is negative if the Nth species is a reactant;
C              NU is positive if the Nth species is a product.
C                   Data type - integer array
C                   Dimension NU(*) at least MAXSPR, the number of
C                   species allowed in a surface reaction.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), KI(*), NU(*)
      INCLUDE 'skstrt.inc'
C
      NSPEC = 0
      DO 50 N = 1, NDIM
         KI(N) = 0
         NU(N) = 0
   50 CONTINUE
C
C
c......ashish modification for small surf mechs
c      IF (ISKWRK(IiNIIS).LE.0 .OR. IS.LE.0 .OR.
c     1    IS.GT.ISKWRK(IiNIIS).OR. MAXSPR.GT.NDIM) RETURN
      IF (ISKWRK(IiNIIS).LE.0 .OR. IS.LE.0 .OR.
     1    IS.GT.ISKWRK(IiNIIS)) RETURN
C
      DO 200 N = 1, MAXSPR
         K = ISKWRK(ISKWRK(IiNUNK) + (IS-1)*MAXSPR + N - 1)
         IF (K .NE. 0) THEN
            NSPEC = NSPEC + 1
            KI(NSPEC) = K
            NU(NSPEC) = ISKWRK(ISKWRK(IiNU) + (IS-1)*MAXSPR + N - 1)
         ENDIF
200   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIORD (IDIM, KDIM, ISKWRK, RSKWRK, NFORD, IFORD, FORD,
     1                   NRORD, IRORD, RORD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIORD (IDIM, KDIM, ISKWRK, RSKWRK, NFORD, IFORD, FORD,
C                     NRORD, IRORD, RORD)
C     Returns the number and indices of surface reactions with modified
C     species orders, and the order values for the species in the
C     surface mechanism.
C
C  INPUT
C     IDIM   - Dimension of arrays IFORD and IRORD, and the second
C              dimension of the matrices FORD and RORD; IDIM must be
C              at least NORD, the number of surface reactions with
C              modified species orders.
C                   Data type - integer scalar
C     KDIM   - First dimension of the arrays FORD and RORD; KDIM must
C              be at least NKKTOT, the number of species in the
C              reaction mechanism.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C
C  OUTPUT
C     NFORD  - Number of surface reactions with modified forward
C              species orders.
C                   Data type - integer scalar
C     IFORD  - Array of indices of surface reactions with modified
C              forward species orders.
C                   Data type - integer array
C                   Dimension IFORD(*) at least NORD, the number
C                   of surface reactions with modified orders.
C     FORD   - Matrix of the modified forward species orders for the
C              NFORD surface reactions;  FORD(K,N) is the forward order
C              of the Kth species for the Nth reaction with modified
C              forward species orders.
C                   Data type - real matrix
C                   Dimension FORD(*,*) at least NKKTOT for the first
C                   dimension, the number of species in the
C                   mechanism, and at least NORD for the second, the
C                   number of surface reactions with modified orders.
C     NRORD  - Number of surface reactions with modified reverse
C              species orders.
C                   Data type - integer scalar
C     IRORD  - Array of indices of surface reactions with modified
C              reverse species orders.
C                   Data type - integer array
C                   Dimension IRORD(*) at least NORD, the number
C                   of surface reactions with modified orders.
C     RORD   - Matrix of the modified reverse species orders for the
C              NRORD surface reactions;  RORD(K,N) is the reverse order
C              of the Kth species for the Nth reaction with modified
C              reverse species orders.
C                   Data type - real matrix
C                   Dimension RORD(*,*) at least NKKTOT for the first
C                   dimension, the number of species in the
C                   mechanism, and at least NORD for the second, the
C                   number of surface reactions with modified orders.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), IFORD(*), FORD(KDIM,*),
     1          IRORD(*), RORD(KDIM,*)
      LOGICAL LFORD, LRORD
      INCLUDE 'skstrt.inc'
C
      NFORD = 0
      NRORD = 0
      DO 100 N = 1, IDIM
         IFORD(N) = 0
         IRORD(N) = 0
         DO 100 K = 1, KDIM
            FORD(K,N) = 0.0
            RORD(K,N) = 0.0
  100 CONTINUE
C
      IF (IDIM.LE.0 .OR. IDIM.LT.ISKWRK(IiNORD) .OR.
     1    KDIM.LT.ISKWRK(IiKTOT)) RETURN
C
      DO 200 N = 1, ISKWRK(IiNORD)
         LFORD = .FALSE.
         LRORD = .FALSE.
         I  = ISKWRK(ISKWRK(IiIORD) + N - 1)
         DO 150 K = 1, MAXORD
            KSPEC = ISKWRK(ISKWRK(IiKORD) + MAXORD*(N-1) + K - 1)
            IF (KSPEC .LT. 0) LFORD = .TRUE.
            IF (KSPEC .GT. 0) LRORD = .TRUE.
  150    CONTINUE
         IF (LFORD) THEN
            NFORD = NFORD + 1
            IFORD(NFORD) = I
            DO 160 K = 1, MAXORD
               KSPEC = ISKWRK(ISKWRK(IiKORD) + MAXORD*(N-1) + K - 1)
               ORD   = RSKWRK(ISKWRK(IrKORD) + MAXORD*(N-1) + K - 1)
               IF (KSPEC .LT. 0) FORD(KSPEC,NFORD) = ORD
  160       CONTINUE
         ENDIF
         IF (LRORD) THEN
            NRORD = NRORD + 1
            IRORD(NRORD) = I
            DO 170 K = 1, MAXORD
               KSPEC = ISKWRK(ISKWRK(IiKORD) + MAXORD*(N-1) + K - 1)
               ORD   = RSKWRK(ISKWRK(IrKORD) + MAXORD*(N-1) + K - 1)
               IF (KSPEC .GT. 0) RORD(KSPEC,NRORD) = ORD
  170       CONTINUE
         ENDIF
  200 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIREV (IR, ISKWRK, RSKWRK, IREV, RAR, RBR, RER)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIREV (IR, ISKWRK, RSKWRK, IREV, RAR, RBR, RER)
C     Returns an integer flag to indicate whether reaction IR
C     has an explicitly assigned reverse rate constant.
C     It also returns the reverse Arrhenius expression values for
C     surface reaction, IR, if it was explicitly assigned
C     in the surface chemkin interpretor.  If the reverse
C     Arrhenius values were not explicity assigned, RAR,
C     RBR, and RER are returned with zero values.
C
C  INPUT
C     IR     - Surface reaction index number.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     IREV   - 1 if reaction IR has an explicitly assigned reverse
C                rate constant
C              0 if reaction IR does not have an explicitly assigned
C                reverse rate constant
C                   Data type - integer scalar
C     RAR    - Pre-exponential constants for the reverse reaction.
C                   cgs units - mole-cm-sec-K
C                   Data type - real scalar
C     RBR    - Temperature dependence exponents for the reverse
C              reaction.
C                   cgs units - none
C                   Data type - real scalar
C     RER    - Activation energy for the reverse reaction.
C                   cgs units - Kelvins
C                   Data type - real scalar
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      RAR = 0.0
      RBR = 0.0
      RER = 0.0
      IREV = 0
C
      DO 10 N = 1, ISKWRK(IiNREV)
         IF (IR .EQ. ISKWRK(ISKWRK(IiIREV) + N - 1)) THEN
            IREV = 1
            IND = ISKWRK(IrRPAR) + (N-1)*(NSPAR+1)
            RAR = RSKWRK(IND)
            RBR = RSKWRK(IND + 1)
            RER = RSKWRK(IND + 2)
         ENDIF
10    CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKIRNU (IDIM, NDIM, ISKWRK, RSKWRK, NIRNU, IRNU, NSPEC,
     1                   KI, RNU)
C
C  START PROLOGUE
C
C  SUBROUTINE SKIRNU (IDIM, NDIM, ISKWRK, RSKWRK, NIRNU, IRNU, NSPEC,
C 1                   KI, RNU)
C     Returns the number and indices of surface reactions with real
C     stoichiometric coefficients, number of species in the reactions,
C     and the species indices and coefficients;
C
C  INPUT
C     IDIM   - Dimension of the arrays IRNU and NSPEC, and the second
C              dimension of the arrays KI, and RNU;  IDIM must be
C              at least NIRNU, the number of surface reactions with
C              real stoichiometric coefficients.
C                   Data type - integer scalar
C     NDIM   - First dimension of the matrices KI and RNU; NDIM must
C              be at least MAXSPR, the number of species allowed in
C              a surface reaction.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C
C  OUTPUT
C     NIRNU  - Number of reactions with real stoichiometric
C              coefficients.
C                   Data type - integer scalar
C     IRNU   - Array of indices of reactions with real
C              stoichiometric coefficients; IRNU(N) is the index of
C              the Nth reaction with real stoichiometric coefficients.
C                   Data type - integer array
C                   Dimension IRNU(*) at least NIRNU, the number of
C                   reactions with real stoichiometric
C                   coefficients.
C     NSPEC  - Array of the number of species in the surface
C              reactions with real stoichiometric coefficients;
C              NSPEC(N) is the number of species in the Nth surface
C              reaction with real stoichiometric coefficients.
C                   Data type - integer array
C                   Dimension NSPEC(*) at least NIRNU, the number
C                   of reactions with real stoichiometric
C                   coefficients.
C     KI     - Matrix of species indices for the species in the
C              NIRNU reactions; KI(M,N) is the species index of
C              the Mth species in the Nth reaction with real
C              stoichiometric coefficients.
C                   Data type - integer matrix
C                   First dimension of KI(NDIM,*) must be at least
C                   MAXSPR, the number of species allowed in a
C                   surface reaction; second dimension of KI must
C                   be at least NIRNU, the number of reactions
C                   with real stoichiometric coefficients.
C     RNU    - Matrix of stoichiometric coefficients for the species
C              in the NIRNU reactions.  RNU(M,N) is the
C              stoichiometric coefficient of the Mth species in
C              the Nth reaction with real stoichiometric
C              coefficients;
C              RNU(M,*) is negative if the Mth species is a reactant;
C              RNU(M,*) is positive if the Mth species is a product.
C                   Data type - real matrix
C                   First dimension RNU(NDIM,*) must be at least
C                   MAXSPR, the number of species allowed in a
C                   surface reaction; second dimension of KI must
C                   be at least NIRNU, the number of reactions
C                   with real stoichiometric coefficients.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), IRNU(*), NSPEC(*),
     1          KI(NDIM,*), RNU(NDIM,*)
      INCLUDE 'skstrt.inc'
C
      DO 100 N = 1, IDIM
         NSPEC(N) = 0
         IRNU(N) = 0
         DO 100 M = 1, NDIM
            KI(M,N) = 0
            RNU(M,N) = 0.0
  100 CONTINUE
C
      IF (ISKWRK(IiNRNU).LE.0 .OR. IDIM.LT.ISKWRK(IiNRNU) .OR.
     1    NDIM.LT.MAXSPR) RETURN
C
      NIRNU = ISKWRK(IiNRNU)
      DO 200 N = 1, ISKWRK(IiNRNU)
         I = ISKWRK(ISKWRK(IiIRNU) + N - 1)
         IRNU(N) = I
         NSPEC(N) = 0
C
         DO 200 M = 1, MAXSPR
            K = ISKWRK(ISKWRK(IiNUNK) + (I-1)*MAXSPR + M - 1)
            IF (K .NE. 0) THEN
               NSPEC(N) = NSPEC(N) + 1
               KI(NSPEC(N),N) = K
               RNU(NSPEC(N),N) =
     1            RSKWRK(ISKWRK(IrRNU) + (N-1)*MAXSPR + M - 1)
            ENDIF
  200 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKISTK (IR, ISKWRK, ISTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKISTK (IR, ISKWRK, ISTFL)
C     Returns an integer flag to indicate whether reaction IR uses
C     sticking coefficients.
C
C  INPUT
C     IR     - Surface reaction index number.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     ISTFL  - 0 if reaction IR does not use sticking coefficients
C              1 if reaction IR does use sticking coefficients
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
      DIMENSION ISKWRK(*)
      INCLUDE 'skstrt.inc'
C
      ISTFL = 0
      DO 10 N = 1, ISKWRK(IiNSTK)
         IF (IR .EQ. ISKWRK(ISKWRK(IiISTK) + N - 1) ) ISTFL = 1
   10 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKLEN  (LINKSK, LOUT, LENI, LENR, LENC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKLEN  (LINKSK, LOUT, LENI, LENR, LENC, IFLAG)
C     Reads the first record of the binary file to return the lengths
C     required for the integer, real, and character work arrays.
C
C  INPUT
C     LINKSK  - Unit number assigned to binary file
C                   Data type - integer scalar
C     LOUT   - Unit number assigned for output
C                   Data type - integer scalar
C
C  OUTPUT
C     LENI   - Dimension required for ISKWRK, the Surface Chemkin
C              integer work array.
C                   Data type - integer scalar
C     LENR   - Dimension required for RSKWRK, the Surface Chemkin
C              real work array.
C                   Data type - integer scalar
C     LENC   - Dimension required for CSKWRK, the Surface Chemkin
C              character work array.
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
      PARAMETER (NLIST=5)
      LOGICAL KERR, VOK, POK
      CHARACTER*16 LIST(NLIST), VERS, PREC
      COMMON /SKCONS/ VERS, PREC, KERR
C      DATA LIST/'4.00','4.01','4.02','4.03','4.04','4.05'/
C      DATA LIST/'4.06','4.07','4.08'/
C      DATA LIST/'4.08b','4.08c','4.09'/
      DATA LIST /'4.2', '4.3', '4.31','4.4','4.5'/
C
      VERS = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0
C
c      IFLAG = 0
      KERR = .FALSE.
      MORE = 0
      READ (LINKSK,err=999) VERS, PREC, KERR, MORE
C
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (VERS .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE
C
       POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the surface binary file...',
     2      ' Check SURFACE INTERPRETER output for error conditions.'
         ENDIF
         IF (.NOT. VOK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Surface binary file is incompatible with Surface',
     3      ' Library Version 5.0'
         ENDIF
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Precision of Surface binary file does not agree with',
     2      ' precision of Surface library'
         ENDIF
c         IFLAG = 20
         BACKSPACE (LINKSK)
         RETURN
      ENDIF
C
      READ (LINKSK,err=1111) LENI, LENR, LENC
C
      BACKSPACE (LINKSK)
      BACKSPACE (LINKSK)
      RETURN
C
  999 CONTINUE
      WRITE (LOUT,*) ' Error reading binary Surface file...'
c      IFLAG = 1
      BACKSPACE (LINKSK)
      RETURN
 1111 CONTINUE
      WRITE (LOUT,*) ' Error reading binary Surface file...'
c      IFLAG = 2
      BACKSPACE (LINKSK)
      BACKSPACE (LINKSK)
      RETURN
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKMXTP (ISKWRK, MXTP)
C
C  START PROLOGUE
C
C  SUBROUTINE SKMXTP (ISKWRK, MXTP)
C     Returns the maximum number of temperatures used in
C     fitting the thermodynamic properties of the species.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     MXTP    - Maximum number of temperatures used in
C              fitting the thermodynamic properties of
C              the species.
C                   Date type - integer scalar
C                   cgs units:  none
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*)
      INCLUDE 'skstrt.inc'
C
      MXTP = MAXTP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNCF  (NELDIM, ISKWRK, NEL)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNCF (NELDIM, ISKWRK, NEL)
C     Returns the elemental composition of the species.
C
C  INPUT
C     NELDIM - First dimension of the matrix NEL.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     NEL    - Elemental compositions of the species.
C              NEL(M,K) is the number of atoms of element M in
C                       species K.
C                   Data type - integer array
C                   The first dimension of NEL(*,*) must be exactly
C                   NELDIM, which is at least NELEM,  the  number of
C                   elements in the problem, and the second dimension at
C                   least KKTOT, the total number of species.
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
      DIMENSION ISKWRK(*), NEL(NELDIM,*)
      INCLUDE 'skstrt.inc'
C
      IF (NELDIM .LT. NELEM) RETURN
      DO 100 M = 1, NELEM
         DO 100 N = 1, ISKWRK(IiKTOT)
            NEL(M,N) = ISKWRK(ISKWRK(IiKCMP) + (N-1)*NELEM + M - 1)
  100 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNCON (ISKWRK, RSKWRK, NCON)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNCON (ISKWRK, RSKWRK, NCON)
C     Returns the total number of surface reactions which do not
C     conserve sites for each phase.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     NCON   - Number of surface reactions which do not conserve sites
C              for each phase.
C                   Date type - integer array
C                   Dimension NCON(*) at least NPHASE,
C                   the total number of phases
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), NCON(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      NIISUR = ISKWRK(IiNIIS)
      DO 50 N = 1, ISKWRK(IiNPHA)
         NCON(N) = 0
         DO 25 I = 1, NIISUR
            RNCF = RSKWRK (ISKWRK(IrNCF) + (N-1)*NIISUR + I - 1)
            IF (ABS(RNCF) .GT. RSKWRK(ISKWRK(IrSKMN)))
     1         NCON(N) = NCON(N) + 1
   25    CONTINUE
   50 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNU   (IDIM, ISKWRK, RSKWRK, KSTOIC, NSTOIC)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNU   (IDIM, ISKWRK, RSKWRK, KSTOIC, NSTOIC)
C     Returns the stoichiometric coefficients of the species and the
C     net change in phases for all of the surface reactions in a
C     mechanism.
C
C  INPUT
C     IDIM   - First dimension of the array NSTOIC.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     KSTOIC - Stoichiometric coefficients for the surface reactions.
C                   cgs units - none
C                   Data type - integer array
C                   The first dimension of NSTOIC(IDIM,*) must be
C                   exactly IDIM, which is at least IISUR, the total
C                   number of surface reactions, and at least KKTOT
C                   for the second dimension, the total number of
C                   species.
C     NSTOIC - Net change of the phases for the surface reactions.
C                   cgs units - none
C                   Data type - integer array
C                   The first dimension of NSTOIC(IDIM,*) must be
C                   exactly IDIM, which is at least IISUR, the total
C                   number of surface reactions, and at least NPHASE
C                   for the second dimension, the total number of
C                   phases.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), KSTOIC(IDIM,*), NSTOIC(IDIM,*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C             PROCESS EACH REACTION
C
      NIISUR = ISKWRK(IiNIIS)
      DO 100 I = 1, NIISUR
C
         DO 10 K = 1, ISKWRK(IiKTOT)
            KSTOIC(I,K) = 0
   10    CONTINUE
C
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR
         INU = ISKWRK(IiNU)   + (I-1)*MAXSPR
         DO 20 N = 1, MAXSPR
            NUNK = ISKWRK (INK + N - 1)
            IF (NUNK .GT. 0) THEN
               NU = ISKWRK(INU + N - 1)
               KSTOIC(I,NUNK) = KSTOIC(I,NUNK)+NU
            ENDIF
   20    CONTINUE
C
         DO 30 N = 1, ISKWRK(IiNPHA)
            NSTOIC(I,N) = RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR + I - 1)
     1                  + RSKWRK(ISKWRK(IrSKMN))
   30    CONTINUE
C
100   CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKNUF   (IDIM, ISKWRK, KSTOIF)
C
C  START PROLOGUE
C
C  SUBROUTINE SKNUF   (IDIM, ISKWRK, KSTOIF)
C     Returns the stoichiometric coefficients of the species
C     for all reactants in all surface reactions in a mechanism.
C     (note - reactants only! - they will all be negative)
C
C  INPUT
C     IDIM   - First dimension of the array NSTOIC.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     KSTOIF - Stoichiometric coefficients for the reactants
C               in all surface reactions.
C                   cgs units - none
C                   Data type - integer array
C                   The first dimension of KSTOIF(IDIM,*) must be
C                   exactly IDIM, which is at least IISUR, the total
C                   number of surface reactions, and at least KKTOT
C                   for the second dimension, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), KSTOIF(IDIM,*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C             PROCESS EACH REACTION
C
      DO 100 I = 1, ISKWRK(IiNIIS)
C
         DO 10 K = 1, ISKWRK(IiKTOT)
            KSTOIF(I,K) = 0
   10    CONTINUE
C
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR
         INU = ISKWRK(IiNU)   + (I-1)*MAXSPR
         DO 20 N = 1, MAXSPR/2
            NUNK = ISKWRK (INK + N - 1)
            IF (NUNK .GT. 0) THEN
               NU = ISKWRK(INU + N - 1)
               KSTOIF(I,NUNK) = KSTOIF(I,NUNK)+NU
            ENDIF
   20    CONTINUE
C
100   CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPCMP (ISTR, IRAY, NN, SETS, NSETS, ISET, IND, NT)
C     This subroutine can do everything that the subroutine SKCOMP can
C     do, and additionally, has the capabilities of separating the
C     elements of IRAY into categories and then search IRAY by element
C     and category.  The categories that each element of IRAY will be
C     assigned to are specified by the input character string vector
C     SETS of vector length NSETS.  Elements of each category in IRAY
C     must be grouped congruously.  The number of elements in each
C     category within IRAY is specified by the input integer vector
C     ISET.  To search for the existence of an element within a
C     category ISTR may additionally be composed of two substrings,
C     ISTR="ELEMENT_NAME/CATEGORY_NAME/", where CATEGORY_NAME is one
C     of the categories specified in SETS.  In this case, IND will
C     return the first position in IRAY where ELEMENT_NAME occurred
C     within the category CATEGORY_NAME.  NT will return the total
C     number of times ELEMENT_NAME occurred within the category
C     CATEGORY_NAME.  If ELEMENT_NAME is not found within the
C     specified category, IND and NT are returned with a value of
C     zero.  If no category is specified within ISTR, IND and NT
C     return with the same values as they would from subroutine SKCOMP.
C
C     Consider the following example,
C        IRAY = {"RED", "BLUE", "JADE", "RUBY", "TOPAZ", "JADE"}
C        NN = 6
C        SETS = {"COLORS", "STONES"},
C        NSETS = 2
C        ISET = {4, 2}.
C     This assumes that the elements of IRAY were grouped into two
C     sets, consisting of 4 and 2 elements, respectively, and the
C     following names
C        "COLORS"  = {"RED", "BLUE", "JADE", "RUBY"}, and
C        "STONES"  = {"TOPAZ", "JADE"}.
C
C     If ISTR="BLUE" then IND=2 and NT=1;
C     if ISTR="PINK" then IND=0 and NT=0; and
C     if ISTR="JADE",then IND=3 and NT=2.
C
C     If ISTR="BLUE/COLORS/" then IND=2 and NT=1;
C     if ISTR="BLUE/STONES/" then IND=0 and NT=0;
C     if ISTR="JADE/GEMS/"   then IND=0 and NT=0; and
C     if ISTR="JADE/STONES/",then IND=6 and NT=1.
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
C     NSETS   - Number of entries in SETS(*)
C                    Data type - integer scalar
C     ISET(*) - Integer total number of entries in a subset of IRAY.
C                    Data type - integer array
C                    Dimension ISET(*) at least NSETS
C
C  OUTPUT
C     IND      - Index of the position in IRAY(*) containing ISTR.
C                If ISTR is not in IRAY(*), IND = 0.
C                If the slash-delimited substring of ISTR is not
C                in SETS(*), IND = 0.
C                If the slash-delimited substring of ISTR is in
C                SETS(N), but the substring before the slash is
C                not a member of the subset associated with SETS(N),
C                IND = 0, whether or not the substring is in IRAY(*).
C                   Data type - integer scalar
C     NT       - Total number of times ISTR occurs in IRAY(*),
C                or total number of times ISTR occurs in a subset
C                of IRAY(*).
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
      NT  = 0
C
      IF (ISTR .EQ. ' ') RETURN
C
      I2 = ILASCH(ISTR)
      IF (ISTR(I2:I2) .NE. '/') THEN
         CALL SKCOMP (ISTR, IRAY, NN, IND, NT)
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
                NT = NT + 1
            ENDIF
   50 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
C
C  START PROLOG
C
C  SUBROUTINE SKPKK  (ISKWRK, KKPHAS, KFIRST, KLAST)
C     Returns arrays of species pointers for the phases.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C
C  OUTPUT
C     KKPHAS - The total number of species in each phase.
C                   Data type - integer array
C                   Dimension KKPHAS(*) at least NPHASE,
C                   the total number of phases
C     KFIRST - The index of the first species in each phase.
C                   Data type - integer array
C                   Dimension KFIRST(*) at least NPHASE,
C                   the total number of phases
C     KLAST  - The index of the last species in each phase.
C                   Data type - integer array
C                   Dimension KLAST(*) at least NPHASE,
C                   the total number of phases
C
C  END PROLOG
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), KFIRST(*), KLAST(*), KKPHAS(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNPHA) .LE. 0) RETURN
C
      DO 50 N = 1, ISKWRK(IiNPHA)
         KKPHAS(N) = ISKWRK(ISKWRK(IiPTOT) + N - 1)
         KFIRST(N) = ISKWRK(ISKWRK(IiPKST) + N - 1)
         KLAST(N)  = ISKWRK(ISKWRK(IiPKND) + N - 1)
   50 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKPNT  (LSAVE, LOUT, V, P, LENI, LENR, LENC, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKPNT  (LSAVE, LOUT, VERS, PREC, LENI, LENR, LENC, IERR)
C     Reads from a binary file information about a Surface Chemkin
C     binary file, pointers for the Surface Chemkin Library, and
C     returns lengths of work arrays.
C
C  INPUT
C     LSAVE  - Integer input unit for binary data file.
C                   Data type - integer scalar
C     LOUT   - Integer output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     VERS   - Version number of the Surface Chemkin binary file.
C                   Data type - real scalar
C     PREC   - Machine precision of the Surface Chemkin binary file.
C                   Data type - character string
C     LENI   - Minimum length required for the integer work array.
C                   Data type - integer scalar
C     LENR   - Minimum length required for the real work array.
C                   Data type - integer scalar
C     LENC   - Minimum length required for the character work array.
C                   Data type - integer scalar
C     KERR   - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'skstrt.inc'
C
      COMMON /MACHN/ SMALL, BIG, EXPARG
      LOGICAL KERR, IERR
      CHARACTER*16 VERS, PREC, V, P
      COMMON /SKCONS/ VERS, PREC, KERR
C----------------------------------------------------------------------C
C     Data about the machine dependent constants is carried in         C
C     COMMON/MACHN/SMALL,BIG,EXPARG
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      EXPARG = LOG(BIG)
C
C----------------------------------------------------------------------C
C
      KERR = .FALSE.
      IERR = KERR
C
      READ (LSAVE, ERR=100) VERS, PREC, LENI, LENR, LENC,
C
C     Integer constants
C
     1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, MAXORD, MAXTP, NCP,
     2   NCP1, NCP2, NCP2T,
C
C     ISKWRK pointers to integer variables
C
     3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
     4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
     5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiMOTZ,
C
C     ISKWRK pointers to integer arrays
C
     6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
     7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
     8   IiISTK, IiIBHM, IiKBHM, IiTBHM, IiIRNU, IiIORD, IiKORD,
C
C     ISKWRK pointers to real variables
C
     *   IrSKMN, IrPATM, IrRU,   IrRUC,
C
C     ISKWRK pointers to real arrays
C
     1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
     2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
     3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,
C
C     ISKWRK pointers to character arrays
C
     4   IcENAM, IcKNAM, IcMNAM, IcPNAM
C
C     END include file for sklib.f
C
      V  = VERS
      P  = PREC
      IERR = KERR
      RETURN
C
  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Surface binary file data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      VERS   = ' '
      V      = VERS
      PREC   = ' '
      P      = PREC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRAEX (IR, ISKWRK, RSKWRK, RA)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRAEX (IR, ISKWRK, RSKWRK, RA)
C
C     Returns the Pre-exponential rate constant
C     (or sticking coefficient) of the IRth reaction, or changes its
C     value, depending on the sign of IR.
C
C  INPUT
C      IR     - Integer reaction number; IR>0 gets RA(I) from RSKWRK,
C                                        IR<0 puts RA(I) into RSKWRK.
C                    Data type - integer scalar
C
C      ISKWRK - Array of integer internal work space.
C                    Data type - integer array
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      RA     - Pre-exponential coefficient, or sticking coefficient,
C               for IRth reaction
C                    Data type - real scalar
C                    cgs units:
C                       mole-cm-sec-K for normal rate expressions
C                       none for sticking coefficient reactions
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1)
      IF (IR .GT. 0) THEN
         RA = RSKWRK(NI)
      ELSE
         RSKWRK(NI) = RA
      ENDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRBEX (IR, ISKWRK, RSKWRK, RB)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRBEX (IR, ISKWRK, RSKWRK, RB)
C.......Created by ashish.
C
C     Returns the temperature exponent (beta in T**beta)
C     of the IRth reaction, or changes its
C     value, depending on the sign of IR.
C
C  INPUT
C      IR     - Integer reaction number; IR>0 gets RB(I) from RSKWRK,
C                                        IR<0 puts RB(I) into RSKWRK.
C                    Data type - integer scalar
C
C      ISKWRK - Array of integer internal work space.
C                    Data type - integer array
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      RB     - Temperature exponent beta,
C               for IRth reaction
C                    Data type - real scalar
C                    cgs units:
C                       unitless
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1)
      IF (IR .GT. 0) THEN
         RB = RSKWRK(NI+1)
      ELSE
         RSKWRK(NI+1) = RB
      ENDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKREEX (IR, ISKWRK, RSKWRK, RE)
C
C  START PROLOGUE
C
C  SUBROUTINE SKREEX (IR, ISKWRK, RSKWRK, RE)
C.......Created by J.E. Sutton, 2012/07/27
C
C     Returns the activation energy
C     of the IRth reaction, or changes its
C     value, depending on the sign of IR.
C
C  INPUT
C      RE     - Activation energy
C               Data type - real scalar
C               cgs units - K
C      IR     - Integer reaction number; IR>0 gets RE(I) from RSKWRK,
C                                        IR<0 puts RE(I) into RSKWRK.
C                    Data type - integer scalar
C
C      ISKWRK - Array of integer internal work space.
C                    Data type - integer array
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      RE     - Activation energy
C               for IRth reaction
C                    Data type - real scalar
C                    cgs units:
C                       Kelvin
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1)
      IF (IR .GT. 0) THEN
         RE = RSKWRK(NI+2)
      ELSE
         RSKWRK(NI+2) = RE
      ENDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRAT  (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT,
     1                   Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRAT  (P, T, ACT, SDEN, ISKWRK, RSKWRK, SDOT, SITDOT)
C     Returns production rates for the species and sites.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     SDOT   - Production rates of the species.
C              1) For K=1,KKGAS, SDOT(K) is the production rate of
C                 gas-phase species K in (moles/cm**2-sec).
C              2) For K=KKGAS+1,KKGAS+KKSUR, SDOT(K) is the production
C                 rate of surface species K in (moles/cm**2-sec).
C              3) For K=KKGAS+KKSUR+1,KKTOT, SDOT(K) is the production
C                 rate of bulk species K in (moles/cm**2-sec).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SDOT(*) at least KKTOT, the total
C                   number of species.
C    SITDOT  - Production rates of the surface phases (subroutine
C              calculates entries for the surface site phases only).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SITDOT(*) at least NPHASE, the total
C                   number of phases.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), SDOT(*),
     1          SITDOT(*)
        character*80 crapa(100)
      INCLUDE 'skstrt.inc'
C
      DO 10 K = 1, ISKWRK(IiKTOT)
         SDOT(K) = 0.0
   10 CONTINUE
C
      SDTOT = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 20 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = 0.0
            SDTOT = SDTOT + SDEN(N)
   20    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
C        PROCESS EACH REACTION
C
      DO 100 I = 1, NIISUR
C
         ROP = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
         INK = ISKWRK(IiNUNK) + (I-1) * MAXSPR
         INU = ISKWRK(IiNU)   + (I-1) * MAXSPR
C
         DO 50 N = 1, MAXSPR
            NUNK = ISKWRK(INK + N - 1)
            IF (NUNK .NE. 0) SDOT(NUNK) = SDOT(NUNK) + ROP *
C*****precision > double
     1         DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C     1         REAL (ISKWRK(INU + N - 1))
C*****END precision > single
50       CONTINUE
C
         IF (NKKSUR .GT. 0) THEN
            DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
               RNCF = RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR + I - 1)
               SITDOT(N) = SITDOT(N) + ROP * RNCF
   60       CONTINUE
         ENDIF
100   CONTINUE
C
      IF (NIIRNU .LE. 0) RETURN
C
      DO 200 L = 1, NIIRNU
C
         I = ISKWRK(ISKWRK(IiIRNU) + L - 1)
         ROP = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
         INK = ISKWRK(IiNUNK) + (I-1) * MAXSPR
         INU = ISKWRK(IrRNU)  + (L-1) * MAXSPR
C
         DO 200 N = 1, MAXSPR
            NUNK = ISKWRK(INK + N - 1)
            IF (NUNK .NE. 0) SDOT(NUNK) = SDOT(NUNK) + ROP *
     1          RSKWRK(INU + N - 1)
  200 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRATI (IR, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRATI (IR, ROP, ISKWRK, RSKWRK, SDOTI, SITDTI)
C     Returns rates of production for each of the species by
C     surface reaction IR.
C
C  INPUT
C     IR     - Reaction index
C                   Data type - integer scalar
C     ROP    - Rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec).
C                   Data type - real array
C                   Dimension ROP(*) at least IISUR, the total number
C                   of surface reactions.
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SDOTI  - Production rates of the species in reaction IR.
C              1) For K=1,KKGAS, SDOTI(K) is the production rate of
C                 gas-phase species K in (moles/cm**2-sec).
C              2) For K=KKGAS+1,KKGAS+KKSUR, SDOTI(K) is the production
C                 rate of surface species K in (moles/cm**2-sec).
C              3) For K=KKGAS+KKSUR+1,KKTOT, SDOTI(K) is the production
C                 rate of bulk species K in (moles/cm**2-sec).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SDOTI(*) at least KKTOT, the total
C                   number of species.
C    SITDTI  - Production rates of the surface phases due to reaction
C              IR (subroutine calculates entries for the surface site
C              phases only).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SITDTI(*) at least NPHASE, the total
C                   number of phases.
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
      DIMENSION ROP(*), ISKWRK(*), RSKWRK(*), SDOTI(*), SITDTI(*)
      INCLUDE 'skstrt.inc'
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         SDOTI(K) = 0.0
   50 CONTINUE
C
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDTI(N) = 0.0
   60    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
C     PROCESS REACTION IR
C
      IRNU = 0
      DO 75 N = 1, ISKWRK(IiNRNU)
         IF (ISKWRK(ISKWRK(IiIRNU)+N-1) .EQ. IR)
     1       IRNU = ISKWRK(IrRNU)+MAXSPR*(N-1)
   75 CONTINUE
C
      IND = (IR-1) * MAXSPR
      DO 100 N = 1, MAXSPR
        NK = ISKWRK(ISKWRK(IiNUNK) + IND + N - 1)
        IF (NK .NE. 0) THEN
C*****precision > double
           RNU = DBLE (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > double
C*****precision > single
C           RNU = REAL (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > single
            IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
            SDOTI(NK) = SDOTI(NK) + ROP(IR)*RNU
         ENDIF
  100 CONTINUE
C
      NIISUR = ISKWRK(IiNIIS)
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 125 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            RNCF = RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR + IR - 1)
            SITDTI(N) = SITDTI(N) + ROP(IR)*RNCF
  125    CONTINUE
      ENDIF
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRDEX (IR, ISKWRK, RSKWRK, RD)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRDEX (IR, ISKWRK, RSKWRK, RD)
C
C     Returns the perturbation factor of the IRth reaction,
C     or changes its value, depending on the sign of IR.
C
C  INPUT
C      IR     - Integer reaction number; IR>0 gets RD(I) from RSKWRK,
C                                        IR<0 puts RD(I) into RSKWRK.
C                    Data type - integer scalar
C
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      RD     - Perturbation factor for the IRth reaction
C                    Data type - real scalar
C                    cgs units:  none
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      NI = ISKWRK(IrPAR) + (ABS(IR)-1)*(NSPAR+1) + NSPAR
      IF (IR .GT. 0) THEN
         RD = RSKWRK(NI)
      ELSE
         RSKWRK(NI) = RD
      ENDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKREXS (K, ISKWRK, RSKWRK, A7)
C
C  START PROLOGUE
C
C  SUBROUTINE SKREXS (K, ISKWRK, RSKWRK, A7)
C  Created by Matthew A. Christiansen, U. Delaware 2012
C     Returns an array of the seventh thermodynamic polynomial
C     coefficients for a species, or changes their value,
C     depending on the sign of K.
C
C  INPUT
C      K      - Integer species number; K>0 gets A7(*) from RSKWRK,
C                                       K<0 puts A7(*) into RSKWRK.
C                    Data type - integer scalar
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      A7     - The array of the 7th thermodynamic polynomial
C               coefficients for the Kth species, over the number
C               of temperature ranges used in fitting thermodynamic
C               properties.
C               Dimension A7(*) at least (MAXTP-1), where MAXTP is
C               the maximum number of temperatures used for fitting
C               the thermodynamic properties of the species.
C                    Data type - real array
C                    cgs units:  none
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
      DIMENSION ISKWRK(*), RSKWRK(*), A7(*)
      INCLUDE 'skstrt.inc'
C
      DO 100 L = 1, MAXTP-1
         NA7 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (ABS(K)-1)*NCP2T + NCP1
         IF (K .GT. 0) THEN
            A7(L) = RSKWRK(NA7)
         ELSE
            RSKWRK(NA7) = A7(L)
         ENDIF
  100 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRHEX (K, ISKWRK, RSKWRK, A6)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRHEX (K, ISKWRK, RSKWRK, A6)
C
C     Returns an array of the sixth thermodynamic polynomial
C     coefficients for a species, or changes their value,
C     depending on the sign of K.
C
C  INPUT
C      K      - Integer species number; K>0 gets A6(*) from RSKWRK,
C                                       K<0 puts A6(*) into RSKWRK.
C                    Data type - integer scalar
C      RSKWRK - Array of real internal work space.
C                    Data type - real array
C
C  OUTPUT
C      A6     - The array of the 6th thermodynamic polynomial
C               coefficients for the Kth species, over the number
C               of temperature ranges used in fitting thermodynamic
C               properties.
C               Dimension A6(*) at least (MAXTP-1), where MAXTP is
C               the maximum number of temperatures used for fitting
C               the thermodynamic properties of the species.
C                    Data type - real array
C                    cgs units:  none
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
      DIMENSION ISKWRK(*), RSKWRK(*), A6(*)
      INCLUDE 'skstrt.inc'
C
      DO 100 L = 1, MAXTP-1
         NA6 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (ABS(K)-1)*NCP2T + NCP
         IF (K .GT. 0) THEN
            A6(L) = RSKWRK(NA6)
         ELSE
            RSKWRK(NA6) = A6(L)
         ENDIF
  100 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
      SUBROUTINE SKROPFWD  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROPFWD,
     1                      Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKROP  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROP)
C     Returns forward rates of progress for the surface reactions.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     ROPFWD    - Forward Rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec).
C                   Data type - real array
C                   Dimension ROP(*) at least IISUR, the total number of
C                   surface reactions.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), ROPFWD(*), T(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
C      DO 100 I = 1, NIISUR
C         ROP(I) = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
C  100 CONTINUE
      DO 100 I = 1, NIISUR
         ROPFWD(I) = RSKWRK(ISKWRK(IrIT1)+I-1)
  100 CONTINUE
C
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C                                                                      C
      SUBROUTINE SKROP  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROP,Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKROP  (P, T, ACT, SDEN, ISKWRK, RSKWRK, ROP)
C     Returns rates of progress for the surface reactions.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     ROP    - Rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec).
C                   Data type - real array
C                   Dimension ROP(*) at least IISUR, the total number of
C                   surface reactions.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), ROP(*), T(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL SKATCZ (P, T, ACT, SDEN, ISKWRK, RSKWRK,
     1             RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      DO 100 I = 1, NIISUR
         ROP(I) = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
  100 CONTINUE
C      DO 100 I = 1, NIISUR
C         ROP(I) = RSKWRK(ISKWRK(IrIT1)+I-1)
C  100 CONTINUE
C
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRP   (ISKWRK, RSKWRK, RU, RUC, PATM)
C     Returns universal gas constants and the pressure of one standard
C     atmosphere.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real work space.
C                   Data type - real array
C
C  OUTPUT
C     RU     - Universal gas constant.
C                   cgs units - 8.314510E7 ergs/(mole*K)
C                   Data type - real scalar
C     RUC    - Universal gas constant used only in conjuction with
C              activation energy.
C                   preferred units - RU / 4.184 cal/(mole*K)
C                   Data type - real scalar
C     PATM   - Pressure of one standard atmosphere.
C                   cgs units - 1.01325E6 dynes/cm**2
C                   Data type - real scalar
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      RU   = RSKWRK(ISKWRK(IrRU))
      RUC  = RSKWRK(ISKWRK(IrRUC))
      PATM = RSKWRK(ISKWRK(IrPATM))
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKRROP (ISKWRK, RSKWRK, IISUR, SMH, MAXSPR, RU, PATM,
     1                   KKGAS, KKSUR, T, CZ, ACT, WT, NREAC, NRPP, NU,
     2                   NUNK, NRNU, IRNU, RNU, NUSUMK, NSPAR, PAR,
     3                   RPAR, NREV, IREV, NCOV, ICOV, KCOV, NDIM, CPAR,
     4                   NSTK, ISTK, NBOHM, IBOHM, IBK, IBT, SDTOT,
     5                   NORD, MAXORD, IORD, KORD, RORD, RKFT, RKRT,
     6                   RKF, RKR, EQKC, EQFAC, MOTZWS, SKMIN,Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE SKRROP (ISKWRK, RSKWRK, IISUR, SMH, MAXSPR, RU, PATM,
C                     KKGAS, KKSUR, T, CZ, ACT, WT, NREAC, NRPP, NU,
C                     NUNK, NRNU,  IRNU, RNU, NUSUMK, NSPAR, PAR,
C                     RPAR, NREV, IREV, NCOV, ICOV, KCOV, NDIM, CPAR,
C                     NSTK, ISTK, NBOHM, IBOHM, IBK, IBT, SDTOT,
C                     NORD, MAXORD, IORD, KORD, RORD, RKFT, RKRT,
C                     RKF, RKR, EQKC, EQFAC, MOTZWS, SKMIN)
C     Returns forward and reverse rates of progress and equilibrium
C     constants for the surface reactions.
C     It is not normally called by the user application code.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     IISUR  - Total number of surface reactions in the mechanism.
C                   Data type - integer scalar
C     SMH    - Entropies minus enthalpies for the species;
C              SMH(K) = S(K)/R - H(K)/RT.
C                   cgs units none
C                   Data type - real array
C                   Dimension SMH(*) at least KKTOT, the total
C                   number of species.
C     MAXSPR - Maximum number of species in a surface reaction.
C                   Data type - integer scalar
C     RU     - Universal gas constant.
C                   cgs units - 8.314510E7 ergs/(mole*K)
C                   Data type - real scalar
C     PATM   - Pressure of one standard atmosphere.
C                   cgs units - 1.01325E6 dynes/cm**2
C                   Data type - real scalar
C     KKGAS  - Total number of gas-phase species.
C                   Data type - integer scalar
C     KKSUR  - Total number of site-phase species.
C                   Data type - integer scalar
C
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CZ     - Matrix of the concentrations of the gas-phase and
C              surface species in the problem, and the activities
C              of the bulk species.  The first KKGAS entries of CZ
C              are the gas-phase molar concentrations (moles/cm**3).
C              The next KKSURF entries are the surface species molar
C              concentrations (moles/cm**2).  The final KKBULK entries
C              are the activities of the bulk species.
C                   Data type - real array
C                   Dimension CZ(*) at least KKTOT, the total
C                   number of species.
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C     WT     - Array of molecular weights for the species.
C                   Data type - real array
C                   Dimension WT(*) at least KKTOT, the total
C                   number of species
C     NREAC  - Array of the number of reactants in the surface
C                   reactions.
C                   Data type - integer array
C                   Dimension NREAC(*) at least IISUR, the total
C                   number of surface reactions.
C     NRPP   - Integer vector that indicates the reversibility and the
C              total number of species (reactants plus products) of the
C              IISUR surface reactions;
C                   +NRPP = reversible surface reaction IR has NRPP
C                           reactants and products
C                   -NRPP = irreversible reaction IR has ABS(NRPP)
C                           reactants and products
C                   Data type - integer array
C                   Dimension NRPP(*) at least IISUR, the total number
C                   of surface reactions.
C     NU     - Matrix of stoichiometric coefficients for the KKTOT
C              species in the IISUR surface reactions.  NU(M,IR) is
C              the stoichiometric coefficient of the Mth species in
C              reaction IR.  These coefficients are negative for
C              reactants and positive for products.  The species
C              number for the Mth species is stored in NUNK.
C                   Data type - integer array
C                   Dimension NU(*,*) exactly MAXSPR for the first
C                   dimension and at least IISUR for the second, the
C                   total number of surface reactions.
C     NUNK   - Matrix of species numbers corresponding to the
C              stoichiometric coefficients in NU.
C                   Data type - integer array
C                   Dimension NUNK(*,*) exactly MAXSPR for the first
C                   dimension and at least IISUR for the second, the
C                   total number of surface reactions.
C     NUSUMK - The total of the coefficients of the gas-phase species
C              in each surface reaction.
C                   Data type - integer array
C                   Dimension NUSUMK(*) at least IISUR, the total
C                   number of surface reactions.
C     NSPAR  - Number of parameters in the rate expression.  In the
C              current formulation NSPAR=3 (see PAR(N,I) below).
C                   Data type - integer scalar
C     PAR    - Matrix of reaction rate parameters in the form:
C              K = A * T**b * EXP(-E/R*T)
C              1) PAR(1,I) are the pre-exponential factors, A.
C              2) PAR(2,I) are the temperature exponents, b.
C              3) PAR(3,I) are the activation energies, E.
C                   Data type - real array
C                   Dimension PAR(*,*) exactly NSPAR for the first
C                   dimension and at least IISUR for the second, the
C                   total number of surface reactions.
C     RPAR   - Matrix of reverse Arrhenius parameters for the NREV
C              reactions
C                   Data type - real array
C                   Dimension RPAR(*,*) exactly NSPAR for the first
C                   dimension and at least NREV for the second,
C                   the total number of reactions with reverse
C                   Arrhenius parameters defined.
C     NREV   - Number of reactions which have reverse Arrhenius
C              parameters defined.
C                   Data type - integer scalar
C     IREV   - Array of reaction numbers which have reverse Arrhenius
C              parameters defined.
C                   Data type - integer array
C                   Dimension IREV(*) at least NREV, the total number
C                   of reactions with reverse Arrhenius parameters
C                   defined.
C     NCOV   - Total number of site species coverage declarations.
C                   Data type - integer scalar
C     ICOV   - Reaction index numbers for the NCOV coverage
C              declarations.
C                   Data type - integer array
C                   Dimension ICOV(*) at least NCOV, the total number
C                   of coverage declarations.
C     KCOV   - Species index numbers for the NCOV coverage declarations.
C                   Data type - integer array
C                   Dimension KCOV(*) at least NCOV, the total number
C                   of coverage declarations.
C     NDIM   - The first dimension of the matrix of coverage parameters
C              for the NCOV coverage declarations.
C                   Data type - integer scalar
C     CPAR   - Matrix of coverage parameters for the NCOV coverage
C              declarations.
C                   Data type - real array
C                   Dimension CPAR(*,*) exactly NSCOV for the first
C                   dimension, the number of coverage parameters
C                   allowed, and at least NCOV for the second
C                   dimension, the total number of coverage
C                   declarations.
C     NSTK   - Number of reactions which have sticking coefficients.
C                   Data type - integer scalar
C     ISTK   - Array of reaction numbers for the NSTK reactions.
C                   Data type - integer array
C                   Dimension ISTK(*) at least NSTK, the total number
C                   of reactions with sticking coefficients.
C     NBOHM  - Number of reactions which have Bohm velocity conditions.
C     IBOHM  - Array of reaction numbers for the NBOHM reactions
C                   Data type - integer array
C                   Dimension IBOHM(*) at least NBOHM.
C     IBK    - Array of species numbers for Bohm velocity condition
C                   Data type - integer array
C                   Dimension IBK(*) at least NBOHM
C     IBT    - Array of species temperature flaggs for the NBOHM reaxs.
C                   Data type - integer array
C                   Dimension IBT(*) at least NBOHM
C     SDTOT  - The sum of the densities of the phases.
C                   Data type - real scalar
C     EQFAC  - Multiplicative constant applied in the calculation of
C              the equilibrium constants for the surface reactions.
C                   Data type - real array
C                   Dimension EQFAC(*) at least IISUR, the total number
C                   of surface reactions.
C     MOTZWS - integer flag indicating if Motz-Wise correction for large
C              sticking coefficients should be applied (1=yes,0=no).
C                   Data type - integer
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C  OUTPUT
C     RKF    - Forward rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension RKF(*) at least IISUR, the total number
C                   of surface reactions.
C     RKR    - Reverse rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension RKR(*) at least IISUR, the total number
C                   of surface reactions.
C     EQKC   - Equilibrium constants in concentration units for
C              the surface reactions.
C                   Data type - real array
C                   Dimension EQKC(*) at least IISUR, the total number
C                   of surface reactions.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SMH(*), CZ(*), ACT(*), WT(*),
     1          NREAC(*), NRPP(*), NU(MAXSPR,*), NUNK(MAXSPR,*),
     2          NUSUMK(*), PAR(NSPAR+1,*), RPAR(NSPAR+1,*), IREV(*),
     3          ICOV(*), KCOV(*), CPAR(NDIM,*), ISTK(*), RKRT(*),
     4          RKFT(*), RKF(*), RKR(*), EQKC(*), EQFAC(*), IRNU(*),
     5          RNU(MAXSPR,*), IORD(*), KORD(MAXORD,*), RORD(MAXORD,*),
     6          T(*), IBOHM(*), IBK(*), IBT(*), EACT(IISUR)
C
C JM insert
      DIMENSION I_CHK_ORD(iisur)
C end insert
C
      COMMON /MACHN/ SMALL,BIG,EXPARG
C
      INCLUDE 'sknew.inc'
      DATA ZERO, ONE/0.0, 1.0/
C
C
      IF (IISUR .LE. 0) RETURN
C
      ALOGT = LOG(T(1))
      PFAC  = PATM / (RU*T(1))
      CALL SKSMH (T(1), ISKWRK, RSKWRK, SMH)
C
C     This loop was moved up from its original location and modified to
C     allow the free energy of reaction to be used in determining the
C     activation energy.
      DO 100 I = 1, IISUR
         SUMSMH = 0.0
C
         DO 90 N = 1, MAXSPR
            IF (NUNK(N,I) .NE. 0) SUMSMH=SUMSMH+NU(N,I)*SMH(NUNK(N,I))
   90    CONTINUE
C
         EQKC(I) = EXP(MIN(SUMSMH,EXPARG)) * PFAC**NUSUMK(I) * EQFAC(I)
C        This sets the activation energy to the /maximum/ of the free energy
C        of reaction, the specified activation energy, and zero. This is to
C        account for the possibility that activation energies estimated with
C        BEP correlations (SKBEP) may be negative. It also accounts for the
C        possibility that any coverage effects in the species binding energies
C        may result in a reaction energy being more positive than the
C        specified activation energy. In this case, the true activation
C        energy must be /at least/ as high as the reaction energy.
         EACT(I) = MAX(-SUMSMH,PAR(3,I),0.0D0)
  100 CONTINUE
C
C     This loop was moved up from its original location and modified to
C     allow the free energy of reaction to be used in determining the
C     activation energy.
      DO 125 N = 1, NRNU
         SUMSMH = 0.0
         I = IRNU(N)
         RNUSUM = 0.0
C
         DO 110 L = 1, MAXSPR
            IF (NUNK(L,I) .NE. 0) THEN
               SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
               IF (NUNK(L,I) .LE. KKGAS) RNUSUM = RNUSUM + RNU(L,N)
            ENDIF
  110    CONTINUE
C
         PFRNU = EXP (RNUSUM * LOG(PFAC))
         EQKC(I) = EXP(MIN(SUMSMH,EXPARG)) * PFRNU * EQFAC(I)
C        This sets the activation energy to the /maximum/ of the free energy
C        of reaction, the specified activation energy, and zero. This is to
C        account for the possibility that activation energies estimated with
C        BEP correlations (SKBEP) may be negative. It also accounts for the
C        possibility that any coverage effects in the species binding energies
C        may result in a reaction energy being more positive than the
C        specified activation energy. In this case, the true activation
C        energy must be /at least/ as high as the reaction energy.
         EACT(I) = MAX(-SUMSMH,PAR(3,I),0.0D0)
  125 CONTINUE
C
C             PROCESS EACH REACTION
C
C      write(*,*)IREV(1),IREV(2),IREV(3),IREV(4)
C      write(*,*)IREV(5),IREV(6),IREV(7),IREV(8)
C       write(*,*)IREV(9),IREV(10),IREV(11),IREV(12)
C       write(*,*)IREV(13),IREV(14),IREV(15),IREV(16)
C        write(*,*)IREV(17),IREV(18),IREV(19),IREV(20)
C       write(*,*)IREV(21),IREV(22),IREV(23),IREV(24)
C      write(*,*)IREV(25),IREV(26),IREV(27),IREV(28)
C       write(*,*)IREV(29),IREV(30),IREV(31),IREV(32)
C       write(*,*)IREV(33),IREV(34),IREV(35)
C      write(*,*)'IREV'
c.......ashish insert
c.......for modification for (T/To)**b instead of T**b
      itnot=Tref_beta
      if (itnot.eq.0) then
         Tnot=300.0
         alogTnot=log(Tnot)
      else
         Tnot=1.0
         alogTnot=0.0
      endif
C
      DO 20 I = 1, IISUR
C
C
C###################################################################
C     Modifications of Vlachos' scheme start here
C  =================================================
C JM: Note:  This version of sklib.f is only for VL's surface scheme
C     with non-linear coverage-dependent activation energies.
C  PAR(3,I) is the activation energy of reaction I in units of Kelvin
C  Multiply coverage dependent factors (which are in units of Kcal/mol)
C  by 4184.0/8.31451 to make them Kelvin.
C
C
C
C
C
C
C     Modifications of Vlachos' scheme end here
C##################################################################
C
C
         RKF(I) = 0.0
         RKR(I) = 0.0
c JM insert
c.......ashish insert (alogtnot)
C        This has been modified to use the minimum activation energy
C        determined above instead of the original, as-specified value
C        (PAR(3,I)).
         RKFT(I) = PAR(1,I) *
     1    EXP(PAR(2,I)*(ALOGT-alogtnot) - EACT(I)/T(1))
C     1    EXP(PAR(2,I)*(ALOGT-alogtnot) - PAR(3,I)/T(1))
c         RKFT(I) = PAR(1,I) *
c     1             EXP(PAR(2,I)*ALOGT - PAR(3,I)/T(1))
         RKRT(I) = 0.0
C
C JM insert: Initialize to zero a new array I_CHK_ORD(i).  This array
C will check to see which surface reactions have an order change
         I_CHK_ORD(I) = 0
C
   20 CONTINUE
C
c      write(*,*)((par(3,20)+add_e(20))/T(1))*(T(1)*8.31451/4184.)
C
C****************************************
C JM print insert
c     open(899,file='adsor.out',status='unknown')
C*************************************************
C
C JM insert: find those surface reactions that have an order change
      DO 21 I = 1,NORD
         I_OR = IORD(I)
         I_CHK_ORD(I_OR) = I
 21   CONTINUE
C end insert
C
      DO 50 N = 1, NSTK
         I = ISTK(N)
         NSUM = 0
         IRSTO = 0
         RSUM = 0.0
C
C JM insert the exponent m in (Gama(tot)**m) need not be integer.
C Define this exponent as real R_NSUM
C
         R_NSUM = 0.0
C
         DO 40 L = 1, NRNU
            IF (I .EQ. IRNU(L)) IRSTO = L
   40    CONTINUE
C
         DO 25 L = 1, NREAC(I)
            IF (NUNK(L,I).GT.0 .AND. NUNK(L,I).LE.KKGAS) THEN
               WT1 = WT(NUNK(L,I))
            ELSEIF (NUNK(L,I).GT.KKGAS .AND. NUNK(L,I).LE.KKGAS+KKSUR)
     1      THEN
C   JM insert: for adsorption reactions with no order change
C   use previous approach
               IF(I_CHK_ORD(I).EQ.0)THEN
C********************************************************
c        write(899,*) 'Adsorb. React with no order change  ',I
C***********************************************************
                  NSUM = NSUM + ABS(NU(L,I))
                  R_NSUM = DFLOAT(NSUM)
               ELSE
C   JM insert: for adsorption reactions with order change
C***********************************************
c        write(899,*) 'Adsorb. React with order change  ',I
C**********************************************
                  R_NSUM = R_NSUM + RORD(L,I_CHK_ORD(I))
               ENDIF
C   end insert
               IF (IRSTO .GT. 0) RSUM = RSUM + ABS(RNU(L,IRSTO))
            ENDIF
   25    CONTINUE
C
C         PAR1 = PAR(1,I) * 3637.6011 / SQRT(WT1)
C         PAR2 = PAR(2,I) + 0.5
C         RKF(I) = PAR1 * EXP(PAR2*ALOGT - PAR(3,I)/T(1)
C
         IF (SDTOT .EQ. 0.0) THEN
            RKFT(I) = 0.0
         ELSE
C JM print:
c        write(6,*) '***************'
c        write(6,*) 'MOTZWS ',MOTZWS
C
            IF (MOTZWS .EQ. 1) THEN
C************************************************
C JM print option (for test purposes)
c        write(899,*) 'reaction, motz ',I,MOTZWS
c        write(899,*) 'reac, gama,T ',I,RKFT(I),T(1)
c        write(899,*) 'KKGAS, KKSUR', KKGAS, KKSUR
c        write(899,*) 'ifr_site, act ',ifr_site, act(ifr_site)
c        do 7444 i_spec =1,kkgas+kksur
c          write(899,*) 'k, act ',I_SPEC,ACT(I_SPEC)
c7444    continue
C************************************************
C
C a) RKFT(I) is the sticking coefficient of reaction I
C b) 3637.6011 = SQRT(Ru/2*PI)
C
C JM : modify the adsorption rate constant to account
C      for partially occupied sites. ACT(ifr_site) is
C      the coverage of the free site (PT or PD usually)
C      and the index ifr_site is imported from the user's
C      program
C
               RKFT(I) = RKFT(I) *3637.6011 * SQRT(T(1)/WT1)
cash     1                   /(1.0-RKFT(I)/2.0)
c    1                   /( 1.0-0.5*ACT(ifr_site)*RKFT(I) )
            ELSE
               RKFT(I) = RKFT(I) * 3637.6011 * SQRT(T(1)/WT1)
            ENDIF
C
c.......Ashish modification for molbeam
c.......Modifies only adsorption As
c      rkft(i)=rkft(i)*sqrt(298.15/T(1))
C
C****************************************************
C JM  NOTE:  THe order of reactions in Vl's chemical scheme should not be changed
C   ===========================================================================
C
C  A) Include  coverage dependence of the sticking coefficients
C  In new scheme the coverage dependence of sticking coefficients has been dropped
cc
cc              RKFT(I)=RKFT(I) * ACT(ifr_site)
ca              write(6,*)'reaction1',PAR(1,1),PAR(2,1),PAR(3,1)
cc            elseif (I.EQ.3) then
cc              RKFT(I)=RKFT(I) * ACT(ifr_site)**2.0
ca              write(6,*)'reaction3',PAR(1,3),PAR(2,3),PAR(3,3)
cc            elseif (I.EQ.11) then
cc              RKFT(I)=RKFT(I) * ACT(ifr_site)
ca              write(6,*)'reaction11',PAR(1,11),PAR(2,11),PAR(3,11)
cc            endif
C****************************************************
            IF (IRSTO.EQ.0 .AND. R_NSUM.GT.0.0) THEN
C              RKFT(I) = RKFT(I) / (SDTOT**NSUM)
               RKFT(I) = RKFT(I) / (SDTOT**R_NSUM)
C****************************************************
C JM print option for test purposes
c        write(899,*) 'Reac, W, R_NSUM  ',i,wt1,r_nsum
c        write(899,*) '   '
C*****************************************************
            ELSEIF (IRSTO.GT.0 .AND. RSUM.GT.SKMIN) THEN
               write(6,*) '******************** real stoich *'
               SDTR = SDTOT**RSUM
               RKFT(I) = RKFT(I) / SDTR
            ENDIF
         ENDIF
C
   50 CONTINUE
C***************************************************
C JM print option for test purposes
c      close(899)
C****************************************************
C
C JM comment
C now calculate the coverage-dependent reaction rates
c      write(*,*)ncov,kcov(1),kcov(2),'ncovkcov'
c      write(*,*)rkft(icov(1)),rkft(icov(2))
c      EE8=((par(3,8))/T(1))*(T(1)*8.31451/4184.)
c      EE10=((par(3,10))/T(1))*(T(1)*8.31451/4184.)
c      write(*,*)EE8,EE10,'Eold'
      DO 60 N = 1, NCOV
         K = KCOV(N)
         RKFT(ICOV(N)) = RKFT(ICOV(N)) *
     1                  10.0**(CPAR(1,N) * MAX(ZERO,ACT(K))) *
     2                  ABS(ACT(K))**CPAR(2, N) *
     3                  EXP(-CPAR(3, N) * ACT(K) / T(1))
c      write(*,*)CPAR(1,N),CPAR(2,N),CPAR(3,N)
c      write(*,*)rkft(icov(n))
   60 CONTINUE
c      Eadd8=
c     #  (-CPAR(3, 1) * ACT(15) / T(1))*(T(1)*8.31451/4184.)
c      Eadd10=
c     #  (-CPAR(3, 2) * ACT(14) / T(1))*(T(1)*8.31451/4184.)
c      EE8=EE8-Eadd8
c      EE10=EE10-Eadd10
c      write(*,*)Eadd8,Eadd10,'Eadd'
c      write(*,*)EE8,EE10,'Enew'
c      stop
C
C
      DO 70 N = 1, NBOHM
         I = IBOHM(N)
         K = IBK(N)
         TEMP = T(IBT(N))
C
C the Bohm velocity is defined as sqrt(kTe/mi) x gamma
C
         NSUM = 0
         IRSTO = 0
         RSUM = 0.0
C
         DO 63 L = 1, NRNU
            IF (I .EQ. IRNU(L)) IRSTO = L
   63    CONTINUE
C
         DO 65 L = 1, NREAC(I)
            IF (NUNK(L,I).GT.KKGAS .AND. NUNK(L,I).LE.KKGAS+KKSUR)
     1      THEN
               NSUM = NSUM + ABS(NU(L,I))
               IF (IRSTO .GT. 0) RSUM = RSUM + ABS(RNU(L,IRSTO))
            ENDIF
   65    CONTINUE
C
         RKFT(I) = PAR(1,I) * 9117.76 * SQRT(TEMP/WT(K))
         IF (IRSTO.EQ.0 .AND. NSUM.GT.0) THEN
            RKFT(I) = RKFT(I) / (SDTOT**NSUM)
         ELSEIF (IRSTO.GT.0 .AND. RSUM.GT.SKMIN) THEN
            SDTR = SDTOT**RSUM
            RKFT(I) = RKFT(I) / SDTR
         ENDIF
C
   70 CONTINUE
C
C     RKR=0.0 for irreversible reactions, else RKR=RKF/MAX(EQKC,SMALL)
C
      DO 150 I = 1, IISUR
         IF (NRPP(I) .GT. 0) RKRT(I) = RKFT(I) / MAX(EQKC(I),SMALL)
  150 CONTINUE
C
      DO 200 N = 1, NREV
         RKRT(IREV(N)) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT-RPAR(3,N)/T(1))
         EQKC(IREV(N)) = RKFT(IREV(N)) / RKRT(IREV(N))
  200 CONTINUE
C
C JM insert
c     OPEN(488,FILE='rorder.out',STATUS= 'UNKNOWN')
c     write(488,*) '  NRNU  ',NRNU
      DO 250 I = 1, IISUR
C
         RKFT(I) = RKFT(I) * PAR(4,I)
         RKRT(I) = RKRT(I) * PAR(4,I)
C
C JM insert to monitor reaction order
c        write(488,*) '  SURFACE REACTION NUMBER  ',I
c        write(488,*)  '   '
c        write(488,*)  '           SPECIES No.    STOICH. COEFF'
         IF (NU(1,I) .EQ. 0) GO TO 250
C
         RKF(I) = RKFT(I) * CZ(NUNK(1,I))**IABS(NU(1,I))
c JM insert
c      write(488,*) '1st species ',nunk(1,i), nu(1,i)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * CZ(NUNK(2,I))**IABS(NU(2,I))
c JM insert
c      write(488,*) '2nd species ',nunk(2,i), nu(2,i)
            IF (NUNK(3,I) .NE. 0) THEN
               RKF(I) = RKF(I) * CZ(NUNK(3,I))**IABS(NU(3,I))
c JM insert
c      write(488,*) '3rd species ',nunk(3,i), nu(3,i)
               IF (NUNK(4,I) .NE. 0) THEN
                  RKF(I) = RKF(I) * CZ(NUNK(4,I))**IABS(NU(4,I))
c JM insert
c      write(488,*) '4th species ',nunk(4,i), nu(4,i)
                  IF (NUNK(5,I) .NE. 0) THEN
                     RKF(I) = RKF(I) * CZ(NUNK(5,I))**IABS(NU(5,I))
                     IF (NUNK(6,I) .NE. 0) THEN
                        RKF(I) = RKF(I) * CZ(NUNK(6,I))**IABS(NU(6,I))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         RKR(I) = RKRT(I) * CZ(NUNK(7,I))**NU(7,I)
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I) * CZ(NUNK(8,I))**NU(8,I)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I) * CZ(NUNK(9,I))**NU(9,I)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I) * CZ(NUNK(10,I))**NU(10,I)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I) * CZ(NUNK(11,I))**NU(11,I)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I) * CZ(NUNK(12,I))**NU(12,I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  250 CONTINUE
C
C  JM insert
c     write(488,*)
c     write(488,*) ' NOW REAL STOICHIOMETRIC COEFF., NRNU ',NRNU
C
      DO 350 N = 1, NRNU
         I = IRNU(N)
         C1 = CZ(NUNK(1,I))**ABS(RNU(1,N))
         C7 = CZ(NUNK(7,I))**RNU(7,N)
         RKF(I) = RKFT(I) * C1
         RKR(I) = RKRT(I) * C7
C
C JM insert
c     write(488,*) ' 1st species Real ',nunk(1,i),rnu(1,i)
         IF (NUNK(2,I) .NE. 0) THEN
            C2 = CZ(NUNK(2,I))**ABS(RNU(2,N))
C JM insert
c     write(488,*) ' 2nd species Real ',nunk(2,i),rnu(2,i)
            RKF(I) = RKF(I) * C2
            IF (NUNK(3,I) .NE. 0) THEN
               C3 = CZ(NUNK(3,I))**ABS(RNU(3,N))
C JM insert
c     write(488,*) ' 3rd species Real ',nunk(3,i),rnu(3,i)
               RKF(I) = RKF(I) * C3
               IF (NUNK(4,I) .NE. 0) THEN
                  C4 = CZ(NUNK(4,I))**ABS(RNU(4,N))
                  RKF(I) = RKF(I) * C4
C JM insert
c     write(488,*) ' 4th species Real ',nunk(4,i),rnu(4,i)
                  IF (NUNK(5,I) .NE. 0) THEN
                     C5 = CZ(NUNK(5,I))**ABS(RNU(5,N))
                     RKF(I) = RKF(I) * C5
                     IF (NUNK(6,I) .NE. 0) THEN
                        C6 = CZ(NUNK(6,I))**ABS(RNU(6,N))
                         RKF(I) = RKF(I) * C6
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
         IF (NUNK(8,I) .NE. 0) THEN
            C8 = CZ(NUNK(8,I))**RNU(8,N)
            RKR(I) = RKR(I) * C8
            IF (NUNK(9,I) .NE. 0) THEN
               C9 = CZ(NUNK(9,I))**RNU(9,N)
               RKR(I) = RKR(I) * C9
               IF (NUNK(10,I) .NE. 0) THEN
                  C10 = CZ(NUNK(10,I))**RNU(10,N)
                  RKR(I) = RKR(I) * C10
                  IF (NUNK(11,I) .NE. 0) THEN
                     C11 = CZ(NUNK(11,I))**RNU(11,N)
                     RKR(I) = RKR(I) * C11
                     IF (NUNK(12,I) .NE. 0) THEN
                        C12 = CZ(NUNK(12,I))**RNU(12,N)
                         RKR(I) = RKR(I) * C12
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  350 CONTINUE
C
C JM insert
c     write(488,*) '  '
c     write(488,*) '    REACTIONS WITH ORDER, NORD = ',NORD
C
      DO 450 N = 1, NORD
         I = IORD(N)
C JM insert
c        write(488,*) ' i, Reaction No. ',N,I
c        write(488,*) '   '
         RKF(I) = RKFT(I)
         RKR(I) = RKRT(I)
C
         DO 450 L = 1, MAXORD
            K = KORD(L,N)
            ORD = RORD(L,N)
            IF (K .LT. 0) THEN
               RKF(I) = RKF(I) * CZ(-K)**ORD
C JM insert
c     write(488,*) 'spec No., order, conc ',K,ORD,CZ(-K)
            ELSEIF (K .GT. 0) THEN
               RKR(I) = RKR(I) * CZ(K)**ORD
            ENDIF
  450 CONTINUE
C
C JM insert
c     close(488)
C
      DO 550 N = 1, NBOHM
         I = IBOHM(N)
         RKF(I) = RKFT(I) * CZ(IBK(N))
         DO 483 NS = 1, 6
            IF (NUNK(NS,I).NE.IBK(N) .AND. NUNK(NS,I).NE.0) THEN
               IF (NUNK(NS,I) .GT. KKGAS) THEN
                  RKF(I) = RKF(I)*CZ(NUNK(NS,I))**IABS(NU(NS,I))
               ENDIF
            ENDIF
 483     CONTINUE
         RKR(I) = 0.0
  550 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSAVE (LOUT, LSAVE, ISKWRK, RSKWRK, CSKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSAVE (LOUT, LSAVE, ISKWRK, RSKWRK, CSKWRK)
C     Writes to a binary file information about a Surface Chemkin
C     binary file, pointers for the Surface Chemkin Library, and
C     Surface Chemkin work arrays.
C
C  INPUT
C     LOUT   - Output file for printed diagnostics.
C                   Data type - integer scalar
C     LSAVE  - Integer output unit.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace containing integer data.
C                   Data type - integer array
C     RSKWRK - Array of real workspace containing real data.
C                   Data type - real array
C     CSKWRK - Array of character workspace containing character data.
C                   Data type - CHARACTER*16 array
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*)
      INCLUDE 'skstrt.inc'
C
      CHARACTER*(*) CSKWRK(*)
      CHARACTER*16 VERS, PREC
      LOGICAL KERR
      COMMON /SKCONS/ VERS, PREC, KERR
C
      LENI = ISKWRK(IiLENI)
      LENR = ISKWRK(IiLENR)
      LENC = ISKWRK(IiLENC)
      WRITE (LSAVE, ERR=999) VERS, PREC, LENI, LENR, LENC,
C
C     Integer constants
C
     1   MAXSPR, NELEM, NKKGAS, NSPAR, NSCOV, MAXORD, MAXTP, NCP,
     2   NCP1, NCP2, NCP2T,
C
C     ISKWRK pointers to integer variables
C
     3   IiLENI, IiLENR, IiLENC, IiKSUR, IiKBLK, IiKTOT, IiNPHA, IiFSUR,
     4   IiLSUR, IiNSUR, IiFBLK, IiLBLK, IiNBLK, IiNIIS, IiNCOV, IiNREV,
     5   IiNSTK, IiNCON, IiNBHM, IiNRNU, IiNORD, IiMOTZ,
C
C     ISKWRK pointers to integer arrays
C
     6   IiPKST, IiPKND, IiPTOT, IiKPHS, IiKCHG, IiKCMP, IiNSCV, IiKNT,
     7   IiNRPP, IiNREA, IiNUNK, IiNU,   IiNSUM, IiICOV, IiKCOV, IiIREV,
     8   IiISTK, IiIBHM, IiKBHM, IiTBHM, IiIRNU, IiIORD, IiKORD,
C
C     ISKWRK pointers to real variables
C
     *   IrSKMN, IrPATM, IrRU,   IrRUC,
C
C     ISKWRK pointers to real arrays
C
     1   IrSDEN, IrKTMP, IrKTHM, IrKDEN, IrAWT,  IrKWT,  IrPAR, IrKCOV,
     2   IrRPAR, IrEQ,   IrRNU,  IrNCF,  IrKORD, IrKFT,  IrKRT, IrKT1,
     3   IrKT2,  IrPT1,  IrIT1,  IrIT2,  IrIT3,
C
C     ISKWRK pointers to character arrays
C
     4   IcENAM, IcKNAM, IcMNAM, IcPNAM
C
C     END include file for sklib.f
C
      WRITE (LSAVE, ERR=999) (ISKWRK(L), L = 1, LENI)
      WRITE (LSAVE, ERR=999) (RSKWRK(L), L = 1, LENR)
      WRITE (LSAVE, ERR=999) (CSKWRK(L), L = 1, LENC)
C
      RETURN
  999 CONTINUE
      WRITE (LOUT, *)
     1 ' Error writing Surface binary file information...'
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSDEN (ISKWRK, RSKWRK, SDEN0)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSDEN (ISKWRK, RSKWRK, SDEN0)
C     Returns a real array of standard-state phase densities as given
C     on input to the interpreter.
C
C  INPUT
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SDEN0  - Standard-state densities for the surface site types,
C              AS READ BY THE INTERPRETER.  The SDEN0 vector has
C              an entry for each phase, including the gas phase,
C              but this subroutine only writes into those entries
C              related to the surface phases,
C              i.e. NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN0(*) at least NPHASE, the total
C                   number of phases.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SDEN0(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNSUR) .LE. 0) RETURN
C
      DO 100 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
         SDEN0(N) = RSKWRK(ISKWRK(IrSDEN) + N - 1)
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSMH  (T, ISKWRK, RSKWRK, SMH)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSMH  (T, ISKWRK, RSKWRK, SMH)
C     Returns the array of dimensionless entropies minus enthalpies
C     for the species.  It is normally not called directly by the user.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SMH    - Dimensionless entropies minus enthalpies for the
C              species;  SMH(K) = S(K)/R - H(K)/RT.
C                   cgs units none
C                   Data type - real array
C                   Dimension SMH(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SMH(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = LOG(T) - 1.0
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)/((N-1)*N)
 100  CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         SMH(K) = SUM + RSKWRK(NA1 + NCP2 - 1)
     1                - RSKWRK(NA1 + NCP1 - 1)/T
 250  CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSML  (T, ISKWRK, RSKWRK, SML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSML  (T, ISKWRK, RSKWRK, SML)
C     Returns an array of the standard state entropies in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SML    - Standard state entropies in molar units for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension SML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SML(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = LOG(T)
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)/(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 140 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  140    CONTINUE
         SML(K) = RSKWRK(ISKWRK(IrRU)) *
     1            (SUM + RSKWRK(NA1 + NCP2 - 1))
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSMS  (T, ISKWRK, RSKWRK, SMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSMS  (T, ISKWRK, RSKWRK, SMS)
C     Returns an array of the standard state entropies in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SMS    - Standard state entropies in mass units for the species.
C                   cgs units - ergs/(gm*K)
C                   Data type - real array
C                   Dimension SMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SMS(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = LOG(T)
      DO 100 N = 2, NCP
         TN(N)=T**(N-1)/(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         SMS(K) = RSKWRK(ISKWRK(IrRU)) *
     1            (SUM + RSKWRK(NA1 + NCP2 - 1)) /
     2            RSKWRK(ISKWRK(IrKWT) + K - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSNUM (LINE, NEXP, LOUT, KNAM, KKTOT, PNAM, NNPHAS,
     1                   KKPHAS, KNUM, NT, NVAL, RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSNUM (LINE, NEXP, LOUT, KNAM, KKTOT, PNAM, NNPHAS,
C                     KKPHAS, KNUM, NT, NVAL, RVAL, KERR)
C     This subroutine is used to read a format-free input line of
C     combined alphanumerical data.  It can be used to parse an input
C     character string, LINE, which may be composed of several blank-
C     delimited substrings.  This subroutine assumes that the first
C     substring in LINE is the name of a species in the Surface Chemkin
C     mechanism.  If the species name is not unique within the Surface
C     Chemkin mechanism, the phase of the species should be input
C     immediately after the species name, delimited by slashes.
C     Upon return from the subroutine, KNUM returns the index position
C     of the species within the Surface Chemkin binary file.  If the
C     species name is not unique, KNUM returns the first position and
C     NT returns the number of the times the species occurs within the
C     binary file.  If the species name is not found, or there is a
C     syntax error, on return, KNUM=0, NT=0, and KERR=.TRUE.
C     The substrings in LINE following the first are expected to
C     represent numbers.  They are converted into floating point
C     values and stored in the output vector, RVAL(*).  Upon input,
C     NEXP is equal to the number of values expected to be found.
C     If NEXP numbers are not found, KERR will be set to .TRUE. on
C     return from the subroutine.
C
C     Example input:
C             LINE     = GA(S)/BULK1/ 1.2
C             NEXP     = 1, the number of values expected
C             LOUT     = 6, a logical unit number on which to write
C                        diagnostic messages
C             KNAM(*)  = Array of character species names
C             KKTOT    = Total number of species
C             PNAM(*)  = Array of character phase names
C             NNPHAS   = Total number of phases
C             KKPHAS(*)= Index array of the number of species in the
C                        phases
C     Output:
C             KNUM     = The index number of the species which
C                        has the name "GA(S)" and resides in phase
C                        "BULK1"
C             NT       = 1, if there is only one species GA(S)
C                        in phase BULK1
C             NVAL     = 1, the number of values found in LINE
C                        following the species name
C             RVAL(1)  = 1.200E+00, the substring converted to a
C                        real number
C             KERR     = .FALSE.
C
C  INPUT
C     LINE   - A character string.
C                   Data type - CHARACTER*80
C     NEXP   - Number of real values expected to be found in
C              character string;  if NEXP < 0, no error message
C              printed to output.
C                   Data type - integer scalar
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C     KNAM   - Array of species names.
C                   Data type - character array
C                   Dimension KNAM(*) at least KKTOT, the total
C                   number of species.
C     KKTOT  - Total number of species.
C                   Data type - integer scalar
C     PNAM   - Array of phase names.
C                   Data type - character array
C                   Dimension PNAM(*) at least NNPHAS, the total
C                   number of phases.
C     NNPHAS - Total number of phases.
C                   Data type - integer scalar
C     KKPHAS - Array of the number of species in the phases.
C                   Data type - integer array
C                   Dimension KKPHAS(*) at least NNPHAS, the
C                   total number of phases.
C  OUTPUT
C     KNUM   - Index number of species which corresponds to the
C              species name in LINE.
C                   Data type - integer scalar
C     NT     - Number of times the species name occurs in
C              the binary file.
C                   Data type - integer scalar
C     NVAL   - Number of real values found in LINE
C                   Data type - integer scalar
C     RVAL   - Array of real values found in LINE
C                   Data type - real array
C                   Dimension RVAL(*) at least NEXP
C     KERR   - Error flag; syntax or dimensioning error,
C              corresponding species not found, or total number of
C              values found is not the number of values expected,
C              will result in KERR = .TRUE.
C                   Data type - logical
C
C  END PROLOG
C     A !' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RVAL(*), KKPHAS(*)
      CHARACTER*(*) KNAM(*), PNAM(*), LINE
      CHARACTER*80 ISTR
      LOGICAL KERR
C
      KNUM = 0
      NVAL = 0
      KERR = .FALSE.
      NT   = 0
C
      ILEN = MIN (IPPLEN(LINE), ILASCH(LINE))
      IF (ILEN .LE. 0) RETURN
C
      IF (INDEX(LINE,'/') .GT. 0) THEN
         I1 = INDEX (LINE, '/')
         I2 = I1 + INDEX(LINE(I1+1:), '/')
      ELSE
         I2 = IFIRCH(LINE) + INDEX(LINE(IFIRCH(LINE):),' ') - 2
      ENDIF
C
      I1 = IFIRCH(LINE)
      IF (I2 .GE. I1) CALL SKPCMP (LINE(I1:I2), KNAM, KKTOT, PNAM,
     1                             NNPHAS, KKPHAS, KNUM, NT)
C
      IF (KNUM .EQ. 0) THEN
         ISTR = ' '
         ISTR = LINE(:ILEN)
         IF (NEXP .GT. 0) THEN
            WRITE (LOUT, '(A)')
     1      ' Error in SKSNUM...'//ISTR(:ILEN)//' not found...'
            KERR = .TRUE.
         ENDIF
      ENDIF
      IF (NEXP .NE. 0)
     1   CALL CKXNUM (LINE(I2+1:), NEXP, LOUT, NVAL, RVAL, KERR)
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSOR  (T, ISKWRK, RSKWRK, SOR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSOR  (T, ISKWRK, RSKWRK, SOR)
C     Returns an array of the nondimensional entropies.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     SOR    - Nondimensional entropies for the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension SOR(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), SOR(*), TN(10)
      INCLUDE 'skstrt.inc'
C
      TN(1) = LOG(T)
      DO 100 N = 2, NCP
         TN(N) = T**(N-1)/(N-1)
100   CONTINUE
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMP = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         SOR(K) = SUM + RSKWRK(NA1 + NCP2 - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYME (ISKWRK, CSKWRK, LOUT, ENAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYME (ISKWRK, CSKWRK, LOUT, ENAM, KERR)
C     Returns a character array of element names.
C
C  INPUT
C     CSKWRK - Array of character workspace.
C                   Data type - CHARACTER*16 array
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     ENAM   - Element names.
C                   Data type - CHARACTER*(*) array
C                   Dimension ENAM(*) at least NELEM,  the total
C                   number of elements in the problem.
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      DIMENSION ISKWRK(*)
      CHARACTER*(*) ENAM(*), CSKWRK(*)
      LOGICAL KERR
      INCLUDE 'skstrt.inc'
C
      KERR = .FALSE.
      ILEN = LEN(ENAM(1))
      DO 100 N = 1, NELEM
         ENAM(N) = ' '
         LT = ILASCH(CSKWRK(ISKWRK(IcENAM) + N - 1))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
         ELSE
            ENAM(N) = CSKWRK(ISKWRK(IcENAM) + N - 1)
         ENDIF
  100 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMM (ISKWRK, CSKWRK, LOUT, MATNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMM (ISKWRK, CSKWRK, LOUT, PNAM, KERR)
C     Returns the character name of a material.
C
C  INPUT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     CSKWRK - Array of character workspace.
C                   Data type - CHARACTER*16 array
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     MATNAM - Material name.
C                   Data type - CHARACTER*(*)
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      DIMENSION ISKWRK(*)
      CHARACTER*(*) CSKWRK(*), MATNAM
      LOGICAL KERR
      INCLUDE 'skstrt.inc'
C
      KERR = .FALSE.
      ILEN = LEN(MATNAM)
      MATNAM = ' '
      LT = ILASCH(CSKWRK(ISKWRK(IcMNAM)))
      IF (LT .GT. ILEN) THEN
         KERR = .TRUE.
         WRITE (LOUT,*)
     1   ' Error in SKSYMM...character string length too small '
      ELSE
         MATNAM = CSKWRK(ISKWRK(IcMNAM))
      ENDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMP (ISKWRK, CSKWRK, LOUT, PNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMP (ISKWRK, CSKWRK, LOUT, PNAM, KERR)
C     Returns a character array of phase names.
C
C  INPUT
C     CSKWRK - Array of character workspace.
C                   Data type - CHARACTER*16 array
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     PNAM   - Phase names.
C                   Data type - CHARACTER*(*) array
C                   Dimension PNAM(*) at least NNPHAS, the total
C                   number of sites.
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      DIMENSION ISKWRK(*)
      CHARACTER*(*) PNAM(*), CSKWRK(*)
      LOGICAL KERR
      INCLUDE 'skstrt.inc'
C
      KERR = .FALSE.
      ILEN = LEN(PNAM(1))
      DO 100 N = 1, ISKWRK(IiNPHA)
         PNAM(N) = ' '
         LT = ILASCH(CSKWRK(ISKWRK(IcPNAM) + N - 1))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
            WRITE (LOUT,*)
     1      ' Error in SKSYMP...character string length too small '
         ELSE
            PNAM(N) = CSKWRK(ISKWRK(IcPNAM) + N - 1)
         ENDIF
  100 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMR (IR, LOUT, ISKWRK, RSKWRK, CSKWRK, LT, RNAM,
     1                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMR (IR, LOUT, ISKWRK, RSKWRK, CSKWRK, LT, RNAM,
C                     KERR)
C     Returns the character string representation of reaction IR.
C
C  INPUT
C     IR     - Reaction index
C                   Data type - integer scalar
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     CSKWRK - Array of character workspace.
C                   Data type - CHARACTER*16 array
C
C  OUTPUT
C     LT     - Total number of non-blank characters in the
C              reaction string.
C                   Data type - integer scalar
C     RNAM   - Reaction string.
C                   Data type - CHARACTER*(*)
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      DIMENSION ISKWRK(*), RSKWRK(*)
      CHARACTER*(*) RNAM, CSKWRK(*)
      CHARACTER*1 IST(10)
      LOGICAL KERR
      INCLUDE 'skstrt.inc'
C
      DATA IST/'0','1','2','3','4','5','6','7','8','9'/
C
      RNAM = ' '
      KERR = .FALSE.
      LT = 0
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      ILEN = LEN(RNAM)
C
      IRNU = 0
      DO 10 N = 1, ISKWRK(IiNRNU)
         IF (IR .EQ. ISKWRK(ISKWRK(IiIRNU) + N - 1)) IRNU = N
   10 CONTINUE
C
      SKMIN = RSKWRK(ISKWRK(IrSKMN))
      DO 100 J = 1, 2
         NS = 0
         DO 50 N = 1, MAXSPR
            NU = ISKWRK(ISKWRK(IiNU)   + (IR-1)*MAXSPR + N - 1)
            K  = ISKWRK(ISKWRK(IiNUNK) + (IR-1)*MAXSPR + N - 1)
C
            IF (IRNU .GT. 0) THEN
               RNU = RSKWRK(ISKWRK(IrRNU) + (IRNU-1)*MAXSPR + N - 1)
               IF (ABS(RNU) .GT. SKMIN) THEN
                  IF (J .EQ. 1) THEN
                     NU = -1
                  ELSE
                     NU = 1
                  ENDIF
               ENDIF
            ENDIF
            IF (J.EQ.1.AND.NU.LT.0 .OR. J.EQ.2.AND.NU.GT.0) THEN
               NS = NS + 1
C
               IF (NS .GT. 1) THEN
                  LL = ILASCH(RNAM)
                  IF (LL+1 .GT. ILEN) THEN
                     KERR = .TRUE.
                     RNAM = ' '
                     WRITE (LOUT, 500)
                     RETURN
                  ENDIF
                  RNAM(LL+1:) = '+'
               ENDIF
               IF (ABS(NU) .GT. 1) THEN
                  LL = ILASCH(RNAM)
                  IF (LL+1 .GT. ILEN) THEN
                     KERR = .TRUE.
                     RNAM = ' '
                     WRITE (LOUT, 500)
                     RETURN
                  ENDIF
                  RNAM(LL+1:) = IST(ABS(NU)+1)
               ENDIF
               LK = ILASCH(CSKWRK(ISKWRK(IcKNAM) + K - 1))
               LL = ILASCH(RNAM)
               IF (LL+LK .GT. ILEN) THEN
                  KERR = .TRUE.
                  RNAM = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               RNAM(LL+1:) = CSKWRK(ISKWRK(IcKNAM)+K-1)(:LK)
C
C              DOES SPECIES NAME OCCUR MORE THAN ONCE?
C
               CALL SKCOMP (CSKWRK(ISKWRK(IcKNAM)+K-1),
     1                      CSKWRK(ISKWRK(IcKNAM)),
     1                      ISKWRK(IiKTOT), KNUM, NK)
               IF (NK .GT. 1) THEN
C
C                 IN WHICH PHASE DOES SPECIES K OCCUR?
C
                  DO 40 NP = 1, ISKWRK(IiNPHA)
                     KFIRST = ISKWRK(ISKWRK(IiPKST) + NP - 1)
                     KLAST  = ISKWRK(ISKWRK(IiPKND) + NP - 1)
                     IF (K.GE.KFIRST .AND. K.LE.KLAST) THEN
                        LK = ILASCH(CSKWRK(ISKWRK(IcPNAM) + NP - 1))
                        LL = ILASCH(RNAM)
                        IF (LL+LK+2 .GT. ILEN) THEN
                           KERR = .TRUE.
                           RNAM = ' '
                           WRITE (LOUT, 500)
                           RETURN
                        ENDIF
                        RNAM(LL+1:) =
     1                  '/'//CSKWRK(ISKWRK(IcPNAM) + NP - 1)(:LK)//'/'
                     ENDIF
   40             CONTINUE
               ENDIF
            ENDIF
   50    CONTINUE
C
         IF (J .EQ. 1) THEN
            LL = ILASCH(RNAM)
            IF (ISKWRK(ISKWRK(IiNRPP)+IR-1) .LT. 0) THEN
               IF (LL+2 .GT. ILEN) THEN
                  KERR = .TRUE.
                  RNAM = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               RNAM(LL+1:) = '=>'
            ELSE
               IF (LL+3 .GT. ILEN) THEN
                  KERR = .TRUE.
                  RNAM = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               RNAM(LL+1:) = '<=>'
            ENDIF
         ENDIF
  100 CONTINUE
C
      LT = ILASCH(RNAM)
  500 FORMAT (' Error in SKSYMR...character string length too small ')
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKSYMS (ISKWRK, CSKWRK, LOUT, KNAM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE SKSYMS (ISKWRK, CSKWRK, LOUT, KNAM, KERR)
C     Returns a character array of species names.
C
C  INPUT
C     CSKWRK - Array of character workspace.
C                   Data type - CHARACTER*16 array
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     KNAM   - Species names.
C                   Data type - CHARACTER*(*) array
C                   Dimension KNAM(*) at least KKTOT, the total
C                   number of species.
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      DIMENSION ISKWRK(*)
      CHARACTER*(*) CSKWRK(*), KNAM(*)
      LOGICAL KERR
      INCLUDE 'skstrt.inc'
C
      KERR = .FALSE.
      ILEN = LEN(KNAM(1))
      DO 100 N = 1, ISKWRK(IiKTOT)
         KNAM(N) = ' '
         LT = ILASCH(CSKWRK(ISKWRK(IcKNAM) + N - 1))
         IF (LT .GT. ILEN) THEN
            KERR = .TRUE.
            WRITE (LOUT,500) CSKWRK(ISKWRK(IcKNAM) + N - 1)
         ELSE
            KNAM(N) = CSKWRK(ISKWRK(IcKNAM) + N - 1)
         ENDIF
  100 CONTINUE
C
  500 FORMAT (' Error in SKSYMS...character string length too small ')
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKUML  (T, ISKWRK, RSKWRK, UML)
C
C  START PROLOGUE
C
C  SUBROUTINE SKUML  (T, ISKWRK, RSKWRK, UML)
C     Returns an array of the internal energies in molar units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     UML    - Internal energies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C                   Dimension UML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), UML(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKHML (T, ISKWRK, RSKWRK, UML)
      RUT = T * RSKWRK(ISKWRK(IrRU))
      DO 100 K = 1, ISKWRK(IiKTOT)
         UML(K) = UML(K) - RUT
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKUMS  (T, ISKWRK, RSKWRK, UMS)
C
C  START PROLOGUE
C
C  SUBROUTINE SKUMS  (T, ISKWRK, RSKWRK, UMS)
C     Returns an array of the internal energies in mass units.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     UMS    - Internal energies in mass units for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension UMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), UMS(*)
      INCLUDE 'skstrt.inc'
C
      CALL SKHMS (T, ISKWRK, RSKWRK, UMS)
      RUT = T * RSKWRK(ISKWRK(IrRU))
      DO 100 K = 1, ISKWRK(IiKTOT)
         UMS(K) = UMS(K) - RUT/RSKWRK(ISKWRK(IrKWT) + K - 1)
100   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE SKWT   (ISKWRK, RSKWRK, WT)
C
C  START PROLOGUE
C
C  SUBROUTINE SKWT   (ISKWRK, RSKWRK, WT)
C     Returns the molecular weights of the species.
C
C  INPUT
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     WT     - Molecular masses for the species.
C                   cgs units - gm/mole
C                   Data type - real array
C                   Dimension WT(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), WT(*)
      INCLUDE 'skstrt.inc'
C
      DO 100 N = 1, ISKWRK(IiKTOT)
         WT(N) = RSKWRK(ISKWRK(IrKWT) + N - 1)
  100 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, CZ)
C     Returns the concentrations of the species, given the pressure,
C     temperature and activities.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array.
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array.
C                   Data type - integer array.
C                   Dimension KTFL(*) at least KKGAS
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     CZ     - Matrix of the concentrations of the gas-phase and
C              surface species in the problem, and the activities
C              of the bulk species.  The first KKGAS entries of CZ
C              are the gas-phase molar concentrations (moles/cm**3).
C              The next KKSURF entries are the surface species molar
C              concentrations (moles/cm**2).  The final KKBULK entries
C              are the activities of the bulk species.
C                   Data type - real array
C                   Dimension CZ(*) at least KKTOT, the total
C                   number of gas-phase + surface + bulk species.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), CZ(*), T(*),
     1          KTFL(*)
      INCLUDE 'skstrt.inc'
C
      DATA ONE/1.0/
C
C     COMPUTE THE MOLAR CONCENTRATIONS (C) FROM MOLE FRACTIONS
C
      SUMXT = 0.0
      DO 50 K = 1, NKKGAS
         SUMXT = SUMXT + ACT(K)*T(KTFL(K))
 50   CONTINUE
      PRUT = P/(RSKWRK(ISKWRK(IrRU))*SUMXT)
      DO 100 K = 1, NKKGAS
         CZ(K) = ACT(K)*PRUT
  100 CONTINUE
C
      IF (ISKWRK(IiKSUR) .GT. 0) THEN
         KCOV = ISKWRK(IiNSCV)
         DO 210 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            KFIRST = ISKWRK(ISKWRK(IiPKST) + N - 1)
            KLAST  = ISKWRK(ISKWRK(IiPKND) + N - 1)
            DO 210 K = KFIRST, KLAST
C*****precision > double
               CZ(K) = ACT(K) * SDEN(N) / DBLE(ISKWRK(KCOV + K - 1))
C*****END precision > double
C*****precision > single
C               CZ(K) = ACT(K) * SDEN(N) / REAL(ISKWRK(KCOV + K - 1))
C*****END precision > single
  210    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiKBLK) .LE. 0) RETURN
      KFIRST = ISKWRK(ISKWRK(IiPKST) + ISKWRK(IiFBLK) - 1)
      KLAST  = ISKWRK(ISKWRK(IiPKND) + ISKWRK(IiLBLK) - 1)
      DO 220 K = KFIRST, KLAST
         CZ(K) = ACT(K)
  220 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKDRDA (IR, P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1                    DKDAI, Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKDRDA (IR, P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
C 1                    DKDAI)
C     Returns the partial of the rates of production for each of the
C     species with respect to the pre-exponential constant of surface
C     reaction IR.
C
C  INPUT
C     IR     - Reaction index
C                   Data type - integer scalar
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array.
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array.
C                   Data type - integer array.
C                   Dimension KTFL(*) at least NKKGAS.
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     DKDAI  - Array of the partial of the production rates of the
C              species with respect to the pre-exponential
C              constant for reaction IR.
C                   cgs units - moles/(cm**2*sec) / (units of A)
C                   Data type - real array
C                   Dimension DKDAI(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), DKDAI(*), T(*),
     1          KTFL(*)
      LOGICAL IFLAG
      INCLUDE 'skstrt.inc'
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDAI(K) = 0.0
   50 CONTINUE
C
      CALL PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1              RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
C     PROCESS REACTION IR
C
      CALL SKRAEX  (IR, ISKWRK, RSKWRK, RA)
C
C     FIND OUT IF REVERSE ARRHENIUS PARAMETERS WERE SPECFIED FOR
C     THIS REACTION
      IFLAG = .TRUE.
      DO 70 N = 1, NIIREV
         IF ( ISKWRK(ISKWRK(IiIREV)+N-1). EQ. IR) IFLAG = .FALSE.
   70 CONTINUE
C
      IRNU = 0
      DO 75 N = 1, NIIRNU
         IF (ISKWRK(ISKWRK(IiIRNU)+N-1) .EQ. IR)
     1       IRNU = ISKWRK(IrRNU)+MAXSPR*(N-1)
   75 CONTINUE
C
      RKF = RSKWRK(ISKWRK(IrIT1) + IR - 1)
      RKR = RSKWRK(ISKWRK(IrIT2) + IR - 1)
      IND = (IR-1) * MAXSPR
      DO 80 N = 1, MAXSPR
         NK = ISKWRK(ISKWRK(IiNUNK) + IND + N - 1)
         IF (NK .NE. 0) THEN
C*****precision > double
            RNU = DBLE (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > double
C*****precision > single
C           RNU = REAL (ISKWRK(ISKWRK(IiNU) + IND + N - 1))
C*****END precision > single
            IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
            IF (IFLAG) THEN
C              REVERSE PARAMETERS WERE NOT GIVEN FOR THIS REACTION
               DKDAI(NK) = DKDAI(NK) - RNU*RKR / RA
            ENDIF
            DKDAI(NK) = DKDAI(NK) + RNU*RKF / RA
         ENDIF
   80 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKDRDC (KSPEC, P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1                    DKDC, Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKDRDC (KSPEC, P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
C 1                   DKDC)
C     Returns the partial derivative of the production rates for
C     each of the species with respect to the concentration of species
C     KSPEC.
C
C  INPUT
C     KSPEC  - Species index
C                   Data type - integer scalar
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array.
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array.
C                   Data type - integer array.
C                   Dimension KTFL(*) at least NKKGAS
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     DKDC   - Array of the partial of the production rates of the
C              species with respect to the concentration
C              of species KSPEC.
C                   cgs units - moles/(cm**2*sec) / (units of KSPEC)
C                   Data type - real array
C                   Dimension DKDC(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ACT(*), ISKWRK(*), RSKWRK(*), SDEN(*), DKDC(*), T(*),
     1          KTFL(*)
      INCLUDE 'skstrt.inc'
      LOGICAL IFLAG
      DATA TEN/10.0/
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      DO 50 K = 1, ISKWRK(IiKTOT)
         DKDC(K) = 0.0
   50 CONTINUE
C
      CALL PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1              RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      CZK = RSKWRK(ISKWRK(IrKT1) + KSPEC - 1)
C
C     PROCESS EACH REACTION
C
      DO 1000 I = 1, NIISUR
C
C        FIND OUT IF REVERSE ARRHENIUS PARAMETERS WERE SPECFIED FOR
C        THIS REACTION
         IFLAG = .TRUE.
         DO 70 N = 1, NIIREV
            IF ( ISKWRK(ISKWRK(IiIREV)+N-1). EQ. I) IFLAG = .FALSE.
   70    CONTINUE
         IRNU = 0
         DO 75 N = 1, NIIRNU
            IF (ISKWRK(ISKWRK(IiIRNU) + N - 1) .EQ. I)
     1          IRNU = ISKWRK(IrRNU)+MAXSPR*(N-1)
   75    CONTINUE
C
         FFACT = 0.0
         RFACT = 0.0
         INK = ISKWRK(IiNUNK) + (I-1)*MAXSPR
         INU = ISKWRK(IiNU)   + (I-1)*MAXSPR
         DO 350 N = 1, MAXSPR
            NK = ISKWRK(INK + N - 1)
            IF (NK .EQ. KSPEC) THEN
C*****precision > double
               RNU = DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C               RNU = REAL (ISKWRK(INU + N - 1))
C*****END precision > single
               IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
               IF (RNU .LT. 0) THEN
                   FFACT = -RNU / CZK
               ELSE
                  RFACT =  RNU / CZK
               ENDIF
            ENDIF
  350    CONTINUE
C
         DO 200 N = 1, NIICOV
            IF ( ISKWRK(ISKWRK(IiICOV)+N-1). EQ. I) THEN
C              COVERAGE PARAMETERS WERE SPECIFIED FOR THIS REACTION
               IF ( ISKWRK(ISKWRK(IiKCOV)+N-1). EQ. KSPEC) THEN
C                 SPECIES NUMBER KSPEC MODIFIES THE RATE EXPRESSION
C                 FOR REACTION I.
                  KI = ISKWRK(IrKCOV) + (N-1)*NSCOV
                  CPAR1 = RSKWRK(KI)
                  CPAR2 = RSKWRK(KI+1)
                  CPAR3 = RSKWRK(KI+2)
C                 ADD THE COV. MODIFICATION FOR FORWARD DIRECTION
                  CMOD = CPAR2/ACT(KSPEC) + LOG(TEN)*CPAR1 - CPAR3/T(1)
C                 SCALE BY THE NUMBER OF SITES COVERED BY THIS SPECIES
C                 DIVIDED BY SITE DENSITY OF THIS SURFACE PHASE.
C                 THIS NUMBER TURNS OUT TO BE THE SITE FRACTION DIVIDED
C                 BY THE MOLAR CONCENTRATION FOR THE SURFACE SPECIES.
                  CMOD = CMOD * ACT(KSPEC) / CZK
                  FFACT = FFACT + CMOD
                  IF (IFLAG) THEN
C                    REVERSE ARRHENIUS PARAMETERS WERE NOT GIVEN FOR
C                    THIS REACTION, SO IT IS OK TO ADD-IN THE
C                    COVERAGE MODIFICATION TO THE REVERSE DIRECTION
                     RFACT = RFACT + CMOD
                  ENDIF
               ENDIF
            ENDIF
  200    CONTINUE
C
         RKF = RSKWRK(ISKWRK(IrIT1) + I - 1)
         RKR = RSKWRK(ISKWRK(IrIT2) + I - 1)
         DO 400 N = 1, MAXSPR
            NK = ISKWRK(INK + N - 1)
            IF (NK .NE. 0) THEN
C*****precision > double
               RNU = DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C               RNU = REAL (ISKWRK(INU + N - 1))
C*****END precision > single
               IF (IRNU .GT. 0) RNU = RSKWRK(IRNU + N - 1)
C
               DKDC(NK) = DKDC(NK) + RNU*(RKF*FFACT - RKR*RFACT)
            ENDIF
  400    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKEEM (TW, EW, ECHRG, BOLTZ, PI, EPS0, AVAG,
     1                   PHIM, WTE, SDOTE)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  this subroutine returns the surface emission rate of electrons
C  given surface temperature and electric field, and the work function
C
      PARAMETER (PLANCK = 6.6256E-27)
C
      EMASS = WTE / AVAG
C
C   INCLUDE SHOTTKY REDUCTION OF EMISSION DUE TO E-FIELD
C
      DELPHI = SQRT (ECHRG * ABS(EW) / (4*PI*EPS0))
      IF (EW .LT. 0.0) THEN
         PHIW = PHIM - DELPHI
      ELSE
         PHIW = PHIM
      ENDIF
C
      ARGEX = ECHRG * PHIW / (BOLTZ * TW)
      XNUE = (4 * PI * EMASS * (BOLTZ*TW)**2 / PLANCK**3)
     1        * EXP(-ARGEX)
      SDOTE = XNUE / AVAG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKEQ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, EQKC,
     1                  Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKEQ (P, T, KTFL,  ACT, SDEN, ISKWRK, RSKWRK, EQKC)
C     Returns the equilibrium constants for the surface reactions
C     given pressure, temperature, species activities, and the site
C     densities.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array.
C                   Data type - integer array
C                   Dimension KTFL(*) at least NKKGAS.
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ISKWRK(*) at least LENISK.
C     RSKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RSKWRK(*) at least LENRSK.
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     EQKC   - Equilibrium constants in concentration units
C              for the reactions.
C                   cgs units - (moles,cm), depends on reaction
C                   Data type - real array
C                   Dimension EQKC(*) at least IISUR, the total
C                   number of surface reactions.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*), RSKWRK(*), ACT(*), SDEN(*), EQKC(*), T(*),
     1          KTFL(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1              RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNIIS) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
c
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      DO 10 I = 1, NIISUR
         EQKC(I) = RSKWRK(ISKWRK(IrIT3) + I - 1)
   10 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKHML  (T, ITFL, ISKWRK, RSKWRK, HML)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKHML  (T, ITFL, ISKWRK, RSKWRK, HML)
C     Returns an array of the enthalpies in molar units.
C
C  INPUT
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real scalar
C     ITFL   - Temperature flag array
C                   Data type -integer array
C                   Dimension ITFL(*) at least KKTOT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C
C  OUTPUT
C     HML    - Enthalpies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C                   Dimension HML(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), HML(*), TN(10), T(*), ITFL(*)
      INCLUDE 'skstrt.inc'
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         IF (K .LE. NKKGAS) THEN
            TEMP = T(ITFL(K))
         ELSE
            TEMP = T(1)
         ENDIF
         RUT = TEMP*RSKWRK(ISKWRK(IrRU))
         TN(1) = 1.0
         DO 100 N = 2, NCP
            TN(N) = TEMP**(N-1)/N
100      CONTINUE
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMPL = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (TEMP .GT. TEMPL) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         HML(K) = RUT*(SUM + RSKWRK(NA1 + NCP1 - 1)/TEMP)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKHMS  (T, ITFL, ISKWRK, RSKWRK, HMS)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKHMS  (T, ITFL, ISKWRK, RSKWRK, HMS)
C     Returns an array of the enthalpies in mass units.
C
C  INPUT
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real scalar
C                   Dimension T(*) at least no. of energy eqs.
C     ITFL   - Temperature flag array
C                   Data type -integer array
C                   Dimension ITFL(*) at least KKTOT
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C  OUTPUT
C     HMS    - Enthalpies in mass units for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension HMS(*) at least KKTOT, the total
C                   number of species.
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
      DIMENSION ISKWRK(*), RSKWRK(*), HMS(*), TN(10), T(*), ITFL(*)
      INCLUDE 'skstrt.inc'
C
      DO 250 K = 1, ISKWRK(IiKTOT)
         TEMP = T(ITFL(K))
         RUT = TEMP*RSKWRK(ISKWRK(IrRU))
         TN(1) = 1.0
         DO 100 N = 2, NCP
            TN(N) = TEMP**(N-1)/N
100      CONTINUE
C
         L = 1
         DO 220 N = 2, ISKWRK(ISKWRK(IiKNT) + K - 1)-1
            TEMPL = RSKWRK(ISKWRK(IrKTMP) + (K-1)*MAXTP + N - 1)
            IF (TEMP .GT. TEMPL) L = L+1
 220     CONTINUE
C
         NA1 = ISKWRK(IrKTHM) + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RSKWRK(NA1 + N - 1)
  240    CONTINUE
         HMS(K) = RUT * (SUM + RSKWRK(NA1 + NCP1 - 1)/TEMP)
     1                / RSKWRK(ISKWRK(IrKWT) + K - 1)
250   CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKRAT (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, SDOT,
     1                   SITDOT, Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKRAT  (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, SDOT,
C 1                   SITDOT)
C     Returns production rates for the species and sites.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array
C                   Data type - integer array
C                   Dimension KTFL(*) at least NKK.
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     SDOT   - Production rates of the species.
C              1) For K=1,KKGAS, SDOT(K) is the production rate of
C                 gas-phase species K in (moles/cm**2-sec).
C              2) For K=KKGAS+1,KKGAS+KKSUR, SDOT(K) is the production
C                 rate of surface species K in (moles/cm**2-sec).
C              3) For K=KKGAS+KKSUR+1,KKTOT, SDOT(K) is the production
C                 rate of bulk species K in (moles/cm**2-sec).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SDOT(*) at least KKTOT, the total
C                   number of species.
C    SITDOT  - Production rates of the surface phases (subroutine
C              calculates entries for the surface site phases only).
C                   cgs units - moles/(cm**2*sec)
C                   Data type - real array
C                   Dimension SITDOT(*) at least NPHASE, the total
C                   number of phases.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), SDOT(*),
     1          SITDOT(*), T(*), KTFL(*)
      INCLUDE 'skstrt.inc'
C
      DO 10 K = 1, ISKWRK(IiKTOT)
         SDOT(K) = 0.0
   10 CONTINUE
C
      SDTOT = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 20 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SITDOT(N) = 0.0
            SDTOT = SDTOT + SDEN(N)
   20    CONTINUE
      ENDIF
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1              RSKWRK(ISKWRK(IrKT1)))
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
C        PROCESS EACH REACTION
C
      DO 100 I = 1, NIISUR
C
         ROP = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
         INK = ISKWRK(IiNUNK) + (I-1) * MAXSPR
         INU = ISKWRK(IiNU)   + (I-1) * MAXSPR
C
         DO 50 N = 1, MAXSPR
            NUNK = ISKWRK(INK + N - 1)
            IF (NUNK .NE. 0) SDOT(NUNK) = SDOT(NUNK) + ROP *
C*****precision > double
     1         DBLE (ISKWRK(INU + N - 1))
C*****END precision > double
C*****precision > single
C     1         REAL (ISKWRK(INU + N - 1))
C*****END precision > single
50       CONTINUE
C
         IF (NKKSUR .GT. 0) THEN
            DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
                RNCF = RSKWRK(ISKWRK(IrNCF) + (N-1)*NIISUR + I - 1)
               SITDOT(N) = SITDOT(N) + ROP * RNCF
   60       CONTINUE
         ENDIF
100   CONTINUE
C
      IF (NIIRNU .LE. 0) RETURN
C
      DO 200 L = 1, NIIRNU
C
         I = ISKWRK(ISKWRK(IiIRNU) + L - 1)
         ROP = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
         INK = ISKWRK(IiNUNK) + (I-1) * MAXSPR
         INU = ISKWRK(IrRNU)  + (L-1) * MAXSPR
C
         DO 200 N = 1, MAXSPR
            NUNK = ISKWRK(INK + N - 1)
            IF (NUNK .NE. 0) SDOT(NUNK) = SDOT(NUNK) + ROP *
     1          RSKWRK(INU + N - 1)
  200 CONTINUE
C
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PSKROP (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, ROP,
     1                   Tref_beta)
C
C  START PROLOGUE
C
C  SUBROUTINE PSKROP  (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK, ROP)
C     Returns rates of progress for the surface reactions.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature array.
C                   cgs units - K
C                   Data type - real array
C                   Dimension T(*) for number of energy eqs.
C     KTFL   - Integer species temperature flag array.
C                   Data type - integer array
C                   Dimension KTFL(*) at least NKKGAS.
C     ACT    - Activities of the species, where
C
C              for the first KKGAS species, ACT(*) are mole fractions
C                   cgs units - none
C
C              for the next KKSURF species, ACT(*) are site fractions,
C              (species density normalized by the site density).
C              The surface concentration in moles/cm**2 is:
C              ACT(K)*SITE_DENSITY / # sites per species
C                   cgs units - none
C
C              for the next KKBULK species, ACT(*) are the bulk species
C              activities.  A species in a bulk phase will have an
C              activity from 0 to 1 and the sum of activities for
C              a bulk phase should be 1.
C                   cgs units - none
C                   Data type - real array
C
C                   Dimension ACT(*) at least KKTOT, the total number
C                   of species.
C
C     SDEN   - Site densities for the surface site types.  This
C              vector may have an entry for each phase, including the
C              gas phase, but the subroutine only uses entries for the
C              surface site phases, NFSURF .LE. N .LE. NLSURF.
C                   cgs units - moles/cm**2
C                   Data type - real array
C                   Dimension SDEN(*) at least NPHASE, the total
C                   number of phases.
C
C     ISKWRK - Array of integer workspace.
C                   Data type - integer array
C     RSKWRK - Array of real workspace.
C                   Data type - real array
C     Tref_beta - flag for reference temperature used in rate constant
C                 calculations; 0: Tref=300K; 1: Tref=1K
C
C  OUTPUT
C     ROP    - Rates of progress for the surface reactions.
C                   cgs units - moles/(cm**2*sec).
C                   Data type - real array
C                   Dimension ROP(*) at least IISUR, the total number of
C                   surface reactions.
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
      DIMENSION ACT(*), SDEN(*), ISKWRK(*), RSKWRK(*), ROP(*), T(*),
     1          KTFL(*)
      INCLUDE 'skstrt.inc'
C
      IF (ISKWRK(IiNIIS) .LE. 0) RETURN
C
      CALL PSKATCZ (P, T, KTFL, ACT, SDEN, ISKWRK, RSKWRK,
     1              RSKWRK(ISKWRK(IrKT1)))
C
      SDTOT  = 0.0
      IF (ISKWRK(IiNSUR) .GT. 0) THEN
         DO 60 N = ISKWRK(IiFSUR), ISKWRK(IiLSUR)
            SDTOT = SDTOT + SDEN(N)
   60    CONTINUE
      ENDIF
C
      NIISUR = ISKWRK(IiNIIS)
      RU     = RSKWRK(ISKWRK(IrRU))
      PA     = RSKWRK(ISKWRK(IrPATM))
      NIIRNU = ISKWRK(IiNRNU)
      NKKSUR = ISKWRK(IiKSUR)
      NIIREV = ISKWRK(IiNREV)
      NIICOV = ISKWRK(IiNCOV)
      NIISTK = ISKWRK(IiNSTK)
      NIIBHM = ISKWRK(IiNBHM)
      NIIORD = ISKWRK(IiNORD)
C
      CALL SKRROP
     1(ISKWRK, RSKWRK, NIISUR, RSKWRK(ISKWRK(IrKT2)), MAXSPR, RU, PA,
     2 NKKGAS, NKKSUR, T, RSKWRK(ISKWRK(IrKT1)), ACT,
     3 RSKWRK(ISKWRK(IrKWT)), ISKWRK(ISKWRK(IiNREA)),
     4 ISKWRK(ISKWRK(IiNRPP)), ISKWRK(ISKWRK(IiNU)),
     5 ISKWRK(ISKWRK(IiNUNK)), NIIRNU, ISKWRK(ISKWRK(IiIRNU)),
     6 RSKWRK(ISKWRK(IrRNU)), ISKWRK(ISKWRK(IiNSUM)), NSPAR,
     7 RSKWRK(ISKWRK(IrPAR)), RSKWRK(ISKWRK(IrRPAR)), NIIREV,
     8 ISKWRK(ISKWRK(IiIREV)), NIICOV, ISKWRK(ISKWRK(IiICOV)),
     9 ISKWRK(ISKWRK(IiKCOV)), NSCOV, RSKWRK(ISKWRK(IrKCOV)), NIISTK,
     * ISKWRK(ISKWRK(IiISTK)), NIIBHM, ISKWRK(ISKWRK(IiIBHM)),
     1 ISKWRK(ISKWRK(IiKBHM)), ISKWRK(ISKWRK(IiTBHM)), SDTOT, NIIORD,
     2 MAXORD, ISKWRK(ISKWRK(IiIORD)), ISKWRK(ISKWRK(IiKORD)),
     3 RSKWRK(ISKWRK(IrKORD)), RSKWRK(ISKWRK(IrKFT)),
     4 RSKWRK(ISKWRK(IrKRT)), RSKWRK(ISKWRK(IrIT1)),
     5 RSKWRK(ISKWRK(IrIT2)), RSKWRK(ISKWRK(IrIT3)),
     6 RSKWRK(ISKWRK(IrEQ)), ISKWRK(IiMOTZ), RSKWRK(ISKWRK(IrSKMN)),
     7 Tref_beta )
C
      DO 100 I = 1, NIISUR
         ROP(I) = RSKWRK(ISKWRK(IrIT1)+I-1) - RSKWRK(ISKWRK(IrIT2)+I-1)
  100 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SKCOVERAGE (ISKWRK,RSKWRK,cov_matrix,cov_matrix2,
     1  cov_matrix3, thresh_cov,thresh_cov2,A6,DH,EA,w,omega,COVFACT)
C This routine changes the thermodynamic fits (A6) based on coverage
C also changes Ea by 1/2 of the change to del H
C Mike Salciccioli Jan 2009
C Updated to split cov_matrix into four pieces: (1) is the interaction
C parameters (still named cov_matrix), (2) is the A6 constants from the
C NASA polynomials, (3) is the base heat of reaction (including StatpQ and
C LSR adjustments), and (4) is the base activation energy.
C J.E. Sutton 2013/09/12
C Updated to incorporate two parameter model initially developed by A.
C Mironenko.
C J.E. Sutton 2014/08/06
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*),RSKWRK(*),cov_matrix(ISKWRK(IiKSUR),*),w(*),
     1      DELH(ISKWRK(IiNIIS)),NUMBERSPEC(12),
     2      NSTOICHCOEF(12),HORT(ISKWRK(IiKTOT)),
     3      omega(ISKWRK(IiNIIS)),A6(ISKWRK(IiKTOT),*),DH(*),EA(*),
     4      COVFACT(*), cov_matrix2(ISKWRK(IiKSUR),*),
     5      cov_matrix3(ISKWRK(IiKSUR),*), thresh_cov(ISKWRK(IiKSUR),*),
     6      thresh_cov2(ISKWRK(IiKSUR),*)
      CHARACTER*80 PNAM(3), KNAM(ISKWRK(IiKTOT))
      CHARACTER*80 ISTR
      LOGICAL KERR
C
      INCLUDE 'skstrt.inc'
C
C     Loop over all surface species
      M=ISKWRK(IiKSUR)
      DO 100 K=1,M
         B=0
         DO 110 L=1,M
C           Calculate the effect of coverage -- applicable only if the
C           coverage exceeds the threshold
C          IF (ABS(cov_matrix(L,K)) .GT. 0.1) THEN
C          write(*,*) thresh_cov(L,K), thresh_cov2(L,K)
C          END IF
            IF (w(NKKGAS+L) .LE. thresh_cov(L,K)) THEN
				B=-cov_matrix(L,K)*w(NKKGAS+L)+B
				
            ELSE IF ((w(NKKGAS+L) .GT. thresh_cov(L,K)) .AND.
     1              (w(NKKGAS+L) .LE. thresh_cov2(L,K))) THEN
	            B=B-cov_matrix(L,K)*thresh_cov(L,K)
				B=B-cov_matrix2(L,K)*(w(NKKGAS+L)-thresh_cov(L,K))	
				
            ELSE IF (w(NKKGAS+L) .GT. thresh_cov2(L,K)) THEN
                B=B-cov_matrix(L,K)*thresh_cov(L,K)
				B=B-cov_matrix2(L,K)*(thresh_cov2(L,K)-thresh_cov(L,K))
                B=B-cov_matrix3(L,K)*(w(NKKGAS+L)-thresh_cov2(L,K))		
	        END IF

  110    CONTINUE
C
C        Find the indices of the thermodynamic coefficients
         IL=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+6-1
         IH=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+13-1
C
C        Find the row in the A6 array where the original
C        thermodynamic coefficients are located
         JK=NKKGAS+K
C
C        Calculate and assign the new A6 values using the coverage effects
C        The COVFACT array handles scaling the coverage parameters as a
C        function of binding energy.
         RSKWRK(IL)=A6(JK,1)+B*COVFACT(K)+C*COVFACT(K)
         RSKWRK(IH)=A6(JK,2)+B*COVFACT(K)+C*COVFACT(K)
  100 CONTINUE
C
      T=w(ISKWRK(IiKTOT))
C
C     Calculate the heats of reaction
      Call SKHORT(T,ISKWRK,RSKWRK,HORT)
C     Convert H/RT to H/R
      HORT=HORT*T
      DO 95 Irxn=1,ISKWRK(IiNIIS)
         DELH(Irxn)=0
         Call SKINU(Irxn,12,ISKWRK,RSKWRK,
     2              NUMPART,NUMBERSPEC,NSTOICHCOEF)
         DO 95 Ispec=1,NUMPART
            DELH(Irxn)=DELH(Irxn)
     2                  +HORT(NUMBERSPEC(Ispec))*NSTOICHCOEF(Ispec)
   95 CONTINUE
C
C  change Ea to reflect change in DHrxn due to temp and coverage
C      DO 105 IIrxn=1,ISKWRK(IiNIIS)
C         EaOLD=EA(IIrxn)
C         dhOLD=DH(IIrxn)
C         dhNew=DELH(IIrxn)
C         EaNEW=EaOLD+omega(IIrxn)*(dhNEW-dhOLD)
C         IND=ISKWRK(IrPAR)+(IIrxn-1)*(NSPAR+1)
C         RSKWRK(IND+2)=EaNew
C  105 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE SKSTATPQ(ISKWRK,RSKWRK,T,StatpQ)
C This subroutine is responsible for applying the StatpQ correction.
C J.E. Sutton 2013/09/12
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*),RSKWRK(*),StatpQ(ISKWRK(IiKSUR))
C
      INCLUDE 'skstrt.inc'
C
C---------------------
C
C StatpQ correction for Temp dependent BE
      M=ISKWRK(IiKSUR)
      DO 110 K=1,M
         IL=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+6-1
         RSKWRK(IL)=RSKWRK(IL)+(StatpQ(K)*(T-298))
         IH=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+13-1
         RSKWRK(IH)=RSKWRK(IH)+(StatpQ(K)*(T-298))
  110 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE SKGETA6 (ISKWRK,RSKWRK,A6)
C This routine stores the original temperature fits and originial delH
C Mike Salciccioli Jan 2009
C Updated to remove cov_matrix and use only the A6 coefficients. It also
C no longer handles the StatpQ correction.
C J.E. Sutton 2013/09/12
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*),RSKWRK(*),A6(ISKWRK(IiKTOT),*)
C
      INCLUDE 'skstrt.inc'
C
C     Store the original A6 coefficients
      DO 150 I=1,ISKWRK(IiKTOT)
C
C        Find the indices of the thermodynamic coefficients
         IL=ISKWRK(IrKTHM)+(I-1)*NCP2T+6-1
         IH=ISKWRK(IrKTHM)+(I-1)*NCP2T+13-1
C
C        Assign the A6 values
         A6(I,1)=RSKWRK(IL)
         A6(I,2)=RSKWRK(IH)
  150 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE SKGETA7 (ISKWRK,RSKWRK,A7)
C This routine stores the original A7 coefficients.
C J.E. Sutton 2013/09/12
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ISKWRK(*),RSKWRK(*),A7(ISKWRK(IiKTOT),*)
C
      INCLUDE 'skstrt.inc'
C
C     Store the original A7 coefficients
      DO 150 I=1,ISKWRK(IiKTOT)
C
C        Find the indices of the thermodynamic coefficients
         IL=ISKWRK(IrKTHM)+(I-1)*NCP2T+7-1
         IH=ISKWRK(IrKTHM)+(I-1)*NCP2T+14-1
C
C        Assign the A7 values
         A7(I,1)=RSKWRK(IL)
         A7(I,2)=RSKWRK(IH)
  150 CONTINUE
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE SKBEP(ISKWRK,RSKWRK,T,nBEP,BEP_def,IBEP_rxn,EA_old)
C This routing hijacks the Ea parameters and replaces them with values
C calculated from BEPs
C Mike Salciccioli June 2010
C
C This code has been overhauled to only consider reactions of the
C form EA=m*DH+b. Other forms (e.g., ETS=m*EFS+b) are not handled.
C Code overhaul performed by J.E. Sutton, January 2013
C Code updated to once again consider ETS=m*EFS/EIS+b where energies are
C defined by heats of formation rather than binding energies.
C Code updates by J.E. Sutton, February 2013
C
C INPUT:
C     T:        The temperature (K)
C     nBEP:     The number of correlations
C     BEP_def:  Real array containing slopes & intercepts of
C               correlations of length nBEP
C     IBEP_rxn: Integer array containing reaction directions and
C               pointers to correlation types of length IiNIIS, the
C               number of surface reactions
C     EA_old:   Real array containing the original activation energies.
C               This has been added to allow efficient perturbation of
C               the BEP estimates during the global sensitivity
C               analysis. The original activation energies are stored
C               in cov_matrix, but these values do not have to be used.
C
C OUTPUT:
C     None, but activation energies in the real work array are altered.
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
      DIMENSION ISKWRK(*),RSKWRK(*),BEP_def(nBEP,*),
     1      DELH(ISKWRK(IiNIIS)),
     2      IBEP_rxn(ISKWRK(IiNIIS),*),
     3      NUMBERSPEC(12),NSTOICHCOEF(12),
     4      HORT(ISKWRK(IiKTOT)), HIS(ISKWRK(IiNIIS)),
     5      HFS(ISKWRK(IiNIIS)), EA_old(*)
C
      INCLUDE 'skstrt.inc'
C
C   Get HIS, HFS, Hrxn for all reactions at temperature and convert to K
      Call SKHORT(T,ISKWRK,RSKWRK,HORT)
      HORT = HORT * T
C
      DO 15 Irxn=1,ISKWRK(IiNIIS)
         DELH(Irxn)=0.
         HIS(Irxn)=0.
         HFS(Irxn)=0.
         Call SKINU(Irxn,12,ISKWRK,RSKWRK,
     1              NUMPART,NUMBERSPEC,NSTOICHCOEF)
         DO 20 Ispec=1,NUMPART
            KI=NUMBERSPEC(Ispec)
            NU=NSTOICHCOEF(Ispec)
            IF (NU<0) THEN
               HIS(Irxn)=HIS(Irxn)+HORT(KI)*ABS(NU)
            ELSE
               HFS(Irxn)=HFS(Irxn)+HORT(KI)*NU
            END IF
            DELH(Irxn)=DELH(Irxn)+HORT(KI)*NU
C
   20 CONTINUE
   15 CONTINUE
C
C   This is the main loop for the program
C   Loop over all rxns and calculate the Ea based on the criteria in BEP.inp
      DO 25 Ir=1,ISKWRK(IiNIIS)
C                initialize variables for calculations
         Ea_new=0.
         SLOPE=0.
         B=0.
C
C   If no correlation, skip rxn and use surf.inp
         IF ((IBEP_rxn(Ir,1)) .NE. 0) THEN
C
            SLOPE=BEP_def(IBEP_rxn(Ir,1),2)
            B=BEP_def(IBEP_rxn(Ir,1),3)
C
            IF (INT(BEP_def(IBEP_rxn(Ir,1),1)) .EQ. 0) THEN
C              Standard BEP (EA)=m*(DH)+b
C              Check to see if the BEP reference direction matches the
C              direction of the reaction. If not, then modify the slope
C              to be 1-m instead of m. The intercept remains unchanged.
C              This is proved in Sutton & Vlachos, ACS Catal. 2(8), 1624
               IF (INT(BEP_def(IBEP_rxn(Ir,1),4))
     1            .NE. IBEP_rxn(Ir,2)) THEN
                  SLOPE=1-SLOPE
               END IF
               Ea_new=SLOPE*DELH(Ir)+B
            ELSE
C              TSS type correlation (ETS)=m*(ES)+b
C              Need to define the IS/FS and calculate EA
C              Possible cases:
C                 1. Decomp. ref., FS as ind. var., decomp. direction
C                 2. Decomp. ref., IS as ind. var., decomp. direction
C                 3. Decomp. ref., FS as ind. var., synth. direction
C                 4. Decomp. ref., IS as ind. var., synth. direction
C                 5. Synth. ref., FS as ind. var., synth. direction
C                 6. Synth. ref., IS as ind. var., synth. direction
C                 7. Synth. ref., FS as ind. var., decomp. direction
C                 8. Synth. ref., IS as ind. var., decomp. direction
C              First check for IS/FS as ind. var. Then check if reaction and
C              reference directions match. If they don't then switch the
C              definition of the IS & TS.
               IF (INT(BEP_def(IBEP_rxn(Ir,1),1)) .EQ. 1) THEN
C                 FS TSS correlation
                  IF (INT(BEP_def(IBEP_rxn(Ir,1),4))
     1               .EQ. IBEP_rxn(Ir,2)) THEN
                     RIVAR=HFS(Ir)
                  ELSE
                     RIVAR=HIS(Ir)
                  END IF
               ELSE
C                 IS TSS correlation
                  IF (INT(BEP_def(IBEP_rxn(Ir,1),4))
     1               .EQ. IBEP_rxn(Ir,2)) THEN
                     RIVAR=HIS(Ir)
                  ELSE
                     RIVAR=HFS(Ir)
                  END IF
               END IF
               Ea_new=SLOPE*RIVAR+B-HIS(Ir)
            END IF
C
C           Assign the BEP estimate to the work array
            IND=ISKWRK(IrPAR)+(Ir-1)*(NSPAR+1)
C           Check for negative barriers
            IF (Ea_New+EA_old(Ir) .LT. 0) THEN
C              Reset to small value
               RSKWRK(IND+2)=0.02
            ELSE
C              Use the corrected values
               RSKWRK(IND+2)=EA_old(Ir)+Ea_New
            ENDIF
C
         END IF
   25 CONTINUE
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE SKSCALE (ISKWRK,RSKWRK,NSCALE,NATOM,ILEVEL,MODE,
     1                    SLOPE,BEREF,BETARG,COVFACT,COVADJ)
C This subroutine will calculate the change in heat of formation when
C moving from one metal to another using the linear scaling relations
C and then apply it to the thermodynamic correlations for the surface
C species.
C
C
C Inputs are:
C   ISKWRK: CHEMKIN integer work array.
C   RSKWRK: CHEMKIN real work array.
C   NSCALE: The number of scaling relations (integer).
C   NATOM: Total number of atoms used (integer).
C   ILEVEL: What types of species should be used with the scaling relations.
C         0 = none, 1 = all but atomic, 2 = all, 3 = map
C   MODE: An integer array of binding modes. Dimensions of this array should be
C         (NSCALE+1,ksmax), where ksmax is the total number of surface species.
C   SLOPE: A real array of size (NSCALE,2) containing pointers to the atomic 
C         species ID numbers in BETARG in the first column and the corresponding
C         slopes in the second column.
C         for the correlations.
C   BEREF: A real array of size ksmax containing the binding energies of all
C          species on the reference metal.
C   BETARG: A real array of length (NATOM,2) containing the atomic species ID
C         numbers in the first column and the binding energies of the target
C         metal in the second column. The binding energies are defined to be
C         positive for exothermic adsorption, and the values supplied should be
C         dimensionless (i.e., normalized by RT).
C   COVADJ: Logical flag indicating whether to calculate adjustments to the
C         adsorbate interaction parameters.
C
C Output:
C   No argument is returned, but the A6 coefficients in the thermodynamic
C     correlations are modified to reflect the new heats of formation adjusted
C     via the linear scaling relations.
C   COVFACT: An array of length ksmax containing multiplication factors to
C         approximately account for the variation of coverage interaction
C         parameters as a function of atomic binding energy.
C
C Written by J. Sutton 2012/01/11
C Rewritten by J. Sutton 2013/12/19
C
C  END PROLOGUE
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      LOGICAL COVADJ
      DIMENSION ISKWRK(*),RSKWRK(*),MODE(NSCALE+1,*),SLOPE(NSCALE,*),
     1      BEREF(*), BETARG(NATOM,*), DHORT(ISKWRK(IiKTOT)), COVFACT(*)
C
      INCLUDE 'skstrt.inc'
C
      IF (ILEVEL==0) RETURN
C
C     Initialize the change to zero. This will account for the gas species,
C     vacant surface sites, and bulk species.
      DHORT=0.0
C
C     The mode value for atomic species should be zero. If a species is not
C     to be used with the scaling relations, the first value should be zero.
C     In this calculation and then again for the atomic species, we use the
C     order BEREF-BETARG rather than BETARG-BEREF as given in the LSRs to
C     account for the positive sign associated with the reference energy
C     definitions in Scale.inp, which is opposite the convention used in
C     developing the LSRs.
      DO 100 ISPEC=1,ISKWRK(IiKSUR)
         DO 100 ISCALE=1,NSCALE
            ITMP=INT(SLOPE(ISCALE,1))
            IATOM=INT(BETARG(ITMP,1))-NKKGAS
            DHORT(ISPEC+NKKGAS)=DHORT(ISPEC+NKKGAS)+
     1         (BEREF(IATOM)-BETARG(ITMP,2))*
     2         SLOPE(ISCALE,2)*MODE(ISCALE+1,ISPEC)*MODE(1,ISPEC)
  100 CONTINUE
C
C     Account for atomic species.
      IF (ILEVEL>1) THEN
         DO 110 ISPEC=1,NATOM
            IATOM=INT(BETARG(ISPEC,1))
            DHORT(IATOM)=BEREF(IATOM-NKKGAS)-BETARG(ISPEC,2)
  110 CONTINUE
      END IF
C
C     Now that we have the dimensionless changes in the heats of formation,
C     we need to apply the changes to the real work array. This is done by
C     adding the value of DH/R to the A6 coefficient in the NASA polynomial.
      M=ISKWRK(IiKSUR)
      DO 120 K=1,M
         IL=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+6-1
         IH=ISKWRK(IrKTHM)+(NKKGAS-1+K)*NCP2T+13-1
         RSKWRK(IL)=RSKWRK(IL)+DHORT(K+NKKGAS)
         RSKWRK(IH)=RSKWRK(IH)+DHORT(K+NKKGAS)
  120 CONTINUE
      IF (ILEVEL==2 .OR. .NOT. COVADJ) RETURN
C
C     Next we calculate multiplication factors for the adsorbate interaction
C     parameters. This uses the reference energies in BEREF as a normalization
C     factor together with DHORT. This only applies to binding energy maps.
      DO 130 ISPEC=1,ISKWRK(IiKSUR)
        COVFACT(ISPEC)=DHORT(ISPEC+NKKGAS)/BEREF(ISPEC)+1.0
  130 CONTINUE
      RETURN
      END
