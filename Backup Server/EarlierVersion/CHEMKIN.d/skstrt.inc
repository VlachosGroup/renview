C
C     Include file for sklib.f
C
      COMMON /SKSTRT/
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
