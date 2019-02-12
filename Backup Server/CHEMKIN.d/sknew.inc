C***********************************************************************
C
C  sknew.inc - file
C
C***********************************************************************
C  PURPOSE:  extra information for SHEMKIN
C-----------------------------------------------------------------------
C  >>>>>>>   grouped fields
C
C  INTEGER ifr_site     the species number that corresponds to the surface
C                       free site PT(S) or PD(S)
C***********************************************************************

      INTEGER ifr_site,i_motze

      COMMON/ISKFR/ifr_site
      COMMON/MOTZE/i_motze
     

c***********************************************************************
c  common block cmn.files for program BIGMIX
c***********************************************************************
c CHEMKIN II - file names
c----------------------------------------------------------------------
c  thrmdata   : name of CHEMKIN-DATABASE (with path)
c  mechnam    : path and name of mechanism - file
c  outpnam    : path and name of output - file
c  linkck_nam : path and name of link - file for the
c               gas phase chemistry
c  linksk_nam : path and name of link - file for the
c               surface chemistry
c----------------------------------------------------------------------
 
      CHARACTER*50 thrmdata,mechnam,outpnam,linksk_nam,linkck_nam
             
      COMMON/FILES/thrmdata,mechnam,outpnam,linksk_nam,linkck_nam

c*****************  end of common block cmn.files
