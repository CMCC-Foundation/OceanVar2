!======================================================================
!
! This file is part of Oceanvar.
!
!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version (GPL-3.0-or-later).
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!======================================================================
!-----------------------------------------------------------------------
!                                                                      !
!> Apply  quality  control                                            
!!
!! eve 999 obs directly assimilated as it is                            
!! eve 998 obs directly assimilated after the merging with 997 flag obs 
!! eve 997 obs indirectly assimilated: merged in 998 flag obs           
!! eve  1  obs excluded after reading                                   
!! eve  2  obs outside model domain                                     
!! eve  3  obs not sufficient land point for interpolation              
!! eve  4  obs with grid point negative index (should not happen)       
!! eve  5  obs large residuals                                          
!! eve  6  obs background check                                         
!! eve  7  obs eof background check                                     
!! eve  8  obs climatological check obs - clim                          
!! eve  9  obs costal rejection                                         
!! eve  10 obs sla too close to equator if D.H.                         
!! eve  11 obs vert check                                               
!! eve  12 obs lvl of no motion                                         
!! eve  13 obs rejected because of the subsampling                      
!                                                                      !
!                                                                      !
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qualitycheck

   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   IF ( qck%clm .OR. qck%eofbgr ) CALL rdclim

   IF ( qck%eofbgr ) CALL get_bgerr

   IF ( mpi%nproc.GT.1 ) THEN
      IF ( qck%eofbgr ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, qck%eta)
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%tem)
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%sal)
      ENDIF
      IF ( qck%clm ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%climtem)
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%climsal)
      ENDIF
      IF ( coastrej%any ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                               &
            -1_i4, 1_i4, grd%im+1,                               &
            -1_i4, 1_i4, grd%jm+1,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, coastrej%distc)
      ENDIF
   ENDIF

   IF ( arg%no .GT. 0 .AND. obs%arg .NE. 0 ) CALL qc_arg

   IF ( gld%no .GT. 0 .AND. obs%gld .NE. 0 ) CALL qc_gld

   IF ( gvl%no .GT. 0 .AND. obs%gvl .NE. 0 ) CALL qc_gvl

   IF ( sla%no .GT. 0 .AND. obs%sla .NE. 0 ) CALL qc_sla

   IF ( sst%no .GT. 0 .AND. obs%sst .NE. 0 ) CALL qc_sst

   IF ( tra%no .GT. 0 .AND. obs%tra .NE. 0 ) CALL qc_tra

   IF ( trd%no .GT. 0 .AND. obs%trd .NE. 0 ) CALL qc_trd

   IF ( vdr%no .GT. 0 .AND. obs%vdr .NE. 0 ) CALL qc_vdr

   IF ( xbt%no .GT. 0 .AND. obs%xbt .NE. 0 ) CALL qc_xbt

   IF ( mpi%nproc.GT.1 ) THEN
      IF ( qck%eofbgr ) THEN
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, qck%eta)
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%tem)
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%sal)
      ENDIF
      IF ( qck%clm ) THEN
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%climtem)
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, qck%climsal)
      ENDIF
      IF ( coastrej%any ) THEN
         CALL exo_mpi( 1_i4, 0_i4,                              &
            1_i4, grd%im+1, 1_i4,                               &
            1_i4, grd%jm+1, 1_i4,                               &
            1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, coastrej%distc)
      ENDIF
   ENDIF


END SUBROUTINE qualitycheck
! ----------------------------------------------------------------------
!>  perform background quality check and reject observations accordingly
!!
!!  it works as follows:
!!
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a USEr-defined factor.
!!  when the relation is satisfied, the obs is rejected.
!
!   Version 1: Andrea Storto 2021
!              Mario Adani 2023
! ----------------------------------------------------------------------
SUBROUTINE qc_bgerr(nval,res,err,bgerr,tresh,flag)

   USE set_knd
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4),INTENT(IN)          :: nval
   REAL(r8),INTENT(IN)             :: res(nval)
   REAL(r8),INTENT(IN)             :: err(nval)
   REAL(r8),INTENT(IN)             :: bgerr(nval)
   REAL(r8),INTENT(IN)             :: tresh(nval)
   INTEGER(i8),INTENT(INOUT)       :: flag(nval)
   INTEGER(i4)                     :: k
   REAL(r8)                        :: zomgp2,ztotvar,zratio

   DO k = 1,nval
      IF ( flag(k) .EQ. 1 ) THEN
         zomgp2  = res(k) * res(k)
         ztotvar = bgerr(k)*bgerr(k) + err(k)*err(k)
         zratio  = zomgp2 / ztotvar
         IF ( zratio .GT. tresh(k) ) flag(k) = 0
      ENDIF
   ENDDO
END SUBROUTINE qc_bgerr
! ----------------------------------------------------------------------
!>  Vertical profile check
!!
!!  Reject measurements if some observations above
!!  within the same profile were rejected
!!
!   Version 1: Andrea Storto 2021
!              Mario Adani 2023
! ----------------------------------------------------------------------
SUBROUTINE vertcheck(no,ino,dpt,par,pararr,flc)

   USE set_knd
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4),INTENT(IN)   :: no
   INTEGER(i4),INTENT(IN)   :: ino(no)
   REAL(r8),INTENT(IN)      :: dpt(no)
   REAL(i8),INTENT(INOUT)   :: flc(no)
   INTEGER(i4),OPTIONAL     :: par
   INTEGER(i4),OPTIONAL     :: pararr(no)
   INTEGER(i4)              :: dimres
   INTEGER(i4)              :: res(no)
   INTEGER(i4),ALLOCATABLE  :: unique(:)
   INTEGER(i4)              :: i
   LOGICAL                  :: msk_profile(no)
   LOGICAL                  :: msk_profile_notgoodflag(no)
   LOGICAL                  :: msk_profile_goodflag(no)
   REAL(r8)                 :: min_dpt_profile_goodflag
   REAL(r8)                 :: max_dpt_profile_notgoodflag

! Check is subroutine is called correctly
   IF ( present(par)  .AND. .NOT. present(pararr) .OR.    &
      present(pararr) .AND. .NOT. present(par   ) ) THEN
      WRITE (drv%dia,*)'vertcheck is not called in appropriate way'
      CALL abort
   ENDIF

! Set profile
   CALL remove_dups(no,ino,dimres,res)
   ALLOCATE ( unique(dimres) )
   unique = res(1:dimres)

! Start quality control
   DO i = 1, dimres
      ! Initialize for each profile
      msk_profile(:)              = .False.
      msk_profile_goodflag(:)     = .False.
      msk_profile_notgoodflag(:)  = .False.
      min_dpt_profile_goodflag    = 1.e16_r8
      max_dpt_profile_notgoodflag = 1.e16_r8
      IF ( present(par) ) THEN
         ! Mask for profile
         WHERE ( unique(i) .EQ. ino  .AND. pararr .EQ. par ) msk_profile(:)       = .True.
         ! Mask for profile and good flag
         WHERE ( unique(i) .EQ. ino .AND. flc .EQ. 1 .AND. pararr .EQ. par ) &
            msk_profile_goodflag(:)    = .True.
         ! Mask for profile and not good flag
         WHERE ( unique(i) .EQ. ino .AND. flc .EQ. 0 .AND. pararr .EQ. par ) &
            msk_profile_notgoodflag(:) = .True.
      ELSE
         ! Mask for profile
         WHERE ( unique(i) .EQ. ino                  ) msk_profile(:)             = .True.
         ! Mask for profile and good flag
         WHERE ( unique(i) .EQ. ino .AND. flc .EQ. 1 ) msk_profile_goodflag(:)    = .True.
         ! Mask for profile and not good flag
         WHERE ( unique(i) .EQ. ino .AND. flc .EQ. 0 ) msk_profile_notgoodflag(:) = .True.
      ENDIF

      ! Find MAX depth for each profile with not good flag
      IF ( ANY(msk_profile_notgoodflag )) THEN
         max_dpt_profile_notgoodflag = MAXVAL(dpt,mask=msk_profile_notgoodflag)
      ENDIF
      ! Find MIN depth for each profile with good flag
      IF ( ANY(msk_profile_goodflag )) THEN
         min_dpt_profile_goodflag = MINVAL(dpt,mask=msk_profile_goodflag)
      ENDIF

      ! Finally...
      ! If min depth with good flag > MAX depth wih not good flag  THEN get rid of the profile
      IF ( min_dpt_profile_goodflag .GT. max_dpt_profile_notgoodflag ) THEN
         WHERE ( msk_profile ) flc(:)=0
      ENDIF
   ENDDO

   DEALLOCATE ( unique )

END SUBROUTINE vertcheck
!---------------------------------------
!> Remove duplicates
!
! Version 1: Mario Adani 2023
!---------------------------------------
SUBROUTINE remove_dups(n,ino,k,res)

   IMPLICIT NONE

   INTEGER,INTENT(IN)  :: n
   INTEGER,INTENT(IN)  :: ino(n)
   INTEGER,INTENT(OUT) :: k                   ! The number of unique elements
   INTEGER,INTENT(OUT) :: res(n)              ! The output
   INTEGER :: i

   k = 1
   res(:) = -999
   res(1) = ino(1)
   DO i=2,n
      ! IF the number already exist in res check next
      IF ( ANY(res == ino(i)) ) cycle
      ! No match found so add it to the output
      k = k + 1
      res(k) = ino(i)
   ENDDO

END SUBROUTINE remove_dups
