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
!> XBT quality control                                                 !
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist or computed from EOFs
!! - Climatology qc abs[y-h(xclim)] > threshold,
!! - Too close tp the coast,
!! - If profile is only subsurface ,
!!
! Version 1: Nome Cognome Anno                                         !
!-----------------------------------------------------------------------
SUBROUTINE qc_xbt

   USE set_knd
   USE obs_str, ONLY : xbt, qck, coastrej
   USE drv_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER(i4)              :: kk, k, j, i
   REAL(r8),ALLOCATABLE     :: treshold(:)
   REAL(r8)                 :: val
   INTEGER                  :: ierr
   INTEGER(i8),ALLOCATABLE  :: allflag(:)
   REAL(r8)                 :: weightsum,pq1,pq2,pq3,pq4
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- XBT observations:                ',xbt%nc

!Initialize
   ALLOCATE ( xbt%bgerr(xbt%no) )
   xbt%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( xbt%flc .EQ. 1 .AND. xbt%par .EQ. 1 .AND. ABS(xbt%res) .GT. qck%res_tem)
         xbt%flc(:) = 0
         xbt%eve(:) = 5
      ENDWHERE
      CALL reset_xbt
      WRITE (drv%dia,*)' ---- after residual check:             ',xbt%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  WHERE sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a USEr-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
! with constant background from NAMELIST
      ALLOCATE ( treshold(xbt%no) )
      ALLOCATE ( flc(xbt%no) )
      flc = xbt%flc
      WHERE ( xbt%par .EQ. 1 ) xbt%bgerr(:) = qck%bgr_tem(1)
      WHERE ( xbt%par .EQ. 1 ) treshold (:) = qck%bgr_tem(2)
      CALL qc_bgerr(xbt%no,xbt%res,xbt%err,xbt%bgerr,treshold,xbt%flc)
      CALL reset_xbt
      WHERE ( ( flc .EQ. 1 ) .AND. ( xbt%flc .EQ. 0 ) ) xbt%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',xbt%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
! with computed background from eofs
      ALLOCATE ( flc(xbt%no) )
      flc = xbt%flc
      IF ( .NOT. ALLOCATED(treshold) ) ALLOCATE ( treshold(xbt%no)  )
      CALL int_bgerr_xbt
      WHERE ( xbt%par .EQ. 1 ) treshold (:) = qck%bgr_tem(2)
      CALL qc_bgerr(xbt%no,xbt%res,xbt%err,xbt%bgerr,treshold,xbt%flc)
      CALL reset_xbt
      WHERE ( ( flc .EQ. 1 ) .AND. ( xbt%flc .EQ. 0 ) ) xbt%eve(:) =  7
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after EOFs background check:      ',xbt%nc
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! if greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      DO kk=1,xbt%no
         IF ( xbt%flc(kk) .EQ. 1 ) THEN
            i = xbt%ib(kk)
            j = xbt%jb(kk)
            k = xbt%kb(kk)
            ! for safe interpolation
            IF ( ANY(qck%climtem(i:i+1,j:j+1,k:k+1) .LE. 0) ) cycle
            val = xbt%pq1(kk) * qck%climtem(i  ,j ,k ) +       &
                  xbt%pq2(kk) * qck%climtem(i+1,j ,k ) +       &
                  xbt%pq3(kk) * qck%climtem(i  ,j+1,k) +       &
                  xbt%pq4(kk) * qck%climtem(i+1,j+1,k)
            IF ( ABS(val-xbt%val(kk)) .GT. qck%clm_lim(1) ) THEN
               xbt%flc(kk) = 0
               xbt%eve(kk) = 8
            ENDIF
         ENDIF
      ENDDO
      CALL reset_xbt
      WRITE (drv%dia,*)' ---- after climatological check:       ',xbt%nc
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%xbt ) THEN
      DO kk = 1,xbt%no
         IF ( xbt%flc(kk) .EQ. 1 ) THEN
            i = xbt%ib(kk)
            j = xbt%jb(kk)
            weightsum = xbt%pq1(kk) + xbt%pq2(kk) + xbt%pq3(kk) + xbt%pq4(kk)
            pq1 = xbt%pq1(kk) / weightsum
            pq2 = xbt%pq2(kk) / weightsum
            pq3 = xbt%pq3(kk) / weightsum
            pq4 = xbt%pq4(kk) / weightsum
            val = pq1 * coastrej%distc(i  ,j  ) +       &
                  pq2 * coastrej%distc(i+1,j  ) +       &
                  pq3 * coastrej%distc(i  ,j+1) +       &
                  pq4 * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_xbt*1000._r8 ) THEN
               xbt%flc(kk) = 0
               xbt%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_xbt
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',xbt%nc
   ENDIF

! 5 ) Perform vertical check
!!  - IF min depth with good flag > max depth wih not good flag
!!  then get rid of the profile
   IF ( qck%vert ) THEN
      ALLOCATE ( allflag(xbt%no) )
      ALLOCATE ( flc(xbt%no) )
      flc = xbt%flc
      !Just in case a profile starts in a region and ENDs in another one
      IF ( mpi%nproc.GT.1 ) THEN
         CALL MPI_ALLREDUCE(  xbt%flc,  allflag, xbt%no, mpi%i8  ,   &
            MPI_MAX,  mpi%comm, ierr)
      ELSE
         allflag = xbt%flc
      ENDIF
      CALL vertcheck(xbt%no,xbt%ino,xbt%dpt,1,xbt%par,xbt%flc)
      ! Find minimim
      WHERE ( allflag .LT. xbt%flc ) xbt%flc = allflag
      CALL reset_xbt
      WHERE ( ( flc .EQ. 1 ) .AND. ( xbt%flc .EQ. 0 ) ) xbt%eve(:) =  11
      WRITE (drv%dia,*)' ---- after vertical check:             ',xbt%nc
      DEALLOCATE ( allflag )
      DEALLOCATE ( flc )
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_xbt
!-----------------------------------------------------------------------
!> Interpolation background error
!                                                                      !
! Version 1: Mario Adani 2023                                          !             
!-----------------------------------------------------------------------
SUBROUTINE int_bgerr_xbt

   USE set_knd
   USE obs_str, ONLY : xbt, qck

   IMPLICIT NONE

   INTEGER(i4)   :: i, j, k, kk

   DO kk = 1, xbt%no
      IF ( xbt%flc(kk) .EQ. 1 ) THEN
         i = xbt%ib(kk)
         j = xbt%jb(kk)
         k = xbt%kb(kk)
         xbt%bgerr(kk) = xbt%pq1(kk) * qck%tem(i  ,j  ,k  ) +       &
                         xbt%pq2(kk) * qck%tem(i+1,j  ,k  ) +       &
                         xbt%pq3(kk) * qck%tem(i  ,j+1,k  ) +       &
                         xbt%pq4(kk) * qck%tem(i+1,j+1,k  ) +       &
                         xbt%pq5(kk) * qck%tem(i  ,j  ,k+1) +       &
                         xbt%pq6(kk) * qck%tem(i+1,j  ,k+1) +       &
                         xbt%pq7(kk) * qck%tem(i  ,j+1,k+1) +       &
                         xbt%pq8(kk) * qck%tem(i+1,j+1,k+1)
      ENDIF
   ENDDO

END SUBROUTINE int_bgerr_xbt
!-----------------------------------------------------------------------
!> Reset observation
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_xbt

   USE obs_str, ONLY : xbt

   IMPLICIT NONE

   WHERE ( xbt%flc .EQ. 0 )
      xbt%bia(:) = 0.
      xbt%res(:) = 0.
      xbt%inc(:) = 0.
      xbt%b_a(:) = 0.
      xbt%pq1(:) = 0.
      xbt%pq2(:) = 0.
      xbt%pq3(:) = 0.
      xbt%pq4(:) = 0.
      xbt%pq5(:) = 0.
      xbt%pq6(:) = 0.
      xbt%pq7(:) = 0.
      xbt%pq8(:) = 0.
   END WHERE

   xbt%nc = SUM(xbt%flc)

END SUBROUTINE reset_xbt
