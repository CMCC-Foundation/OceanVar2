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
!> GLIDER quality control
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist or computed from EOFs
!! - Climatology qc abs[y-h(xclim)] > threshold,
!! - Too close tp the coast,
!! - If profile is only subsurface ,
!!
!                                                                      !
! Version 1: Andrea Storto 2022                                        !             
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_gld

   USE set_knd
   USE obs_str, ONLY : gld, qck, coastrej
   USE drv_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER(i4)              :: kk, k, j, i
   REAL(r8)                 :: val
   REAL(r8)                 :: weightsum,pq1,pq2,pq3,pq4
   INTEGER                  :: ierr
   INTEGER(i8),ALLOCATABLE  :: allflag(:)
   INTEGER(i8),ALLOCATABLE  :: flc(:)
   REAL(r8),ALLOCATABLE     :: treshold(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- GLIDER observations:              ',gld%nc

!Initialize
   ALLOCATE ( gld%bgerr(gld%no) )
   gld%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( gld%flc .EQ. 1 .AND. gld%par .EQ. 1 .AND. ABS(gld%res) .GT. qck%res_tem )
         gld%flc = 0
         gld%eve = 5
      ENDWHERE
      WHERE ( gld%flc .EQ. 1 .AND. gld%par .EQ. 2 .AND. ABS(gld%res) .GT. qck%res_sal )
         gld%flc = 0
         gld%eve = 5
      ENDWHERE
      CALL reset_gld
      WRITE (drv%dia,*)' ---- after residual check:             ',gld%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
! with constant background from NAMELIST
      ALLOCATE ( treshold(gld%no) )
      ALLOCATE ( flc(gld%no) )
      flc = gld%flc
      WHERE ( gld%par .EQ. 1 )
         gld%bgerr(:) = qck%bgr_tem(1)
         treshold (:) = qck%bgr_tem(2)
      END WHERE
      WHERE ( gld%par .EQ. 2 )
         gld%bgerr(:) = qck%bgr_sal(1)
         treshold (:) = qck%bgr_sal(2)
      END WHERE
      CALL qc_bgerr(gld%no,gld%res,gld%err,gld%bgerr,treshold,gld%flc)
      CALL reset_gld
      WHERE ( ( flc .EQ. 1 ) .AND. ( gld%flc .EQ. 0 ) ) gld%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',gld%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
! with computed background from eofs
      IF ( .NOT. ALLOCATED(treshold) ) ALLOCATE ( treshold(gld%no) )
      ALLOCATE ( flc(gld%no) )
      flc = gld%flc
      CALL int_bgerr_gld
      WHERE ( gld%par .EQ. 1 ) treshold (:) = qck%bgr_tem(2)
      WHERE ( gld%par .EQ. 2 ) treshold (:) = qck%bgr_sal(2)
      CALL qc_bgerr(gld%no,gld%res,gld%err,gld%bgerr,treshold,gld%flc)
      CALL reset_gld
      WHERE ( ( flc .EQ. 1 ) .AND. ( gld%flc .EQ. 0 ) ) gld%eve(:) =  7
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after EOFs background check:      ',gld%nc
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate ABS(clim-obs)
!! IF greater THEN a USEr-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      DO kk = 1,gld%no
         IF ( gld%flc(kk) .EQ. 1 ) THEN
            i = gld%ib(kk)
            j = gld%jb(kk)
            k = gld%kb(kk)
            IF ( gld%par(kk) .EQ. 1 )       THEN
               ! for safe interpolation
               IF ( ANY(qck%climtem(i:i+1,j:j+1,k:k+1) .EQ. 0) ) cycle
               val = gld%pq1(kk) * qck%climtem(i  ,j  ,k  ) +       &
                     gld%pq2(kk) * qck%climtem(i+1,j  ,k  ) +       &
                     gld%pq3(kk) * qck%climtem(i  ,j+1,k  ) +       &
                     gld%pq4(kk) * qck%climtem(i+1,j+1,k  ) +       &
                     gld%pq5(kk) * qck%climtem(i  ,j  ,k+1) +       &
                     gld%pq6(kk) * qck%climtem(i+1,j  ,k+1) +       &
                     gld%pq7(kk) * qck%climtem(i  ,j+1,k+1) +       &
                     gld%pq8(kk) * qck%climtem(i+1,j+1,k+1)
            ELSEIF ( gld%par(kk) .EQ. 2 )    THEN
               ! for safe interpolation
               IF ( ANY(qck%climsal(i:i+1,j:j+1,k:k+1) .EQ. 0) ) cycle
               val = gld%pq1(kk) * qck%climsal(i  ,j  ,k  ) +       &
                     gld%pq2(kk) * qck%climsal(i+1,j  ,k  ) +       &
                     gld%pq3(kk) * qck%climsal(i  ,j+1,k  ) +       &
                     gld%pq4(kk) * qck%climsal(i+1,j+1,k  ) +       &
                     gld%pq5(kk) * qck%climsal(i  ,j  ,k+1) +       &
                     gld%pq6(kk) * qck%climsal(i+1,j  ,k+1) +       &
                     gld%pq7(kk) * qck%climsal(i  ,j+1,k+1) +       &
                     gld%pq8(kk) * qck%climsal(i+1,j+1,k+1)
            ELSE
               WRITE (drv%dia,*)' Parameter not supported for GLIDER',gld%par(kk)
               CALL abort
            ENDIF
            IF ( ABS(val-gld%val(kk)) .GT. qck%clm_lim(gld%par(kk)) )  THEN
               gld%flc(kk) = 0
               gld%eve(kk) = 8
            ENDIF
         ENDIF
      ENDDO
      CALL reset_gld
      WRITE (drv%dia,*)' ---- after climatological check:       ',gld%nc
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%gld ) THEN
      DO kk = 1,gld%no
         IF ( gld%flc(kk) .EQ. 1 ) THEN
            i=gld%ib(kk)
            j=gld%jb(kk)
            weightsum = gld%pq1(kk) + gld%pq2(kk) + gld%pq3(kk) + gld%pq4(kk)
            pq1 = gld%pq1(kk) / weightsum
            pq2 = gld%pq2(kk) / weightsum
            pq3 = gld%pq3(kk) / weightsum
            pq4 = gld%pq4(kk) / weightsum
            val = pq1 * coastrej%distc(i  ,j  ) +       &
                  pq2 * coastrej%distc(i+1,j  ) +       &
                  pq3 * coastrej%distc(i  ,j+1) +       &
                  pq4 * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_gld*1000._r8 ) THEN
               gld%flc(kk) = 0
               gld%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_gld
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',gld%nc
   ENDIF

! 5 ) Perform vertical check
!!  - if min depth with good flag > max depth wih not good flag
!!  then get rid of the profile
   IF ( qck%vert ) THEN
      ALLOCATE ( allflag(gld%no) )
      ALLOCATE ( flc(gld%no) )
      flc = gld%flc
      !Just in case a profile starts in a region and ENDs in another one
      IF ( mpi%nproc .GT. 1 ) THEN
         CALL MPI_ALLREDUCE(  gld%flc,  allflag, gld%no, mpi%i8  ,   &
            MPI_MAX,  mpi%comm, ierr)
      ELSE
         allflag = gld%flc
      ENDIF
      CALL vertcheck(gld%no,gld%ino,gld%dpt,1,gld%par,allflag)
      CALL vertcheck(gld%no,gld%ino,gld%dpt,2,gld%par,allflag)
      ! Find minimim
      WHERE ( allflag .LT. gld%flc ) gld%flc = allflag
      CALL reset_gld
      WHERE ( ( flc .EQ. 1 ) .AND. ( gld%flc .EQ. 0 ) ) gld%eve(:) =  11
      WRITE (drv%dia,*)' ---- after vertical check:             ',gld%nc
      DEALLOCATE ( allflag )
      DEALLOCATE ( flc )
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_gld
!-----------------------------------------------------------------------
!                                                                      !
!> Interpolation background error
!                                                                      !
! Version 1: Mario Adani 2023                                          !             
!-----------------------------------------------------------------------
SUBROUTINE int_bgerr_gld

   USE set_knd
   USE obs_str, ONLY : gld, qck

   IMPLICIT NONE

   INTEGER(i4)   :: i, j, k, kk

   DO kk = 1, gld%no
      IF ( gld%flc(kk) .EQ. 1 ) THEN
         i = gld%ib(kk)
         j = gld%jb(kk)
         k = gld%kb(kk)
         IF     ( gld%par(kk) .EQ. 1 ) THEN
            gld%bgerr(kk) = gld%pq1(kk) * qck%tem(i  ,j  ,k  ) +       &
                            gld%pq2(kk) * qck%tem(i+1,j  ,k  ) +       &
                            gld%pq3(kk) * qck%tem(i  ,j+1,k  ) +       &
                            gld%pq4(kk) * qck%tem(i+1,j+1,k  ) +       &
                            gld%pq5(kk) * qck%tem(i  ,j  ,k+1) +       &
                            gld%pq6(kk) * qck%tem(i+1,j  ,k+1) +       &
                            gld%pq7(kk) * qck%tem(i  ,j+1,k+1) +       &
                            gld%pq8(kk) * qck%tem(i+1,j+1,k+1)
         ELSEIF ( gld%par(kk).EQ.2 ) THEN
            gld%bgerr(kk) = gld%pq1(kk) * qck%sal(i  ,j  ,k  ) +       &
                            gld%pq2(kk) * qck%sal(i+1,j  ,k  ) +       &
                            gld%pq3(kk) * qck%sal(i  ,j+1,k  ) +       &
                            gld%pq4(kk) * qck%sal(i+1,j+1,k  ) +       &
                            gld%pq5(kk) * qck%sal(i  ,j  ,k+1) +       &
                            gld%pq6(kk) * qck%sal(i+1,j  ,k+1) +       &
                            gld%pq7(kk) * qck%sal(i  ,j+1,k+1) +       &
                            gld%pq8(kk) * qck%sal(i+1,j+1,k+1)
         ENDIF
      ENDIF
   ENDDO

END SUBROUTINE int_bgerr_gld
!-----------------------------------------------------------------------
!                                                                      !
!> Reset observation
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_gld

   USE obs_str, ONLY : gld

   IMPLICIT NONE

   WHERE ( gld%flc .EQ. 0 )
      gld%bia(:) = 0.
      gld%res(:) = 0.
      gld%inc(:) = 0.
      gld%b_a(:) = 0.
      gld%pq1(:) = 0.
      gld%pq2(:) = 0.
      gld%pq3(:) = 0.
      gld%pq4(:) = 0.
      gld%pq5(:) = 0.
      gld%pq6(:) = 0.
      gld%pq7(:) = 0.
      gld%pq8(:) = 0.
   END WHERE

   gld%nc = SUM(gld%flc)

END SUBROUTINE reset_gld
