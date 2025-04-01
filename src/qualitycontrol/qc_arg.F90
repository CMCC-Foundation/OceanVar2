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
!> ARGO quality control
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
SUBROUTINE qc_arg

   USE set_knd
   USE obs_str, ONLY : arg, qck, coastrej
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
   WRITE (drv%dia,*)' ---- ARGO observations:                ',arg%nc
   CALL FLUSH(drv%dia)

!Initialize
   ALLOCATE ( arg%bgerr(arg%no) )
   arg%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( arg%flc .EQ. 1 .AND. arg%par.EQ.1 .AND. ABS(arg%res) .GT. qck%res_tem )
         arg%flc = 0
         arg%eve = 5
      ENDWHERE
      WHERE ( arg%flc .EQ. 1 .AND. arg%par.EQ.2 .AND. ABS(arg%res) .GT. qck%res_sal )
         arg%flc = 0
         arg%eve = 5
      ENDWHERE
      CALL reset_arg
      WRITE (drv%dia,*)' ---- after residual check:             ',arg%nc
      CALL FLUSH (drv%dia)
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
      ALLOCATE ( treshold(arg%no) )
      ALLOCATE ( flc(arg%no) )
      flc = arg%flc
! with constant background from namelist
      WHERE ( arg%par .EQ. 1 )
         arg%bgerr(:) = qck%bgr_tem(1)
         treshold (:) = qck%bgr_tem(2)
      END WHERE
      WHERE ( arg%par .EQ. 2 )
         arg%bgerr(:) = qck%bgr_sal(1)
         treshold (:) = qck%bgr_sal(2)
      END WHERE
      CALL qc_bgerr(arg%no,arg%res,arg%err,arg%bgerr,treshold,arg%flc)
      CALL reset_arg
      WHERE ( ( flc .EQ. 1 ) .AND. ( arg%flc .EQ. 0 ) ) arg%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',arg%nc
      CALL FLUSH(drv%dia)
   ENDIF

   IF ( qck%eofbgr ) THEN
! with computed background from eofs
      ALLOCATE ( flc(arg%no) )
      flc = arg%flc
      CALL int_bgerr_arg
      IF ( .NOT. ALLOCATED(treshold) ) ALLOCATE ( treshold(arg%no)  )
      WHERE ( arg%par .EQ. 1 ) treshold (:) = qck%bgr_tem(2)
      WHERE ( arg%par .EQ. 2 ) treshold (:) = qck%bgr_sal(2)
      CALL qc_bgerr(arg%no,arg%res,arg%err,arg%bgerr,treshold,arg%flc)
      CALL reset_arg
      WHERE ( ( flc .EQ. 1 ) .AND. ( arg%flc .EQ. 0 ) ) arg%eve(:) =  7
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after EOFs background check:      ',arg%nc
      CALL FLUSH(drv%dia)
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! if greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      DO kk = 1,arg%no
         IF ( arg%flc(kk) .EQ. 1 ) THEN
            i = arg%ib(kk)
            j = arg%jb(kk)
            k = arg%kb(kk)
            IF ( arg%par(kk) .EQ. 1 )       THEN
               ! for safe interpolation
               IF ( ANY(qck%climtem(i:i+1,j:j+1,k:k+1 ) .LE. 0) ) cycle
               val = arg%pq1(kk) * qck%climtem(i  ,j  ,k  ) +       &
                     arg%pq2(kk) * qck%climtem(i+1,j  ,k  ) +       &
                     arg%pq3(kk) * qck%climtem(i  ,j+1,k  ) +       &
                     arg%pq4(kk) * qck%climtem(i+1,j+1,k  ) +       &
                     arg%pq5(kk) * qck%climtem(i  ,j  ,k+1) +       &
                     arg%pq6(kk) * qck%climtem(i+1,j  ,k+1) +       &
                     arg%pq7(kk) * qck%climtem(i  ,j+1,k+1) +       &
                     arg%pq8(kk) * qck%climtem(i+1,j+1,k+1)
            ELSEIF ( arg%par(kk) .EQ. 2 )    THEN
               ! for safe interpolation
               IF ( ANY(qck%climsal(i:i+1,j:j+1,k:k+1) .LE. 0) ) cycle
               val = arg%pq1(kk) * qck%climsal(i  ,j  ,k  ) +       &
                     arg%pq2(kk) * qck%climsal(i+1,j  ,k  ) +       &
                     arg%pq3(kk) * qck%climsal(i  ,j+1,k  ) +       &
                     arg%pq4(kk) * qck%climsal(i+1,j+1,k  ) +       &
                     arg%pq5(kk) * qck%climsal(i  ,j  ,k+1) +       &
                     arg%pq6(kk) * qck%climsal(i+1,j  ,k+1) +       &
                     arg%pq7(kk) * qck%climsal(i  ,j+1,k+1) +       &
                     arg%pq8(kk) * qck%climsal(i+1,j+1,k+1)
            ELSE
               WRITE (drv%dia,*)' Parameter not supported for ARGO',arg%par(kk)
               CALL FLUSH(drv%dia)
               CALL abort
            ENDIF
            IF ( ABS(val-arg%val(kk)) .GT. qck%clm_lim(arg%par(kk)) ) THEN
               arg%flc(kk) = 0
               arg%eve(kk) = 8
            ENDIF
         ENDIF
      ENDDO
      CALL reset_arg
      WRITE (drv%dia,*)' ---- after climatological check:       ',arg%nc
      CALL FLUSH(drv%dia)
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%arg ) THEN
      DO kk = 1,arg%no
         IF ( arg%flc(kk) .EQ. 1 ) THEN
            i = arg%ib(kk)
            j = arg%jb(kk)
            weightsum = arg%pq1(kk) + arg%pq2(kk) + arg%pq3(kk) + arg%pq4(kk)
            pq1 = arg%pq1(kk) / weightsum
            pq2 = arg%pq2(kk) / weightsum
            pq3 = arg%pq3(kk) / weightsum
            pq4 = arg%pq4(kk) / weightsum
            val = pq1 * coastrej%distc(i  ,j  ) +       &
                  pq2 * coastrej%distc(i+1,j  ) +       &
                  pq3 * coastrej%distc(i  ,j+1) +       &
                  pq4 * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_arg*1000._r8 ) THEN
               arg%flc(kk) = 0
               arg%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_arg
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',arg%nc
      CALL FLUSH(drv%dia)
   ENDIF

! 5 ) Perform vertical check
!!  - if min depth with good flag > max depth wih not good flag
!!  then get rid of the profile
   IF ( qck%vert ) THEN
      ALLOCATE ( allflag(arg%no) )
      ALLOCATE ( flc(arg%no) )
      flc = arg%flc
      !Just in case a profile starts in a region and ENDs in another one
      IF ( mpi%nproc .GT. 1 ) THEN
         CALL MPI_ALLREDUCE(  arg%flc,  allflag, arg%no, mpi%i8  ,   &
            MPI_MAX,  mpi%comm, ierr)
      ELSE
         allflag = arg%flc
      ENDIF
      CALL vertcheck(arg%no,arg%ino,arg%dpt,1,arg%par,allflag)
      CALL vertcheck(arg%no,arg%ino,arg%dpt,2,arg%par,allflag)
      ! Find minimim
      WHERE ( allflag .LT. arg%flc ) arg%flc = allflag
      CALL reset_arg
      WHERE ( ( flc .EQ. 1 ) .AND. ( arg%flc .EQ. 0 ) ) arg%eve(:) =  11
      WRITE (drv%dia,*)' ---- after vertical check:             ',arg%nc
      CALL FLUSH (drv%dia)
      DEALLOCATE ( allflag )
      DEALLOCATE ( flc )
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_arg
!-----------------------------------------------------------------------
!> Interpolation background error
!                                                                      !
! Version 1: Mario Adani 2023                                          !             
!-----------------------------------------------------------------------
SUBROUTINE int_bgerr_arg

   USE set_knd
   USE obs_str, ONLY : arg, qck

   IMPLICIT NONE

   INTEGER(i4)   :: i, j, k, kk

   DO kk = 1, arg%no
      IF ( arg%flc(kk) .EQ. 1 ) THEN
         i = arg%ib(kk)
         j = arg%jb(kk)
         k = arg%kb(kk)
         IF ( arg%par(kk) .EQ. 1 ) THEN
            arg%bgerr(kk) = arg%pq1(kk) * qck%tem(i  ,j  ,k  ) +       &
                            arg%pq2(kk) * qck%tem(i+1,j  ,k  ) +       &
                            arg%pq3(kk) * qck%tem(i  ,j+1,k  ) +       &
                            arg%pq4(kk) * qck%tem(i+1,j+1,k  ) +       &
                            arg%pq5(kk) * qck%tem(i  ,j  ,k+1) +       &
                            arg%pq6(kk) * qck%tem(i+1,j  ,k+1) +       &
                            arg%pq7(kk) * qck%tem(i  ,j+1,k+1) +       &
                            arg%pq8(kk) * qck%tem(i+1,j+1,k+1)
         ELSEIF ( arg%par(kk).EQ.2 ) THEN
            arg%bgerr(kk) = arg%pq1(kk) * qck%sal(i  ,j  ,k  ) +       &
                            arg%pq2(kk) * qck%sal(i+1,j  ,k  ) +       &
                            arg%pq3(kk) * qck%sal(i  ,j+1,k  ) +       &
                            arg%pq4(kk) * qck%sal(i+1,j+1,k  ) +       &
                            arg%pq5(kk) * qck%sal(i  ,j  ,k+1) +       &
                            arg%pq6(kk) * qck%sal(i+1,j  ,k+1) +       &
                            arg%pq7(kk) * qck%sal(i  ,j+1,k+1) +       &
                            arg%pq8(kk) * qck%sal(i+1,j+1,k+1)
         ENDIF
      ENDIF
   ENDDO

END SUBROUTINE int_bgerr_arg
!-----------------------------------------------------------------------
!> Reset observation
!                                                                      !
! Version 1: Mario Adani 2023                                        !             
!-----------------------------------------------------------------------
SUBROUTINE reset_arg

   USE obs_str, ONLY : arg

   IMPLICIT NONE

   WHERE ( arg%flc .EQ. 0 )
      arg%bia(:) = 0.
      arg%res(:) = 0.
      arg%inc(:) = 0.
      arg%b_a(:) = 0.
      arg%pq1(:) = 0.
      arg%pq2(:) = 0.
      arg%pq3(:) = 0.
      arg%pq4(:) = 0.
      arg%pq5(:) = 0.
      arg%pq6(:) = 0.
      arg%pq7(:) = 0.
      arg%pq8(:) = 0.
   END WHERE
   arg%nc = SUM(arg%flc)

END SUBROUTINE reset_arg
