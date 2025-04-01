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
!> The main driver for the OceanVar
!!
!! Main program for OceanVar 3D variational data assimilation
!!
!                                                                      !
! Version 0.1: Srdjan Dobricic 2006                                    !
!   Horizontal covariance with recursive filters, vertical with EOFs,  !
!   assimilation of satellite observations of SLA, in situ observations!
!   by XBT and ARGO floats                                             !
!                                                                      !
! Version 0.2: Srdjan Dobricic 2007                                    !
!   Multigrid method. Internal boundaries for horizontal covariances.  !
!                                                                      !
! Versin 2.0: Mario Adani 2024                                         !
!   Quality control                                                    !
!   Thinning                                                           !
!   Diffusive Filter                                                   !
!   Observational Error                                                !
!   Simple Balance Model based on Dynamic Height                       !
!                                                                      !
!                                                                      !
!-----------------------------------------------------------------------
PROGRAM oceanvar


   USE set_knd
   USE drv_str
   USE mpi_str
   USE obs_str, ONLY : obserr, huberqc, thin

   IMPLICIT NONE

   INTEGER(i4)   ::  kts, ktr, ierr

! ---
! Initialize mpi
   CALL MPI_INIT(ierr)

! ---
! Initialize mpi, diagnostics and read namelist
   CALL def_nml

! ---
! Initialize time
   CALL ini_time

! ---
! Outer loop - SST assimilation
   DO kts = 1,drv%nts

      WRITE (drv%dia,*) ' ---------- '
      WRITE (drv%dia,*) ' Outer loop ', kts, drv%nts
      drv%kts = kts
      CALL FLUSH(drv%dia)

! ---
! Outer loop - multigrid
      DO ktr = 1,drv%ntr

         WRITE (drv%dia,*) ' -------------------- '
         WRITE (drv%dia,*) ' Outer loop-multigrid ', ktr, drv%ntr
         drv%ktr = ktr
         CALL FLUSH(drv%dia)

! ---
! Define grid parameters
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ---------------------- '
            WRITE (drv%dia,*) ' Define grid parameters ', ktr, drv%ratio(ktr)
            CALL def_grd
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Define constants for background covariances
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ------------------------------------------- '
            WRITE (drv%dia,*) ' Define constants for background covariances '
            CALL def_cov
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Get observations
         IF ( ktr .EQ. 1 .AND.  kts .EQ. 1 ) THEN
            WRITE (drv%dia,*) ' --------------------------- '
            WRITE (drv%dia,*) ' Get observations            ', ktr, kts
            CALL get_obs
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Define interpolation parameters
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ------------------------------- '
            WRITE (drv%dia,*) ' Define interpolation parameters '
            CALL int_par
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Define observational errors
         IF ( ktr .EQ. 1 .AND.  .NOT. obserr%rd_err_ff ) THEN
            WRITE (drv%dia,*) ' --------------------------- '
            WRITE (drv%dia,*) ' Define observational errors '
            CALL obserrors
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Apply quality control to observations
         IF ( ktr .EQ. 1 .AND.  kts .EQ. 1 ) THEN
            WRITE (drv%dia,*) ' --------------------- '
            WRITE (drv%dia,*) ' Apply quality control '
            CALL qualitycheck
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Thinning
         IF ( ktr .EQ. 1 .AND.  kts .EQ. 1 .AND. thin%any ) THEN
            WRITE (drv%dia,*) ' -------------- '
            WRITE (drv%dia,*) ' Apply thinning '
            CALL thinning
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! If assimilating SST save misfits
         IF ( ktr .EQ. 1 .AND. kts .EQ. 1 .AND. drv%nts .EQ. 2 ) THEN
            WRITE (drv%dia,*) ' ------------ '
            WRITE (drv%dia,*) ' Save Misfits ', ktr, kts
            CALL sav_msf
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Define observational vector
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ---------------------------- '
            WRITE (drv%dia,*) ' Define observational vector  '
            CALL obs_vec
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Initialize huber
         IF (ktr .EQ. 1 .AND. huberqc%any ) THEN
            WRITE (drv%dia,*) ' ---------------------------------- '
            WRITE (drv%dia,*) ' Initialize Huber Norm Distribution '
            CALL ini_huber
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Localize verticaly by using the mixed layer depth
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ------------------------------------------------- '
            WRITE (drv%dia,*) ' Localize verticaly by using the mixed layer depth '
            CALL rdmxd
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Localize horizontaly around observations
         IF ( ktr .EQ. 1 .OR. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' -----------------------------------------'
            WRITE (drv%dia,*) ' Localize horizontaly around observations '
            CALL def_loc
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Initialise barotropic model
         IF ( drv%bmd(drv%ktr) .EQ. 1 ) THEN
            WRITE (drv%dia,*) ' ----------------------------'
            WRITE (drv%dia,*) ' Initialize barotropic model '
            CALL ini_bmd
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Initialise barotropic model
         IF ( drv%bal(drv%ktr) .EQ. 1 ) THEN
            WRITE (drv%dia,*) ' -------------------------------------- '
            WRITE (drv%dia,*) ' Initialize simplified balance operator '
            CALL ini_bal
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Initialize cost function and its gradient
         WRITE (drv%dia,*) ' ----------------------------------------- '
         WRITE (drv%dia,*) ' Initialize cost function and its gradient '
         CALL ini_cfn
         CALL FLUSH(drv%dia)

! ---
! Calculate the initial norm of the gradient
         IF ( ktr .GT. 1 ) THEN
            WRITE (drv%dia,*) ' ------------------------------------------ '
            WRITE (drv%dia,*) ' Calculate the initial norm of the gradient '
            CALL ini_nrm
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Initialise from old iterration
         IF ( ktr .GT. 1 .AND. drv%ratio(ktr) .NE. 1.0_r8 ) THEN
            WRITE (drv%dia,*) ' ------------------------------ '
            WRITE (drv%dia,*) ' Initialize from old iterration '
            CALL ini_itr
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! Minimize the cost function (inner loop)
         WRITE (drv%dia,*) ' ---------------------------------------- '
         WRITE (drv%dia,*) '  Minimize the cost function (inner loop) '
         CALL MIN_cfn
         CALL FLUSH(drv%dia)

! ---
! Save old iterration
         IF ( ktr .LT. drv%ntr ) THEN
            IF ( drv%ratio(ktr+1) .NE. 1.0_r8 ) THEN
               CALL sav_itr
            ENDIF
         ENDIF

! ---
! Convert to innovations
         IF ( drv%ktr .EQ. drv%ntr .AND. drv%kts .EQ. drv%nts ) THEN
            CALL cnv_inn
            CALL FLUSH(drv%dia)
         ENDIF

! ---
! End of outer loop - multigrid
      ENDDO
      WRITE (drv%dia,*) ' -------------------------------'
      WRITE (drv%dia,*) '  End of outer loop - multigrid '
      CALL FLUSH(drv%dia)

! ---
! If assimilating SST modify misfits and save increments
      IF ( drv%ktr .EQ. drv%ntr .AND. drv%kts .EQ. 1 .AND. drv%nts .EQ. 2 ) THEN
         CALL cnv_inn
         CALL obsop
         CALL mod_msf
         CALL sav_itr
      ENDIF

! ---
! End of outer loop - SST assimilation
   ENDDO

   WRITE (drv%dia,*) ' --------------------------------------'
   WRITE (drv%dia,*) '  End of outer loop - SST assimilation '
   CALL FLUSH(drv%dia)

! ---
! If assimilating SST modify increments
   IF ( drv%nts .EQ. 2 ) THEN
      CALL mod_inc
   ENDIF

! ---
! Write outputs and diagnostics
   CALL wrt_dia

! ---
! End of mpi
   CALL MPI_FINALIZE(ierr)

!-----------------------------------------------------------------

END PROGRAM oceanvar
