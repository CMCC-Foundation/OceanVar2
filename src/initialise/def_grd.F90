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
!> Define the grid variables:                                   
!!
!! 1) Read Grids
!! 2) Read temperature and salinity for computation of the 
!!    expansion contraction coefficient (alpha/beta) for 
!!    tangent linear version of equation of state
!! 3) compute the alpha/beta: routine located in oceantools.F90 file
!! 4) Min/max value of dx/dy
!! 5) Close the domain
!! 6) Compute dx*dy and <dx*dy>
!! 7) Compte Coriolis
!! 8) rearrange longitude
!!
!                                                                      !
! Version 1: Srdjan Dobricic                  2006                     !
! Version 2: Mario Adani and Francesco Carere 2023                     !
!-----------------------------------------------------------------------
SUBROUTINE def_grd

   USE set_knd
   USE drv_str
   USE grd_str
   USE mpi_str
   USE cns_str, only : phy

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(I4)               :: i, j, k , ierr, ich
   REAL(r8)                  :: dxdyl,nwp
   REAL(r8), ALLOCATABLE     :: dxdyl_sumarr(:,:),msk_sumarr(:,:)
   REAL(r8)                  :: dxmin,dxmax,dymin,dymax


! ---
! Define grid
   grd%grd_mod  = drv%grid (drv%ktr)

!Read grid definition
   WRITE (drv%dia,*)' ----'
   WRITE (drv%dia,*)' ---- Reading GRIDS : '
   CALL rdgrd

   IF ( ANY(drv%nneos .GE. 2)  )       THEN
      CALL rdeos
   ENDIF

   IF  ( ANY(drv%nneos .EQ. 2)  )  THEN
      ALLOCATE ( grd%alpha3d(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
         grd%beta3d (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
      DO k = 1,grd%km
         DO j = 1-grd%jas,grd%jm+grd%jae
            DO i = 1-grd%ias,grd%im+grd%iae
               CALL alpha_beta(grd%salb(i,j,k),grd%temb(i,j,k),grd%alpha3d(i,j,k),grd%beta3d(i,j,k))
                               grd%alpha3d(i,j,k) = grd%alpha3d(i,j,k) * grd%msk(i,j,k)
                               grd%beta3d(i,j,k)  = grd%beta3d(i,j,k)  * grd%msk(i,j,k)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   grd%dlt = (grd%lat(1,2) - grd%lat(1,1))
   grd%dln = (grd%lon(2,1) - grd%lon(1,1))
   IF ( mpi%nproc.GT.1 ) THEN
      CALL mpi_bcast(grd%dlt,1,mpi%r8,0,mpi%comm,ierr)
      CALL mpi_bcast(grd%dln,1,mpi%r8,0,mpi%comm,ierr)
   ENDIF
   IF (mpi%ir.EQ.1      ) grd%msk(1,:,:)      = 0.0_r8
   IF (mpi%jr.EQ.1      ) grd%msk(:,1,:)      = 0.0_r8
   IF (mpi%ir.EQ.mpi%irm) grd%msk(grd%im,:,:) = 0.0_r8
   IF (mpi%jr.EQ.mpi%jrm) grd%msk(:,grd%jm,:) = 0.0_r8

! ---
   grd%dxdy(:,:) =  grd%dy(:,:) * grd%dx(:,:)
   grd%npsa = grd%img*grd%jmg

   IF ( mpi%nproc .GT. 1 ) THEN

      ALLOCATE(dxdyl_sumarr(grd%img, grd%jmg),msk_sumarr(grd%img, grd%jmg))
      dxdyl_sumarr(:,:) = 0.0_r8
      msk_sumarr(:,:)   = 0.0_r8
      dxdyl_sumarr( grd%igs:grd%ige, grd%jgs:grd%jge) = grd%dxdy(1:grd%im,1:grd%jm)
      msk_sumarr  (grd%igs:grd%ige, grd%jgs:grd%jge)  = grd%msk(1:grd%im,1:grd%jm,1)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, dxdyl_sumarr, grd%img*grd%jmg, mpi%r8, MPI_SUM, mpi%comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, msk_sumarr, grd%img*grd%jmg, mpi%r8, MPI_SUM, mpi%comm, ierr)
      nwp      = SUM( RESHAPE(msk_sumarr,(/grd%img*grd%jmg/) ) )
      dxdyl    = SUM( RESHAPE(dxdyl_sumarr,(/grd%img*grd%jmg/) ), MASK=RESHAPE(msk_sumarr,(/grd%img*grd%jmg/)) .EQ. 1 )
      DEALLOCATE(dxdyl_sumarr,msk_sumarr)

   ELSE

      nwp      = SUM( RESHAPE(grd%msk,(/grd%img*grd%jmg/) ) )
      dxdyl    = SUM( RESHAPE(grd%dxdy,(/grd%img*grd%jmg/) ), MASK=RESHAPE(grd%msk(:,:,1),(/grd%img*grd%jmg/)) .EQ. 1 )

   ENDIF


   grd%adxdy     = dxdyl/nwp
   grd%adxdy     = DSQRT(grd%adxdy)
   grd%dxdy(:,:) = DSQRT(grd%dxdy(:,:))

   grd%nps = grd%im*grd%jm

   DO k = 1,grd%km
      grd%ums(1:grd%im+grd%iae-1,1:grd%jm+grd%jae,k) =                      &
         grd%msk(1:grd%im+grd%iae-1,1:grd%jm+grd%jae,k) * grd%msk(2:grd%im+grd%iae,1:grd%jm+grd%jae,k)
      grd%vms(1:grd%im+grd%iae,1:grd%jm+grd%jae-1,k) =                      &
         grd%msk(1:grd%im+grd%iae,1:grd%jm+grd%jae-1,k) * grd%msk(1:grd%im+grd%iae,2:grd%jm+grd%jae,k)
   ENDDO

   ! Define grid for horizontal covariances
   IF ( drv%filter .EQ. 1) THEN
      grd%msr(:,:,:) = 0.0_r8
      IF ( drv%mask(drv%ktr) .EQ. 1 ) THEN
         grd%msr(:,:,:) = 1.0
      ELSEIF ( drv%mask(drv%ktr) .EQ. 2 ) THEN
         DO k = 1,grd%km
            grd%msr(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,k) =    &
               grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,1)
         ENDDO
      ELSEIF ( drv%mask(drv%ktr) .EQ. 3 ) THEN
         grd%msr(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,:) =    &
            grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,:)
      ELSE
         WRITE (drv%dia,*)'Wrong mask for horizontal covariances ',  &
            drv%mask(drv%ktr)
         STOP
      ENDIF

      IF ( mpi%nproc .GT. 1 .AND. (drv%mask(drv%ktr) .EQ. 2 .OR. drv%mask(drv%ktr) .EQ. 3) ) THEN
         CALL exa_mpi ( 1_i4, 2_i4, grd%im-1, 1-2, grd%im+2,   &
                              2_i4, grd%jm-1, 1-2, grd%jm+2,   &
                              1-3,grd%im+3,1-3,grd%jm+3, grd%km, grd%msr)
         CALL exa_mpi ( 1_i4, 3_i4, grd%im-2, 1-3, grd%im+3,   &
                              3_i4, grd%jm-2, 1-3, grd%jm+3,   &
                              1-3,grd%im+3,1-3,grd%jm+3, grd%km, grd%msr)
      ENDIF  
   ENDIF

   DO j = 1-grd%jas,grd%jm+grd%jae      ! 1,grd%jm
      DO i = 1-grd%ias,grd%im+grd%iae   ! 1,grd%im
         grd%f(i,j) = phy%omega2 *                &
            SIN(grd%lat(i,j)*phy%d2r)
      ENDDO
   ENDDO

   ich = 0
   DO j = 1,grd%jm+grd%jae
      DO i = 1,grd%im+grd%iae-1
         IF ( grd%lon(i+1,j)-grd%lon(i,j) .LT. -300. ) ich = 1
      ENDDO
   ENDDO
   IF ( ich.EQ.1 ) THEN
      DO j = 1-grd%jas,grd%jm+grd%jae
         DO i = 1-grd%ias,grd%im+grd%iae
            IF ( grd%lon(i,j).LT.0.0_r8 ) grd%lon(i,j) = grd%lon(i,j) + 360.
         ENDDO
      ENDDO
   ENDIF

   grd%bsth =  1.e20_r8
   grd%bwst =  1.e20_r8
   grd%bnrt = -1.e20_r8
   grd%beas = -1.e20_r8
   DO j = 1,grd%jm+grd%jae
      DO i = 1,grd%im+grd%iae
         grd%bsth = MIN(grd%bsth,grd%lat(i,j))
         grd%bnrt = MAX(grd%bnrt,grd%lat(i,j))
         grd%bwst = MIN(grd%bwst,grd%lon(i,j))
         grd%beas = MAX(grd%beas,grd%lon(i,j))
      ENDDO
   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, grd%bsth, 1, mpi%r8, MPI_MIN, mpi%comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, grd%bnrt, 1, mpi%r8, MPI_MAX, mpi%comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, grd%bwst, 1, mpi%r8, MPI_MIN, mpi%comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, grd%beas, 1, mpi%r8, MPI_MAX, mpi%comm, ierr)
   ENDIF

!Broadcast first gridpoints for int_obs_hor.F90
   grd%lon1_1 = grd%lon(1,1)
   grd%lat1_1 = grd%lat(1,1)
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL mpi_bcast(grd%lon1_1,1,mpi%r8,0,mpi%comm,ierr)
      CALL mpi_bcast(grd%lat1_1,1,mpi%r8,0,mpi%comm,ierr)
   ENDIF

END SUBROUTINE def_grd

