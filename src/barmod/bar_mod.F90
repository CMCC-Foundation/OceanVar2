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
!> Barotropic model                                                     
!!                                                                     
!! Tangent linear of the barotropic model. 
!! 
!! Ref: S. Dobricic, N. Pinardi. An oceanographic three-dimensional 
!! variational data assimilation scheme. 
!! Ocean Modelling 22 (2008) 89–105.doi:10.1016/j.ocemod.2008.01.004
!!
!
! Version 1: Srdjan Dobricic 2007                                      !
!-----------------------------------------------------------------------
SUBROUTINE bar_mod

   USE set_knd
   USE drv_str
   USE grd_str
   USE bmd_str
   USE mpi_str
   USE cns_str,  ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)    :: i, j, kstp

   IF ( mpi%nproc .GT. 1 ) THEN

      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                    0_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%bx )
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                    0_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%by )
   ENDIF

!---
! Bouyancy forcing
   bmd%bxby(:,:) = 0._r8
   bmd%rgh(:,:)  = 0._r8

   DO j = 2-grd%jas,grd%jm-1+grd%jae                 ! 2,grd%jm-1
      DO i = 2-grd%ias,grd%im-1+grd%iae                ! 2,grd%im-1
         bmd%bxby(i,j) = - (bmd%dt**2)*((grd%bx(i+1,j)-grd%bx(i,j))/grd%dx(i,j) +               &
                                        (grd%by(i,j+1)-grd%by(i,j))/grd%dy(i,j)) * bmd%mst(i,j)
         bmd%rgh(i,j) = bmd%bxby(i,j) * bmd%alp1**2
      ENDDO
   ENDDO

   bmd%etb(:,:) = 0.0_r8
   bmd%ub(:,:) = 0.0_r8
   bmd%vb(:,:) = 0.0_r8
   bmd%un(:,:) = 0.0_r8
   bmd%vn(:,:) = 0.0_r8
   bmd%eta(:,:) = 0.0_r8
   bmd%ua(:,:) = 0.0_r8
   bmd%va(:,:) = 0.0_r8
   bmd%etx(:,:) = 0.0_r8
   bmd%ety(:,:) = 0.0_r8
   bmd%div(:,:) = 0.0_r8
   bmd%cu(:,:) = 0.0_r8
   bmd%cv(:,:) = 0.0_r8
   bmd%dux(:,:) = 0.0_r8
   bmd%dvx(:,:) = 0.0_r8
   bmd%duy(:,:) = 0.0_r8
   bmd%dvy(:,:) = 0.0_r8
   bmd%etm(:,:) = 0.0_r8
   bmd%um(:,:) = 0.0_r8
   bmd%vm(:,:) = 0.0_r8

!---
! Time loop

   DO kstp = 1,bmd%nstps

      DO j = 2-grd%jas,grd%jm-1+grd%jae                ! 2,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae             ! 2,grd%im-1
            bmd%div(i,j) = ((bmd%ub(i+1,j)-bmd%ub(i,j))/grd%dx(i,j) + (bmd%vb(i,j+1)-bmd%vb(i,j))/grd%dy(i,j)) * bmd%mst(i,j)
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm-1+grd%jae      ! 2,grd%jm-1
         DO i = 2-grd%ias,grd%im             ! 2,grd%im
            bmd%etx(i,j) = bmd%alp2*bmd%dt*phy%g*bmd%hgu(i,j)*(bmd%etb(i,j)-bmd%etb(i-1,j  ))/bmd%dxu(i,j) * bmd%msu(i,j)
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm              ! 2,grd%jm
         DO i=2-grd%ias,grd%im-1+grd%iae   ! 2,grd%im-1
            bmd%ety(i,j) = bmd%alp2*bmd%dt*phy%g*bmd%hgv(i,j)*(bmd%etb(i,j)-bmd%etb(i  ,j-1))/bmd%dyv(i,j) * bmd%msv(i,j)
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                                         &
                      -1_i4, 1_i4, grd%im+1,                               &
                       0_i4, 1_i4, grd%jm+1,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%etx)
         CALL exo_mpi( 1_i4, 1_i4,                                         &
                       0_i4, 1_i4, grd%im+1,                               &
                      -1_i4, 1_i4, grd%jm+1,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ety)
      ENDIF

      DO j = 2-grd%jas,grd%jm-1+grd%jae          ! 2,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae       ! 2,grd%im-1
            bmd%rgh(i,j) = bmd%bxby(i,j) - bmd%etb(i,j) + bmd%dt*bmd%div(i,j)           &
                         - bmd%alp1*bmd%dt*(bmd%etx(i+1,j)-bmd%etx(i,j))/grd%dx(i,j)    &
                         - bmd%alp1*bmd%dt*(bmd%ety(i,j+1)-bmd%ety(i,j))/grd%dy(i,j)
         ENDDO
      ENDDO

!---
! Calculate sea level
      CALL invrt( kstp )

      DO j = 1,grd%jm-1+grd%jae                ! 1,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae     ! 2,grd%im-1
            bmd%cu(i,j) = -((bmd%vn(i,j  )*bmd%dxv(i,j  ) + bmd%vn(i-1,j  )*bmd%dxv(i-1,j  ))*grd%f(i,j  ) +     &
                            (bmd%vn(i,j+1)*bmd%dxv(i,j+1) + bmd%vn(i-1,j+1)*bmd%dxv(i-1,j+1))*grd%f(i,j+1))      &
                            * 0.25_r8 / bmd%dxu(i,j)
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm-1+grd%jae      ! 2,grd%jm-1
         DO i = 1,grd%im-1+grd%iae           ! 1,grd%im-1
            bmd%cv(i,j) =  ((bmd%un(i  ,j)*bmd%dyu(i  ,j) + bmd%un(i  ,j-1)*bmd%dyu(i  ,j-1))*grd%f(i  ,j) +      &
                            (bmd%un(i+1,j)*bmd%dyu(i+1,j) + bmd%un(i+1,j-1)*bmd%dyu(i+1,j-1))*grd%f(i+1,j))       &
                            * 0.25_r8 / bmd%dyv(i,j)
         ENDDO
      ENDDO
      DO j = 1,grd%jm                              ! 1,grd%jm
         DO i = 2-grd%ias,grd%im                   ! 2,grd%im
            bmd%dux(i,j) = (bmd%ub(i,j) - bmd%ub(i-1,j))/bmd%dxu(i,j)
            bmd%dvx(i,j) = (bmd%vb(i,j) - bmd%vb(i-1,j))/bmd%dxu(i,j)
         ENDDO
      ENDDO

      DO j = 2-grd%jas,grd%jm                      ! 2,grd%jm
         DO i = 1,grd%im                           ! 1,grd%im
            bmd%duy(i,j) = (bmd%ub(i,j) - bmd%ub(i,j-1))/bmd%dyv(i,j)
            bmd%dvy(i,j) = (bmd%vb(i,j) - bmd%vb(i,j-1))/bmd%dyv(i,j)
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                      -1_i4, 1_i4, grd%im+1,                             &
                       0_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                      -1_i4, 1_i4, grd%im+1,                             &
                       0_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       0_i4, grd%im, 0_i4,                               &
                      -1_i4, 1_i4, grd%jm+1,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       0_i4, grd%im, 0_i4,                               &
                      -1_i4, 1_i4, grd%jm+1,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
      ENDIF

!---
! Calculate new velocity
      DO j = 1,grd%jm-1+grd%jae                    ! 1,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae         ! 2,grd%im
            bmd%ua(i,j) = bmd%ub(i,j) - bmd%dt*(bmd%cu(i,j) +                                                                 &
                                                bmd%alp1*phy%g*bmd%hgu(i,j)*(bmd%eta(i,j)-bmd%eta(i-1,j  ))/bmd%dxu(i,j) +    &
                                                grd%bx(i,j))*bmd%msu(i,j)  -                                                  &
                                        bmd%etx(i,j) +                                                                        &
                                        bmd%df1 * ( (bmd%dux(i+1,j)-bmd%dux(i,j))/bmd%dxu(i,j) +                              &
                                                    (bmd%duy(i,j+1)-bmd%duy(i,j))/bmd%dyu(i,j) ) * bmd%msu(i,j)
         ENDDO
      ENDDO

      DO j = 2-grd%jas,grd%jm-1+grd%jae                    ! 2,grd%jm
         DO i = 1,grd%im-1+grd%iae                         ! 1,grd%im-1
            bmd%va(i,j) = bmd%vb(i,j) - bmd%dt*(bmd%cv(i,j) +                                                                 &
                                                bmd%alp1*phy%g*bmd%hgv(i,j)*(bmd%eta(i,j)-bmd%eta(i  ,j-1))/bmd%dyv(i,j) +    &
                                                grd%by(i,j))*bmd%msv(i,j) -                                                   &
                                        bmd%ety(i,j) +                                                                        &
                                        bmd%df1 * ( (bmd%dvx(i+1,j)-bmd%dvx(i,j))/bmd%dxv(i,j) +                              &
                                                    (bmd%dvy(i,j+1)-bmd%dvy(i,j))/bmd%dyv(i,j) ) * bmd%msv(i,j)
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       1_i4, grd%im, 0_i4,                               &
                       1_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ua)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       1_i4, grd%im, 0_i4,                               &
                       1_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%va)
      ENDIF

      DO j = 1,grd%jm                          ! 1,grd%jm
         DO i = 2-grd%ias,grd%im               ! 2,grd%im
            bmd%dux(i,j) = (bmd%ua(i,j) - bmd%ua(i-1,j))/bmd%dxu(i,j)
            bmd%dvx(i,j) = (bmd%va(i,j) - bmd%va(i-1,j))/bmd%dxv(i,j)
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm                  ! 2,grd%jm
         DO i = 1,grd%im                       ! 1,grd%im
            bmd%duy(i,j) = (bmd%ua(i,j) - bmd%ua(i,j-1))/bmd%dyu(i,j)
            bmd%dvy(i,j) = (bmd%va(i,j) - bmd%va(i,j-1))/bmd%dyv(i,j)
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                      -1_i4, 1_i4, grd%im+1,                             &
                       0_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                      -1_i4, 1_i4, grd%im+1,                             &
                       0_i4, grd%jm, 0_i4,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       0_i4, grd%im, 0_i4,                               &
                      -1_i4, 1_i4, grd%jm+1,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
         CALL exo_mpi( 1_i4, 1_i4,                                       &
                       0_i4, grd%im, 0_i4,                               &
                      -1_i4, 1_i4, grd%jm+1,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
      ENDIF

      DO j = 2-grd%jas,grd%jm-1+grd%jae        ! 2,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae       ! 2,grd%im-1
            bmd%ua(i,j) = bmd%ua(i,j) + bmd%df2 * ((bmd%dux(i+1,j)-bmd%dux(i,j))/bmd%dxu(i,j) +      &
                                                   (bmd%duy(i,j+1)-bmd%duy(i,j))/bmd%dyu(i,j))*bmd%msu(i,j)
         ENDDO
      ENDDO

      DO j = 2-grd%jas,grd%jm-1+grd%jae        ! 2,grd%jm-1
         DO i = 2-grd%ias,grd%im-1+grd%iae     ! 2,grd%im-1
            bmd%va(i,j) = bmd%va(i,j) + bmd%df2 * ((bmd%dvx(i+1,j)-bmd%dvx(i,j))/bmd%dxv(i,j) +      &
                                                   (bmd%dvy(i,j+1)-bmd%dvy(i,j))/bmd%dyv(i,j))*bmd%msv(i,j)
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
         CALL eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%ua )
         CALL eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%va )
      ENDIF

!---
! Asselin filter
      bmd%un(:,:) =  bmd%un(:,:) + ( bmd%ub(:,:) + bmd%ua(:,:) - 2.0_r8*bmd%un(:,:) ) *0.05_r8
      bmd%vn(:,:) =  bmd%vn(:,:) + ( bmd%vb(:,:) + bmd%va(:,:) - 2.0_r8*bmd%vn(:,:) ) *0.05_r8

      bmd%etb(:,:) = bmd%eta(:,:)
      bmd%ub(:,:) =  bmd%un(:,:)
      bmd%vb(:,:) =  bmd%vn(:,:)

      bmd%un(:,:) =  bmd%ua(:,:)
      bmd%vn(:,:) =  bmd%va(:,:)

      bmd%etm(:,:) = bmd%etm(:,:) + bmd%etb(:,:)
      bmd%um(:,:) =  bmd%um(:,:) +  bmd%ub(:,:)
      bmd%vm(:,:) =  bmd%vm(:,:) +  bmd%vb(:,:)

!---
! Temporal average
      IF ( MOD(kstp,bmd%nstpa) .EQ. 0 ) THEN

         bmd%etm(:,:) = bmd%etm(:,:)/DBLE(bmd%nstpa) * bmd%mst(:,:)
         bmd%um(:,:) =  bmd%um(:,:)/DBLE(bmd%nstpa)
         bmd%vm(:,:) =  bmd%vm(:,:)/DBLE(bmd%nstpa)

         grd%eta(:,:) = bmd%etm(:,:)

         bmd%etm(:,:) = 0.0_r8
         bmd%um(:,:) = 0.0_r8
         bmd%vm(:,:) = 0.0_r8

      ENDIF   !Temporal average

   ENDDO   ! kstp

END SUBROUTINE bar_mod
