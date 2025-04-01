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
!> Initialise the barotropic model                                   
!!
!! It computes masks, geometric parameters and 
!! constants for over-relaxation algorithm
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
! Version 1.1: Paolo Oddo    2014                                      !
!-----------------------------------------------------------------------
SUBROUTINE ini_bmd

   USE set_knd
   USE drv_str
   USE grd_str
   USE bmd_str
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)    :: i, j, is, js, ie, je
   REAL(r8)       :: b1
   INTEGER        :: ierr

   bmd%nstp  = INT( 24._r8 * 3600._r8 / bmd%dt )
   bmd%nstps = INT( bmd%ndy * bmd%nstp )
   bmd%nstpa = INT( bmd%ady * bmd%nstp )
   bmd%alp2  = 1.0_r8 - bmd%alp1


   bmd%df1 = bmd%fc1 * grd%adxdy**2
   bmd%df2 = bmd%fc2 * grd%adxdy**2

   WRITE (drv%dia,*)' ---- Inizialization barotropic model:'
   WRITE (drv%dia,*)' Time step: [sec]                     ',( bmd%dt/2._r8)
   WRITE (drv%dia,*)' Diffusion parameter1: [m**2]         ',bmd%df1
   WRITE (drv%dia,*)' Diffusion parameter2: [m**2]         ',bmd%df2
   WRITE (drv%dia,*)' Diffusion parameter1: [m**2/sec]     ',bmd%df1 /( bmd%dt/2._r8)
   WRITE (drv%dia,*)' Diffusion parameter2: [m**2/sec]     ',bmd%df2 /( bmd%dt/2._r8)
   WRITE (drv%dia,*)' ----'

   is=1-grd%ias; ie=grd%im+grd%iae
   js=1-grd%jas; je=grd%jm+grd%jae

   ALLOCATE ( bmd%itr(bmd%nstps) )
   ALLOCATE ( bmd%mst(is:ie,js:je) )
   ALLOCATE ( bmd%hgt(is:ie,js:je) )
   ALLOCATE ( bmd%msu(1:ie, js:je) )
   ALLOCATE ( bmd%hgu(1:ie, js:je) )
   ALLOCATE ( bmd%dxu(1:ie, js:je) )
   ALLOCATE ( bmd%msv(is:ie, 1:je) )
   ALLOCATE ( bmd%hgv(is:ie, 1:je) )
   ALLOCATE ( bmd%dyv(is:ie, 1:je) )
   ALLOCATE ( bmd%dyu(is:ie,js:je) )
   ALLOCATE ( bmd%dxv(is:ie,js:je) )
   ALLOCATE ( bmd%a1 (is+1:ie-1, js+1:je-1) )
   ALLOCATE ( bmd%a2 (is+1:ie-1, js+1:je-1) )
   ALLOCATE ( bmd%a3 (is+1:ie-1, js+1:je-1) )
   ALLOCATE ( bmd%a4 (is+1:ie-1, js+1:je-1) )
   ALLOCATE ( bmd%a0 (is+1:ie-1, js+1:je-1) )
   ALLOCATE ( bmd%bxby(is:ie, js:je) )
   ALLOCATE ( bmd%rgh(is:ie, js:je) )
   ALLOCATE ( bmd%etb(is:ie, js:je) )
   ALLOCATE ( bmd%ub (is:ie, js:je) )
   ALLOCATE ( bmd%vb (is:ie, js:je) )
   ALLOCATE ( bmd%un (is:ie, js:je) )
   ALLOCATE ( bmd%vn (is:ie, js:je) )
   ALLOCATE ( bmd%eta(is:ie, js:je) )
   ALLOCATE ( bmd%ua (is:ie, js:je) )
   ALLOCATE ( bmd%va (is:ie, js:je) )
   ALLOCATE ( bmd%etm(is:ie, js:je) )
   ALLOCATE ( bmd%um (is:ie, js:je) )
   ALLOCATE ( bmd%vm (is:ie, js:je) )
   ALLOCATE ( bmd%div(is:ie, js:je) )
   ALLOCATE ( bmd%cu (is:ie, js:je) )
   ALLOCATE ( bmd%cv (is:ie, js:je) )
   ALLOCATE ( bmd%dux(is:ie, js:je) )
   ALLOCATE ( bmd%duy(is:ie, js:je) )
   ALLOCATE ( bmd%dvx(is:ie, js:je) )
   ALLOCATE ( bmd%dvy(is:ie, js:je) )
   ALLOCATE ( bmd%etx(is:ie, js:je) )
   ALLOCATE ( bmd%ety(is:ie, js:je) )


   bmd%mst(:,:) = 0.0_r8
   bmd%msu(:,:) = 0.0_r8
   bmd%msv(:,:) = 0.0_r8
   bmd%hgt(:,:) = 0.0_r8
   bmd%hgu(:,:) = 0.0_r8
   bmd%hgv(:,:) = 0.0_r8
   bmd%dxu(:,:) = 0.0_r8
   bmd%dyu(:,:) = 0.0_r8
   bmd%dxv(:,:) = 0.0_r8
   bmd%dyv(:,:) = 0.0_r8
   bmd%a1 (:,:) = 0.0_r8
   bmd%a2 (:,:) = 0.0_r8
   bmd%a3 (:,:) = 0.0_r8
   bmd%a4 (:,:) = 0.0_r8
   bmd%a0 (:,:) = 0.0_r8

   bmd%itr(:) = 0


   DO j=1-grd%jas,grd%jm+grd%jae
      DO i=1-grd%ias,grd%im+grd%iae
         IF (grd%hgt(i,j).GT.0.0_r8) THEN
            bmd%hgt(i,j) = grd%hgt(i,j)
         ENDIF
      ENDDO
   ENDDO

   IF (mpi%ir.EQ.1      ) bmd%hgt(1,:)      = 0.0_r8
   IF (mpi%jr.EQ.1      ) bmd%hgt(:,1)      = 0.0_r8
   IF (mpi%ir.EQ.mpi%irm) bmd%hgt(grd%im,:) = 0.0_r8
   IF (mpi%jr.EQ.mpi%jrm) bmd%hgt(:,grd%jm) = 0.0_r8

!-------------------------------------------------

   DO j=1-grd%jas,grd%jm+grd%jae
      DO i=1-grd%ias,grd%im+grd%iae
         IF (bmd%hgt(i,j).GT.0.0_r8) THEN
            bmd%mst(i,j) = 1.0_r8
         ENDIF
      ENDDO
   ENDDO

   bmd%bnm = 0.0_r8
   DO j=2-grd%jas,grd%jm-1+grd%jae
      DO i=2-grd%ias,grd%im-1+grd%iae
         IF (bmd%mst(i,j).EQ.1.) bmd%bnm = bmd%bnm + 1._r8
      ENDDO
   ENDDO

   IF (mpi%nproc.GT.1) THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, bmd%bnm, 1, mpi%r8, MPI_SUM, mpi%comm, ierr)
   ENDIF

   DO j=1-grd%jas,grd%jm+grd%jae
      DO i=2-grd%ias,grd%im+grd%iae
         bmd%hgu(i,j) = MIN(bmd%hgt(i,j),bmd%hgt(i-1,j))
         bmd%dxu(i,j) = (grd%dx(i,j)+grd%dx(i-1,j))*0.5_r8
         bmd%dyu(i,j) = (grd%dy(i,j)+grd%dy(i-1,j))*0.5_r8
         bmd%msu(i,j) =  bmd%mst(i,j)*bmd%mst(i-1,j)
      ENDDO
   ENDDO
   IF (mpi%ir.EQ.1 ) bmd%dxu(1,:) = bmd%dxu(2,:)
   IF (mpi%ir.EQ.1 ) bmd%dyu(1,:) = bmd%dyu(2,:)
   IF (mpi%ir.EQ.1 ) bmd%msu(1,:) = 0.0_r8

   IF (mpi%nproc.GT.1) THEN
      CALL eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
         1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%dyu)
   ENDIF

   DO j=2-grd%jas,grd%jm+grd%jae
      DO i=1-grd%ias,grd%im+grd%iae
         bmd%hgv(i,j) = MIN(bmd%hgt(i,j),bmd%hgt(i,j-1))
         bmd%dxv(i,j) = (grd%dx(i,j)+grd%dx(i,j-1))*0.5_r8
         bmd%dyv(i,j) = (grd%dy(i,j)+grd%dy(i,j-1))*0.5_r8
         bmd%msv(i,j) =  bmd%mst(i,j)*bmd%mst(i,j-1)
      ENDDO
   ENDDO
   IF (mpi%jr.EQ.1 ) bmd%dxv(:,1) = bmd%dxv(:,2)
   IF (mpi%jr.EQ.1 ) bmd%dyv(:,1) = bmd%dyv(:,2)
   IF (mpi%jr.EQ.1 ) bmd%msv(:,1) = 0.0_r8

   IF (mpi%nproc.GT.1) THEN
      CALL eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
         1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%dxv)
   ENDIF

   b1 =  bmd%alp1**2 *(bmd%dt**2)*phy%g

   DO j= 2-grd%jas,grd%jm-1+grd%jae
      DO i= 2-grd%ias,grd%im-1+grd%iae
         bmd%a1(i,j) = b1*bmd%hgu(i+1,j)/grd%dx(i,j)**2*bmd%msu(i+1,j)
         bmd%a2(i,j) = b1*bmd%hgu(i  ,j)/grd%dx(i,j)**2*bmd%msu(i  ,j)
         bmd%a3(i,j) = b1*bmd%hgv(i,j+1)/grd%dy(i,j)**2*bmd%msv(i,j+1)
         bmd%a4(i,j) = b1*bmd%hgv(i,j  )/grd%dy(i,j)**2*bmd%msv(i,j  )
         bmd%a0(i,j) = (bmd%a1(i,j)+bmd%a2(i,j)+bmd%a3(i,j)+bmd%a4(i,j) +1.0_r8)
      ENDDO
   ENDDO

END SUBROUTINE ini_bmd
