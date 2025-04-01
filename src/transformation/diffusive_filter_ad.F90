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
!> Apply Diffusion Filter adjoint                                             
!                                                                      !
! Version 1: Mario Adani      2023                                     !
!            Francesco Carere 2023                                     !
!-----------------------------------------------------------------------
SUBROUTINE diffusive_filter_ad

   USE set_knd
   USE grd_str
   USE dfl_str

   IMPLICIT NONE

   INTEGER(i4) :: i,j,k
   REAL(r8)    :: msk9

   !Mask first grid point close to the coast
   DO k = grd%km, 1, -1
      DO j = grd%jm-1+grd%jae, 2-grd%jas, -1
         DO i = grd%im-1+grd%iae, 2-grd%ias, -1
            msk9 =  grd%msk(i,j,k) + grd%msk(i+1,j,k)+grd%msk(i,j+1,k)+grd%msk(i+1,j+1,k)     &
                                   + grd%msk(i-1,j,k)+grd%msk(i,j-1,k)+grd%msk(i-1,j-1,k)     &
                                   + grd%msk(i+1,j-1,k)+ grd%msk(i-1,j+1,k)
            IF ( msk9 .NE. 9._r8 ) THEN
               grd%tem_ad(i,j,k) = 0._r8
               grd%sal_ad(i,j,k) = 0._r8
               IF ( msk9 .NE. 9._r8 .AND. k .EQ. 1 ) grd%eta_ad(i,j) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   ENDDO

! Diffusion 3d
   DO k = grd%km,1,-1
      CALL dif_flt_ad(k,grd%tem_ad(:,:,k))
      CALL dif_flt_ad(k,grd%sal_ad(:,:,k))
   ENDDO

! Diffusion 2d
   CALL dif_flt_ad(1,grd%eta_ad)


   DO j = grd%jm+grd%jae,1-grd%jas,-1
      DO i = grd%im+grd%iae,1-grd%ias,-1
         grd%eta_ad(i,j) = grd%eta_ad(i,j) * dfl%wgh(i,j,1)
      ENDDO
   ENDDO

   DO k = grd%km,1,-1
      DO j = grd%jm+grd%jae,1-grd%jas,-1
         DO i = grd%im+grd%iae,1-grd%ias,-1
            grd%tem_ad(i,j,k) = grd%tem_ad(i,j,k) * dfl%wgh(i,j,k)
            grd%sal_ad(i,j,k) = grd%sal_ad(i,j,k) * dfl%wgh(i,j,k)
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE dIFfusive_filter_ad
!-----------------------------------------------------------------------
!                                                                      !
!> Diffusion Filter adjoint                                             
!                                                                      !
! Version 1: Mario Adani      2023                                     !
!            Francesco Carere 2023                                     !
!-----------------------------------------------------------------------
SUBROUTINE dif_flt_ad(k, fldb)

   USE grd_str
   USE dfl_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)       :: i, j, it
   REAL(r8), POINTER :: ax(:), ex(:), fx(:)
   REAL(r8), POINTER :: dxb(:), vxb(:)
   REAL(r8), POINTER :: ay(:), ey(:), fy(:)
   REAL(r8), POINTER :: dyb(:), vyb(:)
   REAL(r8)          :: fldb(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae)
   INTEGER(i4)       :: k, nx, ny


   ALLOCATE ( dxb(1-grd%ias:grd%im+grd%iae),&
              ax(1-grd%ias:grd%im+grd%iae), &
              ex(1-grd%ias:grd%im+grd%iae), &
              fx(1-grd%ias:grd%im+grd%iae), &
              vxb(1-grd%ias:grd%im+grd%iae) )
   ALLOCATE ( dyb(1-grd%jas:grd%jm+grd%jae),&
              ay(1-grd%jas:grd%jm+grd%jae), &
              ey(1-grd%jas:grd%jm+grd%jae), &
              fy(1-grd%jas:grd%jm+grd%jae), &
              vyb(1-grd%jas:grd%jm+grd%jae) )

   vxb = 0.0_r8
   vyb = 0.0_r8
   dxb = 0.0_r8
   dyb = 0.0_r8

   nx = size(dxb)
   ny = size(dyb)

   DO it = dfl%nt,1,-1
      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, fldb)
      ENDIF
!LATITUDE
      DO i = grd%im+grd%iae,1-grd%ias,-1
         fy(:) = dfl%Ly(i,:,k)
         ey(:) = dfl%Uy(i,:,k)
         ay(:) = dfl%Ay(i,:,k)
         vyb = vyb + fldb(i, :)
         CALL tridiagLUSolve_adj(0,ny, dyb(:), ay(:), ey(:), fy(:), &
            vyb(:))
         fldb(i, :) = dyb
         dyb = 0.0_r8
      ENDDO
!LONGITUDE
      DO j = grd%jm+grd%jae,1-grd%jas,-1
         fx(:) = dfl%Lx(:,j,k)
         ex(:) = dfl%Ux(:,j,k)
         ax(:) = dfl%Ax(:,j,k)
         vxb = vxb + fldb(:, j)
         CALL tridiagLUSolve_adj(1,nx, dxb(:), ax(:), ex(:), fx(:), &
            vxb(:))
         fldb(:, j) = dxb
         dxb = 0.0_r8
      ENDDO
   ENDDO
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, fldb)
   ENDIF

   DEALLOCATE( ax,ex,fx,dxb,vxb )
   DEALLOCATE( ay,ey,fy,dyb,vyb )

END SUBROUTINE dif_flt_ad
