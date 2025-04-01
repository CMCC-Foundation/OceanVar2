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
!> Apply diffusion filter                                                     
!                                                                      !
! Version 1: Mario Adani      2023                                     !
!            Francesco Carere 2023                                     !
!-----------------------------------------------------------------------
SUBROUTINE diffusive_filter

   USE set_knd
   USE grd_str
   USE dfl_str

   IMPLICIT NONE

   INTEGER(i4) :: i,j,k
   REAL(r8)    :: msk9

   DO j = 1-grd%jas, grd%jm+grd%jae
      DO i = 1-grd%ias, grd%im+grd%iae
         grd%eta(i,j) = grd%eta(i,j) * dfl%wgh(i,j,1)
      ENDDO
   ENDDO

   DO k = 1,grd%km
      DO j = 1-grd%jas, grd%jm+grd%jae
         DO i = 1-grd%ias, grd%im+grd%iae
            grd%tem(i,j,k) = grd%tem(i,j,k) * dfl%wgh(i,j,k)
            grd%sal(i,j,k) = grd%sal(i,j,k) * dfl%wgh(i,j,k)
         ENDDO
      ENDDO
   ENDDO

! Diffusion 2d
   CALL dif_flt(1,grd%eta)

! Diffusion 3d
   DO k = 1,grd%km
      CALL dif_flt(k,grd%tem(:,:,k))
      CALL dif_flt(k,grd%sal(:,:,k))
   ENDDO

   !Mask first grid point close to the coast
   DO k = 1,grd%km
      DO j = 2-grd%jas, grd%jm-1+grd%jae
         DO i = 2-grd%ias, grd%im-1+grd%iae
            msk9 =  grd%msk(i,j,k) + grd%msk(i+1,j,k)+grd%msk(i,j+1,k)+grd%msk(i+1,j+1,k)     &
                                   + grd%msk(i-1,j,k)+grd%msk(i,j-1,k)+grd%msk(i-1,j-1,k)     &
                                   + grd%msk(i+1,j-1,k)+ grd%msk(i-1,j+1,k)
            IF ( msk9 .NE. 9._r8 ) THEN
               grd%tem(i,j,k) = 0._r8
               grd%sal(i,j,k) = 0._r8
               IF ( msk9 .NE. 9._r8 .AND. k .EQ. 1 ) grd%eta(i,j) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   ENDDO


END SUBROUTINE diffusive_filter
!-----------------------------------------------------------------------
!                                                                      !
!> Diffusion filter                                                     
!                                                                      !
! Version 1: Mario Adani      2023                                     !
!            Francesco Carere 2023                                     !
!-----------------------------------------------------------------------
SUBROUTINE dif_flt(k,fld)

   USE grd_str
   USE dfl_str
   USE mpi_str

   IMPLICIT NONE

   REAL(r8), POINTER :: dx(:),ax(:),ex(:),fx(:),vx(:)
   REAL(r8), POINTER :: dy(:),ay(:),ey(:),fy(:),vy(:)
   REAL(r8)          :: fld(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae)
   INTEGER(i4)       :: i,j,it
   INTEGER(i4)       :: k, nx, ny

   ALLOCATE ( dx(1-grd%ias:grd%im+grd%iae),&
              ax(1-grd%ias:grd%im+grd%iae),&
              ex(1-grd%ias:grd%im+grd%iae),&
              fx(1-grd%ias:grd%im+grd%iae),&
              vx(1-grd%ias:grd%im+grd%iae) )
   ALLOCATE ( dy(1-grd%jas:grd%jm+grd%jae),&
              ay(1-grd%jas:grd%jm+grd%jae),&
              ey(1-grd%jas:grd%jm+grd%jae),&
              fy(1-grd%jas:grd%jm+grd%jae),&
              vy(1-grd%jas:grd%jm+grd%jae) )

   nx = size(dx)
   ny = size(dy)

   DO it = 1, dfl%nt
      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, fld )
      ENDIF
!LONGITUDE
      DO j = 1-grd%jas,grd%jm+grd%jae
         dx(:) = fld(:,j)
         ax(:) = dfl%Ax(:,j,k)
         ex(:) = dfl%Ux(:,j,k)
         fx(:) = dfl%Lx(:,j,k)
         CALL tridiagLUSolve(1,nx,dx(:),ax(:),ex(:),fx(:),vx(:))
         fld(:,j) = vx(:)
      ENDDO
!LATITUDE
      DO i = 1-grd%ias,grd%im+grd%iae
         dy(:) = fld(i,:)
         ay(:) = dfl%Ay(i,:,k)
         ey(:) = dfl%Uy(i,:,k)
         fy(:) = dfl%Ly(i,:,k)
         CALL tridiagLUSolve(0,ny,dy(:),ay(:),ey(:),fy(:),vy(:))
         fld(i,:) = vy(:)
      ENDDO
   ENDDO ! nt
   IF ( mpi%nproc .GT. 1) THEN
      CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, fld )
   ENDIF

   DEALLOCATE ( dx,ax,ex,fx,vx )
   DEALLOCATE ( dy,ay,ey,fy,vy )

END SUBROUTINE dif_flt
