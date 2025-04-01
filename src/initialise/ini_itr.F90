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
!> Initialize Iteration                                                
!!
!! It initialize the iteration for the minimization procedure
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE ini_itr

   USE set_knd
   USE drv_str
   USE grd_str
   USE eof_str
   USE ctl_str

   IMPLICIT NONE

   INTEGER(i4)               :: i, j, k , ii, jj, kk
   REAL(r8)                  :: ri, rj, p, q
   REAL(r8)                  :: div_x, div_y
   REAL(r8),    ALLOCATABLE  :: pq1(:,:), pq2(:,:), pq3(:,:), pq4(:,:)
   INTEGER(i4), ALLOCATABLE  :: i1(:,:), j1(:,:)

   ALLOCATE ( pq1(grd%im,grd%jm) )
   ALLOCATE ( pq2(grd%im,grd%jm) )
   ALLOCATE ( pq3(grd%im,grd%jm) )
   ALLOCATE ( pq4(grd%im,grd%jm) )
   ALLOCATE (  i1(grd%im,grd%jm) )
   ALLOCATE (  j1(grd%im,grd%jm) )

! ---
! Interpolate between grids
   DO jj = 1,grd%jm
      DO ii = 1,grd%im
         ri = MAX(1.,MIN(REAL(drv%im-1),REAL(ii-1)/REAL(drv%ratio(drv%ktr)) + 1.))
         i  = INT(ri)
         p  = ri-i
         rj = MAX(1.,MIN(REAL(drv%jm-1),REAL(jj-1)/REAL(drv%ratio(drv%ktr)) + 1.))
         j  = INT(rj)
         q  = rj-j

         i1(ii,jj) = i
         j1(ii,jj) = j

         div_y =  (1.0_r8-q) * MAX(drv%msk(i,j  ),drv%msk(i+1,j  ))      &
                     +    q  * MAX(drv%msk(i,j+1),drv%msk(i+1,j+1))
         div_x =  (1.0_r8-p) * drv%msk(i  ,j) + p * drv%msk(i+1,j)
         pq1(ii,jj) = drv%msk(i,j)                                   &
                      * MAX(drv%msk(i,j),drv%msk(i+1,j))             &
                      * (1.0_r8-p) * (1.0_r8-q)                      &
                      /( div_x * div_y + 1.d-16 )
         pq2(ii,jj) = drv%msk(i+1,j)                                 &
                      * MAX(drv%msk(i,j),drv%msk(i+1,j))             &
                      *     p  * (1.0_r8-q)                          &
                      /( div_x * div_y + 1.d-16 )
         div_x =  (1.0_r8-p) * drv%msk(i  ,j+1) + p * drv%msk(i+1,j+1)
         pq3(ii,jj) = drv%msk(i,j+1)                                 &
                      * MAX(drv%msk(i,j+1),drv%msk(i+1,j+1))         &
                      * (1.0_r8-p) *     q                           &
                      /( div_x * div_y + 1.d-16 )
         pq4(ii,jj) = drv%msk(i+1,j+1)                               &
                      * MAX(drv%msk(i,j+1),drv%msk(i+1,j+1))         &
                      *     p  *     q                               &
                      /( div_x * div_y + 1.d-16 )

      ENDDO
   ENDDO

   DO k = 1,ros%neof
      DO jj = 1,grd%jm
         DO ii = 1,grd%im
            i = i1(ii,jj)
            j = j1(ii,jj)
            grd%ro(ii,jj,k) = pq1(ii,jj) * drv%ro(i,j  ,k) + pq2(ii,jj) * drv%ro(i+1,j  ,k)  &
                            + pq3(ii,jj) * drv%ro(i,j+1,k) + pq4(ii,jj) * drv%ro(i+1,j+1,k)
         ENDDO
      ENDDO
   ENDDO

! ---
! Reconstruct the control vector
   kk = 0
   DO k = 1,ros%neof
      DO j = 1,grd%jm
         DO i = 1,grd%im
            kk = kk+1
            ctl%x_c(kk) = grd%ro(i,j,k)/drv%ratio(drv%ktr)
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE ( drv%ro, drv%ro_ad, drv%msk)
   DEALLOCATE ( pq1 )
   DEALLOCATE ( pq2 )
   DEALLOCATE ( pq3 )
   DEALLOCATE ( pq4 )
   DEALLOCATE (  i1 )
   DEALLOCATE (  j1 )

! Calculate the cost function and its gradient

   CALL costf

END SUBROUTINE ini_itr
