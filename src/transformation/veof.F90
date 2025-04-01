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
!> Vertical transformation                                              
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE veof

   USE set_knd
   USE drv_str
   USE grd_str
   USE eof_str

   IMPLICIT NONE

   INTEGER(i4)            :: i, j, k, l, n
   LOGICAL                :: eta_condition
   REAL(r8), ALLOCATABLE  :: egm(:,:)

   ALLOCATE ( egm( grd%im, grd%jm) )

   grd%eta(:,:  ) = 0.0_r8
   grd%tem(:,:,:) = 0.0_r8
   grd%sal(:,:,:) = 0.0_r8

   eta_condition = ( (drv%bmd(drv%ktr)+drv%bal(drv%ktr) .EQ. 0)                      .OR.  &
                     ((drv%bal(drv%ktr) .EQ. 1) .AND. (drv%ssh_unbalanced(drv%ktr))) )


!cdir noconcur
   DO n=1,ros%neof

      egm(:,:) = 0.0_r8

! egm
      DO l = 1,n
         DO j = 1,grd%jm
            DO i = 1,grd%im
#ifdef opt_huge_memory
               egm(i,j) = egm(i,j) + ros%cor( i, j, l, n) * grd%ro(i,j,l)
#else
               egm(i,j) = egm(i,j) + ros%cor(grd%reg(i,j), l, n) * grd%ro(i,j,l)
#endif
            ENDDO
         ENDDO
      ENDDO

! Eta
      IF (eta_condition) THEN
         DO j = 1,grd%jm
            DO i = 1,grd%im
#ifdef opt_huge_memory
               grd%eta(i,j) = grd%eta(i,j) + ros%evc( i, j,1,n) * egm(i,j)
#else
               grd%eta(i,j) = grd%eta(i,j) + ros%evc(grd%reg(i,j),1,n) * egm(i,j)
#endif
            ENDDO
         ENDDO
      ENDIF

! 3D variables
      DO k = 1,grd%km
         DO j = 1,grd%jm
            DO i = 1,grd%im
#ifdef opt_huge_memory
               grd%tem(i,j,k) = grd%tem(i,j,k) + ros%evc( i, j, k+1       , n)  * egm(i,j) * grd%lcl(i,j,k)
               grd%sal(i,j,k) = grd%sal(i,j,k) + ros%evc( i, j, k+grd%km+1, n)  * egm(i,j) * grd%lcl(i,j,k)
#else
               grd%tem(i,j,k) = grd%tem(i,j,k) + ros%evc(grd%reg(i,j), k+1      ,n) * egm(i,j) * grd%lcl(i,j,k)
               grd%sal(i,j,k) = grd%sal(i,j,k) + ros%evc(grd%reg(i,j),k+grd%km+1,n) * egm(i,j) * grd%lcl(i,j,k)
#endif
            ENDDO
         ENDDO
      ENDDO

   ENDDO !n

   DEALLOCATE ( egm )

END SUBROUTINE veof
