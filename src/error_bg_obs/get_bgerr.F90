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
!> Compute Background error from EOFs 
!!
!! It computes the background error from EOFs.
!! BKG_error = const * sqrt(U**2 * S**2 UT**2)
!!
!                                                                      !
! Version 1:  Andrea Storto 2022                                       !
!             Mario Adani   2024                                       !
!-----------------------------------------------------------------------
SUBROUTINE get_bgerr

   USE set_knd
   USE obs_str
   USE grd_str
   USE eof_str
   USE drv_str
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)          :: i, j, k, n, kbot
   REAL(r8)             :: t, s, tb, sb, rhtl, zrho
!FUNCTION
   REAL(r8)             :: rho_unescotl
   REAL(r8),ALLOCATABLE :: slaeof(:,:)

   ALLOCATE( qck%tem(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
             qck%sal(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
             qck%eta(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae       ) )

   qck%tem(:,:,:) = 0.0_r8
   qck%sal(:,:,:) = 0.0_r8
   qck%eta(:,:  ) = 0.0_r8

! TEMPERATURE / SALINITY
   DO n = 1,ros%neof
      DO k = 1,grd%km
         DO j = 1,grd%jm
            DO i = 1,grd%im
#ifdef opt_huge_memory
               qck%tem(i,j,k) = qck%tem(i,j,k) + ros%evc( i, j, k+1       , n)**2  * ros%eva(i,j,n)**2
               qck%sal(i,j,k) = qck%sal(i,j,k) + ros%evc( i, j, k+grd%km+1, n)**2  * ros%eva(i,j,n)**2
#else
               qck%tem(i,j,k) = qck%tem(i,j,k) + ros%evc(grd%reg(i,j), k+1      ,n)**2 &
                                               * ros%eva(grd%reg(i,j),           n)**2
               qck%sal(i,j,k) = qck%sal(i,j,k) + ros%evc(grd%reg(i,j),k+grd%km+1,n)**2 &
                                               * ros%eva(grd%reg(i,j),           n)**2
#endif
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   qck%tem(:,:,:) = DSQRT(qck%tem(:,:,:))
   qck%sal(:,:,:) = DSQRT(qck%sal(:,:,:))

! SLA 
   IF ( drv%bal(drv%ktr) .EQ. 1 )  THEN
#ifdef opt_huge_memory
      ALLOCATE( grd%im,grd%jm,ros%neof )
      slaeof(:,:,:) = 0._r8
#else
      ALLOCATE( slaeof(ros%nreg,ros%neof) )
      slaeof(:,:) = 0._r8
#endif
      DO j = 1,grd%jm
         DO i = 1,grd%im
            DO n = 1,ros%neof
               kbot = grd%km
               rhtl = 0._r8
               DO k = kbot,1,-1
                  IF ( grd%msk(i,j,k) .LT. 0.9_r8 ) cycle
                  tb      = qck%climtem(i,j,k)
                  sb      = MAX(1.0_r8,qck%climsal(i,j,k))
#ifdef opt_huge_memory
                  t       = ros%evc(i, j,  k+1,n)
                  s       = ros%evc(i, j,  k+grd%km+1, n)
#else
                  t       = ros%evc(grd%reg(i,j), k+1, n)
                  s       = ros%evc(grd%reg(i,j), k+grd%km+1, n)
#endif
                  zrho    = rho_unescotl(sb,tb,s,t)*grd%msk(i,j,k)
                  rhtl    = rhtl + zrho*grd%dz(k)
               ENDDO
#ifdef opt_huge_memory
               slaeof(i,j,n) = -rhtl / phy%rho0
#else
               slaeof(grd%reg(i,j),n) = -rhtl / phy%rho0
#endif
            ENDDO
#ifdef opt_huge_memory
            qck%eta( i, j ) = DSQRT(SUM(slaeof(i, j, :)**2        &
                                     * ros%eva(i, j, :)**2 ))
#else
            qck%eta( i, j ) = DSQRT(SUM(slaeof(grd%reg(i,j),:)**2 &
                                     * ros%eva(grd%reg(i,j),:)**2 ))
#endif
         ENDDO
      ENDDO
   ENDIF

   IF  ( ALLOCATED(slaeof)  ) DEALLOCATE( slaeof )

   IF ( drv%ssh_unbalanced(drv%ktr) .OR. &
        ((drv%bal(drv%ktr) .EQ. 0) .AND. ( drv%bmd(drv%ktr) .EQ. 0)) ) THEN

      DO n = 1,ros%neof
         DO j = 1,grd%jm
            DO i = 1,grd%im
#ifdef opt_huge_memory
               qck%eta(i,j) = qck%eta(i,j) + DSQRT(ros%evc( i, j,        1, n)**2  &
                                                     * ros%eva( i , j ,     n)**2)
#else
               qck%eta(i,j) = qck%eta(i,j) + DSQRT(ros%evc(grd%reg(i,j), 1, n)**2  &
                                                     * ros%eva(grd%reg(i,j),n)**2)
#endif
            ENDDO
         ENDDO
      ENDDO

   ENDIF

END SUBROUTINE get_bgerr
