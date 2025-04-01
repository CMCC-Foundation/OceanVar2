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
!> Define localization parameter                                       
!!
!! It localizes the correction if rcf%loc > rcf%L in namelist. 
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !  
!-----------------------------------------------------------------------
SUBROUTINE def_loc

   USE set_knd
   USE cns_str
   USE grd_str
   USE obs_str
   USE mpi_str
   USE drv_str
   USE netcdf

   IMPLICIT NONE

   INTEGER               :: ierr, i, j, k, iter, ii, jj, iext
   INTEGER(i4)           :: img, jmg
   REAL(r8), ALLOCATABLE :: loct(:,:)
   REAL(r8), ALLOCATABLE :: locs(:,:)
   REAL(r8), ALLOCATABLE ::  lon(:,:)
   REAL(r8), ALLOCATABLE ::  lat(:,:)
   REAL(r8), ALLOCATABLE ::  msk(:,:)
   REAL(r8)              :: dxx, dyy, dst

   IF ( rcf%loc .GT. rcf%L ) THEN

      IF ( mpi%myrank .EQ. 0 ) THEN
         img = grd%img
         jmg = grd%jmg
      ELSE
         img = 1
         jmg = 1
      ENDIF

      ALLOCATE ( loct(img,jmg) )
      ALLOCATE ( locs(img,jmg) )
      ALLOCATE (  lon(img,jmg) )
      ALLOCATE (  lat(img,jmg) )
      ALLOCATE (  msk(img,jmg) )

! --------
! Ajoint of bservational operators
      grd%eta_ad(:,:  ) = 0.0_r8
      grd%tem_ad(:,:,:) = 0.0_r8
      grd%sal_ad(:,:,:) = 0.0_r8
      grd%uvl_ad(:,:,:) = 0.0_r8
      grd%vvl_ad(:,:,:) = 0.0_r8
      obs%gra(:) = obs%res(:)
      CALL obsop_ad

! Localization
      grd%loc(:,:) = 0.0_r8
      DO k = 1,grd%km
         DO j = 1,grd%jm
            DO i = 1,grd%im
               IF ( grd%tem_ad(i,j,k) .NE. 0.0_r8 ) grd%loc(i,j) = 1.0_r8
               IF ( grd%sal_ad(i,j,k) .NE. 0.0_r8 ) grd%loc(i,j) = 1.0_r8
               IF ( grd%uvl_ad(i,j,k) .NE. 0.0_r8 ) grd%loc(i,j) = 1.0_r8
               IF ( grd%vvl_ad(i,j,k) .NE. 0.0_r8 ) grd%loc(i,j) = 1.0_r8
            ENDDO
         ENDDO
      ENDDO
      DO j = 1,grd%jm
         DO i = 1,grd%im
            IF ( grd%eta_ad(i,j) .NE. 0.0_r8 ) grd%loc(i,j) = 1.0_r8
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL gth_mpi( img, jmg, 1_i4, 1_i4, grd%loc, loct)
         CALL gth_mpi( img, jmg, 1_i4, 1_i4, grd%lon,  lon)
         CALL gth_mpi( img, jmg, 1_i4, 1_i4, grd%lat,  lat)
         CALL gth_mpi( img, jmg, 1_i4, 1_i4,                &
                        grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,1),  msk)
      ELSE
         loct(:,:) = grd%loc(:,:)
         lon(:,:)  = grd%lon(:,:)
         lat(:,:)  = grd%lat(:,:)
         msk(:,:)  = grd%msk(:,:,1)
      ENDIF

      IF ( mpi%myrank .EQ. 0 ) THEN

         WRITE (drv%dia,*) 'heavy computation'
         locs(:,:) = loct(:,:)
         DO j = 1,jmg-1
            DO i = 1,img-1
               IF ( locs(i,j) .EQ. 0.0_r8 .AND. msk(i,j) .EQ. 1.0_r8 ) THEN
                  dxx = phy%re*phy%d2r * (lon(i,j)-lon(i+1,j  )) * COS(lat(i,j)*phy%d2r)
                  dyy = phy%re*phy%d2r * (lat(i,j)-lat(i+1,j  ))
                  dst = DSQRT(dxx**2+dyy**2)
                  dxx = phy%re*phy%d2r * (lon(i,j)-lon(i  ,j+1)) * COS(lat(i,j)*phy%d2r)
                  dyy = phy%re*phy%d2r * (lat(i,j)-lat(i  ,j+1))
                  dst = MIN(dst,DSQRT(dxx**2+dyy**2))
                  iext = int(4.*rcf%loc*0.001_r8) / dst
                  DO jj = MAX(j-iext,1),MIN(j+iext,jmg)
                     DO ii = MAX(i-iext,1),MIN(i+iext,img)
                        IF ( loct(ii,jj) .EQ. 1.0_r8 ) THEN
                           dxx = phy%re*phy%d2r * (lon(i,j)-lon(ii,jj)) * COS(lat(i,j)*phy%d2r)
                           dyy = phy%re*phy%d2r * (lat(i,j)-lat(ii,jj))
                           dst = DSQRT((dxx**2+dyy**2)/(2.*(rcf%loc*0.001_r8)**2))
                           locs(i,j) = MAX(locs(i,j),DBLE(EXP(-dst)))
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         loct(:,:) = locs(:,:)

         OPEN  (121,FORM='unformatted',FILE='locs.dat')
         WRITE (121) loct
         CLOSE (121)

      ENDIF

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL gta_mpi( 1_i4, img, jmg, 1_4, 1_i4, grd%loc, loct)
      ELSE
         grd%loc(:,:) = loct(:,:)
      ENDIF

      DEALLOCATE ( loct )
      DEALLOCATE ( locs )
      DEALLOCATE ( lon  )
      DEALLOCATE ( lat  )
      DEALLOCATE ( msk  )

      DO j = 1-grd%jas,grd%jm+grd%jae
         DO i = 1-grd%ias,grd%im+grd%iae
            grd%loc(i,j) = MIN(grd%loc(i,j),1.0_r8)
         ENDDO
      ENDDO

   ELSE

      grd%loc(:,:) = 1.0_r8

   ENDIF

END SUBROUTINE def_loc
