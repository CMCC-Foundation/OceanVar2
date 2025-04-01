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
!> Initialize huber norm distribution for quality control
!!
!! It reads the distribution parameters
!! It fills in the limit in the observational vector to weigth
!! the observational error
!!
!                                                                      !
! Version 1: Andrea Storto 2015                                        !
!            Mario Adani   2024                                        !
!-----------------------------------------------------------------------
SUBROUTINE ini_huber

   USE set_knd
   USE obs_str
   USE drv_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                           :: stat, ncid, idvar
   INTEGER(i4)                           :: Nside, Nnorm, Nparam
   INTEGER(i4)                           :: jstart, jend
   REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: lim_xbt, lim_sla, lim_sst, &
                                            lim_tra, lim_trd, lim_vdr, &
                                            lim_gvl
   REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: lim_arg, lim_gld
   REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: rlim_xbt,rlim_sla, rlim_sst, &
                                            rlim_tra, rlim_trd, rlim_vdr, &
                                            rlim_gvl
   REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: rlim_arg, rlim_gld

! Read coefficient -----------------------------------
   stat = NF90_OPEN(drv%inpdir//'/'//huberqc%flname, NF90_NOWRITE, ncid)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
! Get dimensions
   stat = NF90_INQ_DIMID (ncid, 'Nside', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = Nside)
   stat = NF90_INQ_DIMID (ncid, 'Nnorm', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = Nnorm)
   stat = NF90_INQ_DIMID (ncid, 'Nparam', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = Nparam)

   IF ( Nside .NE. 2 .AND. huberqc%asymm ) THEN
      WRITE (drv%dia,*)' Assymetric Huber Norm defined in the NAMELIST but not two limits found'
      CALL abort
   ENDIF
   IF ( Nside .NE. 1 .AND. .NOT. huberqc%asymm ) THEN
      WRITE (drv%dia,*)' Symmetric Huber Norm defined in the NAMELIST but not one limit found'
      CALL abort
   ENDIF

   IF ( Nnorm .NE. 2 .AND. huberqc%l05 ) THEN
      WRITE (drv%dia,*)' L05 Huber Norm defined in the NAMELIST but not two limits found'
      CALL abort
   ENDIF
   IF ( Nnorm .NE. 1 .AND. .NOT. huberqc%asymm ) THEN
      WRITE (drv%dia,*)' L2-L1 ONLY Huber Norm defined in the NAMELIST but not one limit found '
      CALL abort
   ENDIF

! Get variables
   IF ( huberqc%xbt) THEN
      ALLOCATE ( lim_xbt(Nside,Nnorm), rlim_xbt(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_xbt', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_xbt)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_xbt(:,:) = 1.e+16_r8                            ! Default
      rlim_xbt(:,1) = lim_xbt(1,1)                         ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_xbt(2,1) = lim_xbt(2,1)    ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                             ! If Norm L2-L1 + L05
         rlim_xbt(:,2) = lim_xbt(1,2)                      ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_xbt(2,2) = lim_xbt(2,2) ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%arg ) THEN
      ALLOCATE ( lim_arg(Nside,Nnorm,Nparam), rlim_arg(2,2,Nparam) )
      stat = NF90_INQ_VARID (ncid, 'limits_argo', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_arg)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_arg(:,:,:) = 1.e+16_r8                              ! Default
      rlim_arg(1,1,:) = lim_arg(1,1,:)                         ! First symmetric Norm
      rlim_arg(2,1,:) = lim_arg(1,1,:)                         ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_arg(2,1,:) = lim_arg(2,1,:)    ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                                 ! If Norm L2-L1 + L05
         rlim_arg(1,2,:) = lim_arg(1,2,:)                      ! First symmetric Norm for L05
         rlim_arg(2,2,:) = lim_arg(1,2,:)                      ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_arg(2,2,:) = lim_arg(2,2,:) ! IF Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%sla ) THEN
      ALLOCATE ( lim_sla(Nside,Nnorm), rlim_sla(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_sla', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_sla)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_sla(:,:) = 1.e+16_r8                                ! Default
      rlim_sla(:,1) = lim_sla(1,1)                             ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_sla(2,1) = lim_sla(2,1)        ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                                 ! If Norm L2-L1 + L05
         rlim_sla(:,2) = lim_sla(1,2)                          ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_sla(2,2) = lim_sla(2,2)     ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%gld ) THEN
      ALLOCATE ( lim_gld(Nside,Nnorm,Nparam), rlim_gld(2,2,Nparam) )
      stat = NF90_INQ_VARID (ncid, 'limits_gld', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_gld)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_gld(:,:,:) = 1.e+16_r8                              ! Default
      rlim_gld(1,1,:) = lim_gld(1,1,:)                         ! First symmetric Norm
      rlim_gld(2,1,:) = lim_gld(1,1,:)                         ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_gld(2,1,:) = lim_gld(2,1,:)    ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                                 ! If Norm L2-L1 + L05
         rlim_gld(1,2,:) = lim_gld(1,2,:)                      ! First symmetric Norm for L05
         rlim_gld(2,2,:) = lim_gld(1,2,:)                      ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_gld(2,2,:) = lim_gld(2,2,:) ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%tra ) THEN
      ALLOCATE ( lim_tra(Nside,Nnorm), rlim_tra(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_tra', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_tra)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_tra(:,:) = 1.e+16_r8                               ! Default
      rlim_tra(:,1) = lim_tra(1,1)                            ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_tra(2,1) = lim_tra(2,1)       ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                                ! If Norm L2-L1 + L05
         rlim_tra(:,2) = lim_tra(1,2)                         ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_tra(2,2) = lim_tra(2,2)    ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%trd ) THEN
      ALLOCATE ( lim_trd(Nside,Nnorm), rlim_trd(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_trd', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_trd)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_trd(:,:) = 1.e+16_r8                              ! Default
      rlim_trd(:,1) = lim_trd(1,1)                           ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_trd(2,1) = lim_trd(2,1)      ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                               ! If Norm L2-L1 + L05
         rlim_trd(:,2) = lim_trd(1,2)                        ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_trd(2,2) = lim_trd(2,2)   ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%vdr ) THEN
      ALLOCATE ( lim_vdr(Nside,Nnorm), rlim_vdr(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_vdr', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_vdr)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_vdr(:,:) = 1.e+16_r8                              ! Default
      rlim_vdr(:,1) = lim_vdr(1,1)                           ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_vdr(2,1) = lim_vdr(2,1)      ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                               ! If Norm L2-L1 + L05
         rlim_vdr(:,2) = lim_vdr(1,2)                        ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_vdr(2,2) = lim_vdr(2,2)   ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%gvl ) THEN
      ALLOCATE ( lim_gvl(Nside,Nnorm), rlim_gvl(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_gvl', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_gvl)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_gvl(:,:) = 1.e+16_r8                             ! Default
      rlim_gvl(:,1) = lim_gvl(1,1)                          ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_gvl(2,1) = lim_gvl(2,1)     ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                              ! If Norm L2-L1 + L05
         rlim_gvl(:,2) = lim_gvl(1,2)                       ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_gvl(2,2) = lim_gvl(2,2)  ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   IF ( huberqc%sst ) THEN
      ALLOCATE ( lim_sst(Nside,Nnorm), rlim_sst(2,2) )
      stat = NF90_INQ_VARID (ncid, 'limits_sst', idvar)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, lim_sst)
      IF (stat /= NF90_NOERR)  CALL netcdf_err(stat)
      rlim_sst(:,:) = 1.e+16_r8                             ! Default
      rlim_sst(:,1) = lim_sst(1,1)                          ! First symmetric Norm
      IF ( huberqc%asymm ) rlim_sst(2,1) = lim_sst(2,1)     ! If Assymetric Norm
      IF ( huberqc%l05 )  THEN                              ! If Norm L2-L1 + L05
         rlim_sst(:,2) = lim_sst(1,2)                       ! First symmetric Norm for L05
         IF ( huberqc%asymm ) rlim_sst(2,2) = lim_sst(2,2)  ! If Assymetric Norm for L05
      ENDIF
   ENDIF

   stat = NF90_CLOSE(ncid)

   obs%ahub (:,:) = 1.e+16_r8
   obs%ahub2(:,:) = 1.e+16_r8


!SLA
   jstart = 1
   jend = sla%nc
   IF ( huberqc%sla) CALL fill_ahub2d(jstart,jend,rlim_sla)
!ARGO
   jstart = sla%nc + 1
   jend = sla%nc + arg%nc
   IF ( huberqc%arg) CALL fill_ahub3d(arg%nc,arg%par,jstart,jend,rlim_arg)
!XBT
   jstart = sla%nc + arg%nc + 1
   jend = sla%nc + arg%nc + xbt%nc
   IF ( huberqc%xbt) CALL fill_ahub2d(jstart,jend,rlim_xbt)
!GLIDER
   jstart = sla%nc + arg%nc + xbt%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc
   IF ( huberqc%gld) CALL fill_ahub3d(gld%nc,gld%par,jstart,jend,rlim_gld)
!ARGO TRAJECTORIES
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc +     &
          2 * tra%nc
   IF ( huberqc%tra)  CALL fill_ahub2d(jstart,jend,rlim_tra)
!DRIFTERS TRAJECTORIES
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc
   IF ( huberqc%trd) CALL fill_ahub2d(jstart,jend,rlim_trd)
!DRIFTERS VELOCITIES
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc
   IF ( huberqc%vdr) CALL fill_ahub2d(jstart,jend,rlim_vdr)
!GLIDER VELOCITIES
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc
   IF ( huberqc%gvl) CALL fill_ahub2d(jstart,jend,rlim_gvl)
!SST
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + gvl%nc + 1
   jend = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc + sst%nc
   IF ( huberqc%sst) CALL fill_ahub2d(jstart,jend,rlim_sst)

   IF ( ALLOCATED (lim_sst) ) DEALLOCATE (lim_sst, rlim_sst)
   IF ( ALLOCATED (lim_gvl) ) DEALLOCATE (lim_gvl, rlim_gvl)
   IF ( ALLOCATED (lim_vdr) ) DEALLOCATE (lim_vdr, rlim_vdr)
   IF ( ALLOCATED (lim_tra) ) DEALLOCATE (lim_tra, rlim_tra)
   IF ( ALLOCATED (lim_trd) ) DEALLOCATE (lim_trd, rlim_trd)
   IF ( ALLOCATED (lim_gld) ) DEALLOCATE (lim_gld, rlim_gld)
   IF ( ALLOCATED (lim_arg) ) DEALLOCATE (lim_arg, rlim_arg)
   IF ( ALLOCATED (lim_xbt) ) DEALLOCATE (lim_xbt, rlim_xbt)
   IF ( ALLOCATED (lim_sla) ) DEALLOCATE (lim_sla, rlim_sla)

END SUBROUTINE ini_huber
!===============================================================
SUBROUTINE fill_ahub2d(jstart,jend,rlim)

!-----------------------------------------------------------------------
!                                                                      !
!> Fill in the vector for vertical independent observations
!!
!! 
!!
!                                                                      !
! Version 1: Andrea Storto 2015                                        !
!            Mario Adani   2024                                        !
!-----------------------------------------------------------------------

   USE set_knd
   USE obs_str, ONLY :  obs

   IMPLICIT NONE

   INTEGER(i4) :: jstart, jend, jo
   REAL(r8)    :: rlim(2,2)

   DO jo = jstart,jend
      obs%ahub (jo,:) = ABS(rlim(:,1))
      obs%ahub2(jo,:) = ABS(rlim(:,2))
   ENDDO

END SUBROUTINE fill_ahub2d
!===============================================================
SUBROUTINE fill_ahub3d(npar,par,jstart,jend,rlim)

!-----------------------------------------------------------------------
!                                                                      !
!> Fill in the vector for vertical dependent observations
!!
!!
!!
!                                                                      !
! Version 1: Andrea Storto 2015                                        !
!            Mario Adani   2024                                        !
!-----------------------------------------------------------------------

   USE set_knd
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4) :: jstart, jend, jo, kk
   INTEGER(i8) :: npar
   INTEGER(i8) :: par(npar)
   REAL(r8)    :: rlim(2,2,2)

   kk = 0
   DO jo = jstart,jend
      kk = kk + 1
      obs%ahub (jo,:) = ABS(rlim(:,1,par(kk)))
      obs%ahub2(jo,:) = ABS(rlim(:,2,par(kk)))
   ENDDO

END SUBROUTINE fill_ahub3d
