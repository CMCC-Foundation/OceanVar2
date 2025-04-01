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
!> Driver for computation of observational errors
!!
!! It reads values from files and applyed them to different
!! observational type.
!! If all the flags that specify the errors are false the ones read in 
!! xxx_mis.dat files are used.
!                                                                      !
! Version 1:  Andrea Storto 2022                                       !
!             Mario Adani   2024                                       !
!-----------------------------------------------------------------------
SUBROUTINE obserrors

   USE set_knd
   USE grd_str
   USE obs_str
   USE drv_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                :: stat, ncid, idvar
   INTEGER(i4)                :: km,img,jmg
   INTEGER(i4)                :: start(2), count(2)
   REAL(r4), ALLOCATABLE      :: x1(:),x2(:,:)

!Temperature / Salinity
   IF (obserr%ts_ver_dep_err) THEN

      stat = NF90_OPEN(drv%inpdir//'/'//obserr%ts_flname, NF90_NOWRITE, ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
! Get dimensions
      stat = nf90_inq_dimid (ncid, 'im', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = img)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'jm', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = jmg)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'depth', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = km)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      IF ( ( grd%img .NE. img ) .OR. &
         ( grd%jmg .NE. jmg )      ) THEN
         WRITE (drv%dia,*)' ------------------------------------------ '
         WRITE (drv%dia,*)' Wrong dimension of the ts error grid       '
         WRITE (drv%dia,*)' ------------------------------------------ '
         CALL FLUSH(drv%dia)
         CALL abort
      ENDIF

      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae

      ALLOCATE ( obserr%tem_fc(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae),     &
                 obserr%sal_fc(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae),     &
                 obserr%tem(km), obserr%sal(km),obserr%depth(km) )
      ALLOCATE ( x1(km), x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

! Get Variables
      stat = nf90_inq_varid (ncid, 'tem', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_get_var (ncid, idvar, x1)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      obserr%tem(:) = DBLE(x1(:))
      stat = nf90_inq_varid (ncid, 'sal', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_get_var (ncid, idvar, x1)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      obserr%sal(:) = DBLE(x1(:))
      stat = nf90_inq_varid (ncid, 'depth', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_get_var (ncid, idvar, x1)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      obserr%depth(:) = DBLE(x1(:))
      stat = nf90_inq_varid (ncid, 'tem_fc', idvar)
      IF (stat /= NF90_NOERR) THEN
         WRITE (drv%dia,*)' Multiplication factor for t vertical error set to 1'
         obserr%tem_fc(:,:) = 1._r8
      ELSE
         stat = nf90_get_var (ncid,idvar,x2,start(1:2),count(1:2))
         obserr%tem_fc(:,:) = DBLE(x2(:,:))
      ENDIF
      stat = nf90_inq_varid (ncid, 'sal_fc', idvar)
      IF (stat /= NF90_NOERR) THEN
         WRITE (drv%dia,*)' Multiplication factor for s vertical error set to 1'
         obserr%sal_fc(:,:) = 1._r8
      ELSE
         stat = nf90_get_var (ncid,idvar,x2,start(1:2),count(1:2))
         obserr%sal_fc(:,:) = DBLE(x2(:,:))
      ENDIF
      DEALLOCATE( x1,x2 )

   ENDIF

! SLA
   IF ( obserr%sla_hor_dep_err ) THEN

      ALLOCATE ( obserr%stde(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae), &
                          x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

      stat = NF90_OPEN(drv%inpdir//'/'//obserr%sla_flname, NF90_NOWRITE, ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
! Get dimensions
      stat = nf90_inq_dimid (ncid, 'im', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = img)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'jm', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = jmg)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      IF ( ( grd%img .NE. img ) .OR. &
         ( grd%jmg .NE. jmg )      ) THEN
         WRITE (drv%dia,*)' ------------------------------------------ '
         WRITE (drv%dia,*)' Wrong dimension of the ts error grid       '
         WRITE (drv%dia,*)' ------------------------------------------ '
         CALL FLUSH(drv%dia)
         CALL abort
      ENDIF

      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae

      stat = nf90_inq_varid (ncid, 'stde', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x2,start(1:2),count(1:2))
      obserr%stde(:,:) = DBLE(x2(:,:))
      DEALLOCATE( x2 )

   ENDIF

   WRITE (drv%dia,*) ' -----------------------------------  '
   WRITE (drv%dia,*) ' Computing the obeservational errors  '
   CALL FLUSH(drv%dia)


   IF ( sla%no .GT. 0 .AND. obs%sla .NE. 0 ) THEN
      CALL obserr_sla
   ENDIF

   IF ( arg%no .GT. 0 .AND. obs%arg .NE. 0 ) THEN
      CALL obserr_arg
   ENDIF

   IF ( xbt%no .GT. 0 .AND. obs%xbt .NE. 0 ) THEN
      CALL obserr_xbt
   ENDIF

   IF ( gld%no .GT. 0 .AND. obs%gld .NE. 0 ) THEN
      CALL obserr_gld
   ENDIF

   IF ( tra%no .GT. 0 .AND. obs%tra .NE. 0 ) THEN
      CALL obserr_trd
   ENDIF

   IF ( trd%no .GT. 0 .AND. obs%trd .NE. 0 ) THEN
      CALL obserr_tra
   ENDIF

   IF ( vdr%no .GT. 0 .AND. obs%vdr .NE. 0 ) THEN
      CALL obserr_vdr
   ENDIF

   IF ( gvl%no .GT. 0 .AND. obs%gvl .NE. 0 ) THEN
      CALL obserr_gvl
   ENDIF

   IF ( sst%no .GT. 0 .AND. obs%sst .NE. 0 ) THEN
      CALL obserr_sst
   ENDIF

END SUBROUTINE obserrors
!============================================
REAL(KIND=r8) FUNCTION tim_dep_errfact(time)

!-----------------------------------------------------------------------
!> Function to define the time dependent decay of the error
!!
!! Exponential function based obserr%ztime_weigth defined in namelist
!!
!
! Version 1:  Andrea Storto 2022                                       !
!-----------------------------------------------------------------------


   USE set_knd
   USE obs_str, ONLY : obserr

   IMPLICIT NONE

   REAL(r8) :: time
   REAL(r8) :: td, td2, zt2

   td=MIN(obserr%ztime_weigth,ABS(time))
   td2=td*td
   zt2=obserr%ztime_weigth*obserr%ztime_weigth
   tim_dep_errfact = 1._r8 / EXP( -td2 / zt2 )

END FUNCTION tim_dep_errfact
