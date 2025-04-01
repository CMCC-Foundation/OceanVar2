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
!> Read Eofs                                                         
!!
!! It reads empirical orthogonal functions.
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE rdeofs

   USE set_knd
   USE drv_str
   USE eof_str
   USE grd_str
   USE mpi_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                :: stat, ncid, idvar, neofs, nlevs
   INTEGER(i4)                :: nec, k, nrg
   INTEGER(i4)                :: start(4), count(4)
   REAL(r4), ALLOCATABLE      :: x4(:,:,:,:), x3(:,:,:), x2(:,:)
   INTEGER                    :: i, j, ierr, kt, ke
   INTEGER                    :: tmp


   stat = NF90_OPEN(drv%inpdir//'/'//ros%flname, NF90_NOWRITE, ncid)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

! Get dimensions
   stat = NF90_INQ_DIMID (ncid, 'nreg', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = ros%nreg)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'nlev', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = nlevs)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'neof', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = neofs)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

   WRITE (drv%dia,*)' ---- Eof dimensions are: ',ros%nreg, ros%kmt, neofs
   WRITE (drv%dia,*)' ---- Uses ',ros%neof,' eofs.'

   IF ( ros%neof .GT. neofs ) THEN
      WRITE (drv%dia,*)'Error: Requires more Eofs than available in the input file.'
      STOP
   ENDIF

! Read Regions if necessary
   stat = NF90_INQ_VARID (ncid, 'regs', idvar)
   IF (stat == NF90_NOERR) THEN
      IF ( ros%already_read ) THEN
         WRITE (drv%dia,*)'-- WARNING: EOF region already read in grid file             --'
         WRITE (drv%dia,*)'-- We read them again. You may want to check the consictency --'
      ELSE
         ALLOCATE ( grd%reg(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
      ENDIF
      ALLOCATE ( x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae
      stat = NF90_GET_VAR (ncid, idvar, x2,start(1:2),count(1:2))
      grd%reg(:,:) = INT(x2(:,:))
      WHERE( grd%reg(:,:) .EQ. 0 ) grd%reg(:,:) = 1
      DEALLOCATE ( x2 )
   ELSE
      IF ( .NOT. ros%already_read ) THEN
         WRITE (drv%dia,*)'-- EOF regions  not found neither in the grid file nor in the EOF file --'
         CALL abort
      ENDIF
   ENDIF

! Check the size of nreg. If it is equal to the size of the global grid split the
! matrix on processors. Otherwise READ the whole matrix on each processor.
   IF ( ros%nreg .EQ. grd%img*grd%jmg ) THEN

      ros%nreg = grd%im * grd%jm
      WRITE (drv%dia,*)' ---- NREG equal to ', ros%nreg

      ALLOCATE ( ros%evc( ros%nreg, ros%kmt, neofs), ros%eva( ros%nreg, neofs) )
      ALLOCATE ( ros%cor( ros%nreg, ros%neof, neofs) )
      ALLOCATE ( x3( grd%im, grd%jm, neofs), x4( grd%im, grd%jm, ros%kmt, neofs) )

      k = 0
      DO j = 1,grd%jm
         DO i = 1,grd%im
            k = k + 1
            grd%reg(i,j) = k
         ENDDO
      ENDDO

      start(1) = grd%igs
      start(2) = grd%jgs
      start(3) = 1
      start(4) = 0

      count(1) = grd%im
      count(2) = grd%jm
      count(3) = neofs
      count(4) = 0

      stat = NF90_INQ_VARID (ncid, 'eva', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid,idvar,x3,start(1:3),count(1:3))
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      DO ke = 1,neofs
         k = 0
         DO j = 1,grd%jm
            DO i = 1,grd%im
               k = k + 1
               ros%eva(k,ke) = DBLE(x3(i,j,ke))
            ENDDO
         ENDDO
      ENDDO

      start(1) = grd%igs
      start(2) = grd%jgs
      start(3) = 1
      start(4) = 1
      count(1) = grd%im
      count(2) = grd%jm
      count(3) = ros%kmt
      count(4) = neofs

      stat = NF90_INQ_VARID (ncid, 'evc', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid,idvar,x4,start,count)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      DO ke = 1,neofs
         DO kt = 1,ros%kmt
            k = 0
            DO j = 1,grd%jm
               DO i = 1,grd%im
                  k = k + 1
                  ros%evc(k,kt,ke)  = DBLE(x4(i,j,kt,ke))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE ( x3, x4 )

   ELSE
!  Allocate eof arrays
      ALLOCATE ( ros%evc( ros%nreg, ros%kmt, neofs), ros%eva( ros%nreg, neofs) )
      ALLOCATE ( ros%cor( ros%nreg, ros%neof, neofs) )

      stat = NF90_INQ_VARID (ncid, 'eva', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid,idvar,ros%eva)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_INQ_VARID (ncid, 'evc', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid,idvar,ros%evc)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

   ENDIF

   stat = NF90_CLOSE(ncid)

   DO nec = 1,ros%neof
      DO k = 1,ros%neof
         DO nrg = 1,ros%nreg
            IF ( k .EQ. nec ) THEN
               ros%cor(nrg,k,nec) = ros%eva(nrg,nec)
            ELSE
               ros%cor(nrg,k,nec) = 0.0_r8
            ENDIF
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE rdeofs


