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
!> Read Mix Layer Depth                                                 
!!                                                                      
!!                                                                      
!!                                                                      
! Version 1: Srdjan Dobricic 2006                                       !
!-----------------------------------------------------------------------
SUBROUTINE rdmxd

   USE set_knd
   USE drv_str
   USE grd_str
   USE eof_str
   USE mpi_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                :: stat, ncid, idvar,k, ierr, i, j
   INTEGER(i4)                :: start(3), count(3)
   INTEGER(i4), ALLOCATABLE   :: x2(:,:)

   IF ( drv%kts.EQ.1 .AND. drv%nts.EQ.2 .AND. ros%mld.EQ.1 ) THEN

      stat = NF90_OPEN(drv%inpdir//'/'//'mldn.nc', NF90_NOWRITE, ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      ALLOCATE ( x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      start(3) = 1
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae
      count(3) = 1 

      stat = NF90_INQ_VARID (ncid, 'mldd', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid,idvar,x2,start,count)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      grd%mxd(:,:) =  x2(:,:)

      DEALLOCATE ( x2 )

      stat = NF90_CLOSE(ncid)


      DO j = 1-grd%jas,grd%jm+grd%jae
         DO i = 1-grd%ias,grd%im+grd%iae

            grd%mxd(i,j) = MAX(2,MIN(grd%mxd(i,j),grd%km-1))

            DO k = 1,grd%mxd(i,j)-1
               grd%lcl(i,j,k) = 1.0_r8
            ENDDO
            k = grd%mxd(i,j)
            grd%lcl(i,j,k) = 0.5_r8
            DO k = grd%mxd(i,j)+1,grd%km
               grd%lcl(i,j,k) = 0.0_r8
            ENDDO

         ENDDO
      ENDDO

   ELSE

      grd%lcl(:,:,:) = 1.0_r8
      WRITE (drv%dia,*)'Mixlayer loc=1'

   ENDIF

END SUBROUTINE rdmxd


