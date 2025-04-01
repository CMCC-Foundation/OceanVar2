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
!> Adjoint of the implicit solver - overrelaxation
!!
!!                                                                    
!!
! Version 1: Srdjan Dobricic  2007                                     !
! Version 2: Francesco Carere 2023                                     !
!-----------------------------------------------------------------------
SUBROUTINE invrt_ad(kstp)

   USE set_knd
   USE grd_str
   USE bmd_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)    :: kstp
   REAL(r8), ALLOCATABLE  :: res(:,:)
   INTEGER(i4)    :: i, j, icnt, ierr
   REAL(r8), ALLOCATABLE  :: bufst(:), bufsb(:), bufsr(:), bufsl(:)
   REAL(r8), ALLOCATABLE  :: bufrt(:), bufrb(:), bufrr(:), bufrl(:)

   ALLOCATE ( res(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

   res(:,:) = 0.0_r8

   DO icnt=bmd%itr(kstp),1,-1

      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(3_i4,2_i4,res,1_i4)
      DO j = grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
         DO i = grd%is3,grd%im-1+grd%iae,2       ! 3,im-1,2
            res(i,j) = 0.0_r8
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
         ENDDO
      ENDDO
      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(3_i4,2_i4,res,0_i4)


      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(2_i4,3_i4,res,1_i4)
      DO j = grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
         DO i = grd%is2,grd%im-1+grd%iae,2       ! 2,im-1,2
            res(i,j) = 0.0_r8
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
         ENDDO
      ENDDO
      IF ( mpi%nproc .GT. 1) CALL reproupdate_borders(2_i4,3_i4,res,0_i4)


      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(3_i4,3_i4,res,1_i4)
      DO j = grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
         DO i = grd%is3,grd%im-1+grd%iae,2       ! 3,im-1,2
            res(i,j) = 0.0_r8
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
         ENDDO
      ENDDO
      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(3_i4,3_i4,res,0_i4)

      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(2_i4,2_i4,res,1_i4)
      DO j = grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
         DO i = grd%is2,grd%im-1+grd%iae,2       ! 2,jm-1,2
            res(i,j) = 0.0_r8
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
         ENDDO
      ENDDO
      IF ( mpi%nproc .GT. 1 ) CALL reproupdate_borders(2_i4,2_i4,res,0_i4)
      bmd%eta(:,:) = bmd%eta(:,:) * bmd%mst(:,:)

   ENDDO

   DEALLOCATE (res)

END SUBROUTINE invrt_ad
!=====================================================================================================
SUBROUTINE reproupdate_borders(ist,jst,res,imod)

!-----------------------------------------------------------------------
!> Exchange borders values in adjoint overrelaxation algorithm
!!
!!                                                                    
!!
! Version 1: Francesco Carere 2023                                     !
!-----------------------------------------------------------------------

   USE set_knd
   USE bmd_str
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4), INTENT(IN)                                                    :: ist, jst, imod
   REAL(r8),DIMENSION(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae)      :: res
   REAL(r8),ALLOCATABLE                                                       :: buffsl(:), buffst(:), buffsr(:), buffsb(:)
   REAL(r8),ALLOCATABLE                                                       :: buffrl(:), buffrt(:), buffrr(:), buffrb(:)
   LOGICAL                                                                    :: lsendl, lsendr, lsendt, lsendb
   LOGICAL                                                                    :: lrecvl, lrecvr, lrecvt, lrecvb
   INTEGER(i4)                                                                :: isendl, isendb
   INTEGER(i4)                                                                :: isendr, isendt
   INTEGER                                                                    :: ierr, mpistat(mpi_status_size)
   INTEGER(i8)                                                                :: idx, i, j
   INTEGER(i4)                                                                :: is2_3, js2_3, mpiszbx, mpiszby

   IF ( ist .EQ. 2 ) THEN
      is2_3 = grd%is2
   ELSE
      is2_3 = grd%is3
   ENDIF
   IF ( jst .EQ. 2 ) THEN
      js2_3 = grd%js2
   ELSE
      js2_3 = grd%js3
   ENDIF

   IF ( mpi%jr .EQ. 1 ) THEN
      mpiszbx = (grd%jm-1+grd%jae)/2
   ELSE
      mpiszbx = (grd%jm+grd%jae)/2
   ENDIF
   IF ( mpi%ir .EQ. 1 ) THEN
      mpiszby = (grd%im-1+grd%iae)/2
   ELSE
      mpiszby = (grd%im+grd%iae)/2
   ENDIF

   ALLOCATE( buffrt(mpiszby), buffrb(mpiszby), buffrr(mpiszbx), buffrl(mpiszbx) )
   ALLOCATE( buffst(mpiszby), buffsb(mpiszby), buffsr(mpiszbx), buffsl(mpiszbx) )

   IF ( imod .EQ. 0 ) THEN
      lsendl = is2_3 .NE. 2                         .AND. mpi%ir .NE. 1       !send    left
      lsendb = js2_3 .NE. 2                         .AND. mpi%jr .NE. 1       !send    bot
      lrecvr = MODULO(grd%ige,2) .NE. MODULO(ist,2) .AND. mpi%ir .NE. mpi%irm !receive right
      lrecvt = MODULO(grd%jge,2) .NE. MODULO(jst,2) .AND. mpi%jr .NE. mpi%jrm !receive rop
      IF ( lsendl ) THEN
         idx = 1
         DO j = js2_3,grd%jm-1+grd%jae,2
            buffsl(idx) = res(1,j) *bmd%a2(1,j)
            idx = idx+1
         ENDDO
         CALL MPI_ISEND( buffsl,mpiszbx,mpi%r8,mpi%lft,3,mpi%comm,isendl,ierr)
      ENDIF
      IF ( lsendb ) THEN
         idx = 1
         DO i = is2_3,grd%im-1+grd%iae,2
            buffsb(idx) =  res(i,1)*bmd%a4(i,1)
            idx = idx+1
         ENDDO
         CALL MPI_ISEND( buffsb,mpiszby,mpi%r8,mpi%bot,4,mpi%comm,isendb,ierr)
      ENDIF

      IF ( lrecvr )  THEN
         CALL MPI_RECV( buffrr,mpiszbx,mpi%r8,mpi%rgh,3,mpi%comm,mpistat,ierr)
         idx = 1
         i = grd%im
         DO j = js2_3,grd%jm-1+grd%jae,2            ! 2,jm-1,2
            bmd%eta(i,j) = bmd%eta(i,j) + buffrr(idx)
            idx = idx+1
         ENDDO
      ENDIF
      IF ( lrecvt ) THEN
         CALL MPI_RECV( buffrt,mpiszby,mpi%r8,mpi%top,4,mpi%comm,mpistat,ierr)
         idx = 1
         j = grd%jm
         DO i = is2_3,grd%im-1+grd%iae,2           ! 3,im-1,2
            bmd%eta(i,j) = bmd%eta(i,j) + buffrt(idx)
            idx = idx+1
         ENDDO
      ENDIF
      IF ( lsendl ) CALL MPI_WAIT( isendl, mpistat, ierr)
      IF ( lsendb ) CALL MPI_WAIT( isendb, mpistat, ierr)

   ELSE
      !send l/b, receive r/t
      lsendr = MODULO(grd%ige,2) .EQ. MODULO(ist,2) .AND. mpi%ir .NE. mpi%irm !send    right
      lsendt = MODULO(grd%jge,2) .EQ. MODULO(jst,2) .AND. mpi%jr .NE. mpi%jrm !send    top
      lrecvl = is2_3 .NE. 1                         .AND. mpi%ir .NE. 1       !receive left
      lrecvb = js2_3 .NE. 1                         .AND. mpi%jr .NE. 1       !receive bot
      IF ( lsendr ) THEN
         idx = 1
         i = grd%im
         DO j = js2_3,grd%jm-1+grd%jae,2
            buffsr(idx) = (bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)) *bmd%a1(i,j)
            idx = idx+1
         ENDDO
         CALL MPI_ISEND( buffsr,mpiszbx,mpi%r8,mpi%rgh,1,mpi%comm,isendr,ierr)
      ENDIF
      IF ( lsendt ) THEN
         idx = 1
         j = grd%jm
         DO i = is2_3,grd%im-1+grd%iae,2
            buffst(idx) = ( bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)) *bmd%a3(i,j)
            idx = idx+1
         ENDDO
         CALL MPI_ISEND( buffst,mpiszby,mpi%r8,mpi%top,2,mpi%comm,isendt,ierr)
      ENDIF

      IF ( lrecvl ) THEN
         CALL mpi_recv( buffrl,mpiszbx,mpi%r8,mpi%lft,1,mpi%comm,mpistat,ierr)
         idx = 1
         DO j = js2_3,grd%jm-1+grd%jae,2            ! 2,jm-1,2
            bmd%eta(1,j) = bmd%eta(1,j) + buffrl(idx)
            idx = idx+1
         ENDDO
      ENDIF
      IF ( lrecvb )  THEN
         CALL MPI_RECV( buffrb,mpiszby,mpi%r8,mpi%bot,2,mpi%comm,mpistat,ierr)
         idx = 1
         DO i = is2_3,grd%im-1+grd%iae,2           ! 3,im-1,2
            bmd%eta(i,1) = bmd%eta(i,1) + buffrb(idx)
            idx = idx+1
         ENDDO
      ENDIF
      IF ( lsendr ) CALL MPI_WAIT( isendr, mpistat, ierr)
      IF ( lsendt ) CALL MPI_WAIT( isendt, mpistat, ierr)
   ENDIF

   DEALLOCATE( buffrt, buffrb, buffrr, buffrl )
   DEALLOCATE( buffst, buffsb, buffsr, buffsl )

END SUBROUTINE reproupdate_borders
