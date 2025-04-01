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
!> Tiles definition                                                    
!!
!! It defines the tiles for MPI execution
!!
!                                                                      !
! Version 1: Srdjan Dobricic   2011                                    !
!            Francesco Carere  2023                                    !
!-----------------------------------------------------------------------
SUBROUTINE def_mpid

   USE set_knd
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)               :: iw, jw, kproc, ierr, lr
   INTEGER(i4), ALLOCATABLE  :: ibff (:)
   INTEGER(i4), ALLOCATABLE  :: ibffa(:,:)
   INTEGER, ALLOCATABLE      :: lbf(:,:), lbfa(:,:,:)

   IF ( mpi%nproc .EQ. 1 ) THEN
      mpi%jr  = 1
      mpi%ir  = 1

      grd%igs = 1
      grd%ige = grd%img
      grd%im  = grd%img

      grd%ias = 0
      grd%iae = 0

      grd%jgs = 1
      grd%jge = grd%jmg
      grd%jm  = grd%jmg

      grd%jas = 0
      grd%jae = 0

      mpi%top = MPI_PROC_NULL
      mpi%bot = MPI_PROC_NULL
      mpi%rgh = MPI_PROC_NULL
      mpi%lft = MPI_PROC_NULL

      grd%is2 = 2
      grd%js2 = 2
      grd%is3 = 3
      grd%js3 = 3
   ELSE
! ---
! Find the position of the processor on the 2D grid
      mpi%jr  = mpi%myrank / mpi%irm + 1
      mpi%ir  = mpi%myrank + 1 - (mpi%jr - 1) * mpi%irm

! ---
! Starting and end points on the global grid, dimensions of the tile
      CALL para_reg( grd%img, mpi%irm, mpi%ir, grd%igs, grd%ige, grd%im)

      grd%ias = 1
      grd%iae = 1
      IF ( mpi%ir .EQ. 1       ) grd%ias = 0
      IF ( mpi%ir .EQ. mpi%irm ) grd%iae = 0

      IF (mpi%ir.EQ.1      ) THEN
         grd%is2 = 2
         grd%is3 = 3
      ELSE
         grd%is2 = 1+mod(grd%igs,2)
         grd%is3 = 2-mod(grd%igs,2)
      ENDIF

      CALL para_reg( grd%jmg, mpi%jrm, mpi%jr, grd%jgs, grd%jge, grd%jm)

      grd%jas = 1
      grd%jae = 1
      IF ( mpi%jr.EQ.mpi%jrm ) grd%jae = 0
      IF ( mpi%jr.EQ.1       ) grd%jas = 0

      IF ( mpi%jr.EQ.1       ) THEN
         grd%js2 = 2
         grd%js3 = 3
      ELSE
         grd%js2 = 1+MOD(grd%jgs,2)
         grd%js3 = 2-MOD(grd%jgs,2)
      ENDIF

! ---
! Define surrounding processors
      mpi%top = mpi%myrank + mpi%irm
      mpi%bot = mpi%myrank - mpi%irm
      mpi%lft = mpi%myrank - 1
      mpi%rgh = mpi%myrank + 1
      IF ( mpi%jr .EQ. mpi%jrm ) mpi%top = MPI_PROC_NULL
      IF ( mpi%jr .EQ. 1       ) mpi%bot = MPI_PROC_NULL
      IF ( mpi%ir .EQ. mpi%irm ) mpi%rgh = MPI_PROC_NULL
      IF ( mpi%ir .EQ. 1       ) mpi%lft = MPI_PROC_NULL
   ENDIF

   ALLOCATE ( mpi%thi(2), mpi%thj(2) )
   mpi%thi(1) = MIN( grd%im         , mpi%irm * mpi%thx )
   mpi%thi(2) = MIN( grd%im * grd%km, mpi%irm * mpi%thx )
   mpi%thj(1) = MIN( grd%jm         , mpi%jrm * mpi%thy )
   mpi%thj(2) = MIN( grd%jm * grd%km, mpi%jrm * mpi%thy )

   ALLOCATE ( grd%irs(mpi%thi(2),2), grd%ire(mpi%thi(2),2), grd%imr(mpi%thi(2),2) )
   ALLOCATE ( grd%jrs(mpi%thj(2),2), grd%jre(mpi%thj(2),2), grd%jmr(mpi%thj(2),2) )

   DO lr=1,mpi%thi(1)
      CALL para_reg( grd%im       , mpi%thi(1), lr, grd%irs(lr,1), grd%ire(lr,1), grd%imr(lr,1))
   ENDDO
   DO lr=1,mpi%thi(2)
      CALL para_reg( grd%im*grd%km, mpi%thi(2), lr, grd%irs(lr,2), grd%ire(lr,2), grd%imr(lr,2))
   ENDDO
   DO lr=1,mpi%thj(1)
      CALL para_reg( grd%jm       , mpi%thj(1), lr, grd%jrs(lr,1), grd%jre(lr,1), grd%jmr(lr,1))
   ENDDO
   DO lr=1,mpi%thj(2)
      CALL para_reg( grd%jm*grd%km, mpi%thj(2), lr, grd%jrs(lr,2), grd%jre(lr,2), grd%jmr(lr,2))
   ENDDO

! ---
! Get position of tiles to WRITE the output

   ALLOCATE ( grd%jgsp(mpi%nproc) )
   ALLOCATE ( grd%jgep(mpi%nproc) )
   ALLOCATE ( grd%igsp(mpi%nproc) )
   ALLOCATE ( grd%igep(mpi%nproc) )

   ALLOCATE ( ibff (4) )
   ALLOCATE ( ibffa(4,mpi%nproc) )

   ibff(1) = grd%jgs
   ibff(2) = grd%jge
   ibff(3) = grd%igs
   ibff(4) = grd%ige

   CALL MPI_GATHER(ibff, 4, mpi%i4, ibffa, 4, mpi%i4, 0, mpi%comm, ierr)

   IF ( mpi%myrank .EQ. 0 ) THEN
      grd%npsm = 0
      DO kproc = 1,mpi%nproc
         grd%jgsp(kproc) = ibffa(1,kproc)
         grd%jgep(kproc) = ibffa(2,kproc)
         grd%igsp(kproc) = ibffa(3,kproc)
         grd%igep(kproc) = ibffa(4,kproc)
         grd%npsm = MAX( grd%npsm, (grd%jgep(kproc)-grd%jgsp(kproc)+1)*(grd%igep(kproc)-grd%igsp(kproc)+1))
      ENDDO
   ENDIF

   CALL MPI_BCAST( grd%npsm, 1, mpi%i4, 0, mpi%comm, ierr)

   DEALLOCATE ( ibff )
   DEALLOCATE ( ibffa )

END SUBROUTINE def_mpid
!-----------------------------------------------------------------------------------
SUBROUTINE para_reg( img, irm, ir, igs, ige, im)

   USE set_knd

   IMPLICIT NONE

   INTEGER(i4)    ::  img
   INTEGER        ::  irm, ir
   INTEGER(i4)    ::  igs, ige, im
   INTEGER(i4)    ::  iw1, iw2

   iw1     = img / irm
   iw2     = MOD(img,irm)
   igs = (ir - 1) * iw1 + 1 + MIN((ir-1),iw2)
   ige = igs + iw1 - 1
   IF ( iw2 .GT. (ir-1) ) ige = ige + 1
   im  = ige - igs + 1

END SUBROUTINE para_reg
