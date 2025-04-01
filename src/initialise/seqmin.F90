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
!> Inizialization sequential minimizer                                
!!
!! For reproducibility purpose the minimizer can be call sequentially 
!! It is activated from namelist if flg_min is different from 0
!! To be completely reproducible the element of the vectors entering 
!! in the minimizer must be shuffled. This can be activated using 
!! the compilation key -DREPRO in Makefile
!!
!                                                                      !
! Version 1: Francesco Carere 2024                                     !
!-----------------------------------------------------------------------
SUBROUTINE ini_seqmin

   USE grd_str
   USE eof_str
   USE ctl_str
   USE mpi_str
   USE drv_str

   IMPLICIT NONE

   INTEGER                       :: i, ierr, j
   DOUBLE PRECISION, ALLOCATABLE :: buffr(:)

   INCLUDE 'mpif.h'

! ---
! Allocate memory for global optimization arrays
   ctl_glob%n=grd%npsa * ros%neof
   ALLOCATE( ctl_glob%n_g(mpi%nproc), ctl_glob%n_cum(mpi%nproc) )
   ALLOCATE( ctl_glob%nbd( ctl_glob%n ), ctl_glob%iwa( 3*ctl_glob%n ) )
   ALLOCATE( ctl_glob%x_c( ctl_glob%n ), ctl_glob%g_c( ctl_glob%n ) )
   ALLOCATE( ctl_glob%l_c( ctl_glob%n ), ctl_glob%u_c( ctl_glob%n ) )
   ALLOCATE( ctl_glob%ws(ctl_glob%n,ctl%m), ctl_glob%wy(ctl_glob%n,ctl%m) )
   ALLOCATE( ctl_glob%z_c(ctl_glob%n), ctl_glob%r_c(ctl_glob%n), ctl_glob%d_c(ctl_glob%n), ctl_glob%t_c(ctl_glob%n) )
   ALLOCATE( buffr(ctl_glob%n) )


! ---
! Initialise the arrays
   CALL MPI_GATHER(ctl%n, 1, mpi%i4, ctl_glob%n_g, 1, mpi%i4, 0, mpi%comm, ierr)
   ctl_glob%n_cum(1)=0_i4
   ctl_glob%n_cum(2:) = [(sum(ctl_glob%n_g(1:i)), i = 1, (size(ctl_glob%n_cum)-1))]

   CALL MPI_GATHERV(ctl%nbd, ctl%n, mpi%i4, ctl_glob%nbd, ctl_glob%n_g, ctl_glob%n_cum, mpi%i4, 0, mpi%comm, ierr)
   CALL MPI_GATHERV(ctl%x_c, ctl%n, mpi%r8, ctl_glob%x_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)
   CALL MPI_GATHERV(ctl%g_c, ctl%n, mpi%r8, ctl_glob%g_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)
   CALL MPI_GATHERV(ctl%l_c, ctl%n, mpi%r8, ctl_glob%l_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)
   CALL MPI_GATHERV(ctl%u_c, ctl%n, mpi%r8, ctl_glob%u_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)

   ctl_glob%nbd=0.0_i4
   ctl_glob%iwa=0.0_i4
   ctl_glob%ws=0.0_r8;  ctl_glob%wy=0.0_r8
   ctl_glob%z_c=0.0_r8; ctl_glob%r_c=0.0_r8
   ctl_glob%d_c=0.0_r8; ctl_glob%t_c=0.0_r8

END SUBROUTINE ini_seqmin
!======================================================================
SUBROUTINE gather_seqmin

!-----------------------------------------------------------------------
!                                                                      !
! Gathering vectors for sequential minimizer                           !
!                                                                      !
! Version 1: Francesco Carere 2024                                     !
!-----------------------------------------------------------------------

   USE ctl_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER     :: ierr

   IF ( .NOT. ALLOCATED(ctl_glob%g_c ) ) THEN
      CALL ini_seqMIN
   ELSE
      !gather what has changed: g_c, x_c
      CALL MPI_GATHERV(ctl%g_c, ctl%n, mpi%r8, ctl_glob%g_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)
      CALL MPI_GATHERV(ctl%x_c, ctl%n, mpi%r8, ctl_glob%x_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, 0, mpi%comm, ierr)
   ENDIF

END SUBROUTINE gather_seqmin
!======================================================================
SUBROUTINE scatter_seqmin

!-----------------------------------------------------------------------
!                                                                      !
! Scattering vectors for sequential minimizer                          !
!                                                                      !
! Version 1: Francesco Carere 2024                                     !
!-----------------------------------------------------------------------

   USE ctl_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER     :: ierr

   !Broadcast task, f_c, isave(30)
   CALL mpi_bcast(ctl%f_c, 1, mpi%r8, 0, mpi%comm,ierr)
   CALL mpi_bcast(ctl%task,60,MPI_CHAR,0,mpi%comm,ierr)
   CALL mpi_bcast(ctl%isave(30),1, mpi%i4,0, mpi%comm, ierr)

   !Scatter the updated arrays: x_c, g_c
   CALL MPI_SCATTERv(ctl_glob%x_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, ctl%x_c, ctl%n, mpi%r8, 0, mpi%comm, ierr )
   CALL MPI_SCATTERv(ctl_glob%g_c, ctl_glob%n_g, ctl_glob%n_cum, mpi%r8, ctl%g_c, ctl%n, mpi%r8, 0, mpi%comm, ierr )

END SUBROUTINE scatter_seqmin
!======================================================================
SUBROUTINE get_indexes

!-----------------------------------------------------------------------
!                                                                      !
! mapping indexes between parallel and sequential version              !
!                                                                      !
! Version 1: Mario Adani 2024                                          !
!-----------------------------------------------------------------------

   USE set_knd
   USE grd_str
   USE mpi_str
   USE rpr_str
   USE eof_str

   IMPLICIT NONE

   INTEGER                                  :: ierr
   INTEGER(i4)                              :: proc,i,j,k,kk,kkk
   INTEGER(i4)                              :: im(mpi%nproc),jm(mpi%nproc)
   INTEGER(i4)                              :: igs(mpi%nproc),jgs(mpi%nproc)
   INTEGER(i4)                              :: ige(mpi%nproc),jge(mpi%nproc)
   INTEGER(i4)                              :: tot_size,proc_size
   INTEGER(i4),ALLOCATABLE,dimension(:,:,:) :: map, submap
   INTEGER(i4),ALLOCATABLE,dimension(:)     :: subvec

   INCLUDE 'mpif.h'

   CALL MPI_GATHER(grd%igs, 1, mpi%i4, igs, 1, mpi%i4, 0, mpi%comm, ierr)
   CALL MPI_GATHER(grd%jgs, 1, mpi%i4, jgs, 1, mpi%i4, 0, mpi%comm, ierr)
   CALL MPI_GATHER(grd%ige, 1, mpi%i4, ige, 1, mpi%i4, 0, mpi%comm, ierr)
   CALL MPI_GATHER(grd%jge, 1, mpi%i4, jge, 1, mpi%i4, 0, mpi%comm, ierr)

   IF (mpi%myrank==0) THEN

! Compute dimension of each tile
      DO proc = 1,mpi%nproc
         im(proc)=ige(proc)-igs(proc)+1
         jm(proc)=jge(proc)-jgs(proc)+1
      ENDDO

      tot_size=grd%img*grd%jmg*ros%neof

      ALLOCATE(map(grd%img,grd%jmg,ros%neof))
! 1) create a map with indexes
      kk = 0
      DO k = 1,ros%neof
         DO j = 1,grd%jmg
            DO i = 1,grd%img
               kk = kk+1
               map(i,j,k) = kk
            ENDDO
         ENDDO
      ENDDO

      ALLOCATE(rpr%idxmap(tot_size))
      kkk = 0
      DO proc = 1,mpi%nproc
! 2) subset the global map for each prcessor
         ALLOCATE (submap(im(proc),jm(proc),ros%neof) )
         submap(1:im(proc),1:jm(proc),1:ros%neof) = map( igs(proc):ige(proc), &
                                                         jgs(proc):jge(proc), &
                                                         1        :ros%neof    )
         proc_size = im(proc)*jm(proc)*ros%neof
! 3) from map to vector
         ALLOCATE( subvec(proc_size) )
         kk = 0
         DO k = 1,ros%neof
            DO j = 1,jm(proc)
               DO i = 1,im(proc)
                  kk = kk + 1
                  subvec(kk) = submap(i,j,k)
               ENDDO
            ENDDO
         ENDDO

         DEALLOCATE(submap)

! 3) from each processor vector  to total vector
         DO kk = 1,proc_size
            kkk = kkk+1
            rpr%idxmap(kkk) = subvec(kk)
         ENDDO
         DEALLOCATE(subvec)

      ENDDO
   ENDIF

END SUBROUTINE get_indexes
!======================================================================
SUBROUTINE gth_to_seq(arrayin,arrayou)

!-----------------------------------------------------------------------
!                                                                      !
! remapping from parallel to sequential                                !
!                                                                      !
! Version 1: Mario Adani 2024                                          !
!-----------------------------------------------------------------------

   USE set_knd
   USE rpr_str
   USE eof_str
   USE grd_str

   IMPLICIT NONE

   DOUBLE PRECISION :: arrayin(grd%img*grd%jmg*ros%neof)
   DOUBLE PRECISION :: arrayou(grd%img*grd%jmg*ros%neof)
   INTEGER(i4)      :: tot_size,i

   tot_size = grd%img*grd%jmg*ros%neof

   DO i = 1,tot_size
      arrayou(rpr%idxmap(i)) = arrayin(i)
   ENDDO

END SUBROUTINE gth_to_seq
!======================================================================
SUBROUTINE seq_to_gth(arrayin,arrayou)

!-----------------------------------------------------------------------
!                                                                      !
! remapping from sequential to parallel                                !
!                                                                      !
! Version 1: Mario Adani 2024                                          !
!-----------------------------------------------------------------------

   USE set_knd
   USE rpr_str
   USE eof_str
   USE grd_str

   IMPLICIT NONE

   DOUBLE PRECISION :: arrayin(grd%img*grd%jmg*ros%neof)
   DOUBLE PRECISION :: arrayou(grd%img*grd%jmg*ros%neof)
   INTEGER(i4)      :: tot_size,i

   tot_size = grd%img*grd%jmg*ros%neof

   DO i = 1,tot_size
      arrayou(i) = arrayin(rpr%idxmap(i))
   ENDDO

END SUBROUTINE seq_to_gth

