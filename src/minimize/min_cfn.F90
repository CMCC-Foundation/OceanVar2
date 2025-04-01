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
!> Minimize the cost function                                          
!!
!! It minimize the cost fuction iteratively until the convergence 
!! criteria is reached
!!
!                                                                      !
! Version 1: Srdjan Dobricic                  2006                     !
! Version 2: Mario Adani and Francesco Carere 2024                     !
!-----------------------------------------------------------------------
SUBROUTINE min_cfn

   USE set_knd
   USE drv_str
   USE ctl_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER                                   :: ierr
   REAL(r8)                                  :: mpism
   INTEGER(i4)                               :: cntr, k, ccnt, ccntmx
   REAL(r8)                                  :: maxpg
#ifdef REPRO
   DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: sort_array
#endif

   cntr   = 0
   ccnt   = 0
   ccntmx = 20

   DO WHILE( (ctl%task(1:5).EQ.'NEW_X' .OR. ctl%task(1:2).EQ.'FG' .OR.                    &
      ctl%task(1:5).EQ.'START') .AND. ccnt.LE.drv%cntm(drv%ktr) )

      IF ( mpi%flg_min .EQ. 0 ) THEN !run the minimiser in parallel

         CALL setulb(ctl%n, ctl%m, ctl%x_c, ctl%l_c, ctl%u_c, ctl%nbd, ctl%f_c, ctl%g_c,    &
            ctl%factr, ctl%pgtol, ctl%ws, ctl%wy, ctl%sy, ctl%ss, ctl%yy, ctl%wt,  &
            ctl%wn, ctl%snd, ctl%z_c, ctl%r_c, ctl%d_c, ctl%t_c, ctl%wa, ctl%sg,   &
            ctl%sgo, ctl%yg, ctl%ygo, ctl%iwa,                                     &
            ctl%task, ctl%iprint,  ctl%csave, ctl%lsave, ctl%isave, ctl%dsave,     &
            mpi%comm, mpi%myrank, mpi%nproc)

      ELSE !run the sequential minimiser while in parallel environment


         CALL gather_seqmin !Gather info on processor 0

         IF ( mpi%myrank .EQ. 0 ) THEN
#ifdef REPRO
            ALLOCATE ( sort_array(SIZE(ctl_glob%x_c)) )
            sort_array = ctl_glob%x_c
            CALL gth_to_seq(sort_array,ctl_glob%x_c)
            sort_array = ctl_glob%g_c
            CALL gth_to_seq(sort_array,ctl_glob%g_c)
#endif

            CALL setulb(ctl_glob%n, ctl%m, ctl_glob%x_c, ctl_glob%l_c, ctl_glob%u_c,     &
               ctl_glob%nbd, ctl%f_c, ctl_glob%g_c, ctl%factr, ctl%pgtol, ctl_glob%ws, &
               ctl_glob%wy, ctl%sy, ctl%ss, ctl%yy, ctl%wt, ctl%wn, ctl%snd,           &
               ctl_glob%z_c, ctl_glob%r_c, ctl_glob%d_c, ctl_glob%t_c, ctl%wa, ctl%sg, &
               ctl%sgo, ctl%yg, ctl%ygo, ctl_glob%iwa,                                 &
               ctl%task, ctl%iprint,  ctl%csave, ctl%lsave, ctl%isave, ctl%dsave,      &
               MPI_COMM_SELF, mpi%myrank, 1)

#ifdef REPRO
            sort_array = ctl_glob%x_c
            CALL seq_to_gth(sort_array,ctl_glob%x_c)
            sort_array = ctl_glob%g_c
            CALL seq_to_gth(sort_array,ctl_glob%g_c)
            DEALLOCATE ( sort_array )
#endif
         ENDIF

         CALL scatter_seqmin !scatter/bcast info from processor 0

      ENDIF

      IF ( ctl%task(1:2) .EQ. 'FG' ) THEN

! Calculate the cost function and its gradient
         CALL costf

! Calculate norm of he gradient
         maxpg = 0.0_r8
         IF ( mpi%flg_min .EQ. 0 ) THEN !run the minimiser in parallel

            maxpg = SUM(ctl%g_c**2)
            IF ( mpi%nproc .GT. 1 ) THEN
               CALL MPI_REDUCE( maxpg, mpism, 1, mpi%r8, MPI_SUM, 0, mpi%comm, ierr)
               IF ( mpi%myrank==0 ) THEN
                  maxpg = mpism
               ENDIF
               CALL mpi_bcast(maxpg, 1, mpi%r8, 0, mpi%comm, ierr)
            ENDIF

         ELSE                          !run the sequential minimiser while in parallel environment

            CALL gather_seqmin !Gather info on processor 0
            IF ( mpi%myrank .EQ. 0 ) THEN
#ifdef REPRO
               ALLOCATE ( sort_array(ctl_glob%n) )
               CALL gth_to_seq(ctl_glob%g_c,sort_array)
               maxpg = SUM(sort_array**2)
               DEALLOCATE ( sort_array )
#else
               maxpg = SUM(ctl_glob%g_c**2)
#endif
            ENDIF
            CALL mpi_bcast(maxpg, 1, mpi%r8,  0, mpi%comm, ierr)

         ENDIF

         maxpg = DSQRT(maxpg)

! Modify the stopping criteria
         IF ( cntr .EQ. 0 .AND. ctl%pgper .NE. 0.0_r8 .AND. drv%ktr .EQ. 1 ) THEN
            ctl%pgtol = maxpg * ctl%pgper
            cntr = 1
         ENDIF

         ccnt = ccnt + 1

      ENDIF

      IF ( ccnt .GT. 0 .AND. ctl%f_c .LT. 1.e-32 ) ccnt = drv%cntm(drv%ktr) + 1

      IF ( maxpg .LT. ctl%pgtol ) THEN
         ccnt = drv%cntm(drv%ktr) + 1
         WRITE (drv%dia,*) ' ---- ---- ---- ---- ---- ---- ---- ---- ---- '
         WRITE (drv%dia,*) ' ---- The minimization converged: '
         WRITE (drv%dia,*) ' ---- The cost function:          ',ctl%f_c
         WRITE (drv%dia,*) ' ---- The norm of the gradient:   ',maxpg
         WRITE (drv%dia,*) ' ---- The stopping criteria:      ',ctl%pgtol
         WRITE (drv%dia,*) ' ---- ---- ---- ---- ---- ---- ---- ---- ---- '
      ENDIF

   ENDDO !while

END SUBROUTINE min_cfn
