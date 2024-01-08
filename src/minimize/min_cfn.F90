subroutine min_cfn


!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Minimise the cost function                                           !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str
 use mpi_str

 implicit none

 include 'mpif.h'
 INTEGER               :: ierr
 REAL(r8)              :: mpism

 INTEGER(i4)     :: cntr, k, ccnt, ccntmx 
 REAL(r8)        :: maxpg

  cntr = 0
  ccnt = 0
  ccntmx = 20


  do while( (ctl%task(1:5).eq.'NEW_X' .or. ctl%task(1:2).eq.'FG' .or.                    &
             ctl%task(1:5).eq.'START') .and. ccnt.le.drv%cntm(drv%ktr) )


    call setulb(ctl%n, ctl%m, ctl%x_c, ctl%l_c, ctl%u_c, ctl%nbd, ctl%f_c, ctl%g_c,    &
                ctl%factr, ctl%pgtol, ctl%ws, ctl%wy, ctl%sy, ctl%ss, ctl%yy, ctl%wt,  &
                ctl%wn, ctl%snd, ctl%z_c, ctl%r_c, ctl%d_c, ctl%t_c, ctl%wa, ctl%sg,   &
                ctl%sgo, ctl%yg, ctl%ygo, ctl%iwa,                                     &
                ctl%task, ctl%iprint,  ctl%csave, ctl%lsave, ctl%isave, ctl%dsave,     &
                mpi%comm, mpi%myrank, mpi%nproc) 



      if (ctl%task(1:2) .eq. 'FG') then


! Calculate the cost function and its gradient
          call costf

! Calculate norm of he gradient
            maxpg = 0.0
           do k=1,ctl%n
            maxpg = maxpg+ctl%g_c(k)**2
           enddo

         if(mpi%nproc.gt.1)then
           call mpi_reduce( maxpg, mpism, 1, mpi%r8, mpi_sum, 0, mpi%comm, ierr)
           if(mpi%myrank==0)then
                maxpg = mpism
           endif
           call mpi_bcast(maxpg, 1, mpi%r8, 0, mpi%comm, ierr)
         endif
           maxpg = sqrt(maxpg)

! Modify the stopping criteria
       if(cntr.eq.0 .and. ctl%pgper.ne.0.0 .and. drv%ktr.eq.1 )then
            ctl%pgtol = maxpg * ctl%pgper
            cntr = 1
       endif

            ccnt = ccnt + 1

      endif

      if( ccnt.gt.0 .and. ctl%f_c.lt.1.e-32 ) ccnt = drv%cntm(drv%ktr) + 1

      if( maxpg.lt.ctl%pgtol ) then
          ccnt = drv%cntm(drv%ktr) + 1
          write(drv%dia,*) ' ---- ---- ---- ---- ---- ---- ---- ---- ---- '
          write(drv%dia,*) ' ---- The minimization converged: '
          write(drv%dia,*) ' ---- The cost function:          ',ctl%f_c
          write(drv%dia,*) ' ---- The norm of the gradient:   ',maxpg
          write(drv%dia,*) ' ---- The stopping criteria:      ',ctl%pgtol
          write(drv%dia,*) ' ---- ---- ---- ---- ---- ---- ---- ---- ---- '
      endif
       
  enddo !while



end subroutine min_cfn
