program oceanvar

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006, 2007 Srdjan Dobricic, CMCC, Bologna                   !
!                                                                          !
!    This file is part of OceanVar.                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! The main driver for the OceanVar                                     !
!                                                                      !
! Version 0.1: S.Dobricic 2006                                         !
!   Horizontal covariance with recursive filters, vertical with EOFs,  !
!   assimilation of satellite observations of SLA, in situ observations!
!   by XBT and ARGO floats                                             !
!                                                                      !
! Version 0.2: S.Dobricic 2007                                         !
!   Multigrid method. Internal boundaries for horizontal covariances.  !
!                                                                      !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use mpi_str

 implicit none

 INTEGER(i4)   ::  kts, ktr, ierr

! ---
! Initialize mpi
      call mpi_init(ierr)

! ---
! Initialize mpi, diagnostics and read namelists
      call def_nml

! ---
! Outer loop - SST assimilation
  do kts = 1,drv%nts

   write(drv%dia,*) ' Outer loop', kts, drv%nts
   drv%kts = kts

! ---
! Outer loop - multigrid
  do ktr = 1,drv%ntr

   write(drv%dia,*) ' Outer loop-multigrid', ktr, drv%ntr
   drv%ktr = ktr

! ---
! Define grid parameters
      if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 )then
       write(drv%dia,*) ' Define grid parameters', ktr, drv%ratio(ktr)
       call def_grd
      endif

! ---
! Get observations
      if( ktr.eq.1 .and.  kts.eq.1 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Get observations', ktr, kts
       call get_obs
      endif

! ---
! Define interpolation parameters
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 )then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Define interpolation parameters'
      call int_par
     endif

! ---
! If assimilating SST save misfits
     if( ktr.eq.1 .and. kts.eq.1 .and. drv%nts.eq.2 ) then
      call sav_msf
     endif

! ---
! Define observational vector
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 )then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Define observational vector'
      call obs_vec
     endif

! ---
! Define constants for background covariances
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Define constants for background covariances'
      call def_cov
     endif

! ---
! Localize verticaly by using the mixed layer depth
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Localize verticaly by using the mixed layer depth'
      call rdmxd
     endif

! ---
! Localize horizontaly around observations
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Localize horizontaly around observations'
      call def_loc
     endif

! ---
! Initialise barotropic model
      if(drv%bmd(drv%ktr) .eq. 1) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Initialise barotropic model'
         call ini_bmd
      endif

! ---
! Initialize cost function and its gradient
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Initialize cost function and its gradient'
      call ini_cfn

! ---
! Calculate the initial norm of the gradient
     if( ktr.gt.1 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Calculate the initial norm of the gradient'
      call ini_nrm
     endif

! ---
! Initialise from old iterration
     if( ktr.gt.1 .and. drv%ratio(ktr).ne.1.0 ) then
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) ' Initialise from old iterration'
      call ini_itr
     endif

! ---
! Minimize the cost function (inner loop)
       write(drv%dia,*) ' ------------------------------------'
       write(drv%dia,*) '  Minimize the cost function (inner loop)'
     call min_cfn

! ---
! Save old iterration
   if( ktr.lt.drv%ntr)then
    if(drv%ratio(ktr+1).ne.1.0 ) then
     call sav_itr
    endif
   endif

! ---
! Convert to innovations
   if( drv%ktr.eq.drv%ntr .and. drv%kts.eq.drv%nts ) then
     call cnv_inn
   endif

! ---
! End of outer loop - multigrid
  enddo
     write(drv%dia,*) ' ------------------------------------'
     write(drv%dia,*) '  End of outer loop - multigrid'

! ---
! If assimilating SST modify misfits and save increments
 if( drv%ktr.eq.drv%ntr .and. drv%kts.eq.1 .and. drv%nts.eq.2 ) then
     call cnv_inn
     call obsop
     call mod_msf
     call sav_itr
 endif

! ---
! End of outer loop - SST assimilation
  enddo

     write(drv%dia,*) ' ------------------------------------'
     write(drv%dia,*) '  End of outer loop - SST assimilation'

! ---
! If assimilating SST modify increments
 if( drv%nts.eq.2 ) then
     call mod_inc
 endif

! ---
! Write outputs and diagnostics
     call wrt_dia

 
! ---
! End of mpi
         call mpi_finalize(ierr)

!-----------------------------------------------------------------

end program oceanvar
