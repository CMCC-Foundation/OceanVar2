subroutine obs_vec


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
! Create the observational vector                                      !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str

 implicit none

  INTEGER(i4)    ::  k, i

! -------
! Define observational vector
   write(drv%dia,*) ' ---- Defining the Observational Vector'
   write(drv%dia,*) ' ---- number of good observations (sla%nc): ', sla%nc
   write(drv%dia,*) ' ---- number of good observations (arg%nc): ', arg%nc
   write(drv%dia,*) ' ---- number of good observations (xbt%nc): ', xbt%nc
    obs%no = sla%nc + arg%nc + xbt%nc + gld%nc + 2 * tra%nc + 2 * trd%nc  &
           + vdr%nc + gvl%nc + sst%nc

   write(drv%dia,*) ' ---- Total number of good observations: ', obs%no
   write(drv%dia,*) ' -------------------------------------- '
   ALLOCATE ( obs%inc(obs%no), obs%amo(obs%no), obs%res(obs%no))
   ALLOCATE ( obs%err(obs%no), obs%gra(obs%no))

   k=0

! SLA observations
   do i=1,sla%no
    if(sla%flc(i).eq.1)then
     k=k+1
     obs%res(k) = sla%res(i)
     obs%err(k) = sla%err(i)
    endif
   enddo

! ARGO observations
   do i=1,arg%no
    if(arg%flc(i).eq.1)then
     k=k+1
     obs%res(k) = arg%res(i)
     obs%err(k) = arg%err(i)

    endif
   enddo

! XBT observations
   do i=1,xbt%no
    if(xbt%flc(i).eq.1)then
     k=k+1
     obs%res(k) = xbt%res(i)
     obs%err(k) = xbt%err(i)
    endif
   enddo

! Glider observations
   do i=1,gld%no
    if(gld%flc(i).eq.1)then
     k=k+1
     obs%res(k) = gld%res(i)
     obs%err(k) = gld%err(i)
    endif
   enddo

! Argo trajectory observations
   do i=1,tra%no
    if(tra%flc(i).eq.1)then
     k=k+1
     obs%res(k) = tra%rex(i)
     obs%err(k) = tra%erx(i)
    endif
   enddo
   do i=1,tra%no
    if(tra%flc(i).eq.1)then
     k=k+1
     obs%res(k) = tra%rey(i)
     obs%err(k) = tra%ery(i)
    endif
   enddo

! Trajectory observations of surface drifters
   do i=1,trd%no
    if(trd%flc(i).eq.1)then
     k=k+1
     obs%res(k) = trd%rex(i)
     obs%err(k) = trd%erx(i)
    endif
   enddo
   do i=1,trd%no
    if(trd%flc(i).eq.1)then
     k=k+1
     obs%res(k) = trd%rey(i)
     obs%err(k) = trd%ery(i)
    endif
   enddo

! Observations of drifter velocities
   do i=1,vdr%no
    if(vdr%flc(i).eq.1)then
     k=k+1
     obs%res(k) = vdr%res(i)
     obs%err(k) = vdr%err(i)
    endif
   enddo

! Observations of glider velocities
   do i=1,gvl%no
    if(gvl%flc(i).eq.1)then
     k=k+1
     obs%res(k) = gvl%res(i)
     obs%err(k) = gvl%err(i)
    endif
   enddo

! SST observations
   do i=1,sst%no
    if(sst%flc(i).eq.1)then
     k=k+1
     obs%res(k) = sst%res(i)
     obs%err(k) = sst%err(i)
    endif
   enddo
! Oddo
   !do i=1,k
   !write(*,*)obs%res(i),obs%err(i),i,'Paolo obs_vec'
   !enddo


end subroutine obs_vec
