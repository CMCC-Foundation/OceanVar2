subroutine exo_mpi ( iord, imod, iex, isnd, irec, jex, jsnd, jrec, is, ie, js, je, km, fld)

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
! Exchange borders for MPI        
!                                                                      !
! Version 1: S.Dobricic 2011                                           !
!-----------------------------------------------------------------------


 use set_knd
 use mpi_str

 implicit none

 include 'mpif.h'

 INTEGER(i4)   :: iord, imod, iex, isnd, irec, jex, jsnd, jrec
 INTEGER(i4)   :: is, ie, js, je, km
 REAL(r8)      :: fld(is:ie,js:je,km)
 INTEGER(i4)   :: ierr, mpiszx, mpiszy
 INTEGER       :: isendi, isendj
 INTEGER       :: irecvi, irecvj
 INTEGER       :: mpistat(mpi_status_size)
 REAL(r8), allocatable :: bufrb(:,:)
 REAL(r8), allocatable :: bufrt(:,:)
 REAL(r8), allocatable :: bufrl(:,:)
 REAL(r8), allocatable :: bufrr(:,:)


! ---

   ALLOCATE ( bufrb(ie-is+1,km), bufrt(ie-is+1,km) )
   ALLOCATE ( bufrl(je-js+1,km), bufrr(je-js+1,km) )

     mpiszx = (ie-is+1) * km
     mpiszy = (je-js+1) * km

  if( iord.eq.1 ) then

!-- Ordered exchange (first j, then i)
    if (jex.eq.-1) then
      if (mpi%jr.ne.1      ) bufrb(1:ie-is+1,1:km) = fld(is:ie,jsnd,1:km)
      if (mpi%jr.ne.1      ) call mpi_isend( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendj, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvj, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( isendj, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( irecvj, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm .and. imod.eq.1 ) then
                             fld(is:ie,jrec,1:km) = bufrt(1:ie-is+1,1:km)
      else if (mpi%jr.ne.mpi%jrm) then
                             fld(is:ie,jrec,1:km) = fld(is:ie,jrec,1:km) + bufrt(1:ie-is+1,1:km)
      endif
    else if (jex.eq.1) then
      if (mpi%jr.ne.mpi%jrm) bufrb(1:ie-is+1,1:km) = fld(is:ie,jsnd,1:km)
      if (mpi%jr.ne.mpi%jrm) call mpi_isend( bufrb, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendj, ierr)
      if (mpi%jr.ne.1      ) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvj, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( isendj, mpistat, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( irecvj, mpistat, ierr)
      if (mpi%jr.ne.1      .and. imod.eq.1 ) then
                             fld(is:ie,jrec,1:km) = bufrt(1:ie-is+1,1:km)
      else if (mpi%jr.ne.1      ) then
                             fld(is:ie,jrec,1:km) = fld(is:ie,jrec,1:km) + bufrt(1:ie-is+1,1:km)
      endif
    endif

    if (iex.eq.-1) then
      if (mpi%ir.ne.1      ) bufrl(1:je-js+1,1:km) = fld(isnd,js:je,1:km)
      if (mpi%ir.ne.1      ) call mpi_isend( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
      if (mpi%ir.ne.1      ) call mpi_wait( isendi, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( irecvi, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm .and. imod.eq.1) then 
                             fld(irec,js:je,1:km) = bufrr(1:je-js+1,1:km)
      else if (mpi%ir.ne.mpi%irm) then 
                             fld(irec,js:je,1:km) = fld(irec,js:je,1:km) + bufrr(1:je-js+1,1:km)
      endif
    else if (iex.eq.1) then
      if (mpi%ir.ne.mpi%irm) bufrl(1:je-js+1,1:km) = fld(isnd,js:je,1:km)
      if (mpi%ir.ne.mpi%irm) call mpi_isend( bufrl, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
      if (mpi%ir.ne.1      ) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( isendi, mpistat, ierr)
      if (mpi%ir.ne.1      ) call mpi_wait( irecvi, mpistat, ierr)
      if (mpi%ir.ne.1       .and. imod.eq.1) then 
                             fld(irec,js:je,1:km) = bufrr(1:je-js+1,1:km)
      else if (mpi%ir.ne.1      ) then 
                             fld(irec,js:je,1:km) = fld(irec,js:je,1:km) + bufrr(1:je-js+1,1:km)
      endif
    endif

  else

!-- Unordered exchange
    if (jex.eq.-1) then
      if (mpi%jr.ne.1      ) bufrb(1:ie-is+1,1:km) = fld(is:ie,jsnd,1:km)
      if (mpi%jr.ne.1      ) call mpi_isend( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendj, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvj, ierr)
    else if (jex.eq.1) then
      if (mpi%jr.ne.mpi%jrm) bufrb(1:ie-is+1,1:km) = fld(is:ie,jsnd,1:km)
      if (mpi%jr.ne.mpi%jrm) call mpi_isend( bufrb, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendj, ierr)
      if (mpi%jr.ne.1      ) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvj, ierr)
    endif
    if (iex.eq.-1) then
      if (mpi%ir.ne.1      ) bufrl(1:je-js+1,1:km) = fld(isnd,js:je,1:km)
      if (mpi%ir.ne.1      ) call mpi_isend( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
    else if (iex.eq.1) then
      if (mpi%ir.ne.mpi%irm) bufrl(1:je-js+1,1:km) = fld(isnd,js:je,1:km)
      if (mpi%ir.ne.mpi%irm) call mpi_isend( bufrl, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
      if (mpi%ir.ne.1      ) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
    endif

    if (jex.eq.-1) then
      if (mpi%jr.ne.1      ) call mpi_wait( isendj, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( irecvj, mpistat, ierr)
    else if (jex.eq.1) then
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( isendj, mpistat, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( irecvj, mpistat, ierr)
    endif
    if (iex.eq.-1) then
      if (mpi%ir.ne.1      ) call mpi_wait( isendi, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( irecvi, mpistat, ierr)
    else if (iex.eq.1) then
      if (mpi%ir.ne.mpi%irm) call mpi_wait( isendi, mpistat, ierr)
      if (mpi%ir.ne.1      ) call mpi_wait( irecvi, mpistat, ierr)
    endif

   if(imod.eq.1)then
    if (jex.eq.-1) then
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrec,1:km) = bufrt(1:ie-is+1,1:km)
    else if (jex.eq.1) then
      if (mpi%jr.ne.1      ) fld(is:ie,jrec,1:km) = bufrt(1:ie-is+1,1:km)
    endif
    if (iex.eq.-1) then
      if (mpi%ir.ne.mpi%irm) fld(irec,js:je,1:km) = bufrr(1:je-js+1,1:km)
    else if (iex.eq.1) then
      if (mpi%ir.ne.1      ) fld(irec,js:je,1:km) = bufrr(1:je-js+1,1:km)
    endif
   else
    if (jex.eq.-1) then
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrec,1:km) = fld(is:ie,jrec,1:km) + bufrt(1:ie-is+1,1:km)
    else if (jex.eq.1) then
      if (mpi%jr.ne.1      ) fld(is:ie,jrec,1:km) = fld(is:ie,jrec,1:km) + bufrt(1:ie-is+1,1:km)
    endif
    if (iex.eq.-1) then
      if (mpi%ir.ne.mpi%irm) fld(irec,js:je,1:km) = fld(irec,js:je,1:km) + bufrr(1:je-js+1,1:km)
    else if (iex.eq.1) then
      if (mpi%ir.ne.1      ) fld(irec,js:je,1:km) = fld(irec,js:je,1:km) + bufrr(1:je-js+1,1:km)
    endif
   endif

  endif

!-- 

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )

end subroutine exo_mpi
!------------------------------------------------------------------------------------
subroutine exa_mpi ( imod, isnl, isnr, ircl, ircr, jsnb, jsnt, jrcb, jrct, is, ie, js, je, km, fld)

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
! Exchange borders for MPI        
!                                                                      !
! Version 1: S.Dobricic 2011                                           !
!-----------------------------------------------------------------------


 use set_knd
 use mpi_str

 implicit none

 include 'mpif.h'

 INTEGER(i4)   :: imod
 INTEGER(i4)   :: isnl, isnr, ircl, ircr
 INTEGER(i4)   :: jsnb, jsnt, jrcb, jrct
 INTEGER(i4)   :: is, ie, js, je, km
 REAL(r8)      :: fld(is:ie,js:je,km)
 INTEGER(i4)   :: ierr, mpiszx, mpiszy
 INTEGER       :: isendl, isendr
 INTEGER       :: isendb, isendt
 INTEGER       :: irecvl, irecvr
 INTEGER       :: irecvb, irecvt
 INTEGER       :: mpistat(mpi_status_size)
 REAL(r8), allocatable :: bufrb(:,:)
 REAL(r8), allocatable :: bufrt(:,:)
 REAL(r8), allocatable :: bufrl(:,:)
 REAL(r8), allocatable :: bufrr(:,:)
 REAL(r8), allocatable :: bufsb(:,:)
 REAL(r8), allocatable :: bufst(:,:)
 REAL(r8), allocatable :: bufsl(:,:)
 REAL(r8), allocatable :: bufsr(:,:)


! ---

   ALLOCATE ( bufrb(ie-is+1,km), bufrt(ie-is+1,km) )
   ALLOCATE ( bufrl(je-js+1,km), bufrr(je-js+1,km) )
   ALLOCATE ( bufsb(ie-is+1,km), bufst(ie-is+1,km) )
   ALLOCATE ( bufsl(je-js+1,km), bufsr(je-js+1,km) )

     mpiszx = (ie-is+1) * km
     mpiszy = (je-js+1) * km

      if (mpi%ir.ne.1      ) bufsl(1:je-js+1,1:km) = fld(isnl,js:je,1:km)
      if (mpi%ir.ne.mpi%irm) bufsr(1:je-js+1,1:km) = fld(isnr,js:je,1:km)
      if (mpi%jr.ne.1      ) bufsb(1:ie-is+1,1:km) = fld(is:ie,jsnb,1:km)
      if (mpi%jr.ne.mpi%jrm) bufst(1:ie-is+1,1:km) = fld(is:ie,jsnt,1:km)

      if (mpi%ir.ne.1      ) call mpi_isend( bufsl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendl, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_isend( bufsr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendr, ierr)
      if (mpi%jr.ne.1      ) call mpi_isend( bufsb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendb, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_isend( bufst, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendt, ierr)

      if (mpi%ir.ne.1      ) call mpi_irecv( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvl, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvr, ierr)
      if (mpi%jr.ne.1      ) call mpi_irecv( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvb, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvt, ierr)

      if (mpi%ir.ne.1      ) call mpi_wait( isendl, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( isendr, mpistat, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( isendb, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( isendt, mpistat, ierr)

      if (mpi%ir.ne.1      ) call mpi_wait( irecvl, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( irecvr, mpistat, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( irecvb, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( irecvt, mpistat, ierr)

     if(imod.eq.1)then
      if (mpi%ir.ne.1      ) fld(ircl,js:je,1:km) = bufrl(1:je-js+1,1:km)
      if (mpi%ir.ne.mpi%irm) fld(ircr,js:je,1:km) = bufrr(1:je-js+1,1:km)
      if (mpi%jr.ne.1      ) fld(is:ie,jrcb,1:km) = bufrb(1:ie-is+1,1:km)
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrct,1:km) = bufrt(1:ie-is+1,1:km)
     else
      if (mpi%ir.ne.1      ) fld(ircl,js:je,1:km) = fld(ircl,js:je,1:km) + bufrl(1:je-js+1,1:km)
      if (mpi%ir.ne.mpi%irm) fld(ircr,js:je,1:km) = fld(ircr,js:je,1:km) + bufrr(1:je-js+1,1:km)
      if (mpi%jr.ne.1      ) fld(is:ie,jrcb,1:km) = fld(is:ie,jrcb,1:km) + bufrb(1:ie-is+1,1:km)
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrct,1:km) = fld(is:ie,jrct,1:km) + bufrt(1:ie-is+1,1:km)
     endif

!-- 

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )
   DEALLOCATE ( bufsb, bufst )
   DEALLOCATE ( bufsl, bufsr )

end subroutine exa_mpi
!------------------------------------------------------------------------------------
subroutine eao_mpi ( imod, isnl, isnr, ircl, ircr, jsnb, jsnt, jrcb, jrct, is, ie, js, je, km, fld)

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
! Ordered exchange borders for MPI     
!                                                                      !
! Version 1: S.Dobricic 2011                                           !
!-----------------------------------------------------------------------


 use set_knd
 use mpi_str

 implicit none

 include 'mpif.h'

 INTEGER(i4)   :: imod
 INTEGER(i4)   :: isnl, isnr, ircl, ircr
 INTEGER(i4)   :: jsnb, jsnt, jrcb, jrct
 INTEGER(i4)   :: is, ie, js, je, km
 REAL(r8)      :: fld(is:ie,js:je,km)
 INTEGER(i4)   :: ierr, mpiszx, mpiszy
 INTEGER       :: isendl, isendr
 INTEGER       :: isendb, isendt
 INTEGER       :: irecvl, irecvr
 INTEGER       :: irecvb, irecvt
 INTEGER       :: mpistat(mpi_status_size)
 REAL(r8), allocatable :: bufrb(:,:)
 REAL(r8), allocatable :: bufrt(:,:)
 REAL(r8), allocatable :: bufrl(:,:)
 REAL(r8), allocatable :: bufrr(:,:)
 REAL(r8), allocatable :: bufsb(:,:)
 REAL(r8), allocatable :: bufst(:,:)
 REAL(r8), allocatable :: bufsl(:,:)
 REAL(r8), allocatable :: bufsr(:,:)


! ---

   ALLOCATE ( bufrb(ie-is+1,km), bufrt(ie-is+1,km) )
   ALLOCATE ( bufrl(je-js+1,km), bufrr(je-js+1,km) )
   ALLOCATE ( bufsb(ie-is+1,km), bufst(ie-is+1,km) )
   ALLOCATE ( bufsl(je-js+1,km), bufsr(je-js+1,km) )

     mpiszx = (ie-is+1) * km
     mpiszy = (je-js+1) * km

      if (mpi%ir.ne.1      ) bufsl(1:je-js+1,1:km) = fld(isnl,js:je,1:km)
      if (mpi%ir.ne.mpi%irm) bufsr(1:je-js+1,1:km) = fld(isnr,js:je,1:km)
      if (mpi%ir.ne.1      ) call mpi_isend( bufsl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendl, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_isend( bufsr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendr, ierr)
      if (mpi%ir.ne.1      ) call mpi_irecv( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvl, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_irecv( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvr, ierr)
      if (mpi%ir.ne.1      ) call mpi_wait( isendl, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( isendr, mpistat, ierr)
      if (mpi%ir.ne.1      ) call mpi_wait( irecvl, mpistat, ierr)
      if (mpi%ir.ne.mpi%irm) call mpi_wait( irecvr, mpistat, ierr)
     if(imod.eq.1)then
      if (mpi%ir.ne.1      ) fld(ircl,js:je,1:km) = bufrl(1:je-js+1,1:km)
      if (mpi%ir.ne.mpi%irm) fld(ircr,js:je,1:km) = bufrr(1:je-js+1,1:km)
     else
      if (mpi%ir.ne.1      ) fld(ircl,js:je,1:km) = fld(ircl,js:je,1:km) + bufrl(1:je-js+1,1:km)
      if (mpi%ir.ne.mpi%irm) fld(ircr,js:je,1:km) = fld(ircr,js:je,1:km) + bufrr(1:je-js+1,1:km)
     endif

      if (mpi%jr.ne.1      ) bufsb(1:ie-is+1,1:km) = fld(is:ie,jsnb,1:km)
      if (mpi%jr.ne.mpi%jrm) bufst(1:ie-is+1,1:km) = fld(is:ie,jsnt,1:km)
      if (mpi%jr.ne.1      ) call mpi_isend( bufsb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendb, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_isend( bufst, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendt, ierr)
      if (mpi%jr.ne.1      ) call mpi_irecv( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvb, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_irecv( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvt, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( isendb, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( isendt, mpistat, ierr)
      if (mpi%jr.ne.1      ) call mpi_wait( irecvb, mpistat, ierr)
      if (mpi%jr.ne.mpi%jrm) call mpi_wait( irecvt, mpistat, ierr)
     if(imod.eq.1)then
      if (mpi%jr.ne.1      ) fld(is:ie,jrcb,1:km) = bufrb(1:ie-is+1,1:km)
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrct,1:km) = bufrt(1:ie-is+1,1:km)
     else
      if (mpi%jr.ne.1      ) fld(is:ie,jrcb,1:km) = fld(is:ie,jrcb,1:km) + bufrb(1:ie-is+1,1:km)
      if (mpi%jr.ne.mpi%jrm) fld(is:ie,jrct,1:km) = fld(is:ie,jrct,1:km) + bufrt(1:ie-is+1,1:km)
     endif

!-- 

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )
   DEALLOCATE ( bufsb, bufst )
   DEALLOCATE ( bufsl, bufsr )

end subroutine eao_mpi


subroutine gth_mpi ( img, jmg, k, km, fldin, fld )

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
! Gether a single level on the first processor                         !
!                                                                      !
! Version 1: S.Dobricic 2011                                           !
!-----------------------------------------------------------------------


 use set_knd
 use mpi_str
 use grd_str

 implicit none

 include 'mpif.h'

 INTEGER(i4)   :: img, jmg, k, km
 INTEGER(i4)   :: i, j, kk, kproc
 REAL(r4)      :: fld(img,jmg)
 REAL(r8)      :: fldin(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km)
 INTEGER       :: grd_npsm, ierr, izer

  REAL(r4), allocatable     ::  bff(:)
  REAL(r4), allocatable     ::  bffa(:,:)

! ---

   izer = 0
   grd_npsm = grd%npsm

  ALLOCATE ( bff(grd%npsm) )

  if(mpi%myrank.eq.0) then
     ALLOCATE ( bffa(grd%npsm,mpi%nproc) )
  else
     ALLOCATE ( bffa(1,1) )
  endif


     bff(:) = 0.0
  kk = 0
    do j=1,grd%jm
     do i=1,grd%im
      kk = kk + 1
      bff(kk) = fldin(i,j,k)
     enddo
    enddo


  call mpi_gather(bff, grd_npsm, mpi%r4, bffa, grd_npsm, mpi%r4, izer, mpi%comm, ierr)



  if(mpi%myrank.eq.0)then

   do kproc = 1,mpi%nproc
     kk = 0
     do j=grd%aj1(kproc),grd%ajm(kproc)
      do i=grd%ai1(kproc),grd%aim(kproc)
       kk = kk + 1
       fld(i,j) = bffa(kk,kproc)
      enddo
     enddo
    enddo

  endif

  deallocate ( bff )
  deallocate ( bffa )


end subroutine gth_mpi

subroutine gta_mpi ( imod, img, jmg, k, km, fldou, fld )

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
! Gether a single level on the first processor                         !
!                                                                      !
! Version 1: S.Dobricic 2011                                           !
!-----------------------------------------------------------------------


 use set_knd
 use mpi_str
 use grd_str

 implicit none

 include 'mpif.h'

 INTEGER(i4)   :: imod
 INTEGER(i4)   :: img, jmg, k, km
 INTEGER(i4)   :: i, j, kk, kproc
 REAL(r4)      :: fld(img,jmg)
 REAL(r8)      :: fldou(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km)
 INTEGER       :: grd_npsm, ierr

  REAL(r4), allocatable     ::  bff(:)
  REAL(r4), allocatable     ::  bffa(:,:)

! ---

  grd_npsm = grd%npsm

  ALLOCATE ( bff(grd%npsm) )

  if(mpi%myrank.eq.0) then
     ALLOCATE ( bffa(grd%npsm,mpi%nproc) )
  else
     ALLOCATE ( bffa(1,1) )
  endif

  if(mpi%myrank.eq.0)then

   do kproc = 1,mpi%nproc
     kk = 0
     do j=grd%aj1(kproc),grd%ajm(kproc)
      do i=grd%ai1(kproc),grd%aim(kproc)
       kk = kk + 1
       bffa(kk,kproc) = fld(i,j)
      enddo
     enddo
    enddo

  endif

  call mpi_scatter( bffa, grd_npsm, mpi%r4, bff, grd_npsm, mpi%r4, 0, mpi%comm, ierr)

 if(imod.eq.0)then

  kk = 0
    do j=1,grd%jm
     do i=1,grd%im
      kk = kk + 1
      fldou(i,j,k) = fldou(i,j,k) + bff(kk)
     enddo
    enddo

 else if(imod.eq.1)then

  kk = 0
    do j=1,grd%jm
     do i=1,grd%im
      kk = kk + 1
      fldou(i,j,k) = bff(kk)
     enddo
    enddo

 endif

  deallocate ( bff )
  deallocate ( bffa )


end subroutine gta_mpi
