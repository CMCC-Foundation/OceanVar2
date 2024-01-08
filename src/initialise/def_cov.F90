 subroutine def_cov


!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Define filter constants, EOFs, etc.                                  !
!                                                                      !
! Version 1: S.Dobricic 2006 
! Version 2: S.Dobricic and R.Farina 2013                              !
!-----------------------------------------------------------------------

  use set_knd
  use drv_str
  use grd_str
  use eof_str
  use cns_str
  use mpi_str

  implicit none


  include 'mpif.h'

  REAL(r8)                    :: mpimx, mpimn

  INTEGER(i4)                 :: k, nspl, i, j, kk, ik, jk, iter,my_i
  INTEGER(i4)                 :: iw, kst, ken, kln, kln0, ierr
  REAL(r8)                    :: E, dst
  REAL(r8)    , ALLOCATABLE   :: sfct(:), al(:), bt(:), sc(:)
  REAL(r8)    , ALLOCATABLE   :: sc2(:,:),sc2a(:,:)
  INTEGER(i4) , ALLOCATABLE   :: jnxx(:)
  REAL(r8)                    :: rone, L(0:2)
  REAL(r8)                    :: sigma,q,b0,b1,b2,b3,b
  REAL(r8)                    :: mat_bc_vect(9)
  INTEGER                     :: irecvi, isendi 
  INTEGER,        ALLOCATABLE :: istatus(:)

 ALLOCATE ( istatus(mpi_status_size) )

  rone = 1.0

  if(rcf%loc.gt.rcf%L) then
    rcf%L = rcf%L * rcf%loc / sqrt( rcf%loc**2 - rcf%L**2 )
  endif

    L(0) = rcf%L
    L(1) = rcf%L *0.9
    L(2) = rcf%L*0.1

    mat_bc_vect(:)=0.0

! ---
! Recursive filter constants
!---------
! Create table

    do iter=1,2

    rcf%L = L(iter)
  
       nspl = max(grd%jmg,grd%img)
       ALLOCATE ( sfct(nspl), jnxx(nspl), al(nspl), bt(nspl) ) 

       sfct(:) = 0.0
       jnxx(:) = 0.0
       al(:) = 0.0
       bt(:) = 0.0

       rcf%ntb = min(20,min(grd%jmg,grd%img))

       ALLOCATE ( rcf%al(rcf%ntb))
       ALLOCATE ( rcf%sc(rcf%ntb))

       rcf%dsmn =  1.e20
       rcf%dsmx = -1.e20
       do j=1,grd%jm
        do i=1,grd%im
         rcf%dsmn = min(rcf%dsmn,min(grd%dx(i,j),grd%dy(i,j)))
         rcf%dsmx = max(rcf%dsmx,max(grd%dx(i,j),grd%dy(i,j)))
        enddo
       enddo

   if(mpi%nproc.gt.1) then
      call mpi_reduce( rcf%dsmx, mpimx, 1, mpi%r8, mpi_max, 0, mpi%comm, ierr)
      if(mpi%myrank==0)then
           rcf%dsmx = mpimx
      endif
           call mpi_bcast( rcf%dsmx, 1, mpi%r8, 0, mpi%comm, ierr)
      call mpi_reduce( rcf%dsmn, mpimn, 1, mpi%r8, mpi_min, 0, mpi%comm, ierr)
      if(mpi%myrank==0)then
           rcf%dsmn = mpimn
      endif
           call mpi_bcast( rcf%dsmn, 1, mpi%r8, 0, mpi%comm, ierr)
   endif


       rcf%dsmx = rcf%dsmx + max(rone,(rcf%dsmx-rcf%dsmn)/(rcf%ntb-2.))

       rcf%dsl = (rcf%dsmx-rcf%dsmn) / (rcf%ntb-1.)


       iw = rcf%ntb/mpi%nproc + 1
       kst = min( 1, rcf%ntb + 1)
       ken = min( kst + iw - 1, rcf%ntb)
        kln0 = (ken-kst+1)
       kst = min( mpi%myrank * iw + 1, rcf%ntb + 1)
       ken = min( kst + iw - 1, rcf%ntb)
       kln = ken - kst + 1
        ALLOCATE ( sc(kln))
        ALLOCATE ( sc2(kln0,mpi%nproc))
        ALLOCATE ( sc2a(kln0,mpi%nproc))
        kk = 0
       
        do k=kst,ken
           dst = rcf%dsmn + (k-1.) * rcf%dsl
           sfct(:) = 0.
           sfct(nspl/2+1) = 1.


           sigma=rcf%L/dst
           call def_coef(sigma,b1,b2,b3,b)
           call rcfl_2(nspl,sfct,b1,b2,b3,b)
           call rcfl_2_ad(nspl,sfct,b1,b2,b3,b)

           kk = kk + 1
           sc(kk) = sfct(nspl/2+1)
        enddo


     do k=1,mpi%nproc
      sc2(1:kln,k) = sc(1:kln)
     enddo
      call mpi_alltoall( sc2, kln0, mpi%r8, sc2a, kln0, mpi%r8, mpi%comm, ierr)
     do k=1,mpi%nproc
       kst = min( (k-1) * iw + 1, rcf%ntb + 1)
       ken = min( kst + iw - 1, rcf%ntb)
       kln = ken - kst + 1
       rcf%sc(kst:ken) = sc2a(1:kln,k)
     enddo

        DEALLOCATE ( sc, sc2, sc2a )


       DEALLOCATE ( sfct, jnxx, al, bt )


       do j=1,grd%jm
        do i=1,grd%im
         dst = ( grd%dx(i,j) - rcf%dsmn )/rcf%dsl
         k = int(dst) + 1
         dst = dst - real(k-1)
         grd%scx(i,j,iter) = sqrt( 1./ (rcf%sc(k)*(1.-dst) + rcf%sc(k+1)*dst) ) 
         dst = ( grd%dy(i,j) - rcf%dsmn )/rcf%dsl
         k = int(dst) + 1
         dst = dst - real(k-1)
         grd%scy(i,j,iter) = sqrt( 1./ (rcf%sc(k)*(1.-dst) + rcf%sc(k+1)*dst) ) 
        enddo
       enddo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do j=1,grd%jm
       do i=1,grd%im
        dst =grd%dx(i,j)      
        sigma= rcf%L/dst
        call  def_coef(sigma,grd%alx(i,j,iter),grd%btx(i,j,iter),grd%gmx(i,j,iter),grd%dlx(i,j,iter))
        call  def_coef_bc(mat_bc_vect,grd%alx(i,j,iter), grd%btx(i,j,iter), grd%gmx(i,j,iter))
        grd%mat_bc_x(:,i,j,iter)=mat_bc_vect(:)
     end do
     end do


     
     
      do j=1,grd%im
        do i=1,grd%jm
        dst =grd%dy(j,i)
        sigma= rcf%L/dst
        call  def_coef(sigma,grd%aly(i,j,iter),grd%bty(i,j,iter),grd%gmy(i,j,iter),grd%dly(i,j,iter))
        call  def_coef_bc(mat_bc_vect,grd%aly(i,j,iter), grd%bty(i,j,iter), grd%gmy(i,j,iter))
        grd%mat_bc_y(:,i,j,iter)=mat_bc_vect(:)
        end do
      end do  
   
  
  enddo 
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Mask Selection in x direction!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
        grd%tmx(:) = 0



      do k=1,grd%km
       do j=1,grd%jm

        do i=1-3,grd%im+3
         grd%tmx((k-1)*grd%jm+j) = grd%tmx((k-1)*grd%jm+j) + int(grd%msr(i,j,k))

        end do
        if(grd%tmx((k-1)*grd%jm+j).eq.(grd%im+6))then

            grd%tmx((k-1)*grd%jm+j)=1  !!!! all Ocean

        else if(grd%tmx((k-1)*grd%jm+j).gt.0 .and. grd%tmx((k-1)*grd%jm+j).lt.(grd%im+6))then

            grd%tmx((k-1)*grd%jm+j)=2   !!!Ocean and Costs 

        endif
       end do
      end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Mask Selection in y direction!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        grd%tmy(:) = 0



      do k=1,grd%km
       do j=1,grd%im

        do i=1-3,grd%jm+3

         grd%tmy((k-1)*grd%im+j) = grd%tmy((k-1)*grd%im+j) + int(grd%msr(j,i,k))

        end do
        if(grd%tmy((k-1)*grd%im+j).eq.(grd%jm+6))then

            grd%tmy((k-1)*grd%im+j)=1  !!!! all Ocean

        else if(grd%tmy((k-1)*grd%im+j).gt.0 .and. grd%tmy((k-1)*grd%im+j).lt.(grd%jm+6))then

            grd%tmy((k-1)*grd%im+j)=2   !!!Ocean and Costs

        endif
       end do
      end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     grd%scx(:,:,1) = grd%scx(:,:,1)*0.9 ! * 0.7 !grd%scx(i,j,2)/(grd%scx(i,j,1) + grd%scx(i,j,2)) 
     grd%scx(:,:,2) = 1.0
     grd%scy(:,:,1) = grd%scy(:,:,1)*0.9 ! * 0.7 !grd%scy(i,j,2)/(grd%scy(i,j,1) + grd%scy(i,j,2)) 
     grd%scy(:,:,2) = 1.0


       do k=1,grd%km
        do j=1,grd%jm
         do i=1,grd%im
          if(grd%msr(i,j,k).eq.1.0)then
           grd%fct(i,j,k) = 1.0  
          else
           grd%fct(i,j,k) = 0.0
          endif
         enddo
        enddo
       enddo


! ---
! Vertical EOFs

     ros%kmt = grd%km * 2 + 1

     call rdeofs

     ALLOCATE ( grd%ro( grd%im, grd%jm, ros%neof))
     ALLOCATE ( grd%ro_ad( grd%im, grd%jm, ros%neof))


   rcf%L = L(0)

  DEALLOCATE ( istatus )

 
end subroutine def_cov
