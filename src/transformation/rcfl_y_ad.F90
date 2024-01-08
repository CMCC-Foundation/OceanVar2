subroutine rcfl_y_ad( im, jm, km, fld, alp, bta, gam, del, mat_bc_y)

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
! Recursive filter in y direction - adjoint
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use cns_str
 use mpi_str
 use grd_str

 implicit none

 include 'mpif.h'

  INTEGER(i4)    :: im, jm, km, jas, jae
  REAL(r8)       :: fld(im,jm,km)
  REAL(r8)       :: alp(jm,im), bta(jm,im)
  REAL(r8)       :: gam(jm,im), del(jm,im)
  REAL(r8)       :: mat_bc_y(9,jm,im)
  REAL(r8)       :: coef

 INTEGER(i4)    :: i,j,k, ktr, llr, lr, js, je, is, ik,ii
 INTEGER        :: irecvi, isendi, ierr, npnt


 INTEGER(i4),  ALLOCATABLE :: jp(:)
 REAL(r8),     ALLOCATABLE :: a(:,:), b(:,:), c(:,:),c_r(:,:), msk(:,:)
 REAL(r8),     ALLOCATABLE :: bffr(:), bffs(:)
 INTEGER,      ALLOCATABLE :: istatus(:)
 REAL(r8),     ALLOCATABLE :: v(:),u(:),y(:)


    npnt=3
    jas=3
    jae=3
    is=0

     ALLOCATE ( a(jm,im*km), b(-2:jm,im*km), c(1:jm+jae,im*km), c_r(3,im*km), msk(1-jas:jm+jae,im*km))
     ALLOCATE ( bffr(im*km*npnt), bffs(im*km*npnt) )
     ALLOCATE ( istatus(mpi_status_size) )
     ALLOCATE ( jp(im*km))
     ALLOCATE (v(3),u(3),y(3))





     ik = min(2_i4,km)

        a(:,:) = 0.0
        b(:,:) = 0.0
        c(:,:) = 0.0
        c_r(:,:)=0.0
        msk(:,:)=0.0
        v(:)=0.0
        u(:)=0.0
        is=0  


       do k=1,km
         do j=1,im
          jp((k-1)*im+j) = j
         enddo
      enddo


      do k=1,km
        do j=1,im
        if( grd%tmy((k-1)*grd%im+j).eq.1 .or. grd%tmy((k-1)*grd%im+j).eq.2)then
         do i=1,jm
            c(i,(k-1)*im+j) = fld(j,i,k)
         enddo
        endif
         if(grd%tmy((k-1)*grd%im+j).eq.2)then
           do i=1-jas,jm+jae
              msk(i,(k-1)*im+j) =grd%msr(j,i,k)
            enddo
         endif
        enddo
       enddo

  ! negative direction
  
     do lr = 1,mpi%thi(ik)

        js = grd%irs(lr,ik)
        je = grd%ire(lr,ik)


! MPI receive
         if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) then
            call mpi_irecv( bffr,npnt*grd%imr(lr,ik), mpi%r8, mpi%bot, 1, mpi%comm, irecvi, ierr)
            call mpi_wait( irecvi, istatus, ierr)
               do j = js,je
                c_r(1,j) =  bffr(j-js+1)
                c_r(2,j) =  bffr(j-js+1+grd%imr(lr,ik))
                c_r(3,j) =  bffr(j-js+1+2*grd%imr(lr,ik))
              enddo
          endif



       do j=js,je


        if( grd%tmy(j).eq.1) then
 
           c(1:3,j)= c(1:3,j)+c_r(1:3,j)

             do i=1,jm
              c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
              c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
              c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
              b(i, j)  = b(i  ,j) + del(i,jp(j))*c(i,j)
            end do   
     
        elseif( grd%tmy(j).eq.2) then

           v(:)=0.0
          if(msk(1,j).eq.1 .and.  msk(2,j).eq.1.0 .and. msk(3,j).eq.0.0 )then
           c(1,j)= c(1,j)+c_r(1,j)
           c(2,j)= c(2,j)+c_r(2,j)  
           v(2)= c_r(3,j)
          elseif( msk(1,j).eq.1.0 .and. msk(2,j).eq.0.0 ) then
           c(1,j)=c(1,j)+c_r(1,j)
           v(2)=c_r(2,j)
           v(3)=c_r(3,j)
          else 
           c(1:3,j)= c(1:3,j)+ msk(1:3,j)*c_r(1:3,j)
          endif

          do i=1,jm

           if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.0.0 )then

                 v(1) = v(1) + c(i,j)
                 u(:) = 0.0
                 u(1) =  mat_bc_y(1,i,jp(j))*v(1) + mat_bc_y(4,i,jp(j))*v(2) + mat_bc_y(7,i,jp(j))*v(3)
                 u(2) =  mat_bc_y(2,i,jp(j))*v(1) + mat_bc_y(5,i,jp(j))*v(2) + mat_bc_y(8,i,jp(j))*v(3)
                 u(3) =  mat_bc_y(3,i,jp(j))*v(1) + mat_bc_y(6,i,jp(j))*v(2) + mat_bc_y(9,i,jp(j))*v(3)
                      b(i,    j) = b(i,    j) + del(i,jp(j))*u(1)
                      b(i-1,  j) = b(i-1,  j) + del(i,jp(j))*u(2)
                      b(i-2,  j) = b(i-2,  j) + del(i,jp(j))*u(3)
                 v(:)=0.0  

          elseif(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.0.0) then

               c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
               v(2) = v(2) + bta(i,jp(j))*c(i,j)
               v(3) = v(3) + gam(i,jp(j))*c(i,j)
               b(i, j) = b(i  ,j) + del(i,jp(j))*c(i,j)

          elseif(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.1.0.and.msk(i+3,j).eq.0.0) then
               
               c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
               c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
               v(2) =  v(2)     + gam(i,jp(j))*c(i,j)
               b(i, j) = b(i,j)    + del(i,jp(j))*c(i,j)

          else if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.1.0 .and. msk(i+2,j).eq. 1.0 .and. msk(i+3,j).eq.1.0 )then
           c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
           c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
           c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
           b(i, j) =  b(i  ,j) + del(i,jp(j))*c(i,j)
          end if 

        end do

          if( msk(jm+1,j).eq.1.0 .and. msk(jm+2,j).eq.1.0 .and. msk(jm+3,j).eq.0.0 )then
            c(jm+3,j)=v(2)
          elseif(msk(jm+1,j).eq.1.0 .and. msk(jm+2,j).eq.0.0 )then
            c(jm+2,j)=v(2)
            c(jm+3,j)=v(3)
          end if                

      endif

  enddo




! MPI send
         if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) then
            do j = js,je
               bffs(j-js+1)                   =  c(jm+1,j)   
               bffs(j-js+1+grd%imr(lr,ik))    =  c(jm+2,j)   
               bffs(j-js+1+2* grd%imr(lr,ik)) =  c(jm+3,j)   
            end do
            call mpi_isend( bffs, npnt*grd%imr(lr,ik), mpi%r8, mpi%top, 1, mpi%comm, isendi, ierr)
            call mpi_wait( isendi, istatus, ierr)
         endif


     enddo


 ! go to 1000

     do lr = 1,mpi%thi(ik)

        js = grd%irs(lr,ik)
        je = grd%ire(lr,ik)


! MPI receive
         if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) then
            call mpi_irecv( bffr,npnt*grd%imr(lr,ik), mpi%r8, mpi%top, 1, mpi%comm, irecvi, ierr)
            call mpi_wait( irecvi, istatus, ierr)
               do j = js,je
                 b(jm,j)   =  b(jm,j)  +     bffr(j-js+1)
                 b(jm-1,j) =  b(jm-1,j)+   bffr(j-js+1+grd%imr(lr,ik))
                 b(jm-2,j) =  b(jm-2,j)+   bffr(j-js+1+2*grd%imr(lr,ik))
              enddo
          endif
          
      do j=js,je
      if( grd%tmy(j).eq.1)then

           do i=jm,1+is,-1

             b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
             b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
             b(i-3,j) = b(i-3,j) + gam(i,jp(j))*b(i,j)
             a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
           enddo

      elseif( grd%tmy(j).eq.2)then 


        do i=jm,1+is,-1

          if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.0.0 )then
               a(i,j) = a(i,j)+del(i,jp(j))*b(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j) .eq. 1.0 .and. msk(i-2,j).eq.0)then
                 b(i-1,j) =  b(i-1,j)+alp(i,jp(j))*b(i,j)
                 a(i  ,j) =  a(i  ,j) + del(i,jp(j))*b(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.1.0 .and. msk(i-2,j).eq.1.0 .and. msk(i-3,j).eq.0 )then
                 b(i-1,j) = b(i-1,j) +alp(i,jp(j))*b(i,j)
                 b(i-2,j) = b(i-2,j) +bta(i,jp(j))*b(i,j)
                 a(i  ,j) = a(i  ,j) +del(i,jp(j))*b(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.1.0 .and. msk(i-2,j).eq. 1.0 .and. msk(i-3,j).eq.1 )then
             b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
             b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
             b(i-3,j) = b(i-3,j) + gam(i,jp(j))*b(i,j)
             a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
           endif
           enddo
        endif
     end do


! MPI send
      if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) then
           do j = js,je
                 bffs(j-js+1)                  = b(0,j)  
                 bffs(j-js+1+grd%imr(lr,ik))   = b(-1,j) 
                 bffs(j-js+1+2*grd%imr(lr,ik)) = b(-2,j) 
           enddo
            call mpi_isend( bffs, npnt*grd%imr(lr,ik), mpi%r8, mpi%bot, 1, mpi%comm, isendi, ierr)
            call mpi_wait( isendi, istatus, ierr)
       endif

     enddo


 !    1000 continue


     do k=1,km
        do j=1,im
        if( grd%tmy((k-1)*grd%im+j).eq.1 .or. grd%tmy((k-1)*grd%im+j).eq.2)then
         do i=1,jm
          fld(j,i,k) = a(i,(k-1)*im+j) 
         enddo
        endif
        enddo
       enddo


    DEALLOCATE ( jp )
    DEALLOCATE ( istatus )
    DEALLOCATE ( bffr, bffs )
    DEALLOCATE ( v, u, y )
    DEALLOCATE ( a, b, c, c_r, msk )



end subroutine rcfl_y_ad
