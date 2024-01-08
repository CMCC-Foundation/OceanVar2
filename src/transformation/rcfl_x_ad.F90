subroutine rcfl_x_ad( im, jm, km, fld, alp, bta, gam, del, mat_bc_x)

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
! Recursive filter in x direction - adjoint
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use cns_str
 use mpi_str
 use grd_str

 implicit none

 include 'mpif.h'
 INTEGER(i4)    :: im, jm, km, ias, iae, ias_new, iae_new


 REAL(r8)       :: fld(im,jm,km)
 REAL(r8)       :: alp(im,jm), bta(im,jm)
 REAL(r8)       :: gam(im,jm), del(im,jm)!,del_a(im,jm)
 REAL(r8)       :: mat_bc_x(9,im,jm)
 REAL(r8)       :: coef
 INTEGER(i4)    :: i,j,k, ktr, llr, lr, js, je, is, ik,ii
 INTEGER        :: irecvi, isendi, ierr, npnt


 REAL(r8),     ALLOCATABLE :: v(:), u(:), y(:)
 INTEGER(i4),  ALLOCATABLE :: jp(:)
 REAL(r8),     ALLOCATABLE :: a(:,:), b(:,:), c(:,:), c_r(:,:),   msk(:,:)
 REAL(r8),     ALLOCATABLE :: bffr(:), bffs(:)
 INTEGER,      ALLOCATABLE :: istatus(:)


  npnt = 3
  ias  = 3
  iae  = 3
  is   = 0


    ALLOCATE ( a(1:im,jm*km), b(1-ias:im,jm*km), c(1:im+iae,jm*km), c_r(3,jm*km), msk(1-ias:im+iae,jm*km) )
    ALLOCATE ( v(3), u(3), y(3) )
    aLLOCATE ( bffr(jm*km*npnt), bffs(jm*km*npnt) )
    ALLOCATE ( istatus(mpi_status_size) )
    ALLOCATE ( jp(jm*km))

        ik = min(2_i4,km)

        a(:,:) = 0.0
        b(:,:) = 0.0
        c(:,:) = 0.0
        c_r(:,:) = 0.0
        msk(:,:)=0.0
        v(:)=0.0
        u(:)=0.0

        do k=1,km
         do j=1,jm
          jp((k-1)*jm+j) = j
         enddo
        enddo


       do k=1,km
        do j=1,jm
         if(grd%tmx((k-1)*grd%jm+j).eq.1 .or. grd%tmx((k-1)*grd%jm+j).eq.2)then
           do i=1,im
            c(i,(k-1)*jm+j) = fld(i,j,k)
           enddo
         endif
         if(grd%tmx((k-1)*grd%jm+j).eq.2)then
           do i=1-ias,im+iae
             msk(i,(k-1)*jm+j) =grd%msr(i,j,k)
           enddo
         end if
        enddo
       enddo




! negative direction

 
  do lr = 1,mpi%thj(ik)

        js = grd%jrs(lr,ik)
        je = grd%jre(lr,ik)

 ! MPI receive
          if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) then
              call mpi_irecv( bffr,npnt*grd%jmr(lr,ik),mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
              call mpi_wait( irecvi, istatus, ierr)
                do j = js,je
                  c_r(1,j) = bffr(j-js+1)
                  c_r(2,j) = bffr(j-js+1+grd%jmr(lr,ik))
                  c_r(3,j) = bffr(j-js+1+ 2*grd%jmr(lr,ik))
              enddo
          endif

     do j=js,je

       if(grd%tmx(j).eq.1)then
         
          c(1:3,j)= c(1:3,j)+c_r(1:3,j)
           do i=1,im
            c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
            c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
            c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
            b(i, j) =  b(i  ,j) + del(i,jp(j))*c(i,j)   
           enddo

       else if(grd%tmx(j).eq.2)then

            v(:)=0.0
          if(msk(1,j).eq.1 .and.  msk(2,j).eq.1.0 .and. msk(3,j).eq.0.0 )then
             c(1,j) = c(1,j)+c_r(1,j)
             c(2,j) = c(2,j)+c_r(2,j)
             v(2)   = c_r(3,j)
          elseif( msk(1,j).eq.1.0 .and. msk(2,j).eq.0.0 ) then
             c(1,j) = c(1,j)+c_r(1,j)
             v(2)   = c_r(2,j)
             v(3)   = c_r(3,j)
          else
             c(1:3,j)= c(1:3,j)+ msk(1:3,j)*c_r(1:3,j)
          endif

        do i=1,im
 
         if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.0.0 )then
 
                 v(1) = v(1) + c(i,j)
                 u(:) = 0.0
                 u(1) =  mat_bc_x(1,i,jp(j))*v(1) + mat_bc_x(4,i,jp(j))*v(2) + mat_bc_x(7,i,jp(j))*v(3)
                 u(2) =  mat_bc_x(2,i,jp(j))*v(1) + mat_bc_x(5,i,jp(j))*v(2) + mat_bc_x(8,i,jp(j))*v(3)
                 u(3) =  mat_bc_x(3,i,jp(j))*v(1) + mat_bc_x(6,i,jp(j))*v(2) + mat_bc_x(9,i,jp(j))*v(3)           
                 b(i,    j) = b(i,    j) + del(i,jp(j))*u(1)
                 b(i-1,  j) = b(i-1,  j) + del(i,jp(j))*u(2)
                 b(i-2,  j) = b(i-2,  j) + del(i,jp(j))*u(3)
                 v(:)=0.0

        elseif(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.0.0) then

               c(i+1,j)    =   c(i+1,j)     + alp(i,jp(j))*c(i,j)
               v(2)    =       v(2)         + bta(i,jp(j))*c(i,j)
               v(3)        =   v(3)         + gam(i,jp(j))*c(i,j)
               b(i, j) = b(i  ,j) + del(i,jp(j))*c(i,j)
 
        elseif(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.1.0.and.msk(i+3,j).eq.0.0) then
            
               c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
               c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
               v(2)     = v(2)     + gam(i,jp(j))*c(i,j)
               b(i, j) = b(i,j)    + del(i,jp(j))*c(i,j)

        else if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.1.0 .and. msk(i+2,j).eq. 1.0 .and. msk(i+3,j).eq.1.0 )then

               c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
               c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
               c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
               b(i, j) =  b(i  ,j) + del(i,jp(j))*c(i,j)

        endif

       end do

         if( msk(im+1,j).eq.1.0 .and. msk(im+2,j).eq.1.0 .and. msk(im+3,j).eq.0.0 )then
             c(im+3,j)=v(2)
         elseif(msk(im+1,j).eq.1.0 .and. msk(im+2,j).eq.0.0 )then
             c(im+2,j)=v(2)
             c(im+3,j)=v(3)
         end if

     end if

   end do


! MPI send
          if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) then
             do j = js,je
               bffs(j-js+1)                   =  c(im+1,j)   
               bffs(j-js+1+grd%jmr(lr,ik))    =  c(im+2,j)     
               bffs(j-js+1+2* grd%jmr(lr,ik)) =  c(im+3,j)    
              end do
              call mpi_isend( bffs, npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
              call mpi_wait( isendi, istatus, ierr)
           endif

   end do   

! go to 1000

! positive direction

     
      
     do lr = 1,mpi%thj(ik)

        js = grd%jrs(lr,ik)
        je = grd%jre(lr,ik)


! MPI receive
         if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) then
            call mpi_irecv( bffr,npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
            call mpi_wait( irecvi, istatus, ierr)
             do j = js,je
                 b(im,j)   =      b(im,j)+   bffr(j-js+1)
                 b(im-1,j) =      b(im-1,j)+ bffr(j-js+1+grd%jmr(lr,ik))
                 b(im-2,j) =      b(im-2,j)+ bffr(j-js+1+2*grd%jmr(lr,ik))
             enddo
        endif

         do j=js,je
           if(grd%tmx(j).eq.1)then
            do i=im,1,-1  
              b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
              b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
              b(i-3,j) = b(i-3,j) + gam(i,jp(j))*b(i,j)
              a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
            enddo
      elseif(grd%tmx(j).eq.2)then
        
        do i=im,1,-1
         
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
         if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) then
           do j = js,je
                 bffs(j-js+1)                  =  b(0,j)  
                 bffs(j-js+1+grd%jmr(lr,ik))   =  b(-1,j) 
                 bffs(j-js+1+2*grd%jmr(lr,ik)) =  b(-2,j) 
           enddo
            call mpi_isend( bffs,npnt*grd%jmr(lr,ik), mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
            call mpi_wait( isendi, istatus, ierr)
         endif

     enddo

    1000 continue    
   
       do k=1,km
        do j=1,jm
        if(grd%tmx((k-1)*grd%jm+j).eq.1 .or. grd%tmx((k-1)*grd%jm+j).eq.2)then
         do i=1,im
          fld(i,j,k) = a(i,(k-1)*jm+j) ! *grd%msr(i,j,k)
         enddo
        endif
        enddo
       enddo
    

    DEALLOCATE ( jp )
    DEALLOCATE ( istatus )
    DEALLOCATE ( bffr, bffs ) 
    DEALLOCATE ( a, b, c, msk)

end subroutine rcfl_x_ad
