subroutine rcfl_x( im, jm, km, fld, alp, bta, gam, del, mat_bc_x)

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
! Recursive filter in x direction                                      !
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
 REAL(r8),     ALLOCATABLE :: a(:,:), b(:,:), c(:,:), c_s(:,:) , msk(:,:)
 REAL(r8),     ALLOCATABLE :: bffr(:), bffs(:)
 INTEGER,      ALLOCATABLE :: istatus(:)


  npnt = 3
  ias  = 3
  iae  = 3
    

    ALLOCATE ( a(im,jm*km), b(1-npnt:im,jm*km), c(im+npnt,jm*km), c_s(3,jm*km), msk(1-ias:im+iae,jm*km) )
    ALLOCATE ( v(3), u(3), y(3) )
    aLLOCATE ( bffr(jm*km*npnt), bffs(jm*km*npnt) )
    ALLOCATE ( istatus(mpi_status_size) )
    ALLOCATE ( jp(jm*km))
      
        ik = min(2_i4,km)
   
        a(:,:) = 0.0
        b(:,:) = 0.0
        c(:,:) = 0.0
        msk(:,:)=0.0
        v(:)=0.0
        u(:)=0.0
        c_s(:,:)=0.0

        do k=1,km
         do j=1,jm
          jp((k-1)*jm+j) = j
         enddo
        enddo


       do k=1,km
        do j=1,jm
      if( grd%tmx((k-1)*grd%jm+j).eq.1 .or. grd%tmx((k-1)*grd%jm+j).eq.2)then
         do i=1,im
            a(i,(k-1)*jm+j) = fld(i,j,k)
         enddo
      end if
      if(grd%tmx((k-1)*grd%jm+j).eq.2)then
         do i=1-ias,im+iae
          msk(i,(k-1)*jm+j) =grd%msr(i,j,k)
         enddo
      else if(grd%tmx((k-1)*grd%jm+j).eq.1)then
         do i=1-ias,im+iae
          msk(i,(k-1)*jm+j) = 1.0 ! grd%msr(i,j,k)
         enddo
      endif
        enddo
       enddo


! positive direction
     do lr = 1,mpi%thj(ik)

        js = grd%jrs(lr,ik)
        je = grd%jre(lr,ik)
    
      
        ! MPI receive
          if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) then

          call mpi_irecv( bffr, npnt*grd%jmr(lr,ik), mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
          call mpi_wait( irecvi, istatus, ierr)


            do j = js,je
             b(-2,j) = bffr(j-js+1                    )
             b(-1,j) = bffr(j-js+1 + grd%jmr(lr,ik)   )
             b( 0,j) = bffr(j-js+1 + grd%jmr(lr,ik)*2 )
            enddo

          end if 


        if(mpi%ir.eq.1) then

           is = 1

             do j=js,je
              b( 1,j) = del(1,jp(j))*a(1,j)*msk(1,j)
             end do
        else

           is = 0

        endif


      do j=js,je

       if(grd%tmx(j).eq.1)then

         do i=1+is,im
                b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ gam(i,jp(j))*b(i-3,j)+del(i,jp(j))*a(i,j)
         enddo
       elseif(grd%tmx(j).eq.2)then
         do i=1+is,im
           if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.0.0 )then
               b(i,j) = del(i,jp(j))*a(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j) .eq. 1.0 .and. msk(i-2,j).eq.0)then
                 b(i,j) = alp(i,jp(j))*b(i-1,j)+ del(i,jp(j))*a(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.1.0 .and. msk(i-2,j).eq.1.0 .and. msk(i-3,j).eq.0 )then
                 b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ del(i,jp(j))*a(i,j)
           else if( msk(i,j).eq.1.0 .and. msk(i-1,j).eq.1.0 .and. msk(i-2,j).eq. 1.0 .and. msk(i-3,j).eq.1 )then
                b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ gam(i,jp(j))*b(i-3,j)+del(i,jp(j))*a(i,j)
         endif

         enddo
       endif    
     enddo

! MPI send
         if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) then
            do j = js,je
             bffs(j-js+1                    ) = b(im-2,j)
             bffs(j-js+1 + grd%jmr(lr,ik)   ) = b(im-1,j)
             bffs(j-js+1 + grd%jmr(lr,ik)*2 ) = b(im  ,j)
            enddo
            call mpi_isend( bffs, npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
            call mpi_wait( isendi, istatus, ierr)
         endif


     enddo


  !negative direction
        do lr = 1,mpi%thj(ik)
  
            js = grd%jrs(lr,ik)
            je = grd%jre(lr,ik)
 
 ! MPI receive
         if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) then
               call mpi_irecv( bffr, npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
               call mpi_wait( irecvi, istatus, ierr)
                   do j=js,je
                      c(im+1,j) = bffr(j-js+1                    )
                      c(im+2,j) = bffr(j-js+1 + grd%jmr(lr,ik)   )
                      c(im+3,j) = bffr(j-js+1 + grd%jmr(lr,ik)*2 )
                   enddo
          endif


     is=0

     
        do j=js,je


          if(grd%tmx(j).eq.1)then

            do i=im-is,1,-1
             c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*c(i+3,j)+del(i,jp(j))*b(i,j)
            enddo

             c_s(1:3,j)= c(1:3,j)

          elseif(grd%tmx(j).eq.2)then

             if( msk(im+1,j).eq.1.0 .and. msk(im+2,j).eq.1.0 .and. msk(im+3,j).eq.0.0 )then
                v(2)=c(im+3,j)
             else if( msk(im+1,j).eq.1.0 .and. msk(im+2,j).eq.0.0 )then
                v(2)=c(im+2,j)
                v(3)=c(im+3,j)
             endif

            do i=im-is,1,-1
              if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.0.0 )then
                  coef=del(i,jp(j)) 
                do ii=1,3
                   u(ii)=b(i+1-ii,  j)*coef 
                   v(ii)=0.0  
                end do
               v(1)=v(1)+mat_bc_x(1,i,jp(j))*u(1)+mat_bc_x(2,i,jp(j))*u(2)+mat_bc_x(3,i,jp(j))*u(3)
               v(2)=v(2)+mat_bc_x(4,i,jp(j))*u(1)+mat_bc_x(5,i,jp(j))*u(2)+mat_bc_x(6,i,jp(j))*u(3)
               v(3)=v(3)+mat_bc_x(7,i,jp(j))*u(1)+mat_bc_x(8,i,jp(j))*u(2)+mat_bc_x(9,i,jp(j))*u(3)
                 c(i,j)   = v(1)
              else if(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.0.0) then
                 c(i,j) =  alp(i,jp(j))*C(i+1,j) + bta(i,jp(j))*v(2)+ gam(i,jp(j))*v(3)+del(i,jp(j))*b(i,j)
              else if(msk(i,j).eq.1.0.and.msk(i+1,j).eq.1.0.and.msk(i+2,j).eq.1.0.and.msk(i+3,j).eq.0.0) then
                 c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*v(2)+del(i,jp(j))*b(i,j)
              else if( msk(i,j).eq.1.0 .and. msk(i+1,j).eq.1.0 .and. msk(i+2,j).eq. 1.0 .and. msk(i+3,j).eq.1.0 )then
                 c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*c(i+3,j)+del(i,jp(j))*b(i,j)
              endif

            enddo

           c_s(1:3,j)= c(1:3,j)
          if( msk(1,j).eq.1.0 .and. msk(2,j).eq.1.0 .and. msk(3,j).eq.0.0 )then
           c_s(3,j)=v(2)
          else if( msk(1,j).eq.1.0 .and. msk(2,j).eq.0.0 )then
           c_s(3,j)=v(3)
           c_s(2,j)=v(2)
          endif

     endif

    enddo


! MPI send
          if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) then
             do j=js,je
              bffs(j-js+1                   ) = c_s(1,j)
              bffs(j-js+1 + grd%jmr(lr,ik)  ) = c_s(2,j)
              bffs(j-js+1 + grd%jmr(lr,ik)*2) = c_s(3,j)
             enddo
            call mpi_isend( bffs, npnt*grd%jmr(lr,ik), mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
            call mpi_wait( isendi, istatus, ierr)
        endif

      enddo

! 1000 continue
   

       do k=1,km
        do j=1,jm
        if( grd%tmx((k-1)*grd%jm+j).eq.1 .or. grd%tmx((k-1)*grd%jm+j).eq.2) then
         do i=1,im
          fld(i,j,k) = c(i,(k-1)*jm+j) ! *grd%msr(i,j,k)
         enddo
         endif
        enddo
       enddo


  DEALLOCATE ( jp )
  DEALLOCATE ( istatus )
  DEALLOCATE ( bffr, bffs )
  DEALLOCATE ( v, u, y )
  DEALLOCATE ( a, b, c, c_s, msk )
  

    
end subroutine rcfl_x
