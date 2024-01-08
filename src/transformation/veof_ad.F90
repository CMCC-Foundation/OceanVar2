subroutine veof_ad

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
! Vertical transformation (adjoint)                                    !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use grd_str
 use eof_str
 use mpi_str

 implicit none

 INTEGER(i4)             :: i, j, k, l, n, ierr
 REAL(r8), ALLOCATABLE  :: egm(:,:)

 ALLOCATE ( egm(grd%im,grd%jm) )

   grd%ro_ad(:,:,:) = 0.0


  do n=1,ros%neof

   egm(:,:) = 0.0

! Eta
  if(drv%bmd(drv%ktr) .ne. 1)then
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = egm(i,j) + ros%evc( i, j, 1,n) * grd%eta_ad(i,j)
#else
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), 1,n) * grd%eta_ad(i,j)
#endif
    enddo
   enddo
  endif



! 3D variables
  do k=1,grd%km
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = egm(i,j) + ros%evc( i, j, k+1       ,n) * grd%tem_ad(i,j,k) * grd%lcl(i,j,k)
     egm(i,j) = egm(i,j) + ros%evc( i, j, k+grd%km+1,n) * grd%sal_ad(i,j,k) * grd%lcl(i,j,k)
#else
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), k+1       ,n) * grd%tem_ad(i,j,k) * grd%lcl(i,j,k)
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), k+grd%km+1,n) * grd%sal_ad(i,j,k) * grd%lcl(i,j,k)
#endif
    enddo
   enddo
  enddo

   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
!     egm(i,j) = ros%eva( i, j, n) * egm(i,j) 
#else
!     egm(i,j) = ros%eva(grd%reg(i,j),n) * egm(i,j) 
#endif
    enddo
   enddo

!cdir serial
! 3D variables
  do l=n,ros%neof
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     grd%ro_ad(i,j,n) = grd%ro_ad(i,j,n) + egm(i,j) * ros%cor( i, j, n, l) 
#else
     grd%ro_ad(i,j,n) = grd%ro_ad(i,j,n) + egm(i,j) * ros%cor( grd%reg(i,j), n, l) 
#endif
    enddo
   enddo
  enddo
!cdir end serial

  enddo !n


 DEALLOCATE ( egm )


end subroutine veof_ad
