subroutine int_obs_hor ( no, olat, olon, flc, ib, jb, pb, qb)

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use mpi_str

 implicit none


  integer(i8)   ::  no
  real(r8)      ::  olat(no), olon(no), pb(no), qb(no)
  integer(i8)   ::  flc(no), ib(no), jb(no)

  integer(i4)   ::  k, ierr, kks, kkp
  integer(i4)   ::  i1, kk, i, j1, j, if1, jf1
  real(r8)      ::  p1, q1, pf1, qf1
  real(r8)      ::  msk4, div_x, div_y, rmn, dst, dstm
  logical       ::  ins
  real(r8)      ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
  real(r8)      ::  lat_lt_rot, lon_rb_rot

  ins(i,i1) = i.ge.1 .and. i.lt.i1

  rmn = 1.e-6

 if(grd%prj.eq.0)then

! Equally distant lat-lon grid 

    do kk = 1,no
     q1 = (olat(kk) - grd%lat(1,1)) / grd%dlt + 1.0
     j1 = int(q1)
     p1 = (olon(kk) - grd%lon(1,1)) / grd%dln + 1.0
     i1 = int(p1)
     if(ins(j1,grd%jm+grd%jae) .and. ins(i1,grd%im+grd%iae)) then
       ib(kk) = i1
       jb(kk) = j1
       pb(kk) = max((p1-i1),rmn)
       qb(kk) = max((q1-j1),rmn)
     else
       flc(kk) = 0
     endif
    enddo

 else

! General rotated grid

    do kk = 1,no
     if(olat(kk).ge.grd%bnrt .or. olat(kk).lt.grd%bsth .or.   &
        olon(kk).ge.grd%beas .or. olon(kk).lt.grd%bwst ) then
       flc(kk) = 0
     endif
    enddo

      kkp = 1
      kks = 0

  do kk = 1,no

   if(flc(kk).eq.1 .and. flc(kkp).eq.1 .and. kks.gt.0 .and.          &
          olat(kk).eq.olat(kkp) .and. olon(kk).eq.olon(kkp)) then

       ib(kk) = ib(kkp)
       jb(kk) = jb(kkp)
       pb(kk) = pb(kkp)
       qb(kk) = qb(kkp)
        
   else if(flc(kk).eq.1) then


        dstm = 1.e20
      do j=1,grd%jm
      do i=1,grd%im
        dst = (olat(kk)-grd%lat(i,j))**2 +     &
              ((olon(kk)-grd%lon(i,j))*cos(olat(kk)*3.14/180.))**2
         if(dst.lt.dstm)then
          i1 = i
          j1 = j
          dstm = dst
         endif
      enddo
      enddo


         if1 = -99
         jf1 = -99


      do j=max(j1-1,1),j1
      do i=max(i1-1,1),i1

        tga = ((grd%lon(i,j+1)-grd%lon(i,j))/(grd%lat(i,j+1)-grd%lat(i,j)))
        ang = atan (tga)

        lon_lb_rot = grd%lon(i,j)*cos(ang) - grd%lat(i,j)*sin(ang)
        lat_lb_rot = grd%lon(i,j)*sin(ang) + grd%lat(i,j)*cos(ang)
        lon_rb_rot = grd%lon(i+1,j)*cos(ang) - grd%lat(i+1,j)*sin(ang)
        lat_lt_rot = grd%lon(i,j+1)*sin(ang) + grd%lat(i,j+1)*cos(ang)
        lon_rot = olon(kk) * cos(ang) - olat(kk) * sin(ang)
        lat_rot = olon(kk) * sin(ang) + olat(kk) * cos(ang)

         q1 = (lat_rot - lat_lb_rot) / (lat_lt_rot - lat_lb_rot)
         p1 = (lon_rot - lon_lb_rot) / (lon_rb_rot - lon_lb_rot)

       if(q1.ge.0. .and. q1.lt. 1.0 .and. p1.ge.0. .and. p1.lt.1.0 )then
          if1 = i
          jf1 = j
          pf1 = p1
          qf1 = q1
       endif
      enddo
      enddo

     if(ins(jf1,grd%jm+grd%jae) .and. ins(if1,grd%im+grd%iae)) then
       ib(kk) = if1
       jb(kk) = jf1
       pb(kk) = max(pf1,rmn)
       qb(kk) = max(qf1,rmn)
     else
       flc(kk) = 0
     endif

   endif

      kkp = kk 
      kks = kk

  enddo

 endif


end subroutine int_obs_hor

subroutine int_obs_pq( i1, j1, k1, p1, q1, pq1, pq2, pq3, pq4)

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str

 implicit none

  integer(i8)   ::  i1, j1, k1
  real(r8)      ::  p1, q1
  real(r8)      ::  pq1, pq2, pq3, pq4
  real(r8)      ::  div_x, div_y

   div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
           +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
   div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
     pq1 = grd%msk(i1,j1,k1)                                      &
         * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
            * (1.-p1) * (1.-q1)                                   &
         /( div_x * div_y + 1.e-16 )
     pq2 = grd%msk(i1+1,j1,k1)                                    &
         * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
           *     p1  * (1.-q1)                                    &
         /( div_x * div_y + 1.e-16 )
   div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
     pq3 = grd%msk(i1,j1+1,k1)                                    &
         * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
           * (1.-p1) *     q1                                     &
         /( div_x * div_y + 1.e-16 )
     pq4 = grd%msk(i1+1,j1+1,k1)                                  &
         * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
           *     p1  *     q1                                     &
         /( div_x * div_y + 1.e-16 )


end subroutine int_obs_pq
