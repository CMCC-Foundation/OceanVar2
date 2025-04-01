!======================================================================
!
! This file is part of Oceanvar.
!
!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version (GPL-3.0-or-later).
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!======================================================================
!-----------------------------------------------------------------------
!                                                                      !
!> Load SLA observations and unbias SLA                                    
!!
!!
!!
!                                                                      !
! Version 1:   Srdjan Dobricic 2006                                    !
! Version 1.1: Paolo Oddo      xxxx Modify the track recognition       !
! Version 1.2: Mario Adani     2023 Add choices on bias removal        !
!                                   Add events                         !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_sla

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : sla, obs

   IMPLICIT NONE

   INTEGER(i4)   ::  k

   sla%no = 0
   sla%nc = 0

   OPEN (511,FILE=drv%inpdir//'/sla_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ(511) sla%no

   WRITE (drv%dia,*) ' --- No of SLA obs: ',  sla%no, obs%sla

   IF ( sla%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( sla%ino(sla%no), sla%flg(sla%no), sla%flc(sla%no) )
   ALLOCATE ( sla%lon(sla%no), sla%lat(sla%no), sla%tim(sla%no) )
   ALLOCATE ( sla%val(sla%no), sla%bac(sla%no), sla%inc(sla%no) )
   ALLOCATE ( sla%bia(sla%no), sla%err(sla%no) )
   ALLOCATE ( sla%res(sla%no), sla%b_a(sla%no) )
   ALLOCATE ( sla%ib(sla%no), sla%jb(sla%no) )
   ALLOCATE ( sla%pb(sla%no), sla%qb(sla%no) )
   ALLOCATE ( sla%pq1(sla%no), sla%pq2(sla%no), sla%pq3(sla%no), sla%pq4(sla%no) )
   ALLOCATE ( sla%dpt(sla%no) )
   ALLOCATE ( sla%dtm(sla%no) )
   ALLOCATE ( sla%rss(sla%no), sla%ins(sla%no) )
   ALLOCATE ( sla%fls(sla%no) )
   ALLOCATE ( sla%ksat(sla%no),sla%eve(sla%no) )

! Initialise 
   sla%flc(:) = 1
   sla%inc(:) = 0.0_r8
   sla%bia(:) = 0.0_r8
   sla%b_a(:) = 0.0_r8
   sla%pq1(:) = 0.0_r8; sla%pq2(:) = 0.0_r8; sla%pq3(:) = 0.0_r8; sla%pq4(:) = 0.0_r8
   sla%dpt(:) = 0.0_r8
   sla%rss(:) = 0.0_r8; sla%ins(:) = 0.0_r8
   sla%fls(:) = 0_i8
   sla%rss(:) = 0.0_r8
   sla%ins(:) = 0.0_r8
   sla%fls(:) = 0.0_r8
   sla%eve(:) = 999_i8

! ---
! Level corresponding to the minimum depth
   sla%kdp = grd%km
   DO k = grd%km, 1, -1
      IF ( grd%dep(k) .GE. sla%dep ) sla%kdp = k
   ENDDO

   READ (511)                                                 &
       sla%ino(1:sla%no), sla%ksat(1:sla%no),sla%flg(1:sla%no)&
      ,sla%lon(1:sla%no), sla%lat(1:sla%no), sla%tim(1:sla%no)&
      ,sla%val(1:sla%no), sla%bac(1:sla%no)                   &
      ,sla%err(1:sla%no), sla%res(1:sla%no)                   &
      ,sla%ib(1:sla%no), sla%jb(1:sla%no)                     &
      ,sla%pb(1:sla%no), sla%qb(1:sla%no), sla%dtm(1:sla%no)
   CLOSE (511)

! ---
! Initialise quality flag
   IF (obs%sla .EQ. 0 ) THEN
      sla%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%sla
   ENDIF

   ! Remove bias based on sla%flg
   IF ( sla%unbias )  CALL unbias

   ! Count good observations
   sla%nc = 0
   DO k = 1,sla%no
      IF ( sla%flg(k) .EQ. 1 ) THEN
         sla%nc = sla%nc + 1
      ELSE
         sla%inc(k) = 0.0_r8
         sla%b_a(k) = 0.0_r8
         sla%pq1(k) = 0.0_r8
         sla%pq2(k) = 0.0_r8
         sla%pq3(k) = 0.0_r8
         sla%pq4(k) = 0.0_r8
      ENDIF
   ENDDO

   sla%flc(:) = sla%flg(:)
   WHERE( sla%flc .EQ. 0 ) sla%eve = 1

   WRITE (drv%dia,*) 'Number of good SLA observations after reading: ',  sla%nc

1111 CONTINUE

END SUBROUTINE get_obs_sla
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                              
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_sla

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k, ierr
   INTEGER(i4)   ::  i1, kk, i, j1, j
   INTEGER(i8)   ::  klev
   REAL(r8)      ::  p1, q1
   REAL(r8)      ::  msk4, div_x, div_y, rmn, dst, dstm
   REAL(r8)      ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
   REAL(r8)      ::  lat_lt_rot, lon_rb_rot

   rmn = 1.e-6

   IF ( sla%no .GT. 0 ) THEN

      sla%flc(:) = sla%flg(:)

      sla%dpt(:) = 0.0_r8

! ---
! Adjust longitudes
      IF ( grd%bwst .GT. 180. ) THEN
         DO k=1,sla%no
            IF ( sla%lon(k) .LT. 0.0_r8 ) THEN
               sla%lon(k) = sla%lon(k) + 360.
            ENDIF
         ENDDO
      ENDIF

! ---
! Interpolation parameters
      CALL int_obs_hor ( sla%no, sla%lat, sla%lon, sla%flc, sla%eve, sla%ib, sla%jb, sla%pb, sla%qb)

      DO kk = 1,sla%no
         IF ( sla%flc(kk) .EQ. 1 ) THEN
            i1 = sla%ib(kk)
            j1 = sla%jb(kk)
            sla%dpt(kk) = MAX(grd%hgt(i1,j1),MAX(grd%hgt(i1+1,j1),MAX(grd%hgt(i1,j1+1),grd%hgt(i1+1,j1+1))))
            msk4 = grd%msk(i1,j1,sla%kdp) + grd%msk(i1+1,j1,sla%kdp) + grd%msk(i1,j1+1,sla%kdp) + grd%msk(i1+1,j1+1,sla%kdp)
            IF ( msk4 .LT. 4.0 ) THEN
               sla%flc(kk) = 0
               sla%eve(kk) = 12
            ENDIF
         ENDIF
      ENDDO

! ---
! Horizontal interpolation parameters for each masked grid
      klev = 1
      DO k = 1,sla%no
         IF (sla%flc(k) .EQ. 1) THEN
            CALL int_obs_pq( sla%ib(k), sla%jb(k), klev, sla%pb(k), sla%qb(k),  &
                             sla%pq1(k), sla%pq2(k), sla%pq3(k), sla%pq4(k))
         ENDIF
      ENDDO

! ---
! Count good observations
      sla%nc = 0
      DO k = 1,sla%no
         IF ( sla%flc(k) .EQ. 1 ) THEN
            sla%nc = sla%nc + 1
         ENDIF
      ENDDO

      WRITE (drv%dia,*) 'Number of good SLA observations after parameter interpolation: ',  sla%nc

      sla%inc(:) = 0.0_r8

   ENDIF

END SUBROUTINE int_par_sla
!======================================================================
SUBROUTINE unbias

!-----------------------------------------------------------------------
!                                                                      !
!>  Remove SLA Bias                             
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Mario Adani     2023                                      !
!-----------------------------------------------------------------------


   USE set_knd
   USE obs_str, ONLY : sla
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)   ::  i1, i, iter, k
   REAL(r8)      ::  sumt, sumi, timp, dxx, dyy, dsm, satino, satk

   IF ( sla%bias_at ) THEN
! Remove along track bias ---------------------------------
      DO iter = 1,3

         sla%bia(:) = 0.0_r8
         timp   = sla%tim(1)
         satino = sla%ino(1)
         satk   = sla%ksat(1)
         i1     = 1
         DO k = 2,sla%no

            dxx = phy%re*phy%d2r * (sla%lon(k)-sla%lon(k-1)) * COS(sla%lat(k)*phy%d2r)
            dyy = phy%re*phy%d2r * (sla%lat(k)-sla%lat(k-1))

            IF ((  sla%ino(k) .NE. satino .OR. SQRT(dxx**2+dyy**2).GT.sla%dsm) .AND. k .NE. sla%no ) THEN
               sumt = 0.0_r8
               sumi = 0.0_r8

               DO i = i1,k-1
                  IF ( sla%flg(i) .EQ. 1 ) THEN
                     sumt = sumt + sla%res(i)
                     sumi = sumi + 1.0_r8
                  ENDIF
               ENDDO

               IF (sumi .GT. 0. ) sumt = sumt/sumi

               IF ( sumi .GT. sla%minobspt) THEN
                  DO i = i1,k-1
                     sla%res(i) = sla%res(i) - sumt
                     sla%bia(i) = sumt
                  ENDDO
               ELSE
                  DO i = i1,k-1
                     sla%res(i) = sla%res(i) - sumt
                     sla%bia(i) = sumt
                     sla%flg(i) = 0
                     sla%eve(i) = 14
                  ENDDO
               ENDIF

               timp   = sla%tim(k)
               satino = sla%ino(k)
               satk   = sla%ksat(k)
               i1     = k

            ELSEIF ( k .EQ. sla%no  ) THEN                      
               sumt = 0.0_r8                                   
               sumi = 0.0_r8                            
               DO i = i1,k                            
                  IF ( sla%flg(i) .EQ. 1 ) THEN        
                     sumt = sumt + sla%res(i)     
                     sumi = sumi + 1.0_r8       
                  ENDIF                        
               ENDDO                          
               IF ( sumi .GT. 0. ) sumt = sumt/sumi   

               IF ( sumi .GT. sla%minobspt) THEN
                  DO i = i1,k                      
                     sla%res(i) = sla%res(i) - sumt    
                     sla%bia(i) = sumt               
                  ENDDO                             
               ELSE
                  DO i = i1,k                        
                     sla%res(i) = sla%res(i) - sumt
                     sla%bia(i) = sumt
                     sla%flg(i) = 0
                     sla%eve(i) = 14
                  ENDDO
               ENDIF
               
            ENDIF   

         ENDDO ! sla%no

      ENDDO  ! iter

! Remove global track bias ---------------------------------
   ELSE   ! sla%bias_at

      sla%bia(:) = 0.0_r8
      sumt = 0.0_r8
      sumi = 0.0_r8
      DO i = 1,sla%no
         IF (sla%flg(i).EQ.1) THEN
            sumt = sumt + sla%res(i)
            sumi = sumi + 1.0_r8
         ENDIF
      ENDDO
      IF ( sumi .GT. 0._r8 ) sumt = sumt/sumi
      DO i = 1,sla%no
         sla%res(i) = sla%res(i) - sumt
         sla%bia(i) = sumt
      ENDDO

   ENDIF !sla%bias_at

END SUBROUTINE unbias
