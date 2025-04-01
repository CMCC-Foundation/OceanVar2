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
!> Write outputs and diagnostics                                       
!!
!!
!!
!                                                                      !
! Version 1:   Srdjan Dobricic 2006                                    ! 
! Version 1.1: Paolo Oddo      2014                                    !     
! Version 1.2: Mario Adani     2024                                    !
!-----------------------------------------------------------------------
SUBROUTINE wrt_dia

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str
   USE eof_str
   USE ctl_str
   USE bmd_str
   USE mpi_str
   USE netcdf

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER                    :: status
   INTEGER                    :: ierr
   INTEGER                    :: fh
   INTEGER                    :: filetype
   INTEGER                    :: one
   INTEGER                    :: bffsz
   INTEGER(mpi_offset_KIND)   :: zero
   INTEGER, ALLOCATABLE       ::  map(:)
   INTEGER(i4)                :: i, j, k, kk
   CHARACTER                  :: fgrd
   INTEGER(i8)                :: maxno
   INTEGER(i8), ALLOCATABLE   ::  allflg(:)
   REAL(r8), ALLOCATABLE      ::  allinc(:), allb_a(:)
   INTEGER                    :: nobs
   INTEGER                    :: stat, ncid
   INTEGER                    :: idlon, idlat, iddep
   INTEGER                    :: ideta, idtem, idsal, iduvl, idvvl, idmsk

   WRITE (drv%dia,*) ' --------------------------------------------------'
   WRITE (drv%dia,*) ' Write corrections and statistics '
   WRITE (drv%dia,*) ' --------------------------------------------------'

   WRITE (fgrd,'(i1)')drv%ktr

! ---
! Innovations
   CALL write_inn(   1_i4, grd%eta, 'eta' )
   CALL write_inn( grd%km, grd%tem, 'tem' )
   CALL write_inn( grd%km, grd%sal, 'sal' )
   CALL write_inn( grd%km, grd%uvl, 'uvl' )
   CALL write_inn( grd%km, grd%vvl, 'vvl' )

! ---
! Observations

!Netcdf
   IF ( sla%no .GT. 0 ) THEN
      CALL wrt_sla
   ENDIF
   IF ( arg%no .GT. 0 ) THEN
      CALL wrt_arg
   ENDIF
   IF ( xbt%no .GT. 0 ) THEN
      CALL wrt_xbt
   ENDIF
   IF ( gld%no .GT. 0 ) THEN
      CALL wrt_gld
   ENDIF
   IF ( sst%no .GT. 0 ) THEN
      CALL wrt_sst
   ENDIF
   IF ( vdr%no .GT. 0 ) THEN
      WRITE (drv%dia,*)'vdr netCDF output not implemented yet'
!      CALL wrt_vdr
   ENDIF
   IF ( gvl%no .GT. 0 ) THEN
      WRITE (drv%dia,*)'gvl netCDF output not implemented yet'
!      CALL wrt_gvl
   ENDIF
   IF ( tra%no .GT. 0 ) THEN
      WRITE (drv%dia,*)'tra netCDF output not implemented yet'
!      CALL wrt_tra
   ENDIF
   IF ( trd%no .GT. 0 ) THEN
      WRITE (drv%dia,*)'trd netCDF output not implemented yet'
!      CALL wrt_trd
   ENDIF

!binary 
   IF ( mpi%myrank .EQ. 0 ) open(215,file='obs.dat',form='unformatted')

!ascii 
   IF ( sla%no .GT. 0 ) THEN
      IF ( mpi%myrank .EQ. 0 ) open(216,file='sla_stat.txt',form='formatted')
   ENDIF
   IF ( arg%no .GT. 0 ) THEN
      IF ( mpi%myrank .EQ. 0) open(217,file='arg_stat.txt',form='formatted')
   ENDIF
   IF ( xbt%no .GT. 0 ) THEN
      IF ( mpi%myrank .EQ. 0) open(218,file='xbt_stat.txt',form='formatted')
   ENDIF
   IF ( gld%no .GT. 0 ) THEN
      IF ( mpi%myrank .EQ. 0) open(219,file='gld_stat.txt',form='formatted')
   ENDIF
   IF (sst%no .GT. 0 ) THEN
      IF (mpi%myrank .EQ. 0) open(220,file='sst_stat.txt',form='formatted')
   ENDIF

   maxno = MAX(sst%no,MAX(sla%no,MAX(arg%no,MAX(xbt%no,gld%no))))
   maxno = MAX(gvl%no,MAX(vdr%no,maxno))

   IF ( maxno .GT. 0 ) ALLOCATE ( allflg(maxno), allinc(maxno), allb_a(maxno) )

   IF (sla%no .GT. 0) THEN
      allflg = 0_i8
      allinc = 0.0_r8
      allb_a = 0.0_r8
      nobs = sla%no
      CALL MPI_REDUCE(  sla%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  sla%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  sla%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) sla%no
      IF (sla%no.GT.0) THEN
         sla%flg(1:sla%no) = allflg(1:sla%no)
         sla%inc(1:sla%no) = allinc(1:sla%no)
         sla%b_a(1:sla%no) = allb_a(1:sla%no)
         WRITE (215)                                 &
            sla%ino(1:sla%no), sla%flg(1:sla%no),sla%lon(1:sla%no) &
            ,sla%lat(1:sla%no), sla%tim(1:sla%no),sla%val(1:sla%no) &
            ,sla%bac(1:sla%no), sla%err(1:sla%no),sla%res(1:sla%no) &
            ,sla%bia(1:sla%no), sla%inc(1:sla%no),sla%b_a(1:sla%no) &
            ,sla%dpt(1:sla%no), sla%dtm(1:sla%no)
         WRITE (216,"(i8)") sla%no
         DO kk=1,sla%no
            WRITE (216,"(2I8,12F12.5)")                 &
               sla%ino(kk),sla%flg(kk),sla%lon(kk),sla%lat(kk),sla%tim(kk) &
               ,sla%val(kk),sla%bac(kk),sla%err(kk),sla%res(kk),sla%bia(kk) &
               ,sla%inc(kk),sla%b_a(kk),sla%dpt(kk),sla%dtm(kk)
         ENDDO
      ENDIF
   ENDIF

   IF (arg%no.GT.0) THEN
      allflg = 0_i8
      allinc = 0.0_r8
      allb_a = 0.0_r8
      nobs = arg%no
      CALL MPI_REDUCE(  arg%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  arg%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  arg%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) arg%no
      IF (arg%no.NE.0)  THEN
         arg%flg(1:arg%no) = allflg(1:arg%no)
         arg%inc(1:arg%no) = allinc(1:arg%no)
         arg%b_a(1:arg%no) = allb_a(1:arg%no)
         WRITE (215)                                  &
            arg%ino(1:arg%no),arg%flg(1:arg%no),arg%par(1:arg%no) &
            ,arg%lon(1:arg%no),arg%lat(1:arg%no),arg%dpt(1:arg%no) &
            ,arg%tim(1:arg%no),arg%val(1:arg%no),arg%bac(1:arg%no) &
            ,arg%err(1:arg%no),arg%res(1:arg%no),arg%bia(1:arg%no) &
            ,arg%inc(1:arg%no),arg%b_a(1:arg%no)
         WRITE (217,"(i8)") arg%no
         DO kk=1,arg%no
            WRITE (217,"(3I8,10F12.5)")                  &
               arg%ino(kk),arg%flg(kk),arg%par(kk)         &
               ,arg%lon(kk),arg%lat(kk),arg%dpt(kk)         &
               ,arg%tim(kk),arg%val(kk),arg%bac(kk)         &
               ,arg%err(kk),arg%res(kk),arg%bia(kk)         &
               ,arg%inc(kk)
         ENDDO
      ENDIF
   ENDIF

   IF (xbt%no.GT.0) THEN
      allflg = 0_i8
      allinc = 0.0_r8
      allb_a = 0.0_r8
      nobs = xbt%no
      CALL MPI_REDUCE(  xbt%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  xbt%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  xbt%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) xbt%no
      IF (xbt%no.NE.0) THEN
         xbt%flg(1:xbt%no) = allflg(1:xbt%no)
         xbt%inc(1:xbt%no) = allinc(1:xbt%no)
         xbt%b_a(1:xbt%no) = allb_a(1:xbt%no)
         WRITE (215)                                  &
            xbt%ino(1:xbt%no),xbt%flg(1:xbt%no),xbt%par(1:xbt%no) &
            ,xbt%lon(1:xbt%no),xbt%lat(1:xbt%no),xbt%dpt(1:xbt%no) &
            ,xbt%tim(1:xbt%no),xbt%val(1:xbt%no),xbt%bac(1:xbt%no) &
            ,xbt%err(1:xbt%no),xbt%res(1:xbt%no),xbt%bia(1:xbt%no) &
            ,xbt%inc(1:xbt%no),xbt%b_a(1:xbt%no)
         WRITE (218,"(i8)") xbt%no
         DO kk=1,xbt%no
            WRITE (218,"(3I8,11F12.5)")                  &
               xbt%ino(kk),xbt%flg(kk),xbt%par(kk) &
               ,xbt%lon(kk),xbt%lat(kk),xbt%dpt(kk) &
               ,xbt%tim(kk),xbt%val(kk),xbt%bac(kk) &
               ,xbt%err(kk),xbt%res(kk),xbt%bia(kk) &
               ,xbt%inc(kk),xbt%b_a(kk)
         ENDDO
      ENDIF
   ENDIF

   IF (gld%no.GT.0) THEN
      allflg = 0_i8
      allinc = 0.0_r8
      allb_a = 0.0_r8
      nobs = gld%no
      CALL MPI_REDUCE(  gld%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  gld%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  gld%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) gld%no
      IF (gld%no.NE.0) THEN
         gld%flg(1:gld%no) = allflg(1:gld%no)
         gld%inc(1:gld%no) = allinc(1:gld%no)
         gld%b_a(1:gld%no) = allb_a(1:gld%no)
         WRITE (215)                                  &
            gld%ino(1:gld%no),gld%flg(1:gld%no),gld%par(1:gld%no) &
            ,gld%lon(1:gld%no),gld%lat(1:gld%no),gld%dpt(1:gld%no) &
            ,gld%tim(1:gld%no),gld%val(1:gld%no),gld%bac(1:gld%no) &
            ,gld%err(1:gld%no),gld%res(1:gld%no),gld%err(1:gld%no) &
            ,gld%inc(1:gld%no),gld%b_a(1:gld%no)
         WRITE (219,"(i8)") gld%no
         DO kk=1,gld%no
            WRITE (219,"(3I8,11F12.5)")                  &
               gld%ino(kk),gld%flg(kk),gld%par(kk) &
               ,gld%lon(kk),gld%lat(kk),gld%dpt(kk) &
               ,gld%tim(kk),gld%val(kk),gld%bac(kk) &
               ,gld%err(kk),gld%res(kk),gld%err(kk) &
               ,gld%inc(kk),gld%b_a(kk)
         ENDDO
      ENDIF
   ENDIF

   IF (sst%no.GT.0) THEN
      allflg(1:sst%no) = 0_i8
      allinc(1:sst%no) = 0.0_r8
      allb_a(1:sst%no) = 0.0_r8
      nobs = sst%no
      CALL MPI_REDUCE(  sst%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  sst%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  sst%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) sst%no
      IF (sst%no.NE.0) THEN
         sst%flg(1:sst%no) = allflg(1:sst%no)
         sst%inc(1:sst%no) = allinc(1:sst%no)
         sst%b_a(1:sst%no) = allb_a(1:sst%no)
         WRITE (215)                                  &
            sst%ino(1:sst%no), sst%flg(1:sst%no), sst%par(1:sst%no) &
            ,sst%lon(1:sst%no), sst%lat(1:sst%no)                    &
            ,sst%dpt(1:sst%no), sst%tim(1:sst%no)                    &
            ,sst%val(1:sst%no), sst%bac(1:sst%no)                    &
            ,sst%err(1:sst%no), sst%res(1:sst%no)                    &
            ,sst%ib(1:sst%no), sst%jb(1:sst%no), sst%kb(1:sst%no)    &
            ,sst%pb(1:sst%no), sst%qb(1:sst%no), sst%rb(1:sst%no)
         WRITE (220,"(i8)") sst%no
         DO kk=1,sst%no
            WRITE (220,"(3I8,11F12.5)") &
               sst%ino(kk),sst%flg(kk),sst%par(kk) &
               ,sst%lon(kk),sst%lat(kk),sst%dpt(kk) &
               ,sst%tim(kk),sst%val(kk),sst%bac(kk) &
               ,sst%err(kk),sst%res(kk),sst%bia(kk) &
               ,sst%inc(kk),sst%b_a(kk)
         ENDDO
      ENDIF
   ENDIF

   IF (vdr%no.NE.0) THEN
      allflg(1:vdr%no) = 0_i8
      allinc(1:vdr%no) = 0.0_r8
      allb_a(1:vdr%no) = 0.0_r8
      nobs = vdr%no
      CALL MPI_REDUCE(  vdr%flg, allflg, nobs, mpi%i8  ,   &
         MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  vdr%inc, allinc, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE(  vdr%b_a, allb_a, nobs, mpi%r8  ,   &
         MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF (mpi%myrank.EQ.0) THEN
      WRITE (215) vdr%no
      IF (vdr%no.NE.0) THEN
         vdr%flg(1:vdr%no) = allflg(1:vdr%no)
         vdr%inc(1:vdr%no) = allinc(1:vdr%no)
         vdr%b_a(1:vdr%no) = allb_a(1:vdr%no)
         WRITE (215)                                  &
            vdr%ino(1:vdr%no), vdr%flg(1:vdr%no), vdr%par(1:vdr%no) &
            ,vdr%lon(1:vdr%no), vdr%lat(1:vdr%no)                    &
            ,vdr%dpt(1:vdr%no), vdr%tim(1:vdr%no)                    &
            ,vdr%tms(1:vdr%no), vdr%tme(1:vdr%no)                    &
            ,vdr%val(1:vdr%no), vdr%bac(1:vdr%no)                    &
            ,vdr%err(1:vdr%no), vdr%res(1:vdr%no)                    &
            ,vdr%bia(1:vdr%no), vdr%inc(1:vdr%no)                    &
            ,vdr%b_a(1:vdr%no)
      ENDIF
   ENDIF

   IF ( gvl%no .NE. 0 ) THEN
      allflg(1:gvl%no) = 0_i8
      allinc(1:gvl%no) = 0.0_r8
      allb_a(1:gvl%no) = 0.0_r8
      nobs = gvl%no
      CALL MPI_REDUCE( gvl%flg, allflg, nobs, mpi%i8  ,   &
                       MPI_MAX, 0, mpi%comm, ierr)
      CALL MPI_REDUCE( gvl%inc, allinc, nobs, mpi%r8  ,   &
                       MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE( gvl%b_a, allb_a, nobs, mpi%r8  ,   &
                       MPI_SUM, 0, mpi%comm, ierr)
   ENDIF
   IF ( mpi%myrank .EQ. 0 ) THEN
      WRITE (215) gvl%no
      IF ( gvl%no .NE. 0) THEN
         gvl%flg(1:gvl%no) = allflg(1:gvl%no)
         gvl%inc(1:gvl%no) = allinc(1:gvl%no)
         gvl%b_a(1:gvl%no) = allb_a(1:gvl%no)
         WRITE (215)                                                 &
            gvl%ino(1:gvl%no), gvl%flg(1:gvl%no), gvl%par(1:gvl%no)  &
            ,gvl%lon(1:gvl%no), gvl%lat(1:gvl%no)                    &
            ,gvl%dpt(1:gvl%no), gvl%tim(1:gvl%no)                    &
            ,gvl%tms(1:gvl%no), gvl%tme(1:gvl%no)                    &
            ,gvl%val(1:gvl%no), gvl%bac(1:gvl%no)                    &
            ,gvl%err(1:gvl%no), gvl%res(1:gvl%no)                    &
            ,gvl%bia(1:gvl%no), gvl%inc(1:gvl%no)                    &
            ,gvl%b_a(1:gvl%no)
      ENDIF
   ENDIF

   IF ( mpi%myrank .EQ. 0 ) THEN
      WRITE (215) tra%no
      IF ( tra%no .NE. 0 ) THEN
         WRITE (215)                                                 &
            tra%dpt                                                  &
            ,tra%ino(1:tra%no), tra%flg(1:tra%no)                    &
            ,tra%loi(1:tra%no), tra%lai(1:tra%no)                    &
            ,tra%lof(1:tra%no), tra%laf(1:tra%no)                    &
            ,tra%lob(tra%nt+1,1:tra%no), tra%lab(tra%nt+1,1:tra%no)  &
            ,tra%rex(1:tra%no), tra%inx(1:tra%no)                    &
            ,tra%rey(1:tra%no), tra%iny(1:tra%no)                    &
            ,tra%loa(1:tra%no), tra%laa(1:tra%no)                    &
            ,tra%erx(1:tra%no), tra%ery(1:tra%no)
      ENDIF
   ENDIF

   IF ( mpi%myrank .EQ. 0 ) THEN
      WRITE (215) trd%no
      IF ( trd%no .NE. 0 ) THEN
         WRITE (215)                                                 &
            trd%dpt                                                  &
            ,trd%ino(1:trd%no), trd%flg(1:trd%no)                    &
            ,trd%loi(1:trd%no), trd%lai(1:trd%no)                    &
            ,trd%lof(1:trd%no), trd%laf(1:trd%no)                    &
            ,trd%lob(trd%nt+1,1:trd%no), trd%lab(trd%nt+1,1:trd%no)  &
            ,trd%rex(1:trd%no), trd%inx(1:trd%no)                    &
            ,trd%rey(1:trd%no), trd%iny(1:trd%no)                    &
            ,trd%loa(1:trd%no), trd%laa(1:trd%no)                    &
            ,trd%erx(1:trd%no), trd%ery(1:trd%no)
      ENDIF
   ENDIF

   IF ( mpi%myrank .EQ. 0 ) CLOSE (215)
   IF ( mpi%myrank .EQ. 0 ) CLOSE (216)
   IF ( mpi%myrank .EQ. 0 ) CLOSE (217)
   IF ( mpi%myrank .EQ. 0 ) CLOSE (218)
   IF ( mpi%myrank .EQ. 0 ) CLOSE (219)
   IF ( mpi%myrank .EQ. 0 ) CLOSE (220)

   IF ( maxno .GT. 0 ) DEALLOCATE ( allflg, allinc, allb_a )

END SUBROUTINE wrt_dia
!-----------------------------------------------------------------------
!
!> Write increments 
!!
!!
!!
! Version 1:   Srdjan Dobricic 2006                                   ! 
! Version 1.1: Paolo Oddo      2014                                   !     
! Version 1.2: Mario Adani     2024                                   !
!----------------------------------------------------------------------
SUBROUTINE write_inn( km, fldin, cpr)

   USE set_knd
   USE grd_str
   USE mpi_str
   USE netcdf

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER(i4)               :: km
   INTEGER(i4)               :: i, j, kk, k
   REAL(r8)                  ::  fldin(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km)
#ifdef REPRO
   REAL(r8), ALLOCATABLE     ::  bff(:)
   REAL(r8), ALLOCATABLE     ::  bffa(:,:)
   REAL(r8), ALLOCATABLE     ::  fld(:,:)
#else
   REAL(r4), ALLOCATABLE     ::  bff(:)
   REAL(r4), ALLOCATABLE     ::  bffa(:,:)
   REAL(r4), ALLOCATABLE     ::  fld(:,:)
#endif
   INTEGER                    :: ierr, kproc
   INTEGER                    :: stat, ncid, idvar, idlon, idlat, iddep
   CHARACTER*3                :: cpr

   IF ( mpi%myrank .EQ. 0 ) ALLOCATE ( fld(grd%img,grd%jmg) )

   ALLOCATE ( bff(grd%npsm), bffa(grd%npsm,mpi%nproc) )

   DO k = 1,km

      IF ( mpi%nproc .GT. 1 ) THEN

#ifdef REPRO
         bff(:) = 0.0_r8
#else
         bff(:) = 0.0_r4
#endif
         kk = 0
         DO j = 1,grd%jm
            DO i = 1,grd%im
               kk = kk + 1
#ifdef REPRO
               bff(kk) = fldin(i,j,k)
#else
               bff(kk) = SNGL(fldin(i,j,k))
#endif
            ENDDO
         ENDDO

#ifdef REPRO
         CALL MPI_GATHER(bff, grd%npsm, MPI_DOUBLE_PRECISION, bffa, grd%npsm, MPI_DOUBLE_PRECISION, 0, mpi%comm, ierr)
#else
         CALL MPI_GATHER(bff, grd%npsm, MPI_REAL, bffa, grd%npsm, MPI_REAL, 0, mpi%comm, ierr)
#endif

         IF ( mpi%myrank .EQ. 0 ) THEN
#ifdef REPRO
            fld(:,:) = 0.0_r8
#else
            fld(:,:) = 0.0_r4
#endif
            DO kproc = 1,mpi%nproc
               kk = 0
               DO j = grd%jgsp(kproc),grd%jgep(kproc)
                  DO i = grd%igsp(kproc),grd%igep(kproc)
                     kk = kk + 1
                     fld(i,j) = bffa(kk,kproc)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ELSE

#ifdef REPRO
         fld(:,:) = fldin(:,:,k)
#else
         fld(:,:) = SNGL(fldin(:,:,k))
#endif
      ENDIF

      IF ( mpi%myrank .EQ. 0 ) THEN

         IF ( k .EQ. 1 ) THEN
            stat = NF90_CREATE('corr_'//cpr//'.nc', NF90_SHARE, ncid)
            stat = NF90_OPEN  ('corr_'//cpr//'.nc', NF90_WRITE, ncid)

            stat = NF90_REDEF(ncid)

            stat = NF90_DEF_DIM(ncid, 'im', grd%img-2*grd%iex, idlon)
            stat = NF90_DEF_DIM(ncid, 'jm', grd%jmg, idlat)
            stat = NF90_DEF_DIM(ncid, 'km', grd%km , iddep)

#ifdef REPRO
            IF ( km .EQ. 1 ) THEN
               stat = NF90_DEF_VAR ( ncid, cpr, NF90_DOUBLE, (/ idlon, idlat /), idvar)
            ELSE
               stat = NF90_DEF_VAR ( ncid, cpr, NF90_DOUBLE, (/ idlon, idlat, iddep /), idvar)
            ENDIF
#else
            IF ( km .EQ. 1 ) THEN
               stat = NF90_DEF_VAR ( ncid, cpr, NF90_FLOAT, (/ idlon, idlat /), idvar)
            ELSE
               stat = NF90_DEF_VAR ( ncid, cpr, NF90_FLOAT, (/ idlon, idlat, iddep /), idvar)
            ENDIF
#endif

            stat = NF90_ENDDEF(ncid)
            IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
         ENDIF

         IF ( km .EQ. 1 ) THEN
            stat = NF90_PUT_VAR ( ncid, idvar,                                 &
                                  fld(grd%iex+1:grd%img-grd%iex,1:grd%jmg),    &
                                  start = (/ 1, 1 /),                          &
                                  count = (/ grd%img-2*grd%iex, grd%jmg/) )
         ELSE
            stat = NF90_PUT_VAR ( ncid, idvar,                                 &
                                  fld(grd%iex+1:grd%img-grd%iex,1:grd%jmg),    &
                                  start = (/ 1, 1, k /),                       &
                                  count = (/ grd%img-2*grd%iex, grd%jmg, 1/) )
         ENDIF

         IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

         IF ( k .EQ. km ) THEN
            stat = NF90_CLOSE (ncid)
         ENDIF

      ENDIF

   ENDDO ! k

   DEALLOCATE ( bff, bffa )

   IF ( mpi%myrank .EQ. 0) DEALLOCATE ( fld )

END SUBROUTINE write_inn

