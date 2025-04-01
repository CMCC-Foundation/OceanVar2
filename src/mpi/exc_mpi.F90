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
!> Exchange borders for MPI
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2011                                      !
!-----------------------------------------------------------------------
SUBROUTINE exo_mpi ( iord, imod, iex, isnd, irec, jex, jsnd, jrec, is, ie, js, je, km, fld)

   USE set_knd
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)   :: iord, imod, iex, isnd, irec, jex, jsnd, jrec
   INTEGER(i4)   :: is, ie, js, je, km
   REAL(r8)      :: fld(is:ie,js:je,km )
   INTEGER(i4)   :: ierr, mpiszx, mpiszy
   INTEGER       :: isendi, isendj
   INTEGER       :: irecvi, irecvj
   INTEGER       :: mpistat( mpi_status_size)
   REAL(r8), ALLOCATABLE :: bufrb(:,:)
   REAL(r8), ALLOCATABLE :: bufrt(:,:)
   REAL(r8), ALLOCATABLE :: bufrl(:,:)
   REAL(r8), ALLOCATABLE :: bufrr(:,:)

! ---

   ALLOCATE ( bufrb(ie-is+1,km ), bufrt(ie-is+1,km ) )
   ALLOCATE ( bufrl(je-js+1,km ), bufrr(je-js+1,km ) )

   mpiszx = (ie-is+1) * km
   mpiszy = (je-js+1) * km

   IF ( iord .EQ. 1 ) THEN

!-- Ordered exchange (first j, THEN i)
      IF ( jex .EQ. -1 ) THEN
         IF ( mpi%jr .NE. 1       ) bufrb(1:ie-is+1,1:km ) = fld(is:ie,jsnd,1:km )
         IF ( mpi%jr .NE. 1       ) CALL MPI_ISEND( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendj, ierr)
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvj, ierr)
         IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( isendj, mpistat, ierr)
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( irecvj, mpistat, ierr)
         IF ( mpi%jr .NE. mpi%jrm .AND. imod .EQ. 1 ) THEN
            fld(is:ie,jrec,1:km ) = bufrt(1:ie-is+1,1:km )
         ELSE IF ( mpi%jr .NE. mpi%jrm ) THEN
            fld(is:ie,jrec,1:km ) = fld(is:ie,jrec,1:km ) + bufrt(1:ie-is+1,1:km )
         ENDIF
      ELSE IF ( jex.EQ.1 ) THEN
         IF ( mpi%jr .NE. mpi%jrm ) bufrb(1:ie-is+1,1:km ) = fld(is:ie,jsnd,1:km )
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_ISEND( bufrb, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendj, ierr)
         IF ( mpi%jr .NE. 1       ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvj, ierr)
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( isendj, mpistat, ierr)
         IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( irecvj, mpistat, ierr)
         IF ( mpi%jr .NE. 1      .AND. imod .EQ. 1 ) THEN
            fld(is:ie,jrec,1:km ) = bufrt(1:ie-is+1,1:km )
         ELSE IF ( mpi%jr .NE. 1      ) THEN
            fld(is:ie,jrec,1:km ) = fld(is:ie,jrec,1:km ) + bufrt(1:ie-is+1,1:km )
         ENDIF
      ENDIF

      IF ( iex.EQ.-1 ) THEN
         IF ( mpi%ir .NE. 1       ) bufrl(1:je-js+1,1:km ) = fld(isnd,js:je,1:km )
         IF ( mpi%ir .NE. 1       ) CALL MPI_ISEND( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
         IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( isendi, mpistat, ierr)
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( irecvi, mpistat, ierr)
         IF ( mpi%ir .NE. mpi%irm .AND. imod .EQ. 1 ) THEN
            fld(irec,js:je,1:km ) = bufrr(1:je-js+1,1:km )
         ELSE IF ( mpi%ir .NE. mpi%irm ) THEN
            fld(irec,js:je,1:km ) = fld(irec,js:je,1:km ) + bufrr(1:je-js+1,1:km )
         ENDIF
      ELSE IF ( iex .EQ. 1 ) THEN
         IF ( mpi%ir .NE. mpi%irm ) bufrl(1:je-js+1,1:km ) = fld(isnd,js:je,1:km )
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_ISEND( bufrl, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
         IF ( mpi%ir .NE. 1       ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( isendi, mpistat, ierr)
         IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( irecvi, mpistat, ierr)
         IF ( mpi%ir .NE. 1       .AND. imod .EQ. 1 ) THEN
            fld(irec,js:je,1:km ) = bufrr(1:je-js+1,1:km )
         ELSE IF ( mpi%ir .NE. 1      ) THEN
            fld(irec,js:je,1:km ) = fld(irec,js:je,1:km ) + bufrr(1:je-js+1,1:km )
         ENDIF
      ENDIF

   ELSE

!-- Unordered exchange
      IF ( jex .EQ. -1 ) THEN
         IF ( mpi%jr .NE. 1       ) bufrb(1:ie-is+1,1:km ) = fld(is:ie,jsnd,1:km )
         IF ( mpi%jr .NE. 1       ) CALL MPI_ISEND( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendj, ierr)
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvj, ierr)
      ELSE IF ( jex .EQ.1 ) THEN
         IF ( mpi%jr .NE. mpi%jrm ) bufrb(1:ie-is+1,1:km ) = fld(is:ie,jsnd,1:km )
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_ISEND( bufrb, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendj, ierr)
         IF ( mpi%jr .NE. 1       ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvj, ierr)
      ENDIF
      IF ( iex .EQ. -1 ) THEN
         IF ( mpi%ir .NE. 1       ) bufrl(1:je-js+1,1:km ) = fld(isnd,js:je,1:km )
         IF ( mpi%ir .NE. 1       ) CALL MPI_ISEND( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendi, ierr)
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
      ELSE IF ( iex .EQ. 1 ) THEN
         IF ( mpi%ir .NE. mpi%irm ) bufrl(1:je-js+1,1:km ) = fld(isnd,js:je,1:km )
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_ISEND( bufrl, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendi, ierr)
         IF ( mpi%ir .NE. 1       ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
      ENDIF

      IF ( jex .EQ. -1 ) THEN
         IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( isendj, mpistat, ierr)
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( irecvj, mpistat, ierr)
      ELSE IF ( jex .EQ. 1 ) THEN
         IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( isendj, mpistat, ierr)
         IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( irecvj, mpistat, ierr)
      ENDIF
      IF ( iex .EQ. -1 ) THEN
         IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( isendi, mpistat, ierr)
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( irecvi, mpistat, ierr)
      ELSE IF ( iex .EQ. 1 ) THEN
         IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( isendi, mpistat, ierr)
         IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( irecvi, mpistat, ierr)
      ENDIF

      IF ( imod .EQ. 1 ) THEN
         IF ( jex .EQ. -1 ) THEN
            IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrec,1:km ) = bufrt(1:ie-is+1,1:km )
         ELSE IF ( jex .EQ. 1 ) THEN
            IF ( mpi%jr .NE. 1       ) fld(is:ie,jrec,1:km ) = bufrt(1:ie-is+1,1:km )
         ENDIF
         IF ( iex .EQ. -1 ) THEN
            IF ( mpi%ir .NE. mpi%irm ) fld(irec,js:je,1:km ) = bufrr(1:je-js+1,1:km )
         ELSE IF ( iex .EQ. 1 ) THEN
            IF ( mpi%ir .NE. 1       ) fld(irec,js:je,1:km ) = bufrr(1:je-js+1,1:km )
         ENDIF
      ELSE
         IF ( jex .EQ. -1 ) THEN
            IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrec,1:km ) = fld(is:ie,jrec,1:km ) + bufrt(1:ie-is+1,1:km )
         ELSE IF ( jex .EQ. 1 ) THEN
            IF ( mpi%jr .NE. 1       ) fld(is:ie,jrec,1:km ) = fld(is:ie,jrec,1:km ) + bufrt(1:ie-is+1,1:km )
         ENDIF
         IF ( iex .EQ. -1 ) THEN
            IF ( mpi%ir .NE. mpi%irm ) fld(irec,js:je,1:km ) = fld(irec,js:je,1:km ) + bufrr(1:je-js+1,1:km )
         ELSE IF (iex .EQ. 1 ) THEN
            IF ( mpi%ir .NE. 1       ) fld(irec,js:je,1:km ) = fld(irec,js:je,1:km ) + bufrr(1:je-js+1,1:km )
         ENDIF
      ENDIF

   ENDIF

!--

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )

END SUBROUTINE exo_mpi
!-----------------------------------------------------------------------
!                                                                      !
!> Exchange borders for MPI
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2011                                      !
!-----------------------------------------------------------------------
SUBROUTINE exa_mpi ( imod, isnl, isnr, ircl, ircr, jsnb, jsnt, jrcb, jrct, is, ie, js, je, km, fld)

   USE set_knd
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)   :: imod
   INTEGER(i4)   :: isnl, isnr, ircl, ircr
   INTEGER(i4)   :: jsnb, jsnt, jrcb, jrct
   INTEGER(i4)   :: is, ie, js, je, km
   REAL(r8)      :: fld(is:ie,js:je,km )
   INTEGER(i4)   :: ierr, mpiszx, mpiszy
   INTEGER       :: isendl, isendr
   INTEGER       :: isendb, isendt
   INTEGER       :: irecvl, irecvr
   INTEGER       :: irecvb, irecvt
   INTEGER       :: mpistat( mpi_status_size)
   REAL(r8), ALLOCATABLE :: bufrb(:,:)
   REAL(r8), ALLOCATABLE :: bufrt(:,:)
   REAL(r8), ALLOCATABLE :: bufrl(:,:)
   REAL(r8), ALLOCATABLE :: bufrr(:,:)
   REAL(r8), ALLOCATABLE :: bufsb(:,:)
   REAL(r8), ALLOCATABLE :: bufst(:,:)
   REAL(r8), ALLOCATABLE :: bufsl(:,:)
   REAL(r8), ALLOCATABLE :: bufsr(:,:)


! ---

   ALLOCATE ( bufrb(ie-is+1,km ), bufrt(ie-is+1,km ) )
   ALLOCATE ( bufrl(je-js+1,km ), bufrr(je-js+1,km ) )
   ALLOCATE ( bufsb(ie-is+1,km ), bufst(ie-is+1,km ) )
   ALLOCATE ( bufsl(je-js+1,km ), bufsr(je-js+1,km ) )

   mpiszx = (ie-is+1) * km
   mpiszy = (je-js+1) * km

   IF ( mpi%ir .NE. 1       ) bufsl(1:je-js+1,1:km ) = fld(isnl,js:je,1:km )
   IF ( mpi%ir .NE. mpi%irm ) bufsr(1:je-js+1,1:km ) = fld(isnr,js:je,1:km )
   IF ( mpi%jr .NE. 1       ) bufsb(1:ie-is+1,1:km ) = fld(is:ie,jsnb,1:km )
   IF ( mpi%jr .NE. mpi%jrm ) bufst(1:ie-is+1,1:km ) = fld(is:ie,jsnt,1:km )

   IF ( mpi%ir .NE. 1       ) CALL MPI_ISEND( bufsl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendl, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_ISEND( bufsr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendr, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_ISEND( bufsb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendb, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_ISEND( bufst, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendt, ierr)

   IF ( mpi%ir .NE. 1       ) CALL MPI_IRECV( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvl, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvr, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_IRECV( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvb, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvt, ierr)

   IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( isendl, mpistat, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( isendr, mpistat, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( isendb, mpistat, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( isendt, mpistat, ierr)

   IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( irecvl, mpistat, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( irecvr, mpistat, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( irecvb, mpistat, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( irecvt, mpistat, ierr)

   IF ( imod .EQ. 1 ) THEN
      IF ( mpi%ir .NE. 1       ) fld(ircl,js:je,1:km ) = bufrl(1:je-js+1,1:km )
      IF ( mpi%ir .NE. mpi%irm ) fld(ircr,js:je,1:km ) = bufrr(1:je-js+1,1:km )
      IF ( mpi%jr .NE. 1       ) fld(is:ie,jrcb,1:km ) = bufrb(1:ie-is+1,1:km )
      IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrct,1:km ) = bufrt(1:ie-is+1,1:km )
   ELSE
      IF ( mpi%ir .NE. 1       ) fld(ircl,js:je,1:km ) = fld(ircl,js:je,1:km ) + bufrl(1:je-js+1,1:km )
      IF ( mpi%ir .NE. mpi%irm ) fld(ircr,js:je,1:km ) = fld(ircr,js:je,1:km ) + bufrr(1:je-js+1,1:km )
      IF ( mpi%jr .NE. 1       ) fld(is:ie,jrcb,1:km ) = fld(is:ie,jrcb,1:km ) + bufrb(1:ie-is+1,1:km )
      IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrct,1:km ) = fld(is:ie,jrct,1:km ) + bufrt(1:ie-is+1,1:km )
   ENDIF

!--

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )
   DEALLOCATE ( bufsb, bufst )
   DEALLOCATE ( bufsl, bufsr )

END SUBROUTINE exa_mpi
!-----------------------------------------------------------------------
!                                                                      !
!> Ordered exchange borders for MPI
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2011                                      !
!-----------------------------------------------------------------------
SUBROUTINE eao_mpi ( imod, isnl, isnr, ircl, ircr, jsnb, jsnt, jrcb, jrct, is, ie, js, je, km, fld)

   USE set_knd
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)   :: imod
   INTEGER(i4)   :: isnl, isnr, ircl, ircr
   INTEGER(i4)   :: jsnb, jsnt, jrcb, jrct
   INTEGER(i4)   :: is, ie, js, je, km
   REAL(r8)      :: fld(is:ie,js:je,km )
   INTEGER(i4)   :: ierr, mpiszx, mpiszy
   INTEGER       :: isendl, isendr
   INTEGER       :: isendb, isendt
   INTEGER       :: irecvl, irecvr
   INTEGER       :: irecvb, irecvt
   INTEGER       :: mpistat( mpi_status_size)
   REAL(r8), ALLOCATABLE :: bufrb(:,:)
   REAL(r8), ALLOCATABLE :: bufrt(:,:)
   REAL(r8), ALLOCATABLE :: bufrl(:,:)
   REAL(r8), ALLOCATABLE :: bufrr(:,:)
   REAL(r8), ALLOCATABLE :: bufsb(:,:)
   REAL(r8), ALLOCATABLE :: bufst(:,:)
   REAL(r8), ALLOCATABLE :: bufsl(:,:)
   REAL(r8), ALLOCATABLE :: bufsr(:,:)


! ---

   ALLOCATE ( bufrb(ie-is+1,km ), bufrt(ie-is+1,km ) )
   ALLOCATE ( bufrl(je-js+1,km ), bufrr(je-js+1,km ) )
   ALLOCATE ( bufsb(ie-is+1,km ), bufst(ie-is+1,km ) )
   ALLOCATE ( bufsl(je-js+1,km ), bufsr(je-js+1,km ) )

   mpiszx = (ie-is+1) * km
   mpiszy = (je-js+1) * km

   IF ( mpi%ir .NE. 1       ) bufsl(1:je-js+1,1:km ) = fld(isnl,js:je,1:km )
   IF ( mpi%ir .NE. mpi%irm ) bufsr(1:je-js+1,1:km ) = fld(isnr,js:je,1:km )
   IF ( mpi%ir .NE. 1       ) CALL MPI_ISEND( bufsl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, isendl, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_ISEND( bufsr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, isendr, ierr)
   IF ( mpi%ir .NE. 1       ) CALL MPI_IRECV( bufrl, mpiszy, mpi%r8, mpi%lft, 1, mpi%comm, irecvl, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_IRECV( bufrr, mpiszy, mpi%r8, mpi%rgh, 1, mpi%comm, irecvr, ierr)
   IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( isendl, mpistat, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( isendr, mpistat, ierr)
   IF ( mpi%ir .NE. 1       ) CALL MPI_WAIT( irecvl, mpistat, ierr)
   IF ( mpi%ir .NE. mpi%irm ) CALL MPI_WAIT( irecvr, mpistat, ierr)
   IF ( imod .EQ. 1 ) THEN
      IF ( mpi%ir .NE. 1       ) fld(ircl,js:je,1:km ) = bufrl(1:je-js+1,1:km )
      IF ( mpi%ir .NE. mpi%irm ) fld(ircr,js:je,1:km ) = bufrr(1:je-js+1,1:km )
   ELSE
      IF ( mpi%ir .NE. 1       ) fld(ircl,js:je,1:km ) = fld(ircl,js:je,1:km ) + bufrl(1:je-js+1,1:km )
      IF ( mpi%ir .NE. mpi%irm ) fld(ircr,js:je,1:km ) = fld(ircr,js:je,1:km ) + bufrr(1:je-js+1,1:km )
   ENDIF

   IF ( mpi%jr .NE. 1       ) bufsb(1:ie-is+1,1:km ) = fld(is:ie,jsnb,1:km )
   IF ( mpi%jr .NE. mpi%jrm ) bufst(1:ie-is+1,1:km ) = fld(is:ie,jsnt,1:km )
   IF ( mpi%jr .NE. 1       ) CALL MPI_ISEND( bufsb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, isendb, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_ISEND( bufst, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, isendt, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_IRECV( bufrb, mpiszx, mpi%r8, mpi%bot, 1, mpi%comm, irecvb, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_IRECV( bufrt, mpiszx, mpi%r8, mpi%top, 1, mpi%comm, irecvt, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( isendb, mpistat, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( isendt, mpistat, ierr)
   IF ( mpi%jr .NE. 1       ) CALL MPI_WAIT( irecvb, mpistat, ierr)
   IF ( mpi%jr .NE. mpi%jrm ) CALL MPI_WAIT( irecvt, mpistat, ierr)
   IF ( imod .EQ. 1 ) THEN
      IF ( mpi%jr .NE. 1       ) fld(is:ie,jrcb,1:km ) = bufrb(1:ie-is+1,1:km )
      IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrct,1:km ) = bufrt(1:ie-is+1,1:km )
   ELSE
      IF ( mpi%jr .NE. 1       ) fld(is:ie,jrcb,1:km ) = fld(is:ie,jrcb,1:km ) + bufrb(1:ie-is+1,1:km )
      IF ( mpi%jr .NE. mpi%jrm ) fld(is:ie,jrct,1:km ) = fld(is:ie,jrct,1:km ) + bufrt(1:ie-is+1,1:km )
   ENDIF

!--

   DEALLOCATE ( bufrb, bufrt )
   DEALLOCATE ( bufrl, bufrr )
   DEALLOCATE ( bufsb, bufst )
   DEALLOCATE ( bufsl, bufsr )

END SUBROUTINE eao_mpi
!-----------------------------------------------------------------------
!                                                                      !
!> Gether a single level on the first processor                        
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic   2011                                    !
!            Francesco Carere  2024                                    !
!-----------------------------------------------------------------------
SUBROUTINE gth_mpi ( img, jmg, k, km, fldin, fld )

   USE set_knd
   USE mpi_str
   USE grd_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)   :: img, jmg, k, km
   INTEGER(i4)   :: i, j, kk, kproc
   REAL(r8)      :: fld(img,jmg)
   REAL(r8)      :: fldin(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km )
   INTEGER       :: grd_npsm, ierr, izer
   REAL(r8), ALLOCATABLE     ::  bff(:)
   REAL(r8), ALLOCATABLE     ::  bffa(:,:)

! ---

   izer = 0
   grd_npsm = grd%npsm

   ALLOCATE ( bff(grd%npsm ) )

   IF ( mpi%myrank .EQ. 0 ) THEN
      ALLOCATE ( bffa(grd%npsm,mpi%nproc) )
   ELSE
      ALLOCATE ( bffa(1,1) )
   ENDIF


   bff(:) = 0.0_r8
   kk = 0
   DO j = 1,grd%jm
      DO i = 1,grd%im
         kk = kk + 1
         bff(kk) = fldin(i,j,k)
      ENDDO
   ENDDO


   CALL MPI_GATHER(bff, grd_npsm, mpi%r8, bffa, grd_npsm, mpi%r8, izer, mpi%comm, ierr)



   IF ( mpi%myrank .EQ. 0 ) THEN

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

   DEALLOCATE ( bff )
   DEALLOCATE ( bffa )


END SUBROUTINE gth_mpi
!-----------------------------------------------------------------------
!                                                                      !
!> Gether a single level on the first processor                         
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2011                                      !    
!-----------------------------------------------------------------------
SUBROUTINE gta_mpi ( imod, img, jmg, k, km, fldou, fld )

   USE set_knd
   USE mpi_str
   USE grd_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)   :: imod
   INTEGER(i4)   :: img, jmg, k, km
   INTEGER(i4)   :: i, j, kk, kproc
   REAL(r8)      :: fld(img,jmg)
   REAL(r8)      :: fldou(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km )
   INTEGER       :: grd_npsm, ierr

   REAL(r8), ALLOCATABLE     ::  bff(:)
   REAL(r8), ALLOCATABLE     ::  bffa(:,:)

! ---

   grd_npsm = grd%npsm

   ALLOCATE ( bff(grd%npsm ) )

   IF ( mpi%myrank .EQ. 0 ) THEN
      ALLOCATE ( bffa(grd%npsm,mpi%nproc) )
   ELSE
      ALLOCATE ( bffa(1,1) )
   ENDIF

   IF ( mpi%myrank .EQ. 0 ) THEN
      DO kproc = 1,mpi%nproc
         kk = 0
         DO j = grd%jgsp(kproc),grd%jgep(kproc)
            DO i = grd%igsp(kproc),grd%igep(kproc)
               kk = kk + 1
               bffa(kk,kproc) = fld(i,j)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   CALL MPI_SCATTER( bffa, grd_npsm, mpi%r8, bff, grd_npsm, mpi%r8, 0, mpi%comm, ierr)

   IF ( imod .EQ. 0 ) THEN
      kk = 0
      DO j = 1,grd%jm
         DO i = 1,grd%im
            kk = kk + 1
            fldou(i,j,k) = fldou(i,j,k) + bff(kk)
         ENDDO
      ENDDO

   ELSE IF ( imod .EQ. 1 ) THEN
      kk = 0
      DO j = 1,grd%jm
         DO i = 1,grd%im
            kk = kk + 1
            fldou(i,j,k) = bff(kk)
         ENDDO
      ENDDO
   ENDIF

   DEALLOCATE ( bff )
   DEALLOCATE ( bffa )


END SUBROUTINE gta_mpi
