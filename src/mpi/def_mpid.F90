subroutine def_mpid

 use set_knd
 use grd_str
 use mpi_str

 implicit none

 include 'mpif.h'

 integer(i4) :: iw, jw, kproc, ierr, lr
 integer(i4), ALLOCATABLE  :: ibff (:)
 integer(i4), ALLOCATABLE  :: ibffa(:,:)

 INTEGER, ALLOCATABLE       :: lbf(:,:), lbfa(:,:,:)


    if(mpi%nproc.eq.1)then

      mpi%jr  = 1
      mpi%ir  = 1

      grd%igs = 1
      grd%ige = grd%img
      grd%im  = grd%img

      grd%ias = 0
      grd%iae = 0

      grd%jgs = 1
      grd%jge = grd%jmg
      grd%jm  = grd%jmg

      grd%jas = 0
      grd%jae = 0

    mpi%top = MPI_PROC_NULL
    mpi%bot = MPI_PROC_NULL
    mpi%rgh = MPI_PROC_NULL
    mpi%lft = MPI_PROC_NULL

   ALLOCATE ( mpi%thi(2), mpi%thj(2) )

   mpi%thi(1) = 1
   mpi%thi(2) = grd%km
   mpi%thj(1) = 1
   mpi%thj(2) = grd%km

      grd%is2 = 2
      grd%js2 = 2
      grd%is3 = 3
      grd%js3 = 3

    else

! ---
! Find the position of the processor on the 2D grid
      mpi%jr  = mpi%myrank / mpi%irm + 1
      mpi%ir  =  mpi%myrank + 1 - (mpi%jr - 1) * mpi%irm
      
! ---
! Starting and end points on the global grid, dimensions of the tile
      call para_reg( grd%img, mpi%irm, mpi%ir, grd%igs, grd%ige, grd%im)

      grd%ias = 1
      grd%iae = 1
    if (mpi%ir.eq.1      ) grd%ias = 0
    if (mpi%ir.eq.mpi%irm) grd%iae = 0

    if (mpi%ir.eq.1      ) then
      grd%is2 = 2
      grd%is3 = 3
    else
      grd%is2 = 1+mod(grd%igs,2)
      grd%is3 = 2-mod(grd%igs,2)
    endif
    if (mpi%ir.eq.mpi%irm) then
      grd%ipe = 1
    else
      grd%ipe = mod(grd%ige,2)
    endif

      call para_reg( grd%jmg, mpi%jrm, mpi%jr, grd%jgs, grd%jge, grd%jm)


      grd%jas = 1
      grd%jae = 1
    if (mpi%jr.eq.mpi%jrm) grd%jae = 0
    if (mpi%jr.eq.1      ) grd%jas = 0

    if (mpi%jr.eq.1      ) then
      grd%js2 = 2
      grd%js3 = 3
    else
      grd%js2 = 1+mod(grd%jgs,2)
      grd%js3 = 2-mod(grd%jgs,2)
    endif
    if (mpi%jr.eq.mpi%jrm) then
      grd%jpe = 1
    else
      grd%jpe = mod(grd%jge,2)
    endif
! ---
! Define surrounding processors
      mpi%top = mpi%myrank + mpi%irm
      mpi%bot = mpi%myrank - mpi%irm
      mpi%lft = mpi%myrank - 1
      mpi%rgh = mpi%myrank + 1
    if (mpi%jr.eq.mpi%jrm) mpi%top = MPI_PROC_NULL
    if (mpi%jr.eq.1      ) mpi%bot = MPI_PROC_NULL
    if (mpi%ir.eq.mpi%irm) mpi%rgh = MPI_PROC_NULL
    if (mpi%ir.eq.1      ) mpi%lft = MPI_PROC_NULL

! ---
! Define sub-tiles for the recursive filter
   ALLOCATE ( mpi%thi(2), mpi%thj(2) )
!   ALLOCATE ( mpi%psi(mpi%jrm+1), mpi%psj(mpi%irm+1) )


    endif

   mpi%thi(1) = min( grd%im         , mpi%irm * mpi%thx ) 
   mpi%thi(2) = min( grd%im * grd%km, mpi%irm * mpi%thx ) 
   mpi%thj(1) = min( grd%jm         , mpi%jrm * mpi%thy ) 
   mpi%thj(2) = min( grd%jm * grd%km, mpi%jrm * mpi%thy ) 

   ALLOCATE ( grd%irs(mpi%thi(2),2), grd%ire(mpi%thi(2),2), grd%imr(mpi%thi(2),2) )
   ALLOCATE ( grd%jrs(mpi%thj(2),2), grd%jre(mpi%thj(2),2), grd%jmr(mpi%thj(2),2) )

   do lr=1,mpi%thi(1)
      call para_reg( grd%im       , mpi%thi(1), lr, grd%irs(lr,1), grd%ire(lr,1), grd%imr(lr,1))
   enddo
   do lr=1,mpi%thi(2)
      call para_reg( grd%im*grd%km, mpi%thi(2), lr, grd%irs(lr,2), grd%ire(lr,2), grd%imr(lr,2))
   enddo
   do lr=1,mpi%thj(1)
      call para_reg( grd%jm       , mpi%thj(1), lr, grd%jrs(lr,1), grd%jre(lr,1), grd%jmr(lr,1))
   enddo
   do lr=1,mpi%thj(2)
      call para_reg( grd%jm*grd%km, mpi%thj(2), lr, grd%jrs(lr,2), grd%jre(lr,2), grd%jmr(lr,2))
   enddo

! ---
! Get position of tiles to write the output

   ALLOCATE ( grd%aj1(mpi%nproc) )
   ALLOCATE ( grd%ajm(mpi%nproc) )
   ALLOCATE ( grd%ai1(mpi%nproc) )
   ALLOCATE ( grd%aim(mpi%nproc) )

   ALLOCATE ( ibff (4) )
   ALLOCATE ( ibffa(4,mpi%nproc) )

   ibff(1) = grd%jgs
   ibff(2) = grd%jge
   ibff(3) = grd%igs
   ibff(4) = grd%ige

   call mpi_gather(ibff, 4, mpi%i4, ibffa, 4, mpi%i4, 0, mpi%comm, ierr)

  if(mpi%myrank.eq.0)then

     grd%npsm = 0
   do kproc = 1,mpi%nproc
     grd%aj1(kproc) = ibffa(1,kproc)
     grd%ajm(kproc) = ibffa(2,kproc)
     grd%ai1(kproc) = ibffa(3,kproc)
     grd%aim(kproc) = ibffa(4,kproc)
     grd%npsm = max( grd%npsm, (grd%ajm(kproc)-grd%aj1(kproc)+1)*(grd%aim(kproc)-grd%ai1(kproc)+1))
   enddo

  endif

  call mpi_bcast( grd%npsm, 1, mpi%i4, 0, mpi%comm, ierr)

   DEALLOCATE ( ibff )
   DEALLOCATE ( ibffa )

end subroutine def_mpid
!-----------------------------------------------------------------------------------
subroutine para_reg( img, irm, ir, igs, ige, im)

 use set_knd

 implicit none

  INTEGER(i4)    ::  img
  INTEGER        ::  irm, ir
  INTEGER(i4)    ::  igs, ige, im
  INTEGER(i4)    ::  iw1, iw2

      iw1     = img / irm
      iw2     = mod(img,irm)
      igs = (ir - 1) * iw1 + 1 + min((ir-1),iw2)
      ige = igs + iw1 - 1
      if(iw2.gt.(ir-1)) ige = ige + 1
      im  = ige - igs + 1

end subroutine para_reg
