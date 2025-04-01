subroutine read_nml

use set_knd
use grd_str

implicit none

! Logical variables
LOGICAL          :: rd_corr            ! Read correlation from file
LOGICAL          :: use_bc             ! use boundary condition
LOGICAL          :: use_cst            ! use coast to reduce correlation
LOGICAL          :: use_wgh            ! Normalize by previously compued weight
! Variables
REAL(r8)         :: rx,ry              ! Correlation radius
REAL(r8)         :: cst_dst            ! Coastal distance from which correlation radious is decreased
INTEGER(i4)      :: nt                 ! number of iteration (should be the same of the one used for the computation of the weights)
CHARACTER(100)   :: bc_type            ! Type of BC



NAMELIST /dflst/  rd_corr,                      &
                  use_bc, use_cst, use_wgh,     &
                  rx, ry, nt, bc_type, cst_dst

open(11,file='weight_dif_flt_nml',form='formatted')
read(11,dflst)

! number of iteration
diflt%nt      = nt
! read correlation from file
diflt%rd_corr = rd_corr
! ... if not use the one from namelist
diflt%rx      = rx
diflt%ry      = ry
! Use boundary conditions
diflt%use_bc  = use_bc 
! if used type of boundary conditions
ALLOCATE( character(len=len(trim(bc_type))) :: diflt%bc_type )
diflt%bc_type      = trim(bc_type)
! use coastal distance to reduce the correlation radius
diflt%use_cst = use_cst
diflt%cst_dst = cst_dst
close(11)

end subroutine read_nml
