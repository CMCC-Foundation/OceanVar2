MODULE grd_str

 use set_knd

implicit none

public

   TYPE grid_t

        REAL(r8),    POINTER     ::  lon(:,:)       ! Longitude
        REAL(r8),    POINTER     ::  lat(:,:)       ! Latitude
        REAL(r8),    POINTER     ::  dx(:,:)        ! Longitude
        REAL(r8),    POINTER     ::  dy(:,:)        ! Latitude
        REAL(r8),    POINTER     ::  msk(:,:,:)     ! tmask
        REAL(r8),    POINTER     ::  distc3d(:,:,:) ! distance from coast
        INTEGER(i4)              ::  im             ! No. points in x direction
        INTEGER(i4)              ::  jm             ! No. points in y direction
        INTEGER(i4)              ::  km             ! No. points in z direction


   END TYPE grid_t

   TYPE (grid_t)                 :: grd

   TYPE diflt_t

       INTEGER(i4)               :: nt      ! number of iteration ( must be the same of the one used for computing weigths)
       LOGICAL                   :: rd_corr ! Read correlation length from file
       REAL(r8)                  :: rx      ! correlation radius in x  ( if not read from file )  
       REAL(r8),POINTER          :: wgh(:,:,:) ! Numerical weight for diffusion filter
       REAL(r8)                  :: ry      ! correlation radius in y  ( if not read from file )
       LOGICAL                   :: use_bc  ! use boundary condition
       CHARACTER(len=:),ALLOCATABLE::  bc_type   ! type of boundary condition (if use_bc = .T. then choose NEUMANN/DIRICHLET)
       LOGICAL                   :: use_cst ! use coastal distance to reduce the correlation radius
       REAL(r8)                  :: cst_dst ! costal distance value ( if use_cst = .T. )
   END TYPE diflt_t
   TYPE (diflt_t)                :: diflt

END MODULE grd_str

