program mk_weights

use set_knd
use grd_str
use netcdf

implicit none

! Constants
REAL(r8) :: pi = 4._r8*datan(1._r8)
REAL(r8) :: r_earth = 6371._r8
! Indices
INTEGER(i4) :: it, i, j, k
INTEGER(i4) :: nx,ny,nz

! Local variables
REAL(r8)         :: r_earth2 
REAL(r8),POINTER :: dlon(:,:),dlat(:,:)
REAL(r8),POINTER :: lon_rad(:,:),lat_rad(:,:)
REAL(r8),POINTER :: rx(:,:,:), ry(:,:,:)
REAL(r8),POINTER :: kx(:,:,:), ky(:,:,:)
REAL(r8),POINTER :: gdx(:,:), gdy(:,:)
REAL(r8),POINTER :: smallF(:,:)      
REAL(r8),POINTER :: smallM(:,:)     
REAL(r8),POINTER :: smallKX(:,:)     
REAL(r8),POINTER :: smallKY(:,:)     
REAL(r8),POINTER :: smallDX(:,:)     
REAL(r8),POINTER :: smallDY(:,:) 
REAL(r8),POINTER :: smallL(:,:) 
REAL(r8),POINTER :: msk_eps(:,:,:)
REAL(r8),POINTER :: msk(:,:,:)
REAL(r8),POINTER :: cf(:,:,:)
REAL(r8),ALLOCATABLE :: ax(:),bx(:),cx(:),dx(:),ex(:),fx(:)
REAL(r8),ALLOCATABLE :: ay(:),by(:),cy(:),dy(:),ey(:),fy(:)
REAL(r8),ALLOCATABLE :: dkx(:),dky(:)
REAL(r8)         :: mu, eps
INTEGER(i4)      :: npx,npy,npx2,npy2
INTEGER(i4)      :: minx,maxx,miny,maxy,midx,midy
INTEGER(i4)      :: xstartS,xendS,ystartS,yendS
INTEGER(i4)      :: xstartL,xendL,ystartL,yendL


!-----------------------------------------

!- GET NUMBER OF PERTURBATIONS
!

!- READ NAMELIST
!
call read_nml

!- CHECK STABILITY 
if ( diflt%nt .le. 2 ) then
   print*,'NT can must be greater then 2'
   stop
endif

!- SQUARE EARTH RADIUS 
!
r_earth  = r_earth*1000._r8 !km -> m
r_earth2 = r_earth**2

!- READ GRID
!
call read_grid

! DEGREE TO RADIANS
!
allocate( lon_rad(grd%im,grd%jm), &
          lat_rad(grd%im,grd%jm) )
lon_rad=grd%lon*pi/180._r8!
lat_rad=grd%lat*pi/180._r8!

!- COMPUTE LON,LAT DIFFERENCE
allocate( dlon(grd%im,grd%jm), &
          dlat(grd%im,grd%jm) )
do j = 1,grd%jm
do i = 1,grd%im-1
   dlon(i,j) = abs(lon_rad(i,j)-lon_rad(i+1,j  ))
enddo
enddo
! Boundary condition in the last grid point
dlon(grd%im,:) = dlon(grd%im-1,:) 
do j = 1,grd%jm-1
do i = 1,grd%im
   dlat(i,j) = abs(lat_rad(i,j)-lat_rad(i  ,j+1))
enddo
enddo
! Boundary condition in the last grid point
dlat(:,grd%jm) = dlat(:,grd%jm-1) 
where( dlat == 0 ) dlat=minval(dlat,mask=(dlat .ne. 0))

!- DX,DY
allocate(gdx(grd%im,grd%jm), &
         gdy(grd%im,grd%jm))
gdx(:,:) = r_earth*cos(lat_rad)*dlat
gdy(:,:) = r_earth*dlon

!-CORRELATION RADIUS
!
allocate(rx (grd%im,grd%jm,grd%km), &
         ry (grd%im,grd%jm,grd%km) )
if ( diflt%rd_corr ) then
   call read_corrad
else
   rx(:,:,:) = diflt%rx ! meters
   ry(:,:,:) = diflt%ry ! meters
endif

!- DECREASE CORRELATION RADIUS CLOSE TO THE COAST
! 
if ( diflt%use_cst ) then
    do i = 1, grd%km
       where (grd%distc3d(:,:,i) < diflt%cst_dst )
         rx(:,:,i)  = max( (rx(:,:,i)*grd%distc3d(:,:,i))/diflt%cst_dst, gdx )
         ry(:,:,i)  = max( (ry(:,:,i)*grd%distc3d(:,:,i))/diflt%cst_dst, gdy )
       end where
    enddo
endif

!-DIFFUSION COEFFICIENT
!
allocate(kx(grd%im,grd%jm,grd%km), &
         ky(grd%im,grd%jm,grd%km) )
kx(:,:,:) = rx(:,:,:)**2/dble(2*(diflt%nt-2)) ! meter**2/s
ky(:,:,:) = ry(:,:,:)**2/dble(2*(diflt%nt-2)) ! meter**2/s

!- DEFINE BOUNDARY CONDITION
!
if ( diflt%use_bc ) then
 allocate(msk_eps(grd%im,grd%jm,grd%km))
 if (diflt%bc_type == 'DIRICHLET') mu = 1
 if (diflt%bc_type == 'NEUMANN')   mu = 0
 eps = (2*mu)/(mu+2*(1-mu))
 msk_eps(:,:,:) = 1._r8
 where ( grd%msk(:,:,:) .eq. 0 ) msk_eps(:,:,:) = eps
 kx(:,:,:) = kx(:,:,:) * msk_eps(:,:,:)
 ky(:,:,:) = ky(:,:,:) * msk_eps(:,:,:)
 deallocate(msk_eps)
endif

!COMPUTE ANALYTICAL SOLUTION
!
allocate(cf(grd%im,grd%jm,grd%km))
do k = 1, grd%km
   cf(:,:,k) = ( 4._r8 * pi * (diflt%nt-1) * sqrt(kx(:,:,k) * ky(:,:,k)) ) / &
               (gdx(:,:)*gdy(:,:) ) *grd%msk(:,:,k)
enddo

! CREATE MASK WHERE NUMERICAL SOLUTION IS NEEDED
!
allocate(msk(grd%im,grd%jm,grd%km))
msk = 0._r8
where ( grd%msk(:,:,:) .eq.  1._r8     .and. &
       (grd%distc3d(:,:,:) .lt. 2._r8*rx(:,:,:) .or.  &
        grd%distc3d(:,:,:) .lt. 2._r8*ry(:,:,:)) ) msk(:,:,:) = 1

!$OMP PARALLEL
!$OMP PARALLEL DO                       &
!$OMP COLLAPSE(3)                       &
!$OMP DEFAULT(NONE)                     &
!$OMP PRIVATE(i,j,it)                   &
!$OMP PRIVATE(dkx,ax,bx,cx,dx,ex,fx)   &
!$OMP PRIVATE(dky,ay,by,cy,dy,ey,fy)   &
!$OMP PRIVATE(npx,npy,npx2,npy2)       &
!$OMP PRIVATE(smallF,smallM)           &
!$OMP PRIVATE(smallKX,smallKY)         &
!$OMP PRIVATE(smallDX,smallDY,smallL)  &
!$OMP PRIVATE(midx,midy,maxx,minx)     &
!$OMP PRIVATE(maxy,miny)               &
!$OMP PRIVATE(xstartL,ystartL)         &
!$OMP PRIVATE(xstartS,ystartS)         &
!$OMP PRIVATE(xendL,yendL)             &
!$OMP PRIVATE(xendS,yendS)             &
!$OMP SHARED(grd,cf,diflt,kx,ky,msk)   &
!$OMP SHARED(dlon,dlat,lat_rad)        &
!$OMP SHARED(gdx,gdy,rx,ry,r_earth2)    
do nz = 1,grd%km
do ny = 1,grd%jm
do nx = 1,grd%im
   if (MSK(nx,ny,nz) .eq. 1) then
      npx = int(ceiling(2*rx(nx,ny,nz)/gdx(nx,ny)))
      npy = int(ceiling(2*ry(nx,ny,nz)/gdy(nx,ny)))
      npx2=npx*2+1
      npy2=npy*2+1
      ALLOCATE( smallF (npx2,npy2) &
               ,smallM (npx2,npy2) &
               ,smallKX(npx2,npy2) &
               ,smallKY(npx2,npy2) &
               ,smallDX(npx2,npy2) &
               ,smallDY(npx2,npy2) &
               ,smallL (npx2,npy2) &
               ,dkx(npx2)          &
               ,ax (npx2)          &
               ,bx (npx2)          &
               ,cx (npx2)          &
               ,dx (npx2)          &
               ,ex (npx2)          &
               ,fx (npx2)          &
               ,dky(npy2)          &
               ,ay (npy2)          &
               ,by (npy2)          &
               ,cy (npy2)          &
               ,dy (npy2)          &
               ,ey (npy2)          &
               ,fy (npy2)          &
              )
      midx = npx+1
      midy = npy+1
      minx = nx-npx
      maxx = nx+npx
      miny = ny-npy
      maxy = ny+npy
      smallF(:,:)     = 0._r8
      smallF(midx,midy) = 1._r8
      smallM  = 0._r8
      smallKX = 0.!-999._r8
      smallKY = 0.!-999._r8
      smalldx = -999._r8
      smallDY = -999._r8
      xstartS = 1
      xendS   = npx2
      ystartS = 1
      yendS   = npy2

      xstartL = max(1,minx)
      xendL   = min(maxx,grd%im)
      ystartL = max(1,miny)
      yendL   = min(maxy,grd%jm)

      xstartL = max(1,minx)
      xendL   = min(maxx,grd%im)
      ystartL = max(1,miny)
      yendL   = min(maxy,grd%jm)

      if (minx .le. 0) then
         xstartS = abs(minx)+2
         xendL   = npx2+minx-1
      endif
      if (miny .le. 0) then
         ystartS = abs(miny)+2
         yendL   = npy2+miny-1
      endif
      if (maxx .gt. grd%im) then
         xendS = grd%im-minx+1
      endif
      if (maxy .gt. grd%jm) then
         yendS = grd%jm-miny+1
      endif

      smallM( xstartS:xendS,ystartS:yendS) = grd%msk(xstartL:xendL,ystartL:yendL,nz)
      smallKX(xstartS:xendS,ystartS:yendS) = kx     (xstartL:xendL,ystartL:yendL,nz) 
      smallKY(xstartS:xendS,ystartS:yendS) = ky     (xstartL:xendL,ystartL:yendL,nz) 
      smallDX(xstartS:xendS,ystartS:yendS) = dlon   (xstartL:xendL,ystartL:yendL)
      smallDY(xstartS:xendS,ystartS:yendS) = dlat   (xstartL:xendL,ystartL:yendL)
      smallL( xstartS:xendS,ystartS:yendS) = lat_rad(xstartL:xendL,ystartL:yendL)
      do it = 1,diflt%nt
         !- LONGITUDE
         do j = 1,npy2
            dkx(:)    =   smallKX(:,j)/(r_earth2*cos(smallL(:,j))**2)
            call cmp_coefficients(npx2,dkx,              &
                             smallDX(:,j),ax,bx,cx)
           call tridiagLU(npx2,ax,bx,cx,ex,fx)
           dx = smallF(:,j)
           call tridiagLUSolve(npx2,dx,ax,ex,fx,smallF(:,j))
         enddo
         !- LATITUDE
         do i = 1,npx2
            dky(:) =   smallKY(i,:)/r_earth2
            call cmp_coefficients(npy2,dky,            &
                             smallDY(i,:),ay,by,cy)
            call tridiagLU(npy2,ay,by,cy,ey,fy)
            dy = smallF(i,:)
            call tridiagLUSolve(npy2,dy,ay,ey,fy,smallF(i,:))
         enddo
      enddo ! diflt%nt
      cf(nx,ny,nz) =  1._r8/smallF(midx,midy)*grd%msk(nx,ny,nz)
      DEALLOCATE( smallF &
                 ,smallM &
                 ,smallKX&
                 ,smallKY &
                 ,smallDX &
                 ,smallDY &
                 ,smallL  &
                 ,dkx     &
                 ,ax      &
                 ,bx      &
                 ,cx      &
                 ,dx      &
                 ,ex      &
                 ,fx      &
                 ,dky     &
                 ,ay      &
                 ,by      &
                 ,cy      &
                 ,dy      &
                 ,ey      &
                 ,fy      &
                )

    endif
enddo !grd%im
enddo !grd%jm
enddo !grd%km
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP END PARALLEL 


call writenc(cf,'weights.nc')

deallocate ( rx, ry, kx, ky, cf, msk, lon_rad, lat_rad )

end program mk_weights
!----------------------------------------
subroutine cmp_coefficients(n,alpha,dx,a,b,c)

use set_knd

implicit none

INTEGER(i4),INTENT(IN) :: n
REAL(r8),   INTENT(IN) :: alpha(n)
REAL(r8),   INTENT(IN) :: dx(n)
REAL(r8),   INTENT(OUT):: a(n),b(n),c(n)

INTEGER(i4) :: i
REAL(r8)    :: coefficient(n)

coefficient(:) = alpha(:)/dx(:)**2
do i=1,n-1
!   a(i) = (-alpha(i)/dx(i)**2) ! subdiagonal a: coefficients of phi(i-1)
!   c(i) = a(i)                 ! superdiagonal c: coefficients of phi(i+1)
!   b(i) = 1 - 2*a(i)           ! diagonal b: coefficients of phi(i)
    a(i) = -coefficient(i)
    c(i) = -coefficient(i+1)
    b(i) = 1+coefficient(i+1)+coefficient(i)
enddo
!Not sure
!a(1) = a(2)
!b(1) = b(2)
!a(1) = a(2)
a(n) = a(n-1)
b(n) = b(n-1)
c(n) = c(n-1)

end subroutine cmp_coefficients
!------------------------------------------
subroutine tridiagLU(n,a,b,c,e,f)
! tridiagLU Obtain the LU factorization of a tridiagonal matrix
!
! Synopsis: [e,f] = tridiag(a,b,c)
!
! Input: a,b,c = vectors defining the tridiagonal matrix. a is the
! subdiagonal, b is the main diagonal, and c is the superdiagonal
!
! Output: e,f = vectors defining the L and U factors of the tridiagonal matrix

use set_knd
use grd_str

implicit none

INTEGER(i4),INTENT(IN) :: n
REAL(r8),INTENT(IN)    :: a(n),b(n),c(n)
REAL(r8),INTENT(OUT)   :: e(n),f(n)
INTEGER(i4)            :: i


e(:) = 0._r8
f(:) = 0._r8
e(1) = b(1)
f(1) = c(1)/b(1)

do i = 2,n
   e(i) = b(i) - a(i)*f(i-1)
   f(i) = c(i) / e(i)
enddo

end subroutine tridiagLU
!--------------------------------
subroutine tridiagLUSolve(n,d,a,e,f,v)
! tridiagLUsolve Solve (LU)*v = d where L and U are LU factors of a tridiagonal
! matrix
!
! Synopsis: v = tridiagLUsolve(d,e,f)
! v = tridiagLUsolve(d,e,f,v)
!
! Input: d = right hand side vector of the system of equatoins
! e,f = vectors defining the L and U factors of the tridiagonal matrix.
! e and f are obtained with the tridiagLU function
! v = solution vector. If v is supplied, the elements of v are over-
! written (thereby saving the memory allocation step). If v is not
! supplied, it is created. v is used as a scratch vector in the
! forward solve.
!
! Output: v = solution vector


use set_knd

implicit none

INTEGER(i4),INTENT(IN)    :: n
REAL(r8),   INTENT(IN)    :: d(n),e(n),f(n),a(n)
REAL(r8),   INTENT(INOUT) :: v(n)
INTEGER(i4)               :: i

! --- Forward substitution to solve L*w = d
v(1) = d(1) / e(1)
do i = 2,n
v(i) = ( d(i) - a(i)*v(i-1) ) / e(i)
enddo
! --- Backward substitution to solve U*v = w
do i = n-1,1,-1
v(i) = v(i) - f(i) * v(i+1)
enddo
end subroutine tridiagLUsolve

