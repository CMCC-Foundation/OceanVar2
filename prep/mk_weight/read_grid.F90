subroutine read_grid

use set_knd
use grd_str
use netcdf

implicit none

integer(i4) :: stat,ncid,idvar
integer(i4) :: start(3),count(3)


stat = nf90_open('grid1.nc', NF90_NOWRITE, ncid)
!stat = nf90_open('GRID_orca.nc', NF90_NOWRITE, ncid)

stat = nf90_inq_dimid (ncid, 'im', idvar)
stat = nf90_inquire_dimension (ncid, idvar, len = grd%im)
stat = nf90_inq_dimid (ncid, 'jm', idvar)
stat = nf90_inquire_dimension (ncid, idvar, len = grd%jm)
stat = nf90_inq_dimid (ncid, 'km', idvar)
stat = nf90_inquire_dimension (ncid, idvar, len = grd%km)

ALLOCATE(grd%lon(grd%im,grd%jm),grd%lat(grd%im,grd%jm),&
         grd%dx(grd%im,grd%jm),grd%dy(grd%im,grd%jm),  &
         grd%distc3d(grd%im,grd%jm,grd%km) )
ALLOCATE(grd%msk(grd%im,grd%jm,grd%km))

start(1) = 1
start(2) = 1
start(3) = 1
count(1) = grd%im
count(2) = grd%jm
count(3) = grd%km

stat = nf90_inq_varid (ncid, 'lon', idvar)
stat = nf90_get_var (ncid,idvar,grd%lon,start(1:2),count(1:2))
stat = nf90_inq_varid (ncid, 'lat', idvar)
stat = nf90_get_var (ncid,idvar,grd%lat,start(1:2),count(1:2))
stat = nf90_inq_varid (ncid, 'dx', idvar)
stat = nf90_get_var (ncid,idvar,grd%dx,start(1:2),count(1:2))
stat = nf90_inq_varid (ncid, 'dy', idvar)
stat = nf90_get_var (ncid,idvar,grd%dy,start(1:2),count(1:2))
stat = nf90_inq_varid (ncid, 'tmsk', idvar)
stat = nf90_get_var (ncid,idvar,grd%msk,start(1:3),count(1:3))
stat = nf90_inq_varid (ncid, 'dcoast', idvar)
stat = nf90_get_var (ncid, idvar, grd%distc3d,start(1:3),count(1:3))
stat = nf90_close(ncid)

end subroutine read_grid
