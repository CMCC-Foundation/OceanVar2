subroutine writenc(coeff,filename)

use set_knd
use grd_str
use netcdf

implicit none

INTEGER(i4)      :: stat,ncid,idlon,idlat,idtim,idvar,idz
REAL(r8)         :: coeff(grd%im,grd%jm,grd%km)
character(len=*) :: filename

stat = nf90_create      ( filename, NF90_NETCDF4, ncid                               )
!Dimensions
stat = nf90_def_dim     ( ncid, "im",   grd%im, idlon                                )
stat = nf90_def_dim     ( ncid, "jm",   grd%jm, idlat                                )
stat = nf90_def_dim     ( ncid, "km",   grd%km, idz                                  )
!stat = nf90_def_dim     ( ncid, "time", nf90_unlimited, idtim                        )
!Variables
stat = nf90_def_var     ( ncid, "lon", nf90_float, (/ idlon, idlat /), idvar         )
stat = nf90_put_var     ( ncid, idvar,  grd%lon                                      )
stat = nf90_put_att     ( ncid, idvar,  "coordinate", "lon lat"                      )
stat = nf90_def_var     ( ncid, "lat",  nf90_float, (/ idlon, idlat /), idvar        )
stat = nf90_put_var     ( ncid, idvar,  grd%lat                                      )
stat = nf90_put_att     ( ncid, idvar,  "coordinate", "lon lat"                      )
stat = nf90_def_var     ( ncid, "dy",  nf90_float, (/ idlon, idlat /), idvar        )
stat = nf90_put_var     ( ncid, idvar,  grd%dy                                      )
stat = nf90_put_att     ( ncid, idvar,  "girspacing y", "lon lat"                   )
stat = nf90_def_var     ( ncid, "dx",  nf90_float, (/ idlon, idlat /), idvar        )
stat = nf90_put_var     ( ncid, idvar,  grd%dx                                      )
stat = nf90_put_att     ( ncid, idvar,  "girspacing x", "lon lat"                   )
stat = nf90_def_var     ( ncid, "coef", nf90_double, (/ idlon, idlat, idz /), idvar,deflate_level = 4 )
stat = nf90_enddef      ( ncid )
stat = nf90_put_var ( ncid, idvar, coeff,                      &
                      start = (/    1,      1,         1  /),    &
                      count = (/ grd%im, grd%jm,    grd%km/) )
stat = nf90_close (ncid)
end subroutine writenc

