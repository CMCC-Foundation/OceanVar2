       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!Subroutine def_Coef!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine def_coef(sigma,b1,b2,b3,b)

       implicit none

       real*8  sigma,b0,b1,b2,b3,b,q

       if (sigma .ge.  2.5) then
        q = 0.98711 * sigma - 0.96330
       else
        q = 3.97156 - 4.14554 * sqrt(1 - 0.26891 * sigma);
       endif



     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     !!Compute Recursive Gaussian Filter Coefficient!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       b0 = 1.57825 + (2.44413 * q) + (1.4281 * q**2) + (0.422205 * q**3)
       b1 =( (2.44413 * q)  + (2.85619 * q**2) + (1.26661 * q**3))/b0
       b2 =( -(1.4281 * q**2) - (1.26661 * q**3))/b0
       b3 =( 0.422205 * q**3)/b0
       b  = 1 -  ( b1 + b2 + b3 );







      end subroutine def_coef




