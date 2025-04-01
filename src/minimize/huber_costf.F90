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
!> Compute the cost function with Huber norm distribution               
!!
!!  
!!
!                                                                      !
! Version 1: Andrea Storto 2021                                        !
!            Mario Adani   2023                                        !
!                                                                      !
! ----------------------------------------------------------------------
SUBROUTINE huber_costf(jo_tot)

   USE set_knd
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)          :: jstart, jend
   REAL(r8)             :: jo_sla, jo_arg, jo_xbt, jo_gld, &
                           jo_tra, jo_trd, jo_vdr, jo_gvl, &
                           jo_sst
   ! Function
   REAL(r8)             :: huberf

   REAL(r8),INTENT(OUT) :: jo_tot

   !Initialize
   jo_sla = 0._r8
   jo_arg = 0._r8
   jo_xbt = 0._r8
   jo_gld = 0._r8
   jo_tra = 0._r8
   jo_trd = 0._r8
   jo_vdr = 0._r8
   jo_gvl = 0._r8
   jo_sst = 0._r8
   jo_tot = 0._r8

! SLA observations
   jstart = 1
   jend   = sla%nc
   IF ( huberqc%sla ) THEN
      jo_sla =  huberf(jstart,jend)
   ELSE
      jo_sla = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! ARGO observations
   jstart = sla%nc + 1
   jend   = sla%nc + arg%nc
   IF ( huberqc%arg ) THEN
      jo_arg = huberf(jstart,jend)
   ELSE
      jo_arg = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! XBT observations
   jstart = sla%nc + arg%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc
   IF ( huberqc%xbt ) THEN
      jo_xbt = huberf(jstart,jend)
   ELSE
      jo_gld = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! Glider observations
   jstart = sla%nc + arg%nc + xbt%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc
   IF ( huberqc%gld ) THEN
      jo_gld = huberf(jstart,jend)
   ELSE
      jo_gld = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! Argo trajectory observations
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc
   IF ( huberqc%tra ) THEN
      jo_tra = huberf(jstart,jend)
   ELSE
      jo_tra = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! Trajectory observations of surface drIFters
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc
   IF ( huberqc%trd ) THEN
      jo_trd = huberf(jstart,jend)
   ELSE
      jo_trd = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! Observations of drIFter velocities
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc
   IF ( huberqc%vdr ) THEN
      jo_vdr = huberf(jstart,jend)
   ELSE
      jo_vdr = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! Observations of glider velocities
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc
   IF ( huberqc%gvl ) THEN
      jo_gvl = huberf(jstart,jend)
   ELSE
      jo_gvl = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

! SST observations
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + gvl%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc + sst%nc
   IF ( huberqc%sst ) THEN
      jo_sst = huberf(jstart,jend)
   ELSE
      jo_sst = 0.5_r8 * DOT_PRODUCT( obs%amo(jstart:jend), obs%amo(jstart:jend) )
   ENDIF

   jo_tot = jo_sla + jo_arg + jo_xbt + jo_gld + jo_tra + &
            jo_trd + jo_vdr + jo_gvl + jo_sst


END SUBROUTINE huber_costf
!======================================================================
REAL (KIND=r8) FUNCTION huberf(jstart,jend)

!-----------------------------------------------------------------------
!                                                                      !
!>  Huber  Function                                                     !
!!
!!
!!
!                                                                      !
! Version 1: Andrea Storto 2021                                        !
!                                                                      !
! ----------------------------------------------------------------------

   USE set_knd
   USE obs_str, ONLY : obs

   IMPLICIT NONE

   INTEGER(i4)  :: jo, jstart, jend
   REAL(r8)     :: zcst, zj, zk, jo_cost

   jo_cost = 0._r8

   DO jo = jstart,jend

      IF ( obs%amo(jo) .GT. ABS(obs%ahub(jo,1)) .AND.                              &
           obs%amo(jo) .LE. ABS(obs%ahub2(jo,1)) ) THEN

         zcst = obs%amo(jo)*obs%ahub(jo,1)-0.5_r8*obs%ahub(jo,1)*obs%ahub(jo,1)

      ELSEIF ( obs%amo(jo) .LT. -ABS(obs%ahub(jo,2)) .AND.                         &
               obs%amo(jo) .GE. -ABS(obs%ahub2(jo,2)) ) THEN

         zcst = -obs%amo(jo)*obs%ahub(jo,2)-0.5_r8*obs%ahub(jo,2)*obs%ahub(jo,2)

      ELSEIF ( obs%amo(jo) .GT. ABS(obs%ahub2(jo,1)) ) THEN

         zj   = 2._r8*obs%ahub(jo,1)*DSQRT(obs%ahub2(jo,1))
         zk   = 0.5_r8*obs%ahub(jo,1)*obs%ahub(jo,1)-obs%ahub(jo,1)*obs%ahub2(jo,1)+zj*DSQRT(obs%ahub2(jo,1))
         zcst = zj*DSQRT(obs%amo(jo))-zk

      ELSEIF ( obs%amo(jo) .LT. -ABS(obs%ahub2(jo,2)) ) THEN

         zj   = -2._r8*obs%ahub(jo,2)*DSQRT(obs%ahub2(jo,2))
         zk   = 0.5_r8*obs%ahub(jo,2)*obs%ahub(jo,2)-obs%ahub(jo,2)*obs%ahub2(jo,2)-zj*DSQRT(obs%ahub2(jo,2))
         zcst = -zj*DSQRT(ABS(obs%amo(jo)))-zk

      ELSE

         zcst = 0.5_r8 * obs%amo(jo) * obs%amo(jo)

      ENDIF

      jo_cost = jo_cost + zcst

   ENDDO

   huberf = jo_cost

END FUNCTION huberf
