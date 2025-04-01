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
! Version 1: Andrea Storto 2021                                        !
!            Adani  Mario  2024                                        !
!-----------------------------------------------------------------------
!==========================================================================
!>           rho_unescotl : tangent-linear of density 
!==========================================================================
REAL(KIND=r8) FUNCTION rho_unescotl(salinity, temperature,&
& salinity_tl,temperature_tl)

   USE set_knd

   IMPLICIT NONE

   ! CALLing argument declaraions:
   REAL(r8) ::  salinity, temperature
   REAL(r8) ::  salinity_tl, temperature_tl

   ! internal variable declarations:
   REAL(r8) ::  s, t, p, roots
   REAL(r8) ::  s_tl, t_tl, roots_tl
   REAL(r8) ::  a, b, c, d, e
   REAL(r8) ::  aw, bw, kw
   REAL(r8) ::  a2, b2, c2, kzero, k
   REAL(r8) ::  a2_tl, b2_tl
   REAL(r8) ::  rhow, rhozero
   REAL(r8) ::  rhow_tl

   ! initialize
   s = salinity
   s_tl = salinity_tl
   t = temperature
   t_tl = temperature_tl
   roots = DSQRT (s)

   rhow_tl = ( 5._r8*6.536332d-09 * t**4._r8 -4._r8*1.120083d-06*t**3._r8 + &
   & 3._r8*1.001685d-04*t**2._r8 - 2._r8*9.095290d-03*t + &
   & 6.793952d-02 )

   a2 = (((5.3875d-09 * t - 8.2467d-07) * t + 7.6438d-05) &
   &          * t - 4.0899d-03) * t + 8.24493d-01

   a2_tl = (4._r8*5.3875d-09 * t**3._r8 - 3._r8*8.2467d-07 * t**2._r8 + &
   &  2._r8*7.6438d-05*t - 4.0899d-03)

   b2 = (-1.6546d-06 * t + 1.0227d-04) * t - 5.72466d-03
   b2_tl = (-2._r8*1.6546d-06 * t + 1.0227d-04)

   c2 = 4.8314d-04

   rho_unescotl = 2._r8*c2*s*s_tl + b2_tl*roots*s*t_tl + &
   & 1.5_r8*b2*roots*s_tl + a2_tl*s*t_tl + a2*s_tl + &
   & rhow_tl*t_tl

END FUNCTION rho_unescotl
!==========================================================================
!>           rho_unescoad : adjoint of density       
!==========================================================================
SUBROUTINE rho_unescoad(rho_ad,salinity, temperature,&
& s_ad,t_ad)

   USE set_knd

   IMPLICIT NONE

   ! CALLing argument declaraions
   REAL(r8),INTENT(IN) ::  rho_ad
   REAL(r8),INTENT(IN) ::  salinity, temperature
   REAL(r8),INTENT(OUT) ::  s_ad, t_ad
   REAL(r8) ::  pressure

   ! internal variable declarations
   REAL(r8) ::  s, t, p, roots
   REAL(r8) ::  s_tl, t_tl, roots_tl
   REAL(r8) ::  a, b, c, d, e
   REAL(r8) ::  aw, bw, kw
   REAL(r8) ::  a2, b2, c2, kzero, k
   REAL(r8) ::  a2_tl, b2_tl
   REAL(r8) ::  rhow, rhozero
   REAL(r8) ::  rhow_tl, rhozero_tl
   logical  :: lcont

   ! initialize
   s = salinity
   t = temperature
   roots = DSQRT (s)

   rhow_tl = ( 5._r8*6.536332d-09 * t**4._r8 -4._r8*1.120083d-06*t**3._r8 + &
   & 3._r8*1.001685d-04*t**2._r8 - 2._r8*9.095290d-03*t + &
   & 6.793952d-02 )

   a2 = (((5.3875d-09 * t - 8.2467d-07) * t + 7.6438d-05) &
   &          * t - 4.0899d-03) * t + 8.24493d-01

   a2_tl = (4._r8*5.3875d-09 * t**3._r8 - 3._r8*8.2467d-07 * t**2._r8 + &
   &  2._r8*7.6438d-05*t - 4.0899d-03)

   b2 = (-1.6546d-06 * t + 1.0227d-04) * t - 5.72466d-03
   b2_tl = (-2._r8*1.6546d-06 * t + 1.0227d-04)

   c2 = 4.8314d-04


   s_ad = (2._r8*c2*s + 1.5_r8*b2*roots + a2 )*rho_ad
   t_ad = (b2_tl*roots*s + a2_tl*s + rhow_tl )*rho_ad

END SUBROUTINE rho_unescoad
!==========================================================================
!> gsw_ct_from_pt : Conservative Temperature                [deg C]
!!
!! Calculates Conservative Temperature from potential temperature of seawater
!!
!! sa      : Absolute Salinity                              [g/kg]
!! pt      : potential temperature with                     [deg C]
!!           reference pressure of 0 dbar
!!
!==========================================================================
REAL(KIND=r8) FUNCTION gsw_ct_from_pt(sa, pt)

   USE set_knd

   IMPLICIT NONE

   !  cp0  =  The "specific heat" for use                         [ J/(kg K) ]
   !          with Conservative Temperature
   REAL (r8), PARAMETER :: gsw_cp0 = 3991.86795711963_r8
   !  sfac  =  1/(40*gsw_ups)
   REAL (r8), PARAMETER :: gsw_sfac = 0.0248826675584615_r8

   REAL (r8), INTENT(IN) :: sa, pt
   REAL (r8) :: pot_enthalpy, x2, x, y

   x2 = gsw_sfac*sa
   x = DSQRT(x2)
   y = pt*0.025_r8        ! normalize for F03 and F08

   pot_enthalpy =  61.01362420681071_r8 + y*(168776.46138048015_r8 + &
      y*(-2735.2785605119625_r8 + y*(2574.2164453821433_r8 + &
      y*(-1536.6644434977543_r8 + y*(545.7340497931629_r8 + &
      (-50.91091728474331_r8 - 18.30489878927802_r8*y)*y))))) + &
      x2*(268.5520265845071_r8 + y*(-12019.028203559312_r8 + &
      y*(3734.858026725145_r8 + y*(-2046.7671145057618_r8 + &
      y*(465.28655623826234_r8 + (-0.6370820302376359_r8 - &
      10.650848542359153_r8*y)*y)))) + &
      x*(937.2099110620707_r8 + y*(588.1802812170108_r8 + &
      y*(248.39476522971285_r8 + (-3.871557904936333_r8 - &
      2.6268019854268356_r8*y)*y)) + &
      x*(-1687.914374187449_r8 + x*(246.9598888781377_r8 + &
      x*(123.59576582457964_r8 - 48.5891069025409_r8*x)) + &
      y*(936.3206544460336_r8 + &
      y*(-942.7827304544439_r8 + y*(369.4389437509002_r8 + &
      (-33.83664947895248_r8 - 9.987880382780322_r8*y)*y))))))

   gsw_ct_from_pt = pot_enthalpy/gsw_cp0

END FUNCTION gsw_ct_from_pt
!==========================================================================
!>  Thermal and Aline expantion contraction coefficients
!!
!!  Calculates thermal expansion and saline contraction coefficients
!!  of seawater from absolute salinity, conservative temperature
!!
!!  sa    =  Absolute Salinity                                        [ g/kg ]
!!  ct    =  Protential Temperature                                  [ deg C ]
!!  alpha =  Thermal expansion (=-drho/dCT)                       [ kg/m^3/K ]
!!  beta  =  Haline contraction (=drho/dSA)                  [ kg/m^3/(g/kg) ]
!
!==========================================================================
SUBROUTINE alpha_beta (sa, pt, alpha, beta)

   USE set_knd

   IMPLICIT NONE

   REAL (r8), INTENT(IN)  :: sa, pt
   REAL (r8), INTENT(OUT) :: alpha, beta

   REAL (r8) :: ss,tt,ct, gsw_ct_from_pt
!  SSO  =  Standard Ocean Reference Salinity.                      [ g/kg ]
   REAL (r8), PARAMETER :: gsw_sso = 35.16504_r8
!  UPS  =  unit conversion factor for salinities                   [ g/kg ]
   REAL (r8), PARAMETER :: gsw_ups = gsw_sso/35.0_r8

   REAL (r8), PARAMETER :: deltas = 32.0_r8
   REAL (r8), PARAMETER :: sau = 40.0_r8*gsw_ups
   REAL (r8), PARAMETER :: ctu = 40.0_r8
   REAL (r8), PARAMETER :: zu = 1d4

   REAL (r8), PARAMETER :: alp050 =  2.8648472338d-02
   REAL (r8), PARAMETER :: alp040 = -6.7823124325d-02
   REAL (r8), PARAMETER :: alp140 = -5.9955123381d-02
   REAL (r8), PARAMETER :: alp030 =  8.3587258634d-01
   REAL (r8), PARAMETER :: alp130 = -1.1301873278d+00
   REAL (r8), PARAMETER :: alp230 =  5.3494903247d-01
   REAL (r8), PARAMETER :: alp020 = -1.6234028145d+00
   REAL (r8), PARAMETER :: alp120 =  2.5061411737d+00
   REAL (r8), PARAMETER :: alp220 = -1.4770589780d+00
   REAL (r8), PARAMETER :: alp010 =  1.8523793418d+00
   REAL (r8), PARAMETER :: alp110 = -3.0734838779d+00
   REAL (r8), PARAMETER :: alp210 =  3.0136782240d+00
   REAL (r8), PARAMETER :: alp310 = -1.4543073694d+00
   REAL (r8), PARAMETER :: alp410 =  2.7320572723d-01
   REAL (r8), PARAMETER :: alp320 =  2.3782870611d-01
   REAL (r8), PARAMETER :: alp000 = -6.5047345980d-01
   REAL (r8), PARAMETER :: alp100 =  1.6337444787d+00
   REAL (r8), PARAMETER :: alp200 = -2.0484575392d+00
   REAL (r8), PARAMETER :: alp300 =  1.4268760685d+00
   REAL (r8), PARAMETER :: alp400 = -4.4447427136d-01
   REAL (r8), PARAMETER :: alp500 =  4.8463173700d-02

   REAL (r8), PARAMETER :: bet000 =  1.0785939671d+01
   REAL (r8), PARAMETER :: bet100 = -4.4465045269d+01
   REAL (r8), PARAMETER :: bet200 =  7.6072094337d+01
   REAL (r8), PARAMETER :: bet300 = -6.3964420131d+01
   REAL (r8), PARAMETER :: bet400 =  2.6898783594d+01
   REAL (r8), PARAMETER :: bet500 = -4.5234968986d+00
   REAL (r8), PARAMETER :: bet010 = -8.1303841476d-01
   REAL (r8), PARAMETER :: bet110 =  2.0388435182d+00
   REAL (r8), PARAMETER :: bet210 = -2.1302689715d+00
   REAL (r8), PARAMETER :: bet310 =  8.8477644261d-01
   REAL (r8), PARAMETER :: bet410 = -1.2058930400d-01
   REAL (r8), PARAMETER :: bet020 =  7.6476477580d-01
   REAL (r8), PARAMETER :: bet120 = -1.4997670675d+00
   REAL (r8), PARAMETER :: bet220 =  1.0856114040d+00
   REAL (r8), PARAMETER :: bet320 = -2.7192349143d-01
   REAL (r8), PARAMETER :: bet030 = -4.1572985119d-01
   REAL (r8), PARAMETER :: bet130 =  4.9004223351d-01
   REAL (r8), PARAMETER :: bet230 = -1.1835625260d-01
   REAL (r8), PARAMETER :: bet040 =  1.4061037779d-01
   REAL (r8), PARAMETER :: bet140 = -1.3310958936d-01
   REAL (r8), PARAMETER :: bet050 =  5.9673736141d-03

   ct=gsw_ct_from_pt(sa, pt)

   ss = DSQRT((sa+deltas)/sau)
   tt = ct/ctu

   alpha = ((((alp050*tt                                                  &
      + alp140*ss+alp040)*tt                                       &
      + (alp230*ss+alp130)*ss+alp030)*tt                           &
      + ((alp320*ss+alp220)*ss+alp120)*ss+alp020)*tt               &
      + (((alp410*ss+alp310)*ss+alp210)*ss+alp110)*ss+alp010)*tt   &
      + ((((alp500*ss+alp400)*ss+alp300)*ss+alp200)*ss+alp100)*ss  &
      + alp000

   beta =  ((((bet050*tt                                                  &
      + bet140*ss+bet040)*tt                                       &
      + (bet230*ss+bet130)*ss+bet030)*tt                           &
      + ((bet320*ss+bet220)*ss+bet120)*ss+bet020)*tt               &
      + (((bet410*ss+bet310)*ss+bet210)*ss+bet110)*ss+bet010)*tt   &
      + ((((bet500*ss+bet400)*ss+bet300)*ss+bet200)*ss+bet100)*ss  &
      + bet000
   beta = beta / ss

END SUBROUTINE alpha_beta
