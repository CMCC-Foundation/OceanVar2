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
!> Define analysis parameters from namelists                        
!!
!!
!!
!                                                                      !
! Version 1:   Srdjan Dobricic 2006                                    !
! Version 2:   Mario Adani and Francesco Carere   2024                 !
!-----------------------------------------------------------------------
SUBROUTINE def_nml

   USE set_knd
   USE drv_str
   USE grd_str
   USE bmd_str
   USE bal_str
   USE obs_str
   USE eof_str
   USE cns_str
   USE ctl_str
   USE dfl_str
   USE mpi_str
   USE adjck_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4), PARAMETER    :: ngrids        = 3
   INTEGER(i4), PARAMETER    :: sla_sat_maxnu = 20

   CHARACTER*128 :: flag_a
   INTEGER(i4)   :: sdat_f, shou_f
   INTEGER(i4)   :: neof, nreg, ntr
   CHARACTER*128 :: eof_flname
   INTEGER(i4)   :: ctl_m, prj
   INTEGER(i4)   :: obs_sla, obs_arg, obs_xbt, obs_gld, obs_tra, obs_trd, obs_gvl, obs_sst
   INTEGER(i4)   :: obs_vdr, bmd_ncnt, mld
   LOGICAL       :: bias_at
   LOGICAL       :: unbias
   REAL(r8)      :: sla_dep, rcf_L, ctl_tol, ctl_per, bmd_fc1, bmd_fc2, rcf_loc
   REAL(r8)      :: sla_dsm, sla_MINobspt
   REAL(r8)      :: bmd_dt, bmd_ndy, bmd_ady, bmd_alp, bmd_ovr, bmd_resem
   LOGICAL       :: bal_ck, bmd_ck, byg_ck, dfl_ck
   CHARACTER*128 :: eos_flname
   INTEGER(i4)   :: nneos(ngrids)
   LOGICAL       :: ssh_unbalanced(ngrids)
   REAL(r8)      :: exp_coef_t, exp_coef_s
   INTEGER(i4)   :: grid (ngrids)
   REAL(r8)      :: ratio(ngrids)
   INTEGER(i4)   :: filter
   INTEGER(i4)   :: mask (ngrids)
   INTEGER(i4)   :: barmd(ngrids)
   INTEGER(i4)   :: balmd(ngrids)
   INTEGER(i4)   :: divda(ngrids)
   INTEGER(i4)   :: divdi(ngrids)
   INTEGER(i4)   :: cntm(ngrids)
   LOGICAL       :: err_from_file, ts_ver_dep_err
   LOGICAL       :: tim_dep_err, sla_hor_dep_err,sla_sat_dep_err
   CHARACTER*128 :: obserr_ts_flname
   LOGICAL       :: qc_res, qc_conbgr, qc_clm, qc_eofbgr, qc_vert
   REAL(r8)      :: qc_res_sla, qc_res_tem, qc_res_sal, &
                    qc_res_vel, qc_res_sst, qc_res_dis
   REAL(r8)      :: qc_bgr_sla(2), qc_bgr_tem(2), qc_bgr_sal(2), &
                    qc_bgr_vel(2), qc_bgr_sst(2), qc_bgr_dis(2)
   REAL(r8)      :: qc_clm_limT, qc_clm_limS
   CHARACTER*128 :: qc_clm_flname
   LOGICAL       :: huberqc_arg, huberqc_gld, huberqc_gvl, huberqc_sla, &
                    huberqc_sst, huberqc_tra, huberqc_trd, huberqc_vdr, &
                    huberqc_xbt, huberqc_asymm, huberqc_L05
   CHARACTER*128 :: huberqc_flname
   INTEGER(i4)   :: huberqc_iter
   LOGICAL       :: coastrej_arg, coastrej_gld, coastrej_gvl, coastrej_sla,   &
                    coastrej_sst, coastrej_tra, coastrej_trd, coastrej_vdr,   &
                    coastrej_xbt
   REAL(r8)      :: coastrej_km_arg, coastrej_km_gld, coastrej_km_gvl, coastrej_km_sla,   &
      coastrej_km_sst, coastrej_km_tra, coastrej_km_trd, coastrej_km_vdr,   &
      coastrej_km_xbt
   LOGICAL       :: thin_arg, thin_gld, thin_gvl, thin_sla,   &
                    thin_sst, thin_tra, thin_trd, thin_vdr,   &
                    thin_xbt
   REAL(r8)      :: thin_tim, thin_spc
   INTEGER(i4)   :: sla_sat_nu
   CHARACTER*13  :: sla_sat_na(sla_sat_maxnu)
   REAL(r8)      :: sla_sat_err(sla_sat_maxnu)
   REAL(r8)      :: tem_err,sal_err,sla_err
   REAL(r8)      :: ztime_weigth
   REAL(r8)      :: sla_hor_dep_ecf
   CHARACTER*128 :: obserr_sla_flname

   CHARACTER*128 :: wgh_flname,bc_type,crl_flname
   INTEGER(i4)   :: nt
   REAL(r8)      :: rx, ry, cst_dst
   LOGICAL       :: rd_corr, USE_bc,  USE_cst, rd_wgh

   INTEGER(i4)   :: ierr
   INTEGER(i4)   :: mpi_irm, mpi_jrm, mpi_thx, mpi_thy!, mpi_flgmin
   CHARACTER*4   :: cproc
   CHARACTER*512 :: inp_dir

   INTEGER(i4)   :: i

   NAMELIST /runlst/   flag_a, sdat_f, shou_f
   NAMELIST /obslst/   obs_sla, obs_arg, obs_xbt, obs_gld, obs_tra,          &
                       obs_trd, obs_vdr, obs_gvl, obs_sst
   NAMELIST /grdlst/   ntr, prj, grid, ratio, mask, barmd,         &
                       divda, divdi, cntm, balmd, nneos, ssh_unbalanced,     &
                       exp_coef_t, exp_coef_s, eos_flname, filter
   NAMELIST /ctllst/   ctl_m, ctl_tol, ctl_per
   NAMELIST /covlst/   neof, nreg, rcf_L, rcf_loc, mld, eof_flname
   NAMELIST /slalst/   sla_dep, bias_at, unbias, sla_dsm, sla_MINobspt
   NAMELIST /bmdlst/   bmd_dt, bmd_ndy, bmd_ady, bmd_alp, bmd_fc1, bmd_fc2,  &
                       bmd_ovr, bmd_resem, bmd_ncnt
   NAMELIST /adjcklst/ bal_ck, bmd_ck, byg_ck, dfl_ck
   NAMELIST /mpilst/   mpi_irm, mpi_jrm, mpi_thx, mpi_thy!, mpi_flgmin
   NAMELIST /iolst/    inp_dir
   NAMELIST /errlst/   err_from_file, tim_dep_err, ztime_weigth,             &
                       sla_hor_dep_err, sla_hor_dep_ecf, sla_sat_dep_err,    &
                       sla_sat_nu, sla_sat_na, sla_sat_err,                  &
                       ts_ver_dep_err, sla_err, tem_err, sal_err,            &
                       obserr_ts_flname, obserr_sla_flname
   NAMELIST /qcklst/   qc_res, qc_conbgr, qc_clm, qc_eofbgr, qc_vert,           &
                       qc_res_sla, qc_res_tem, qc_res_sal,                      &
                       qc_res_vel, qc_res_sst, qc_res_dis,                      &
                       qc_bgr_sla, qc_bgr_tem, qc_bgr_sal,                      &
                       qc_bgr_vel, qc_bgr_sst, qc_bgr_dis,                      &
                       qc_clm_limT, qc_clm_limS, qc_clm_flname
   NAMELIST /hublst/   huberqc_arg, huberqc_gld, huberqc_gvl, huberqc_sla,   &
                       huberqc_sst, huberqc_tra, huberqc_trd, huberqc_vdr,   &
                       huberqc_xbt, huberqc_asymm, huberqc_L05,              &
                       huberqc_flname, huberqc_iter
   NAMELIST /crjlst/   coastrej_arg, coastrej_gld, coastrej_gvl, coastrej_sla,   &
                       coastrej_sst, coastrej_tra, coastrej_trd, coastrej_vdr,   &
                       coastrej_xbt,                                             &
                       coastrej_km_arg, coastrej_km_gld, coastrej_km_gvl, coastrej_km_sla,   &
                       coastrej_km_sst, coastrej_km_tra, coastrej_km_trd, coastrej_km_vdr,   &
                       coastrej_km_xbt
   NAMELIST /thnlst/   thin_arg, thin_gld, thin_gvl, thin_sla,   &
                       thin_sst, thin_tra, thin_trd, thin_vdr,   &
                       thin_xbt, thin_tim, thin_spc
   NAMELIST /diflst/   nt, rd_corr, rx, ry, USE_bc, bc_type, USE_cst, cst_dst, &
                       wgh_flname, crl_flname, rd_wgh

! -------------------------------------------------------------------
! Define some mpi constants
   mpi%comm = mpi_comm_world

   CALL mpi_comm_size(mpi%comm,mpi%nproc ,ierr)
   CALL mpi_comm_rank(mpi%comm,mpi%myrank,ierr)

   IF (i4.EQ.4) THEN
      mpi%i4 = mpi_INTEGER4
   ELSE
      mpi%i4 = mpi_INTEGER8
   ENDIF
   IF (i8.EQ.4) THEN
      mpi%i8 = mpi_INTEGER4
   ELSE
      mpi%i8 = mpi_INTEGER8
   ENDIF
   IF (r4.EQ.4) THEN
      mpi%r4 = MPI_REAL4
   ELSE
      mpi%r4 = MPI_REAL8
   ENDIF
   IF (r8.EQ.4) THEN
      mpi%r8 = MPI_REAL4
   ELSE
      mpi%r8 = MPI_REAL8
   ENDIF

! -------------------------------------------------------------------
! Open a formatted file for the diagnostics
   drv%dia = 12

   WRITE (cproc,'(i4)')mpi%myrank
   IF ( mpi%myrank .LT. 1000) cproc(1:1) = '0'
   IF ( mpi%myrank .LT. 100 ) cproc(2:2) = '0'
   IF ( mpi%myrank .LT. 10  ) cproc(3:3) = '0'

   OPEN ( drv%dia, FILE='OceanVar.diagnostics_'//cproc, FORM='formatted' )

!---------------------------------------------------------------------
! Open the NAMELIST
   OPEN (11,FILE='OceanVar_nml',FORM='formatted')

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) '  '
   WRITE (drv%dia,*) '                      NAMELISTS: '
   WRITE (drv%dia,*) '  '

! ---
   READ (11,runlst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' RUN NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Flag for the analysis:           flag_a   = ', flag_a
   WRITE (drv%dia,*) ' Starting day of the forecast:    sdat_f   = ', sdat_f
   WRITE (drv%dia,*) ' Starting hour of the forecast:   shou_f   = ', shou_f

   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(flag_a))) :: drv%flag )
   drv%flag = TRIM(flag_a)
   drv%sdat = sdat_f
   drv%shou = shou_f

   CALL FLUSH (drv%dia)
! ---
   READ (11,obslst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' OBSERVATIONS NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Use SLA observations:                   obs_sla = ', obs_sla
   WRITE (drv%dia,*) ' Use ARGO observations:                  obs_arg = ', obs_arg
   WRITE (drv%dia,*) ' Use XBT observations:                   obs_xbt = ', obs_xbt
   WRITE (drv%dia,*) ' Use glider observations:                obs_gld = ', obs_gld
   WRITE (drv%dia,*) ' Use Argo trajectory observations:       obs_tra = ', obs_tra
   WRITE (drv%dia,*) ' Use drifter trajectory observations:    obs_trd = ', obs_trd
   WRITE (drv%dia,*) ' Use drifter observations of velocity:   obs_vdr = ', obs_vdr
   WRITE (drv%dia,*) ' Use glider observations of velocity:    obs_gvl = ', obs_gvl
   WRITE (drv%dia,*) ' Use SST observations:                   obs_sst = ', obs_sst

   obs%sla = obs_sla
   obs%arg = obs_arg
   obs%xbt = obs_xbt
   obs%gld = obs_gld
   obs%tra = obs_tra
   obs%trd = obs_trd
   obs%vdr = obs_vdr
   obs%gvl = obs_gvl
   obs%sst = obs_sst

   CALL FLUSH(drv%dia)
! ---
   sla_sat_na(:) = '             '
   READ (11,errlst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' OBSERVTIONAL ERROR NAMELIST INPUT                          '
   WRITE (drv%dia,*) ' Read error from file:                err_from_file   =     ', err_from_file
   IF ( .NOT. err_from_file ) THEN
      WRITE (drv%dia,*) ' Temporal dependency:                 tim_dep_err     =     ', tim_dep_err
      IF ( tim_dep_err ) THEN
         WRITE (drv%dia,*) ' Temporal decay (days) for obs err:   ztime_weigth    =     ', ztime_weigth
      ENDIF
      WRITE (drv%dia,*) ' SLA satellite dependency:            sla_sat_dep_err =     ', sla_sat_dep_err
      WRITE (drv%dia,*) ' SLA horizontal dependency:           sla_hor_dep_err =     ', sla_hor_dep_err
      IF ( sla_hor_dep_err ) THEN
         WRITE (drv%dia,*) ' SLA horizontal correction factor:    sla_hor_dep_ecf   =   ', sla_hor_dep_ecf
      ENDIF
      IF ( sla_hor_dep_err ) THEN
         WRITE (drv%dia,*) ' SLA horizontal dependency file name: obserr_sla_flname =   ', TRIM(obserr_sla_flname)
      ENDIF
      IF ( sla_sat_nu .GT. sla_sat_maxnu ) THEN
         WRITE (drv%dia,*) ' Number of satellites greater than max number of satellites               '
         WRITE (drv%dia,*) ' Maximum number of satellite (hardcoded in def_nml.F90): maxnu_sla_sat  = ', sla_sat_maxnu
         WRITE (drv%dia,*) ' Number of satellite used:                               nu_sla_sat     = ', sla_sat_nu
         WRITE (drv%dia,*) ' Please increase maxnu_sla_sat in def_nml.F90                             '
         CALL abort
      ENDIF
      IF ( ( .NOT.  sla_sat_dep_err ) .AND. ( .NOT.  sla_sat_dep_err ) ) THEN
         WRITE (drv%dia,*) ' SLA constant satellite error:        sla_err         = ',sla_err
      ENDIF
      IF (  sla_sat_dep_err ) THEN
         WRITE (drv%dia,*) ' Satellite name and error associated sla_na_sat / sla_sat_err:'
         DO i = 1, sla_sat_nu
            WRITE (drv%dia,*) '    ', sla_sat_na(i),'/ ',sla_sat_err(i)
         ENDDO
      ENDIF
      WRITE (drv%dia,*) ' Insitu vertical dependency:          ts_ver_dep_err   =    ', ts_ver_dep_err
      IF  ( ts_ver_dep_err ) THEN
         WRITE (drv%dia,*) ' T/S error filename:              obserr_ts_flname =    ', TRIM(obserr_ts_flname)
      ELSE
         WRITE (drv%dia,*) ' Temperature constant   error:        tem_err         = ',tem_err
         WRITE (drv%dia,*) ' Salinity    constant   error:        sal_err         = ',sal_err
      ENDIF
   ENDIF

   obserr%rd_err_ff          = err_from_file
   obserr%tim_dep_err        = tim_dep_err
   obserr%ztime_weigth       = ztime_weigth
   obserr%sla_hor_dep_err    = sla_hor_dep_err
   obserr%sla_hor_dep_ecf    = sla_hor_dep_ecf
   IF ( sla_hor_dep_err ) THEN
      ALLOCATE ( CHARACTER(LEN=LEN(TRIM(obserr_sla_flname))) :: obserr%sla_flname )
      obserr%sla_flname = TRIM(obserr_sla_flname)
   ENDIF
   obserr%sla_sat_dep_err    = sla_sat_dep_err
   obserr%sla_sat_nu         = sla_sat_nu
   IF (  sla_sat_dep_err ) THEN
      ALLOCATE ( obserr%sla_sat_name(sla_sat_nu) )
      ALLOCATE ( obserr%sla_sat_err (sla_sat_nu) )
      DO i = 1,obserr%sla_sat_nu
         obserr%sla_sat_name(i) = sla_sat_na(i)
         obserr%sla_sat_err (i) = sla_sat_err(i)
      ENDDO
   ENDIF
   obserr%ts_ver_dep_err = ts_ver_dep_err
   IF ( ts_ver_dep_err) THEN
      ALLOCATE ( CHARACTER(LEN=LEN(TRIM(obserr_ts_flname))) :: obserr%ts_flname )
      obserr%ts_flname = obserr_ts_flname
   ENDIF

   obserr%sla_con_err = sla_err
   obserr%tem_con_err = tem_err
   obserr%sal_con_err = sal_err

   CALL FLUSH (drv%dia)
! ---
   READ (11,thnlst)
   WRITE (drv%dia,*) '--------------------------------------------------------'
   WRITE (drv%dia,*) ' THINNING NAMELIST INPUT                                '
   WRITE (drv%dia,*) ' Apply vertical thinning to ARGO:                       ', thin_arg
   WRITE (drv%dia,*) ' Apply vertical thinning to GLIDER:                     ', thin_gld
   WRITE (drv%dia,*) ' Apply vertical thinning to GLIDER  velocity:           ', thin_gvl
   WRITE (drv%dia,*) ' Apply vertical thinning to SLA:                        ', thin_sla
   WRITE (drv%dia,*) ' Apply vertical thinning to SST:                        ', thin_sst
   WRITE (drv%dia,*) ' Apply vertical thinning to TRA:                        ', thin_tra
   WRITE (drv%dia,*) ' Apply vertical thinning to TRD:                        ', thin_trd
   WRITE (drv%dia,*) ' Apply vertical thinning to DRIFTER velocity:           ', thin_vdr
   WRITE (drv%dia,*) ' Apply vertical thinning to XBT:                        ', thin_xbt

   thin%any = thin_arg  .OR. thin_gld .OR. &
              thin_gvl  .OR. thin_sla .OR. &
              thin_sst  .OR. thin_tra .OR. &
              thin_trd  .OR. thin_vdr .OR. &
              thin_xbt
   IF ( thin%any )    THEN
      WRITE (drv%dia,*) ' Time  window for thinning observation:                 ', thin_tim
      WRITE (drv%dia,*) ' Space window for thinning observation:                 ', thin_spc
   ENDIF
   thin%arg = thin_arg
   thin%gld = thin_gld
   thin%gvl = thin_gvl
   thin%sla = thin_sla
   thin%sst = thin_sst
   thin%tra = thin_tra
   thin%trd = thin_trd
   thin%vdr = thin_vdr
   thin%xbt = thin_xbt
   thin%tim = thin_tim
   thin%spc = thin_spc

   CALL FLUSH (drv%dia)

! ---
   READ (11,qcklst)

   WRITE (drv%dia,*) '--------------------------------------------------------'
   WRITE (drv%dia,*) ' QUALITY CHECK  NAMELIST INPUT                          '
   WRITE (drv%dia,*) ' Absolute misfits:                                 ', qc_res
   WRITE (drv%dia,*) ' EOF Background  :                                 ', qc_eofbgr
   WRITE (drv%dia,*) ' Background/Observation ratio:                     ', qc_conbgr
   WRITE (drv%dia,*) ' Vertical check  :                                 ', qc_vert
   WRITE (drv%dia,*) ' Climatological  :                                 ', qc_clm
   IF ( qc_res ) THEN
      WRITE (drv%dia,*) ' Max SLA          residual  :                     ', qc_res_sla
      WRITE (drv%dia,*) ' Max Temperature  residual  :                     ', qc_res_tem
      WRITE (drv%dia,*) ' Max Salinity     residual  :                     ', qc_res_sal
      WRITE (drv%dia,*) ' Max Velocity     residual  :                     ', qc_res_vel
      WRITE (drv%dia,*) ' Max SST          residual  :                     ', qc_res_sst
      WRITE (drv%dia,*) ' Max distance     residual  :                     ', qc_res_dis
   ENDIF
   IF ( qc_conbgr ) THEN
      WRITE (drv%dia,*) ' Background / Treshold for SLA            :                     ', qc_bgr_sla(1),' / ',qc_bgr_sla(2)
      WRITE (drv%dia,*) ' Background / Treshold for Temperature    :                     ', qc_bgr_tem(1),' / ',qc_bgr_tem(2)
      WRITE (drv%dia,*) ' Background / Treshold for Salinity       :                     ', qc_bgr_sal(1),' / ',qc_bgr_sal(2)
      WRITE (drv%dia,*) ' Background / Treshold for Velocity       :                     ', qc_bgr_vel(1),' / ',qc_bgr_vel(2)
      WRITE (drv%dia,*) ' Background / Treshold for SST            :                     ', qc_bgr_sst(1),' / ',qc_bgr_sst(2)
      WRITE (drv%dia,*) ' Background / Treshold for distance       :                     ', qc_bgr_dis(1),' / ',qc_bgr_dis(2)
   ENDIF
   IF ( qc_clm ) THEN
      WRITE (drv%dia,*) ' Climatological Treshold for Temperature ABS(clim-obs)  :           ', qc_clm_limT
      WRITE (drv%dia,*) ' Climatological Treshold for Salinity    ABS(clim-obs)  :           ', qc_clm_limS
   ENDIF
   IF ( ( qc_eofbgr  ) .OR. ( qc_clm ) ) THEN
      WRITE (drv%dia,*) ' File name for climatology:                                      ', TRIM(qc_clm_flname)
   ENDIF


   qck%res     = qc_res
   qck%conbgr  = qc_conbgr
   qck%eofbgr  = qc_eofbgr
   qck%vert    = qc_vert
   qck%clm     = qc_clm
   qck%res_sla = qc_res_sla
   qck%res_tem = qc_res_tem
   qck%res_sal = qc_res_sal
   qck%res_vel = qc_res_vel
   qck%res_sst = qc_res_sst
   qck%res_dis = qc_res_dis
   qck%bgr_sla = qc_bgr_sla
   qck%bgr_tem = qc_bgr_tem
   qck%bgr_sal = qc_bgr_sal
   qck%bgr_vel = qc_bgr_vel
   qck%bgr_sst = qc_bgr_sst
   qck%bgr_dis = qc_bgr_dis
   qck%clm_lim(1) = qc_clm_limT
   qck%clm_lim(2) = qc_clm_limS
   IF ( ( qck%eofbgr ) .OR. ( qck%clm ) ) THEN
      ALLOCATE ( CHARACTER(LEN=LEN(TRIM(qc_clm_flname))) :: qck%flname )
      qck%flname = TRIM(qc_clm_flname)
   ENDIF

   CALL FLUSH (drv%dia)
! ---
   READ (11,hublst)

   WRITE (drv%dia,*) '------------------------------------------------'
   WRITE (drv%dia,*) ' HUBER NORM DISTRIBUTION QUALITY CONTROL INPUT: '
   WRITE (drv%dia,*) ' Apply to ARGO:                                 ', huberqc_arg
   WRITE (drv%dia,*) ' Apply to GLD :                                 ', huberqc_gld
   WRITE (drv%dia,*) ' Apply to GLV :                                 ', huberqc_gvl
   WRITE (drv%dia,*) ' Apply to SLA :                                 ', huberqc_sla
   WRITE (drv%dia,*) ' Apply to SST :                                 ', huberqc_sst
   WRITE (drv%dia,*) ' Apply to TRA :                                 ', huberqc_tra
   WRITE (drv%dia,*) ' Apply to TRD :                                 ', huberqc_trd
   WRITE (drv%dia,*) ' Apply to VDR :                                 ', huberqc_vdr
   WRITE (drv%dia,*) ' Apply to XBT :                                 ', huberqc_xbt

   huberqc%arg   = huberqc_arg
   huberqc%gld   = huberqc_gld
   huberqc%gvl   = huberqc_gvl
   huberqc%sla   = huberqc_sla
   huberqc%sst   = huberqc_sst
   huberqc%tra   = huberqc_tra
   huberqc%trd   = huberqc_trd
   huberqc%vdr   = huberqc_vdr
   huberqc%xbt   = huberqc_xbt
   huberqc%any   = huberqc_arg  .OR. huberqc_gld .OR. &
      huberqc_gvl  .OR. huberqc_sla .OR. &
      huberqc_sst  .OR. huberqc_tra .OR. &
      huberqc_trd  .OR. huberqc_vdr .OR. &
      huberqc_xbt
   huberqc%asymm = huberqc_asymm
   huberqc%L05   = huberqc_L05
   huberqc%iter  = huberqc_iter
   IF ( huberqc%any ) THEN
      WRITE (drv%dia,*) ' Coefficent file name:                          ', TRIM(huberqc_flname)
      WRITE (drv%dia,*) ' Non-symmetric coeffients:                      ', huberqc_asymm
      WRITE (drv%dia,*) ' L2-L1 + L05 Norm                               ', huberqc_L05
      WRITE (drv%dia,*) ' Iteration at which the variational QC start:   ', huberqc_iter
      ALLOCATE ( CHARACTER(LEN=LEN(TRIM(huberqc_flname))) :: huberqc%flname )
      huberqc%flname = TRIM(huberqc_flname)
   ENDIF

   CALL FLUSH (drv%dia)
! ---
   READ (11,crjlst)

   WRITE (drv%dia,*) '------------------------------------------------'
   WRITE (drv%dia,*) ' COASTAL REJECTION INPUT: '
   WRITE (drv%dia,*) ' Apply to ARGO:                                 ', coastrej_arg
   WRITE (drv%dia,*) ' Apply to GLD :                                 ', coastrej_gld
   WRITE (drv%dia,*) ' Apply to GLV :                                 ', coastrej_gvl
   WRITE (drv%dia,*) ' Apply to SLA :                                 ', coastrej_sla
   WRITE (drv%dia,*) ' Apply to SST :                                 ', coastrej_sst
   WRITE (drv%dia,*) ' Apply to TRA :                                 ', coastrej_tra
   WRITE (drv%dia,*) ' Apply to TRD :                                 ', coastrej_trd
   WRITE (drv%dia,*) ' Apply to VDR :                                 ', coastrej_vdr
   WRITE (drv%dia,*) ' Apply to XBT :                                 ', coastrej_xbt
   IF ( coastrej_arg ) THEN
      WRITE (drv%dia,*) 'Km from coast for ARGO',coastrej_km_arg
   ENDIF
   IF ( coastrej_gld ) THEN
      WRITE (drv%dia,*) 'Km from coast for GLD ',coastrej_km_gld
   ENDIF
   IF ( coastrej_gvl ) THEN
      WRITE (drv%dia,*) 'Km from coast for GVL ',coastrej_km_gvl
   ENDIF
   IF ( coastrej_sla ) THEN
      WRITE (drv%dia,*) 'Km from coast for SLA ',coastrej_km_sla
   ENDIF
   IF ( coastrej_sst ) THEN
      WRITE (drv%dia,*) 'Km from coast for SST ',coastrej_km_sst
   ENDIF
   IF ( coastrej_tra ) THEN
      WRITE (drv%dia,*) 'Km from coast for TRA ',coastrej_km_tra
   ENDIF
   IF ( coastrej_trd ) THEN
      WRITE (drv%dia,*) 'Km from coast for TRD ',coastrej_km_trd
   ENDIF
   IF ( coastrej_vdr ) THEN
      WRITE (drv%dia,*) 'Km from coast for VDR ',coastrej_km_vdr
   ENDIF
   IF ( coastrej_xbt ) THEN
      WRITE (drv%dia,*) 'Km from coast for XBT ',coastrej_km_xbt
   ENDIF

   coastrej%arg   = coastrej_arg
   coastrej%gld   = coastrej_gld
   coastrej%gvl   = coastrej_gvl
   coastrej%sla   = coastrej_sla
   coastrej%sst   = coastrej_sst
   coastrej%tra   = coastrej_tra
   coastrej%trd   = coastrej_trd
   coastrej%vdr   = coastrej_vdr
   coastrej%xbt   = coastrej_xbt
   coastrej%any   = coastrej_arg  .OR. coastrej_gld .OR. &
      coastrej_gvl  .OR. coastrej_sla .OR. &
      coastrej_sst  .OR. coastrej_tra .OR. &
      coastrej_trd  .OR. coastrej_vdr .OR. &
      coastrej_xbt
   coastrej%km_arg   = coastrej_km_arg
   coastrej%km_gld   = coastrej_km_gld
   coastrej%km_gvl   = coastrej_km_gvl
   coastrej%km_sla   = coastrej_km_sla
   coastrej%km_sst   = coastrej_km_sst
   coastrej%km_tra   = coastrej_km_tra
   coastrej%km_trd   = coastrej_km_trd
   coastrej%km_vdr   = coastrej_km_vdr
   coastrej%km_xbt   = coastrej_km_xbt

   CALL FLUSH (drv%dia)
! ---
   READ (11,grdlst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' GRID NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Multigrid iterrations:                  ntr    = ', ntr
   WRITE (drv%dia,*) ' Projection:                             prj    = ', prj
   WRITE (drv%dia,*) ' Filter 1-Recursive,2-DIFfusive,3-None: filter  = ', filter
   WRITE (drv%dia,*) ' Grids:                                 grid    = ', grid (1:ntr)
   WRITE (drv%dia,*) ' Ratio:                                ratio    = ', ratio(1:ntr)
   WRITE (drv%dia,*) ' Masks:                                 mask    = ',  mask(1:ntr)
   WRITE (drv%dia,*) ' Run barotropic model:                 barmd    = ', barmd(1:ntr)
   WRITE (drv%dia,*) ' Simplified balance model (D.H.):      balmd    = ', balmd(1:ntr)
   WRITE (drv%dia,*) ' Divergence damping in analysis:       divda    = ', divda(1:ntr)
   WRITE (drv%dia,*) ' Divergence damping in initialisation: divdi    = ', divdi(1:ntr)
   WRITE (drv%dia,*) ' Maximum number of cost function calls: cntm    = ', cntm(1:ntr)
   WRITE (drv%dia,*) ' T/S filename to compute EOS:          flname   = ', TRIM(eos_flname)
   WRITE (drv%dia,*) ' Equation of state:                    nneos    = ', nneos(1:ntr)
   WRITE (drv%dia,*) ' Unbalanced component of ssh:   ssh_unbalanced  = ', ssh_unbalanced(1:ntr)
   WRITE (drv%dia,*) ' Thermal expansion   coefficient:   exp_coef_t  = ', exp_coef_t
   WRITE (drv%dia,*) ' Haline  contraction coefficient:   exp_coef_s  = ', exp_coef_s

   drv%ntr = ntr
   grd%prj = prj
   drv%filter = filter
   ALLOCATE ( drv%grid (drv%ntr))
   ALLOCATE ( drv%ratco(drv%ntr))
   ALLOCATE ( drv%ratio(drv%ntr))
   ALLOCATE ( drv%mask (drv%ntr))
   ALLOCATE ( drv%bmd(drv%ntr))
   ALLOCATE ( drv%bal(drv%ntr))
   ALLOCATE ( drv%dda(drv%ntr))
   ALLOCATE ( drv%ddi(drv%ntr))
   ALLOCATE ( drv%cntm(drv%ntr))
   ALLOCATE ( drv%nneos(drv%ntr))
   ALLOCATE ( drv%ssh_unbalanced(drv%ntr))
   drv%grid (1:drv%ntr)           = grid (1:drv%ntr)
   drv%ratco(1:drv%ntr)           = ratio(1:drv%ntr)
   drv%mask (1:drv%ntr)           = mask (1:drv%ntr)
   drv%bmd  (1:drv%ntr)           = barmd(1:drv%ntr)
   drv%bal  (1:drv%ntr)           = balmd(1:drv%ntr)
   drv%dda  (1:drv%ntr)           = divda(1:drv%ntr)
   drv%ddi  (1:drv%ntr)           = divdi(1:drv%ntr)
   drv%cntm (1:drv%ntr)           = cntm(1:drv%ntr)
   drv%nneos(1:drv%ntr)           = nneos(1:drv%ntr)
   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(eos_flname))) :: drv%eosflname )
   drv%eosflname = TRIM(eos_flname)
   drv%ssh_unbalanced(1:drv%ntr)  = ssh_unbalanced(1:drv%ntr)
   grd%alpha                      = exp_coef_t
   grd%beta                       = exp_coef_s

   drv%ratio(        1)    = 1.0
   IF (drv%ntr.GT.1) drv%ratio(2:drv%ntr)    = drv%ratco(1:drv%ntr-1) / drv%ratco(2:drv%ntr)

   CALL FLUSH (drv%dia)
! ---
   READ (11,ctllst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' MINIMIZER NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Number of saved vectors:         ctl_m    = ', ctl_m
   WRITE (drv%dia,*) ' Minimum gradient of J:           ctl_tol  = ', ctl_tol
   WRITE (drv%dia,*) ' Percentage of initial gradient:  ctl_per  = ', ctl_per

   ctl%m     = ctl_m
   ctl%pgtol = ctl_tol
   ctl%pgper = ctl_per

   CALL FLUSH (drv%dia)
! ---
   READ (11,covlst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' COVARIANCE NAMELIST INPUT: '
   WRITE (drv%dia,*) ' EOFs filename :                           = ', TRIM(eof_flname)
   WRITE (drv%dia,*) ' Number of EOFs:                  neof     = ', neof
   WRITE (drv%dia,*) ' Number of regions:               nreg     = ', nreg
   WRITE (drv%dia,*) ' Horizontal correlation radius:   rcf_L    = ', rcf_L
   WRITE (drv%dia,*) ' Horizontal localization radius:  rcf_loc  = ', rcf_loc
   WRITE (drv%dia,*) ' Mixed layer depth information:   ros_mld  = ', mld

   ros%neof     = neof
   ros%nreg     = nreg
   rcf%L        = rcf_L
   rcf%loc      = rcf_loc
   ros%mld      = mld
   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(eof_flname))) :: ros%flname )
   ros%flname = TRIM(eof_flname)


   CALL FLUSH (drv%dia)
! ---
   READ (11,diflst)
   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' DIFFUSIVE FILTER NAMELIST INPUT:                           '
   WRITE (drv%dia,*) ' Read normalization factor:                rd_wgh       =   ',rd_wgh
   WRITE (drv%dia,*) ' Normalization factor filename:            wgh_flname   =   ',TRIM(wgh_flname)
   WRITE (drv%dia,*) ' Read horizontal correlation from file:    rd_corr      =   ',rd_corr
   WRITE (drv%dia,*) ' Horizontal correlation filename:          hcorr_flname =   ',TRIM(crl_flname)
   WRITE (drv%dia,*) ' Longitudinal correlation length:          rx           =   ',rx
   WRITE (drv%dia,*) ' Latitudinal correlation length:           ry           =   ',ry
   WRITE (drv%dia,*) ' Apply boundary condition:                 USE_bc       =   ',USE_bc
   WRITE (drv%dia,*) ' Type of boundary condition:               bc_type      =   ',TRIM(bc_type)
   WRITE (drv%dia,*) ' Use coastal distance:                     USE_cst      =   ',USE_cst
   WRITE (drv%dia,*) ' Coastal distance:                         cst_dst      =   ',cst_dst
   WRITE (drv%dia,*) ' No iteration:                             nt           =   ',nt

   dfl%rd_wgh    = rd_wgh
   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(wgh_flname))) :: dfl%wgh_flname )
   dfl%wgh_flname = TRIM(wgh_flname)
   dfl%rd_corr    = rd_corr
   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(crl_flname))) :: dfl%crl_flname )
   dfl%crl_flname = TRIM(crl_flname)
   dfl%rx         = rx
   dfl%ry         = ry
   dfl%USE_bc     = USE_bc
   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(bc_type))) :: dfl%bc_type )
   dfl%bc_type    = TRIM(bc_type)
   dfl%USE_cst    = USE_cst
   dfl%cst_dst    = cst_dst
   dfl%nt         = nt


   CALL FLUSH (drv%dia)
! ---
   READ (11,slalst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' SLA NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Minimum depth for observations          :  sla_dep       = ', sla_dep
   WRITE (drv%dia,*) ' Remove bias                             :  unbias        = ', unbias
   WRITE (drv%dia,*) ' Remove bias along track                 :  bias_at       = ', bias_at
   WRITE (drv%dia,*) ' Distance to recongnize the track        :  sla_dsm       = ', sla_dsm
   WRITE (drv%dia,*) ' Minimum number of observation per track :  sla_minobspt  = ', sla_minobspt

   sla%dep      = sla_dep
   sla%dsm      = sla_dsm
   sla%unbias   = unbias
   sla%bias_at  = bias_at
   sla%minobspt = sla_minobspt

   CALL FLUSH (drv%dia)
! ---
   READ (11,bmdlst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' BAROTROPIC MODEL NAMELIST INPUT: '
   WRITE (drv%dia,*) ' Time step:           bmd_dt                        = ', bmd_dt
   WRITE (drv%dia,*) ' Simulation days:     bmd_ndy                       = ', bmd_ndy
   WRITE (drv%dia,*) ' Averaged days:       bmd_ady                       = ', bmd_ady
   WRITE (drv%dia,*) ' Implicit weight:     bmd_alp                       = ', bmd_alp
   WRITE (drv%dia,*) ' Friction intensity:  bmd_fc1                       = ', bmd_fc1
   WRITE (drv%dia,*) ' Friction intensity:  bmd_fc2                       = ', bmd_fc2
   WRITE (drv%dia,*) ' Over-relaxation:     bmd_ovr                       = ', bmd_ovr
   WRITE (drv%dia,*) ' Minimum residual:    bmd_resem                     = ', bmd_resem
   WRITE (drv%dia,*) ' Maximum iterations   bmd_ncnt                      = ', bmd_ncnt

   bmd%dt    = bmd_dt 
   bmd%ndy   = bmd_ndy 
   bmd%ady   = bmd_ady 
   bmd%alp1  = bmd_alp 
   bmd%fc1   = bmd_fc1
   bmd%fc2   = bmd_fc2
   bmd%ovr   = bmd_ovr 
   bmd%resem = bmd_resem 
   bmd%ncnt   = bmd_ncnt 

   CALL FLUSH (drv%dia)
! ---
   READ (11,adjcklst)

   WRITE (drv%dia,*) '---------------------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK NAMELIST INPUT:                     '
   WRITE (drv%dia,*) ' Check simplIFied balance operator:   bal_ck     = ', bal_ck
   WRITE (drv%dia,*) ' Check barotropic model           :   bmd_ck     = ', bmd_ck
   WRITE (drv%dia,*) ' Check buoyancy                   :   byg_ck     = ', byg_ck
   WRITE (drv%dia,*) ' Check diffusion filter           :   dfl_ck     = ', dfl_ck

   adjck%bal = bal_ck
   adjck%bmd = bmd_ck
   adjck%byg = byg_ck
   adjck%dfl = dfl_ck


   CALL FLUSH (drv%dia)
! ---
   READ (11,mpilst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' MPI DISTRIBUTION OF PROCESSORS: '
   WRITE (drv%dia,*) ' Number of tiles in x direction              = ', mpi_irm
   WRITE (drv%dia,*) ' Number of tiles in y direction              = ', mpi_jrm
   WRITE (drv%dia,*) ' Mult. factor for RF threads in x direction  = ', mpi_thx
   WRITE (drv%dia,*) ' Mult. factor for RF threads in y direction  = ', mpi_thy

   mpi%irm = mpi_irm
   mpi%jrm = mpi_jrm
   mpi%thx = mpi_thx
   mpi%thy = mpi_thy

   CALL FLUSH (drv%dia)
! ---
   READ (11,iolst)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) ' DIRECTORIES FOR DOING I/O: '
   WRITE (drv%dia,*) ' Directory for reading input                 = ', inp_dir

   ALLOCATE ( CHARACTER(LEN=LEN(TRIM(inp_dir))) :: drv%inpdir )
   drv%inpdir=TRIM(inp_dir)

   WRITE (drv%dia,*) '------------------------------------------------------------'
   WRITE (drv%dia,*) '-----------------END NAMELIST INIT -------------------------'
   WRITE (drv%dia,*) '------------------------------------------------------------'

   CALL FLUSH (drv%dia)

! Define number of outer loop iterations for assimilating SST
! ---
   drv%nts = obs%sst + 1

! Checks inconsistency NAMELIST options
! ---
   DO i = 1,ntr
      IF ( ( barmd(i) .EQ. 1 ) .AND. ( balmd(i) .EQ. 1 ) ) THEN
         WRITE (drv%dia,*)'---------------------------------------------------------------------------------'
         WRITE (drv%dia,*)' Barotropic model and simplified balance operator are both activated for grid:   '
         WRITE (drv%dia,*)' Please choose either balmd or barmd 1, or both  0                               '
         WRITE (drv%dia,*)'---------------------------------------------------------------------------------'
         CALL FLUSH (drv%dia)
         CALL abort
      ENDIF
   ENDDO

! Correct the flags for sequentially running modules
! ---
#ifdef REPRO
mpi%flg_min = 1
#else
mpi%flg_min = 0
#endif

   IF (mpi%nproc.EQ.1) THEN
      mpi%flg_min = 0
   ENDIF

END SUBROUTINE def_nml
