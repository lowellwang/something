SUBROUTINE TDDFT(xatom, fatom, workr_n, Etot, iforce_cal, ido_rho, &
                 ido_vr, tolug, tolE, niter, nline, iCGmth, iscfmth, &
                 FermidE, itypeFermi, mCGbad, E_st, err_st, AL, &
                 nkpt, ntype, convergE, islda, igga, iwg_out, fwg_out, &
                 ivr_out, amx_mth, xgga, iwg_in, icoul, ivext_in, &
                 ilocal, iatom, ityatom, imov_at, totNel, Ealphat, t_start_petot)

  USE fft_data
  USE load_data
  USE data
  USE data_TDDFT
  USE data_pulay_6

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Parameters passed in
  INTEGER :: iforce_cal, ido_rho, ido_vr, niter, nline, mCGbad, nkpt, &
             ntype, islda, igga, iwg_out, iwg_in, ivr_out, &
             icoul, ivext_in, ilocal
  INTEGER, DIMENSION(100) :: iCGmth, iscfmth, itypeFermi
  INTEGER, DIMENSION(matom) :: iatom, ityatom
  INTEGER, DIMENSION(3,matom) :: imov_at

  REAL(8) :: Etot, tolug, tolE, convergE, totNel, Ealphat, xgga, t_start_petot
  REAL(8), DIMENSION(3, 3) :: AL
  REAL(8), DIMENSION(100) :: FermidE, amx_mth
  REAL(8), DIMENSION(3, matom) :: xatom, fatom
  REAL(8), DIMENSION(mst, nkpt, islda) :: E_st, err_st

  COMPLEX(8), DIMENSION(mr_n) :: workr_n

  CHARACTER(1), DIMENSION(2) :: fwg_out

  ! Etotcalc local variables
  INTEGER :: kpt, iislda
  REAL(8) :: TS, E_dDrho, TS0, E_dDrho0
  REAL(8), EXTERNAL :: UxcCA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The TDDFT variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: iiscf, itime, iscale, mmn0
  INTEGER, DIMENSION(100) :: iCGmth_td, iscfmth_td, itypeFermi_td
  INTEGER, DIMENSION(100) :: iCGmth_bo, iscfmth_bo, itypeFermi_bo

  REAL(8) :: EE, time, time_in, dt_step, dt, dt_fac, rkappa(3)
  REAL(8), DIMENSION(100) :: FermidE_td, FermidE_bo
  REAL(8), DIMENSION(mst,nkpt,islda) :: E_st0, E1, E2, occ, dos

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: dxatom_in

  COMPLEX(8), DIMENSION(mst,mst) :: cdot_tmp, cdot_in, cdot_in_tmp
  COMPLEX(8), DIMENSION(mst,mst,nkpt,islda) :: cdot
  COMPLEX(8), DIMENSION(mg_nx,mst,nkpt,islda) :: ug_n_bp_0

  COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: cc_pp_tmp
  COMPLEX(8), DIMENSION(:,:,:,:), ALLOCATABLE :: cc_pp, cc_pp0

  ! Internal time integral
  INTEGER :: itime_in, n_dt_now, n_dt_max
  COMPLEX(8) :: cfact1, cfact2, cfact3
  COMPLEX(8), DIMENSION(mst) :: HT_D, HT_0, cc_a, cc_a_tmp2
  COMPLEX(8), DIMENSION(mst,mst) :: HT_OffD, HT_S, HT_C, MTMP
  COMPLEX(8), DIMENSION(mst,mst,nkpt,islda) :: H_T, cc_a_0, cc_a_tmp

  ! Boltzmann factor (test)
  REAL(8) :: x_ltok
  REAL(8), DIMENSION(mst,mst) :: r_factor
  COMPLEX(8), DIMENSION(mst,mst) :: R_ltok_1, S_boltz, S_boltz_f, S_boltz_st, &
                                     U_boltz_st

  ! Other output
  REAL(8), DIMENSION(2, mst) :: frac
  REAL(8), DIMENSION(mst, mst, 2) :: dipole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Force & MD calc

  INTEGER :: ireplace, ivlct, ntemp
  INTEGER :: itherm, nhchain
  REAL(8), DIMENSION(7) :: qmass, xi, vxi, axi
  REAL(8) :: vxc1, vxc2, uxc1, uxc2
  REAL(8) :: InitTemp, DesiredTemp, TotalEn, Ekin, Box, Delt, rtmp
  REAL(8), DIMENSION(3,matom) :: fatom0, fatom1, DeltR
  REAL(8), DIMENSION(3,matom) :: Vi, V_output

  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: fatom_tmp, fatom_tmp2

  REAL(8), DIMENSION(mr_nL) :: vvion, vvionT, rrhocr
  REAL(8), DIMENSION(mr_nL,islda) :: vvr_nL
  REAL(8), DIMENSION(mr_n, islda) :: vvr_n
  REAL(8) :: vv0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! For temp & storage
  REAL(8), DIMENSION(mr_nL,islda,3) :: rho_storage_nL
  REAL(8), DIMENSION(mr_n) :: rho_tmp, rho_tmp2, rho_tmp3
  REAL(8), DIMENSION(mr_nL) :: rho_tmp_nL, rho_tmp2_nL, rho_tmp3_nL
  REAL(8), DIMENSION(mr_nL,islda) :: rho0, rho1
!  REAL(8) :: W(mst), rwork(10*mst)
  COMPLEX(8), DIMENSION(10*mst) :: work
!  COMPLEX(8), DIMENSION(mst,mst) :: SS1, SS2, SS3

  ! For timing
  INTEGER :: step_num, step_avail
  REAL(8) :: t_timewall, t_left, t_loop(10), t_loop_ave
  REAL(8) :: t_0, t_1, t_loop0, t_loop1, t_step0, t_step1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Misc
  INTEGER :: info, idone
  INTEGER, DIMENSION(mst) :: ipiv
  INTEGER :: imn, mx1, mx2, m, m1, m2, rho_flag, kkk, ii, it, nn1, nn2, itmp, itmp2
  INTEGER :: nnn, i, j, k, l, ierr, i1, i2, j1, j2, im, jj, ix, iy, iz, imax
  INTEGER :: inode_mst
  REAL(8) :: kT, sum0, sum2, sum3, sum4(islda), rsum, sum_max, rho_fac, omega, z1, z2
  REAL(8) :: tot_ma, dipole_x, dipole_y, dipole_z
  COMPLEX(8) :: cc,cc3,cc2,cc1,cmax,csum,csum1,csum2
  COMPLEX(8), DIMENSION(mr_n,mst) :: workr

  REAL(8), PARAMETER :: Hart = 27.211396d0, Bohr = 0.529177d0
  REAL(8), PARAMETER :: hbar2 = 0.65822d0 / Hart   ! hartree
  REAL(8), PARAMETER :: PI = 3.14159265d0
  COMPLEX(8), PARAMETER :: cai = (0.d0,1.d0), one = (1.d0,0.d0), zero = (0.d0,0.d0)

  CHARACTER(50) filename, fileindex, string
  CHARACTER(20), DIMENSION(2) :: frho_out, frho_out1, frho_out2, frho_out3
  CHARACTER(20), DIMENSION(2) :: fwg, fwg_st, fdw, fdR
  CHARACTER(20) focc_st, focc, fcc_st, ftmp_st, fpsi, fpsi_ug
  CHARACTER(20) frep_td, finp_td, fdx_td
  CHARACTER(20) fdate, ftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(mx .ne. nblock_band_mx) THEN
    PRINT *, "nblock_band_mx.ne.mx,stop", nblock_band_mx, mx
    CALL MPI_abort(MPI_COMM_WORLD, ierr)
    STOP
  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  finp_td = "TDDFT.input"
  frep_td = "report_TDDFT"

  focc_st = "occ_st"
  focc = "occ"
  fcc_st = "cc"
  fpsi = "rho_psi"
  fpsi_ug = "ug_psi"
  ftmp_st = "btmp_st"

  frho_out1(1) = "rho_st_1"
  frho_out1(2) = "rho_st_1_2"
  frho_out2(1) = "rho_st_2"
  frho_out2(2) = "rho_st_2_2"
  frho_out3(1) = "rho_st_3"
  frho_out3(2) = "rho_st_3_2"
  fwg_st(1) = "wg_st"
  fwg_st(2) = "wg_st_2"
!      fdw(1) = "dw_st_spin1"
!      fdw(2) = "dw_st_spin2"
!      fdR(1) = "dR_st_spin1"
!      fdR(2) = "dR_st_spin2"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iCGmth_td = 3
  iscfmth_td = 0
  FermidE_td = FermidE
  itypeFermi_td = 1  ! for test

  iCGmth_bo = 3
  iscfmth_bo = 1
  FermidE_bo = FermidE
  itypeFermi_bo = 1  ! for test

  kpt = 1
  ng_n = ngtotnod(inode, kpt)
  nr_nL = n1L * n2L * n3L / nnodes

  npulay_max_6 = 40
  n_dt_max = 1000

  fatom = 0.d0
  DeltR = 0.d0
  occ = 1.d0
  E1 = 0.d0
  E2 = 0.d0
!      E_output=0.d0
  rho_storage_nL = 0.d0
!      cc_phase=one

  qmass(1:7) = 0.d0
  xi(1:7) = 0.d0
  vxi(1:7) = 0.d0
  axi(1:7) = 0.d0

  gz_1 = 0.d0
  gz_2 = 0.d0

  frac = 0.d0

  step_num = 0
  t_loop = 0.d0
  t_loop_ave = 0.d0

  dipole = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL input_TDDFT(totNel, mx, mmn0, iwg_in, imax, rkappa, &
                   ivlct, ntemp, Vi, itherm, nhchain, qmass, &
                   xi, vxi, axi, &
                   InitTemp, DesiredTemp, &
                   finp_td, ftmp_st, fdx_td, t_timewall)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! Try to heat the system
!  Vi(:,:) = Vi(:,:) * sqrt(300.d0/4.37d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Box = 1.d0
  dt_step = dtMD                                           ! unit fs
  Delt = dtMD / 2.418884D-2                                ! a.u.
  qmass = qmass * Hart                                     ! a.u.
  temperature = InitTemp
  kT = temperature * 8.6173324d-5 / Hart
  time = init_time - dt_step

  DO i = 1, natom
    MDatom(i) = MDtype(ityatom(i)) / 5.485799D-4
  END DO

  CALL data_allocate_6(mr_nL, islda, npulay_max_6)
  CALL data_allocate_TDDFT(mr_n, mg_nx, mst, mmn, nkpt, islda)
  cpsi_td = zero

  ALLOCATE(cc_pp0(mst, mmn, nkpt, islda))
  ALLOCATE(cc_pp(mst, mmn, nkpt, islda))
  ALLOCATE(cc_pp_tmp(mst, mmn))
  IF(ifatom.eq.2) ALLOCATE(dxatom_in(3, matom, (ntime_init+1):(ntime_init+ntime)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(iCGmth(1).eq.-1) CALL output_psi()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open file frep_td for output
  ! if it already exists, rename it then open a new one

  IF(inode.eq.1) THEN

    OPEN(17, FILE=frep_td, STATUS='old', ACTION='READ', IOSTAT=ierr)

    IF(ierr.eq.0) THEN
      CLOSE(17)
      CALL DATE_AND_TIME(fdate, ftime)
      filename=trim(adjustl(frep_td))//"."//trim(adjustl(fdate))//trim(adjustl(ftime))
      CALL system('mv '//trim(adjustl(frep_td))//' '//filename)
    ELSE
      CLOSE(17)
    END IF

    OPEN(17, FILE=frep_td)
    REWIND(17)

    WRITE(17,*)       "Input parameters"
    WRITE(17,*)       "*********************************************"
    WRITE(17,FMT1701) "  nbasis                =", mst2
    WRITE(17,FMT1701) "  nbands                =", mmn

    IF(mmn0.gt.1) THEN
    WRITE(17,FMT1701) "  -> nband0             =", mmn0
    END IF

    WRITE(17,FMT1703) "  dt(fs)                =", dtMD
    WRITE(17,FMT1701) "  imd                   =", iMD
    WRITE(17,FMT1703) "  temperature(K)        =", temperature

    IF(itherm.eq.1) THEN
    WRITE(17,'(A)')   "  -> Nose-Hoover Chain:"
    WRITE(17,FMT1701) "  -> chain length       =", nhchain
    WRITE(17,'(A)')   "  -> N-H masses(eV^-1), frequences(eV)"
    WRITE(17,'(A,E15.8,1X,E15.8)') &
                      "       ", &
                      qmass(1)/Hart, &
                      sqrt(dble(ntemp)*3.d0*DesiredTemp*8.6173324d-5*Hart/qmass(1))
    WRITE(17,FMT1703) "  -> Final temp(K)      =",DesiredTemp
    END IF

    WRITE(17,FMT1701) "  nelm                  =", mscf
    WRITE(17,FMT1702) "  rhodiff(e-)           =", tolrho
    WRITE(17,FMT1702) "  intglerr              =", tolintgl
    WRITE(17,FMT1701) "  istep                 =", ntime_init
    WRITE(17,FMT1703) "  starttime(fs)         =", init_time
    WRITE(17,FMT1701) "  nstep                 =", ntime
    WRITE(17,FMT1701) "  iforce                =", ifatom
    WRITE(17,FMT1701) "  ivdt                  =", ivext_dt
    WRITE(17,FMT1701) "  ikappa                =", ikappa

    IF(ikappa.gt.0) THEN
    WRITE(17,'(a,3(1x,e10.3),a)') &
                      "  -> rkappa             ="
    WRITE(17,'(a,3(1x,e10.3),a)') &
                      "     (",(rkappa(j), j=1,3), ")"
    END IF
    WRITE(17,FMT1701) "  iboltz                =", iboltz
    WRITE(17,FMT1701) "  ibandshift            =", ibshift
    WRITE(17,FMT1701) "  iscale                =", i_scale
    WRITE(17,FMT1701) "  iexci                 =", iocc

    IF(iocc.gt.0) THEN
    WRITE(17,FMT1701) "  -> hole state         =", jhole
    WRITE(17,FMT1701) "  -> elect state        =", jelec
    WRITE(17,FMT1701) "  -> nexci              =", jxi
    WRITE(17,FMT1703) "  -> rexci(e-)          =", rxi
    END IF

    WRITE(17,'(a)')   "  Nuclei masses         ="
    WRITE(17,'(a,11(1x,f10.3))') &
                      "  ", (MDtype(j), j=1,imax)

    IF(imp_k.ne.0.d0) THEN
    WRITE(17,FMT1703) "  Particle E_k          =", imp_k
    END IF

    IF(ibo_md.gt.0) THEN
    WRITE(17,'(a)')   "  Do Adiabatic MD       =  TRUE"
    END IF

    WRITE(17,*)       "*********************************************"
    CALL system_flush(17)

  END IF ! IF(inode.eq.1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(inode.eq.1) THEN
    WRITE(6,*)
    WRITE(6,*) '-------------------------------------------------'
    WRITE(6,*) '|                                               |'
    WRITE(6,*) '|                 TDDFT BEGIN                   |'
    WRITE(6,*) '|                                               |'
    WRITE(6,*) '-------------------------------------------------'
    WRITE(17,*)
    WRITE(17,*) "TDDFT start"
  END IF

  IF(ikappa.gt.0) THEN
   inormal = 0
  ELSE
   inormal = 1
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize TDDFT
  CALL initial_TDDFT()
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(nkpt.eq.1) THEN
    kpt = 1
    CALL gen_G_comp(kpt, 0)
    CALL fftprep_comp(n1, n2, n3)
    ng_n = ngtotnod(inode, kpt)
  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! 1st loop level: For every time step
  ! 1) at t=t2, check xatom, v_ext, potential correction, initial rho, etc.
  ! 2) do 2nd loop: leapfrog between t1 and t2
  ! 3) obtain converged time-evolved wavefuntion, goto output
  string = "initial itime timing"
  CALL timing_mpi(string,t_0)
  t_step0 = t_0

Loop_itime: DO itime = ntime_init + 1, ntime_init + ntime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    time=time+dt_step
    idone=-2

    IF(inode.eq.1) THEN
      WRITE(6,*)
      WRITE(6,*) "##################################################"
      WRITE(6,*) " Time step =", itime-1
      WRITE(6,*) " Time      =", time
      WRITE(6,*)
      WRITE(17,*)
      WRITE(17,*) "****************************************"
      WRITE(17,FMT1701) " Time step         =",itime-1
      WRITE(17,*) "------------------------------"
      CALL system_flush(17)
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if ifatom == 2, read xatom for each step from external file
    IF(ifatom.eq.2) xatom(:, :) = dxatom_in(:, :, itime)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iscale = i_scale
    IF(itime.le.3) iscale = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calc the semiconductor surface polarization effect
    if(inode.eq.1.and.ibshift.eq.1) call bandshift0(xatom, iatom, AL)
    !if(inode.eq.1) call bandshift0(xatom, iatom, AL)
    if(ibshift.eq.1) then
      call mpi_bcast(dP, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(gz_1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(gz_2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! check external field
    CALL calc_vext()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! give an initial value to rho(t), if needed

    IF(itime.gt.1.and.(itime-1).le.jxi.and.iocc.gt.0) THEN
    ! If exitation exists

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iscale = 0
      ! Change the original occupation occ0
      occ0(jhole, 1, 1) = occ0(jhole, 1, 1) - rxi / dble(jxi)
      occ0(jelec, 1, 1) = occ0(jelec, 1, 1) + rxi / dble(jxi)

      CALL unfmtIO_r(18, focc_st, occ0, mst*nkpt*islda, 1, 0)

      IF(inode.eq.1) THEN
        WRITE(6,*)
        WRITE(6,*) "------------------------------"
        WRITE(6,*) "| Change Occupation"
        WRITE(6,*) '| ', jhole, '=', jhole, '-', rxi/dble(jxi), '=', occ0(jhole,1,1)
        WRITE(6,*) '| ', jelec, '=', jelec, '+', rxi/dble(jxi), '=', occ0(jelec,1,1)
        WRITE(6,*) "------------------------------"
        WRITE(6,*)
      END IF

      ! use new occ0 to gen new rho
      ! this rho will become the initial input for next step

!      DO iislda=1,islda
!
!        rho_tmp=0.d0
!        DO kpt=1,nkpt
!          IF(nkpt.ne.1) THEN
!            CALL gen_G_comp(kpt,0)
!            CALL fftprep_comp(n1,n2,n3)
!            ng_n = ngtotnod(inode,kpt)
!          END IF
!          workr=zero
!          CALL d3fft_comp_block(cpsi_td(1,1,kpt,iislda), workr, -1, kpt, mst)
!          DO m=1,mst
!            rho_tmp(1:nr_n) = rho_tmp(1:nr_n) + occ0(m,kpt,iislda) * ABS(workr(1:nr_n,m))**2
!          END DO
!        END DO
!
!        CALL convert_SLvr(rho_tmp, rho_tmp_nL, 1)
!        rho_nL(:, iislda) = rho_tmp_nL(:)
!
!      END DO

      IF(itime.gt.1.and.itime.lt.4) THEN 
        rho_nL(:,:) = rho_storage_nL(:,:,3)
      ELSEIF(itime.ge.4) THEN
        rho_nL(:,:) = 3 * rho_storage_nL(:,:,3) - &
                      3 * rho_storage_nL(:,:,2) + &
                      rho_storage_nL(:,:,1)
      END IF

      ! end of IF(itime.gt.1.and.(itime-1).le.jxi.and.iocc.gt.0)

    ELSEIF(itime.gt.1.and.itime.lt.4) THEN 

      rho_nL(:,:) = rho_storage_nL(:,:,3)

    ELSEIF(itime.ge.4) THEN

      rho_nL(:,:) = 3 * rho_storage_nL(:,:,3) - &
                    3 * rho_storage_nL(:,:,2) + &
                    rho_storage_nL(:,:,1)
            
    END IF
     
    rho0=rho_nL
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 2nd loop level, do charge density leapfrog until charge density converged:
    ! 1) rho0 --Etotcalc.f--> adiabatic ug_n_bp
    ! 2) ug_n_bp_0 + ug_n_bp -> H(t)
    ! 3) do 3rd loop level: calc time-evolved wavefunction cpsi_td
    ! 4) cpsi_td + occ0 -> rho1, rho0 + rho1 --pulay mixing--> rho0, back to 1)

Loop_iiscf: DO iiscf = 1, mscf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(inode.eq.1) THEN

        WRITE(6,*)
        WRITE(6,*) "------------------------------"
        WRITE(6,*) "| Leapfrog Step",iiscf
        WRITE(6,*) "------------------------------"
        WRITE(6,*)
        WRITE(6,*) "------------------------------"
        WRITE(6,*) "| Adiabatic States Calc"
        WRITE(6,*) "------------------------------"
        WRITE(6,*)

      END IF

      string = "initial iiscf timing"
      CALL timing_mpi(string,t_0)
      t_loop0 = t_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(itime.eq.1) THEN
      ! t=0, no more work need to do.
      ! calc force and Etot, then exit loop, goto output routine

        idone=1

        IF(ibo_md.le.0) THEN
          CALL recalc()
        ELSEIF(ibo_md.gt.0) THEN
          occ0 = occ
          E_st0 = E_st
          cc_pp0 = cc_pp
          rho0 = rho_nL
          rho_storage_nL(:,:,1) = rho_storage_nL(:,:,2)
          rho_storage_nL(:,:,2) = rho_storage_nL(:,:,3)
          rho_storage_nL(:,:,3) = rho0(:,:)
        END IF

        EXIT Loop_iiscf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSEIF(itime.gt.1) THEN
      ! t>0

        IF(iiscf.gt.1) THEN
        ! not the first leapfrog step, use rho from last step
          rho_nL=rho0
        ENDIF

        string = "before Etotcalc"
        CALL timing_mpi(string,t_0)

        IF(ibo_md.le.0) THEN

          iforce_cal = 0
          ireplace = 0
          nnn = 1
          CALL Etotcalc(xatom, fatom, workr_n, Etot, iforce_cal, 0, 1, tolug, &
                        tolE, nnn, nline, iCGmth_td, iscfmth_td, FermidE_td, &
                        itypeFermi_td, mCGbad, E_st, err_st, AL, nkpt, ntype, &
                        convergE, islda, igga, iwg_out, fwg_out, ivr_out, &
                        amx_mth, xgga, occ, ireplace, TS, E_dDrho)

        ELSEIF(ibo_md.gt.0) THEN
        ! BOMD only needs Etotcalc to get the ground states, so exit 2nd
        ! loop level right after calling Etotcalc.f
        ! inner loops like 3rd and 4nd loop level will be skipped

          fatom = 0.d0
          iforce_cal = 1
          ireplace = 0

          ! itime>1, already has good rho, use iscfmth_bo (no nonsc steps) to save time
          CALL Etotcalc(xatom, fatom, workr_n, Etot, iforce_cal, 0, 1, tolug, &
                        tolE, niter, nline, iCGmth_bo, iscfmth_bo, FermidE_bo, &
                        itypeFermi_bo, mCGbad, E_st, err_st, AL, nkpt, ntype, &
                        convergE, islda, igga, iwg_out, fwg_out, ivr_out, &
                        amx_mth, xgga, occ, ireplace, TS, E_dDrho)

          IF(islda.eq.1.and.nkpt.eq.1) THEN
            ug_n_bp_0(:,:,1,1)=ug_n_bp(:,:)
          ELSE
            !ug_n_bp_0=ug_all
          END IF

          occ0 = occ
          E_st0 = E_st
          cc_pp0 = cc_pp
          rho0 = rho_nL

          rho_storage_nL(:,:,1) = rho_storage_nL(:,:,2)
          rho_storage_nL(:,:,2) = rho_storage_nL(:,:,3)
          rho_storage_nL(:,:,3) = rho0(:,:)

        END IF ! if BOMD

        string = "Etotcalc"
        CALL timing_mpi(string,t_0)

        IF(ibo_md.gt.0) THEN
          idone = 1
          EXIT Loop_iiscf
        END IF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate <phi(t1)|phi(t2)>

        cdot=zero
        cc_pp=cc_pp0

        DO iislda=1,islda
        DO kpt=1,nkpt

          IF(nkpt.ne.1) THEN
            CALL gen_G_comp(kpt,0)
            CALL fftprep_comp(n1,n2,n3)
            ng_n = ngtotnod(inode,kpt)
          END IF

          CALL ugIOBP(ug_n_bp, kpt, 2, 0, iislda, -1, nkpt, islda)

          E1=E_st0
          E2=E_st

          cdot_tmp=zero
          CALL dot_product_BP(ug_n_bp_0(1,1,kpt,iislda),ug_n_bp(1,1),ng_n,cdot_tmp)
          cdot_tmp=cdot_tmp*vol

          string = "<phi_t1|phi_t2>"
          CALL timing_mpi(string,t_0)

          ! Orthonomalization cdot
          call ortho_Gram(cdot_tmp, mst, mst)
          IF(mst2.lt.mst) cdot_tmp(mst2+1:mst,mst2+1:mst) = zero

          string = "orthonormal"
          CALL timing_mpi(string,t_0)

          ! fix the phase difference between phi(t1) and phi(t2)
          DO m2=1,mst2
            !cmax=zero
            !DO m1=1,mst2
            !  IF(abs(cdot_tmp(m1,m2)).gt.abs(cmax)) THEN
            !    cmax = cdot_tmp(m1,m2)
            !  END IF
            !END DO
            cmax = cdot_tmp( maxloc( abs(cdot_tmp(1:mst2, m2)), 1 ) , m2)
            cmax = cmax / abs(cmax)
            cdot_tmp(:,m2) = conjg(cmax) * cdot_tmp(:,m2)
            ug_n_bp(:,m2) = conjg(cmax) * ug_n_bp(:,m2)
          END DO

          cdot(:,:,kpt,iislda)=cdot_tmp(:,:)

        END DO ! kpt
        END DO ! islda

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSE ! errors on itime or iocc

        WRITE(6,*) "--ERROR: itime error",itime
        CALL MPI_barrier(MPI_COMM_WORLD,ierr)
        CALL MPI_abort(MPI_COMM_WORLD,ierr)
        STOP

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END IF

      ! Calc how many states on each processor
      inode_mst = mod(mst2, nnodes_tot)
      IF(inode_mst.eq.0) THEN
        inode_mst = mst2 / nnodes_tot
      ELSE
        inode_mst = (mst2 - inode_mst) / nnodes_tot + 1
      END IF
      IF(inode.eq.1) THEN
        WRITE(6,*) "States num on each core:",inode_mst
        WRITE(6,*) "basis size:",mst2
        WRITE(6,*)
      END IF

      cc_a_0=zero
      cc_a=zero
      H_T=zero

      ! H_T = <phi(t1)| H(t2) |phi(t1)>

      DO iislda=1,islda
      DO kpt=1,nkpt

        DO i1=1,mst2
          cdot_tmp(1:mst2,i1)=conjg(cdot(i1,1:mst2,kpt,iislda))*E2(1:mst2,kpt,iislda)
        END DO

        CALL zgemm('N','N',mst,mst,mst,one, &
                   cdot(1,1,kpt,iislda),mst,cdot_tmp,mst,one, &
                   H_T(1,1,kpt,iislda),mst)

        ! Boltzmann factor (in test)
        IF(iboltz.gt.0) CALL calc_boltz()

      END DO
      END DO

      string = "parallel ready to go"
      CALL timing_mpi(string,t_0)
      t_0=mpi_wtime()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 3rd loop level

Loop_imn: DO imn=(inode-1)*inode_mst+1,inode*inode_mst

        n_dt_now = 0
        IF(imn.gt.mst2) EXIT Loop_imn  ! all states on this proccesser are done

in_lp1: DO iislda = 1,islda
in_lp2: DO kpt = 1,nkpt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! HT_D and HT_OffD is the diagonal and off-diagonal terms of H_T
          ! cc_a_0 stores D_ij == <phi(t1)|phi(t)>
          ! At the beginning of each time step, D_ij is unit diagonal matrix,
          ! then D_ij evaluates with H(t)

          cc_a_0(imn,imn,kpt,iislda)=one
          HT_D = zero
          HT_OffD(:,:) = H_T(:,:,kpt,iislda)
          DO i=1,mst2
            HT_D(i)=HT_OffD(i,i)
            HT_OffD(i,i)=zero
          END DO

          ! define HT_S == HT_OffD^2, HT_C == HT_OffD^3
          HT_S=zero
          CALL zgemm('N','N',mst,mst,mst,one,HT_OffD,mst,HT_OffD,mst,one,HT_S,mst)
          HT_C=zero
          CALL zgemm('N','N',mst,mst,mst,one,HT_S,mst,HT_OffD,mst,one,HT_C,mst)

          ! n_dt_now is the internal time step number
          ! set min(n_dt_now) = 10
          sum_max = 0.d0
          DO i2=1,mst2
            DO i1=1,mst2
              csum = HT_OffD(i1,i2)
              IF(abs(csum).gt.sum_max) sum_max = abs(csum)
            END DO
          END DO
          sum_max = sum_max * dt_step / hbar2
          n_dt_now = INT(sum_max / tolintgl) + 1
          IF(n_dt_now.lt.10) n_dt_now = 10
          IF(n_dt_now.gt.n_dt_max) n_dt_now = n_dt_max

          WRITE(6,'(a,2(i5,1x),i7)') "     inode, imn, steps, ",inode,imn,n_dt_now

          dt = dt_step / (n_dt_now - 0.d0)
          dt_fac = dt / hbar2
          time_in = 0.d0

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 4th loop level, D_ij integral in [t1, t2]

Loop_itime_in: DO itime_in = 1,n_dt_now

            time_in = time_in + dt

            ! HT_0 = exp(-i * H0(t) * dt / 2)

            HT_0(1:mst2) = cai * dt_fac * ( E1(1:mst2,kpt,iislda) &
                                            + (time_in - dt / 2.d0) &
                                            * (HT_D(1:mst2) - E1(1:mst2,kpt,iislda)) &
                                            / dt_step &
                                          )
            HT_0(1:mst2) = exp(-HT_0(1:mst2) / 2.d0)

            ! MTMP = HT_0 * Taylor_term^3(-i * HT_OffD * dt) * HT_0

            cfact1 = cai * (time_in-dt/2.d0) * dt_fac / dt_step
            cfact2 = cfact1**2
            cfact3 = cfact1**3
            cfact1 = -cfact1
            cfact2 = cfact2 / 2.d0
            cfact3 = -cfact3 / 6.d0

            DO i2=1,mst2
            DO i1=1,mst2

              cc2=zero
              IF(i1.eq.i2) cc2=one

              MTMP(i1,i2) = HT_0(i1) * &
                            ( cc2 + cfact1 * HT_OffD(i1,i2) + &
                                    cfact2 * HT_S(i1,i2) + &
                                    cfact3 * HT_C(i1,i2) &
                            ) * HT_0(i2)

            END DO
            END DO

            cc_a_tmp2(:) = cc_a_0(:,imn,kpt,iislda)
            cc_a = zero
            CALL zgemm('N','N',mst,1,mst,one,MTMP,mst,cc_a_tmp2,mst,one,cc_a,mst)
            cc_a_0(:,imn,kpt,iislda) = cc_a(:)

          END DO Loop_itime_in   ! 4th loop level | itime_in=1,n_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       END DO in_lp2
       END DO in_lp1

      END DO Loop_imn    ! 3rd loop level | imn

      t_1 = mpi_wtime()
      WRITE(6,'(a,i5,2x,e10.3)') "   Parallel done, inode, time, ",inode,t_1-t_0

      ! Now D_ij(t2) has been stored in cc_a_0

      cc_a_tmp=zero
      CALL mpi_allreduce(cc_a_0,cc_a_tmp,mst*mst*nkpt*islda, &
                         MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      cc_a_0=cc_a_tmp

      string = "parallel ends, calc psi_td"
      CALL timing_mpi(string,t_0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! enforce orthogonalization
      DO iislda=1,islda
      DO kpt=1,nkpt
        call ortho_Gram(cc_a_0(1,1,kpt,iislda), mst, mst)
      END DO  ! kpt
      END DO  ! iislda

      string = "D_ij orthonormal"
      CALL timing_mpi(string,t_0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! new cc_pp = D_ij * old cc_pp

      cc_pp=cc_pp0

      DO iislda=1,islda
      DO kpt=1,nkpt

        IF(nkpt.ne.1) THEN
          CALL gen_G_comp(kpt,0)
          CALL fftprep_comp(n1,n2,n3)
          ng_n=ngtotnod(inode,kpt)
        END IF
        CALL ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)

        cc_pp_tmp = zero
        CALL zgemm('N','N',mst,mmn,mst,one, &
                   cc_a_0(1,1,kpt,iislda),mst,cc_pp(1,1,kpt,iislda), &
                   mst,one,cc_pp_tmp,mst)
        cc_pp(:,:,kpt,iislda) = cc_pp_tmp(:,:)

        ! transfer basis from phi(t1) to phi(t2)

        cc_pp_tmp = zero
        cdot_tmp = zero

        DO i1=1,mst2
        DO i2=1,mst2
          cdot_tmp(i2,i1) = conjg(cdot(i1,i2,kpt,iislda))
        END DO
        END DO

        CALL zgemm('N','N',mst,mmn,mst,one, &
                   cdot_tmp,mst,cc_pp(1,1,kpt,iislda), &
                   mst,one,cc_pp_tmp,mst)

        ! normalize if needed
        IF(inormal.eq.1) THEN
          IF(mmn0.gt.1) THEN
            cc_pp_tmp(:,1:mmn0-1) = zero
            DO i = 1,mmn0-1
              cc_pp_tmp(i,i) = one
            ENDDO
          END IF
          call ortho_Gram(cc_pp_tmp, mst, mmn)
          IF(mst2.lt.mst) cc_pp_tmp(mst2+1:mst,1:mmn) = zero
        END IF
        cc_pp(:,:,kpt,iislda)=cc_pp_tmp(:,:)

      END DO
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL calc_psi()
      CALL diff_rho()
!      GOTO 5328
      idone = check_converge()

      IF(idone.eq.0) THEN
        CYCLE Loop_iiscf
      ELSEIF(idone.eq.1.or.idone.eq.-1) THEN
        EXIT Loop_iiscf
      ELSE
        PRINT *, "--ERROR: Convergency error"
        STOP
      END IF

    END DO Loop_iiscf  ! iiscf


    IF( iiscf.gt.mscf .or. &
        (iiscf.le.mscf.and.idone.eq.0) .or. &
        (idone.ne.-1.and.idone.ne.0.and.idone.ne.1) &
      ) THEN
      PRINT *, "--ERROR: Loop_iiscf error"
      STOP
    END IF


    IF(idone.eq.-1 .and. mod((itime-1),noutput).eq.0) THEN
    ! if iiscf reaches mscf, we should output the data of last step, then stop
    ! However, if output has already been done for last step, then just skip the output
      GOTO 9998
    END IF


    IF(idone.eq.1) THEN
    ! only do this part while iiscf gets converged
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0. calc states distribution and dipole(i,j)
      IF(ibshift.eq.1) THEN

      string = "Before calc_dipole"
      CALL timing_mpi(string,t_0)

      itmp = 300
      itmp2 = mmn
      call calc_dipole(1, AL, nkpt, islda, frac, dipole, workr, totNel, itmp, itmp2)

      string = "calc_dipole"
      CALL timing_mpi(string,t_0)

      IF(1.eq.1.and.inode.eq.1) THEN

        filename="phi_dipole."
        WRITE(fileindex,'(i)') nkpt
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') islda
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') itime-1
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))

        OPEN(29,FILE=filename)
        REWIND(29)
        DO i = itmp, itmp2
          WRITE(29,180) (dipole(i,j,1), j=itmp, itmp2)
        END DO
        CLOSE(29)

        filename="psi_dipole."
        WRITE(fileindex,'(i)') nkpt
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') islda
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') itime-1
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))

        OPEN(29,FILE=filename)
        REWIND(29)
        DO i = itmp, itmp2
          WRITE(29,180) (dipole(i,j,2), j=itmp, itmp2)
        END DO
        CLOSE(29)

      ENDIF

      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 1. calc the total dipole moment, by Jie Ma
      tot_ma = 0.d0
      dipole_x = 0.d0
      dipole_y = 0.d0
      dipole_z = 0.d0

      DO ii=1,nr_n
        jj=ii+(inode-1)*nr_n
        ix=(jj-1)/(n2*n3)+1
        iy=(jj-1-(ix-1)*n2*n3)/n3+1
        iz=jj-(ix-1)*n2*n3-(iy-1)*n3
        tot_ma=tot_ma+rho_storage_nL(ii,1,3)*vol/nr_n/nnodes
        dipole_x=dipole_x+rho_storage_nL(ii,1,3)*vol/nr_n/nnodes*dble(ix-dble(n1)/2-1)
        dipole_y=dipole_y+rho_storage_nL(ii,1,3)*vol/nr_n/nnodes*dble(iy-dble(n2)/2-1)
        dipole_z=dipole_z+rho_storage_nL(ii,1,3)*vol/nr_n/nnodes*dble(iz-dble(n3)/2-1)
      END DO
      CALL global_sumr(tot_ma)
      CALL global_sumr(dipole_x)
      CALL global_sumr(dipole_y)
      CALL global_sumr(dipole_z)
      IF(inode.eq.1) THEN
        WRITE(17,FMT1704) "  Total charge     =",tot_ma
        WRITE(17,FMT1704) "  Dipole X         =",dipole_x
        WRITE(17,FMT1704) "  Dipole Y         =",dipole_y
        WRITE(17,FMT1704) "  Dipole Z         =",dipole_z
        t_loop1 = mpi_wtime()
        WRITE(17,FMT1703) "  Time             =",t_loop1-t_loop0
        t_loop0 = t_loop1
        CALL system_flush(17)
      END IF
      ! 1. end
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 2.
      IF(inode.eq.1) THEN
        WRITE(6,*) 
        WRITE(6,*) "------------------------------"
        WRITE(6,*) "| Move Atoms"
        WRITE(6,*) "------------------------------"
        WRITE(6,*)
        t_0=mpi_wtime()
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DEBUG Artificial cooling down
!      IF( itime.gt.1 .and. mod((itime-1),10).eq.0 ) THEN
!        Vi=0.d0
!        IF(inode.eq.1) WRITE(6,*) " Artificial cooling down at itime=",itime
!      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(ifatom.ne.2) THEN

        IF(ifatom.eq.1) THEN
          fatom0=fatom
          fatom1=fatom
        ELSE
          fatom0=0.d0
          fatom1=0.d0
        END IF

        CALL MVATOMS(itime, iscale, Delt, Box, xatom, fatom0, fatom1, &
                     AL, Etot, InitTemp, DesiredTemp, TotalEn, Ekin, &
                     ivlct, ntemp, Vi, V_output, &
                     itherm, nhchain, qmass, xi, vxi, axi)

        string = "atoms moving ends"
        CALL timing_mpi(string,t_0)

      ELSE

        IF(itime.ne.(ntime_init+ntime)) THEN
          xatom(:, :) = dxatom_in(:, :, itime+1)
        ELSE
          xatom(:, :) = 0.d0
        ENDIF

      ENDIF
      ! 2. end
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 3.
      IF(inode.eq.1) THEN

        OPEN(29,FILE="dfatom0", access="append")
        WRITE(29,*) "Direct configuration=", itime-1
        DO i=1,natom
          WRITE(29,296) iatom(i), (fatom0(j,i), j=1,3), (imov_at(j,i), j=1,3)
        END DO
        CLOSE(29)

        OPEN(29,FILE="dfatom1", access="append")
        WRITE(29,*) "Direct configuration=", itime-1
        DO i=1,natom
          WRITE(29,296) iatom(i), (fatom1(j,i), j=1,3), (imov_at(j,i), j=1,3)
        END DO
        CLOSE(29)

        OPEN(29,FILE="dV", access="append")
        WRITE(29,*) "Direct configuration=", itime-1
        DO i=1,natom
          WRITE(29,296) iatom(i), (V_output(j,i), j=1,3), (imov_at(j,i), j=1,3)
        END DO
        CLOSE(29)

        OPEN(29,FILE="dxatom", access="append")
        WRITE(29,*) "Direct configuration=", itime
        DO i=1,natom
          WRITE(29,295) iatom(i), (xatom(j,i), j=1,3), (imov_at(j,i), j=1,3)
        END DO
        CLOSE(29)

        IF(itherm.eq.1) THEN
          OPEN(29,FILE="dtherm", access="append")
          WRITE(29,*) "NH THERMOSTAT", itime
          DO i=1,nhchain
            WRITE(29,'(3(1X,E15.8))') xi(i), vxi(i), axi(i)
          ENDDO
          CLOSE(29)
        ENDIF

        IF(ibo_md.gt.0) THEN
          OPEN(29,FILE="update_xatom")
          REWIND(29)
          WRITE(29,'(i8)') natom
          WRITE(29,'(3(f13.7,1x))') AL(:,1)
          WRITE(29,'(3(f13.7,1x))') AL(:,2)
          WRITE(29,'(3(f13.7,1x))') AL(:,3)
          DO i=1,natom
            WRITE(29,295) iatom(i), (xatom(j,i), j=1,3), (imov_at(j,i), j=1,3)
          END DO
          CLOSE(29)
        END IF

295     FORMAT(i4,2x,3(f15.9,1x),2x,3(i2,1x))
296     FORMAT(i4,2x,3(E15.8,1x),2x,3(i2,1x))

!        filename="r_factor."
!        WRITE(fileindex,'(i6)') itime-1
!        filename=trim(adjustl(filename))//trim(adjustl(fileindex))
!        OPEN(18,FILE=filename,FORM="unformatted")
!        REWIND(18)
!        WRITE(18) r_factor
!        CLOSE(18)
!
!        filename="S_boltz."
!        WRITE(fileindex,'(i6)') itime-1
!        filename=trim(adjustl(filename))//trim(adjustl(fileindex))
!        OPEN(18,FILE=filename,FORM="unformatted")
!        REWIND(18)
!        WRITE(18) S_boltz_st
!        CLOSE(18)
!
!        filename="U_boltz."
!        WRITE(fileindex,'(i6)') itime-1
!        filename=trim(adjustl(filename))//trim(adjustl(fileindex))
!        OPEN(18,FILE=filename,FORM="unformatted")
!        REWIND(18)
!        WRITE(18) U_boltz_st
!        CLOSE(18)

      END IF
      ! 3. end
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 4. occupation
      DO iislda=1,islda
      DO kpt=1,nkpt

        filename=trim(adjustl(focc))//"."
        WRITE(fileindex,'(i)') kpt
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') iislda
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex,'(i)') itime-1
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))

        IF(inode.eq.1) THEN
          OPEN(18,FILE=filename)
          REWIND(18)
          DO i=1,mst
            IF(ibo_md.le.0) THEN
              IF(i.le.mmn) THEN
                WRITE(18,180) E_st(i,kpt,iislda) * Hart, dos(i,kpt,iislda), &
                              E_td(i,kpt,iislda) * Hart, occ0(i,kpt,iislda), &
                              frac(1, i), frac(2, i)
              ELSE
                WRITE(18,180) E_st(i,kpt,iislda) * Hart, dos(i,kpt,iislda), &
                              E_st(i,kpt,iislda) * Hart, dos(i,kpt,iislda), &
                              frac(1, i), frac(2, i)
              END IF
            ELSE
              WRITE(18,180) E_st(i,kpt,iislda) * Hart, occ0(i,kpt,iislda), &
                            frac(1, i), frac(2, i)
            END IF
          END DO
          CLOSE(18)
        END IF

      END DO
      END DO
      ! 4. end
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 5.
      IF(ibo_md.le.0) THEN
        !!!!!!!!!!!!!!!!!!
        ! output cc
        filename=trim(adjustl(fcc_st))//"."
        WRITE(fileindex,'(i)') itime-1
        filename=trim(adjustl(filename))//trim(adjustl(fileindex))
        CALL unfmtIO_c(18,filename,cc_pp0,mst*mmn*nkpt*islda,1,0)
        !!!!!!!!!!!!!!!!!!
        ! output H(t2)
        ! filename="H_t2."
        ! WRITE(fileindex,'(i)') itime-1
        ! filename=trim(adjustl(filename))//trim(adjustl(fileindex))
        ! CALL unfmtIO_c(18,filename,H_T,mst*mst*nkpt*islda,1,0)
      END IF
      ! 5. end
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! if idone.eq.1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF


    ! output rho and ug when 1. mod(itime,noutput)=0,
    ! or 2. itime reaches mtime,
    ! or 3. iiscf reaches mscf
    IF( mod((itime-1),noutput).eq.0 .or. &
        itime.eq.ntime_init+ntime .or. &
        idone.eq.-1 &
      ) THEN

      IF(idone.ne.-1) THEN
        WRITE(fileindex,'(i)') itime-1
      ELSE
        WRITE(fileindex,'(i)') itime-2
      ENDIF

      IF(inode.eq.1) THEN

        WRITE(17,*) "------------------------------"
        IF(idone.ne.-1) THEN
          WRITE(17,*) " Output temp data"
        ELSE
          WRITE(17,*) " Output temp data for last step"
        END IF
        WRITE(17,*) "------------------------------"
        CALL system_flush(17)

        filename=trim(adjustl(ftmp_st))//"."//trim(adjustl(fileindex))

        OPEN(18,FILE=filename,FORM='unformatted')
        REWIND(18)
        IF(idone.ne.-1) THEN
          WRITE(18) itime, time
        ELSE
          WRITE(18) itime-1, time-dt_step
        END IF
        WRITE(18) xatom
        WRITE(18) Vi
        WRITE(18) cc_pp0
        WRITE(18) E_st0
        WRITE(18) Ealphat,TS0,E_dDrho0
        WRITE(18) Delt,Box,InitTemp,DesiredTemp,TotalEn,Ekin
        IF(itherm.eq.1) THEN
          WRITE(18) xi, vxi, axi
        ENDIF
        CLOSE(18)

!        OPEN(18,FILE="pulay_st",FORM="unformatted")
!        REWIND(18)
!        WRITE(18) AA_pulay_6
!        CLOSE(18)

      END IF

      DO iislda=1,islda

!      CALL pulayIO(AL,dw_6(1,1,iislda),mr_nL,npulay_max_6,
!     &           1,n1L,n2L,n3L,fdw(iislda))
!      CALL pulayIO(AL,dR_6(1,1,iislda),mr_nL,npulay_max_6,
!     &           1,n1L,n2L,n3L,fdR(iislda))

        filename=trim(adjustl(frho_out1(iislda)))//"."//trim(adjustl(fileindex))
        CALL rhoIO(AL,rho_storage_nL(1,iislda,1),mr_nL,1,n1L,n2L,n3L,filename)
        filename=trim(adjustl(frho_out2(iislda)))//"."//trim(adjustl(fileindex))
        CALL rhoIO(AL,rho_storage_nL(1,iislda,2),mr_nL,1,n1L,n2L,n3L,filename)
        filename=trim(adjustl(frho_out3(iislda)))//"."//trim(adjustl(fileindex))
        CALL rhoIO(AL,rho_storage_nL(1,iislda,3),mr_nL,1,n1L,n2L,n3L,filename)

!        do i = 203, 208
!  
!         rho_tmp = 0.d0
!         call d3fft_comp(ug_n_bp(1,i), workr_n, -1, 1)
!         rho_tmp(1:nr_n) = cdabs(workr_n(1:nr_n))**2
!         rho_tmp_nL = 0.d0
!         call convert_SLvr(rho_tmp, rho_tmp_nL, 1)
!  
!         filename="rho"
!         write(fileindex,'(i6)') i
!         filename=trim(adjustl(filename))//"."//trim(adjustl(fileindex))
!         write(fileindex,'(i6)') itime-1
!         filename=trim(adjustl(filename))//"."//trim(adjustl(fileindex))
!         call rhoIO(AL,rho_tmp_nL,mr_nL,1,n1L,n2L,n3L,filename)
!  
!        enddo

      END DO

      CALL MPI_barrier(MPI_COMM_WORLD,ierr)

      !!!!!  ONLY FOR K=1 !!!!!
      IF(islda.eq.1.and.nkpt.eq.1) THEN
        ug_n_bp(:,:)=ug_n_bp_0(:,:,1,1)
      END IF

      DO iislda=1,islda
          fwg(iislda)=trim(adjustl(fwg_st(iislda)))//"."//trim(adjustl(fileindex))
      END DO
      CALL write_wg_BP(fwg,AL,islda,nkpt)
      CALL MPI_barrier(MPI_COMM_WORLD,ierr)

      ! call write_wg_BP() will destory ug_n_bp on node 0
      ! so read it again from storage
      IF(islda.eq.1.and.nkpt.eq.1) THEN
        ug_n_bp(:,:)=ug_n_bp_0(:,:,1,1)
      END IF

      string = "output temp data"
      CALL timing_mpi(string,t_0)

!    END IF
    END IF

    ! Special: output specific psi charge density
    IF( iocc.gt.0 .and. mod((itime-1),10).eq.0 ) THEN

      filename = trim(adjustl(fpsi))//"."
      WRITE(fileindex, '(i)') jelec
      filename = trim(adjustl(filename))//trim(adjustl(fileindex))//"."
      WRITE(fileindex, '(i)') itime-1
      filename = trim(adjustl(filename))//trim(adjustl(fileindex))
      CALL rhoIO(AL,rho_tmp2_nL,mr_nL,1,n1L,n2L,n3L,filename)

      filename = trim(adjustl(fpsi))//"."
      WRITE(fileindex, '(i)') jhole
      filename = trim(adjustl(filename))//trim(adjustl(fileindex))//"."
      WRITE(fileindex, '(i)') itime-1
      filename = trim(adjustl(filename))//trim(adjustl(fileindex))
      CALL rhoIO(AL,rho_tmp3_nL,mr_nL,1,n1L,n2L,n3L,filename)

      CALL MPI_barrier(MPI_COMM_WORLD,ierr)

    END IF

9998  CONTINUE
    ! end of output

180 FORMAT(300(E15.8,1x))

    IF(idone.eq.-1) EXIT Loop_itime

    CALL MPI_barrier(MPI_COMM_WORLD,ierr)

    t_step1=mpi_wtime()

    IF(inode.eq.1) THEN

      WRITE(6,*)
      WRITE(6,*) "##################################################"
      WRITE(6,*)
      WRITE(17,*) "------------------------------"
      WRITE(17,FMT1703) " This step cost    =",t_step1-t_step0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  check if we still have time to do next time step
      IF(itime.gt.1) THEN

        step_num = step_num + 1
        IF(step_num.le.10) THEN
          t_loop(step_num) = t_step1 - t_step0
          t_loop_ave = sum(t_loop(1:step_num)) / step_num
        ELSE
          t_loop(1:9) = t_loop(2:10)
          t_loop(10) = t_step1 - t_step0
          t_loop_ave = sum(t_loop) / 10
        END IF
        WRITE(17,FMT1703) " Ave. step cost    =",t_loop_ave

        IF(t_timewall.gt.0.d0) THEN
          t_left = t_timewall - (t_step1 - t_start_petot)
          step_avail = floor(t_left / t_loop_ave)
          WRITE(17,'(A,I10)') " Expected to end at time step   :", step_avail + itime - 1
          IF(step_avail.le.4.and.itime.lt.(ntime_init+ntime-1)) THEN
          ! almost reach timewall, so just do one more step
            ntime = itime - ntime_init + 1
            write(17,*) "Do one more time step, then stop"
          END IF
        END IF

      END IF

      WRITE(17,*) "------------------------------"

    END IF

    call mpi_bcast(ntime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    t_step0 = t_step1

    IF(itime.eq.ntime_init+ntime) EXIT Loop_itime

  END DO Loop_itime


  IF(inode.eq.1) THEN
    WRITE(6,*) "test end"
    IF(idone.eq.-1) THEN
      WRITE(6,*) "Too many loops on time step",itime-1
      WRITE(6,*) "Temp data has been output for time step",itime-2
      WRITE(17,*) "Please restart from time step",itime-1
      WRITE(17,*) "or other savepoints"
    ENDIF
    WRITE(17,*) "****************************************"
    CLOSE(17)
  END IF

  DEALLOCATE(cc_pp0)
  DEALLOCATE(cc_pp)
  DEALLOCATE(cc_pp_tmp)
  IF(ifatom.eq.2) DEALLOCATE(dxatom_in)
  CALL data_deallocate_6()
  CALL data_deallocate_TDDFT()


9999 CONTINUE

  CALL MPI_barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_abort(MPI_COMM_WORLD,ierr)

  STOP

  RETURN




  CONTAINS

!cccccccccccccccccccccccccccccccccccc
!     TDDFT initialization
!cccccccccccccccccccccccccccccccccccc
SUBROUTINE initial_TDDFT()

  IF(inode.eq.1) THEN
    WRITE(6,*)
    WRITE(6,*) "| Initialization"
    WRITE(6,*) "------------------------------"
    WRITE(6,*)
    WRITE(17,*)
    WRITE(17,*) "Initialization"
    WRITE(17,*) "****************************************"
    WRITE(17,*) " see <report> FILE for details"
    CALL system_flush(17)
  END IF

  ! if ifatom=2, read dxatom
  IF(ifatom.eq.2) THEN

    IF(inode.eq.1) THEN

      PRINT *, "Checking if "//trim(adjustl(fdx_td))//" is available..."

      OPEN(107, file=trim(adjustl(fdx_td)), status='old', action='read', &
           iostat=ierr)
      IF(ierr.ne.0) THEN
        PRINT *, "--ERROR: File "//trim(adjustl(fdx_td))//" not found, stop"
        CALL MPI_abort(MPI_COMM_WORLD,ierr)
      ENDIF

      REWIND(107)
      DO itime = ntime_init + 1, ntime_init + ntime
        READ(107,'(A23,I)') string, itmp
        IF(itmp.ne.(itime-1)) THEN
          PRINT *, "--ERROR: Cannot read dxatom for istep=",itime-1
          CALL MPI_abort(MPI_COMM_WORLD,ierr)
        ENDIF
        DO i = 1, natom
          READ(107,'(i4,2x,3(f15.9,1x),2x,3(i2,1x))') &
               itmp, (dxatom_in(j, i, itime), j=1,3), itmp, itmp, itmp
        ENDDO
      ENDDO
      CLOSE(107)
      
      PRINT *, "OK!"

    ENDIF

    CALL mpi_bcast(dxatom_in,3*matom*ntime,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    xatom(:, :) = dxatom_in(:, :, ntime_init+1)

  ENDIF

  IF(ntime_init.eq.0) THEN
  ! Initialize everything from ground state calculation (self-cons Etotcalc.f)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! calc the semiconductor surface polarization effect
    if(inode.eq.1.and.ibshift.eq.1) call bandshift0(xatom, iatom, AL)
    !if(inode.eq.1) call bandshift0(xatom, iatom, AL)
    if(ibshift.eq.1) then
      call mpi_bcast(dP, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(gz_1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(gz_2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    endif

    IF(inode.eq.1) THEN
      t_0=mpi_wtime()
    END IF

    itime = ntime_init
    CALL calc_vext()

    string = "before Etotcalc"
    CALL timing_mpi(string,t_0)

    IF(ibo_md.le.0) THEN
      fatom = 0.d0
      iforce_cal = 0
      ireplace = 0
    ELSEIF(ibo_md.gt.0) THEN
      fatom = 0.d0
      iforce_cal = 1
      ireplace = 0
    END IF

    CALL Etotcalc(xatom, fatom, workr_n, Etot, iforce_cal, ido_rho, ido_vr, tolug, &
                  tolE, niter, nline, iCGmth, iscfmth, FermidE, &
                  itypeFermi, mCGbad, E_st, err_st, AL, nkpt, ntype, &
                  convergE, islda, igga, iwg_out, fwg_out, ivr_out, &
                  amx_mth, xgga, occ, ireplace, TS, E_dDrho)
    string = "Etotcalc"
    CALL timing_mpi(string,t_0)

    IF(inode.eq.1) THEN
      OPEN(18, FILE= trim(adjustl(focc))//".t0")
      REWIND(18)
      DO i=1,mst
        WRITE(18,*) E_st(i,1,1) * Hart, occ(i,1,1)
      END DO
      CLOSE(18)
    END IF

    occ0=occ
    IF(inode.eq.1) THEN
      OPEN(18, FILE=focc_st, FORM="unformatted")
      REWIND(18)
      WRITE(18) occ0
      CLOSE(18)
    END IF

    CALL calc_kappa()

    cc_pp0=cc_pp
    IF(islda.eq.1.and.nkpt.eq.1) THEN
      ug_n_bp_0(:,:,1,1)=ug_n_bp(:,:)
    ELSE
!        ug_n_bp_0=ug_all
    END IF
    E_st0=E_st
    rho0=rho_nL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELSE ! ntime_init.gt.0
  ! read everything from external files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(fileindex, '(i)') ntime_init-1

    IF(inode.eq.1) THEN

      t_0=mpi_wtime()
      filename = trim(adjustl(ftmp_st))//"."//trim(adjustl(fileindex))
      WRITE(6,*) "Read temp data from "//filename
      WRITE(17,*) "------------------------------"
      WRITE(17,'(a)') "  Read data from "//filename
      WRITE(17,*) "------------------------------"
      CALL system_flush(17)

      OPEN(18,FILE=filename,FORM='unformatted',STATUS='old',ACTION='READ',IOSTAT=ierr)
      REWIND(18)
      READ(18) itmp,rtmp
      READ(18) xatom
      !!!!!!!!!!!!!!!!!!!
      IF(ifatom.eq.2) xatom(:, :) = dxatom_in(:, :, ntime_init+1)
      !!!!!!!!!!!!!!!!!!!
      READ(18) Vi
      READ(18) cc_pp0
      READ(18) E_st0
      READ(18) Ealphat,TS0,E_dDrho0
      READ(18) rtmp,Box,InitTemp,DesiredTemp,TotalEn,Ekin
      IF(itherm.eq.1) THEN
        READ(18) xi, vxi, axi
      ENDIF
      CLOSE(18)

    END IF

    CALL mpi_bcast(xatom,3*matom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(Vi,3*matom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(cc_pp0,mst*mmn*nkpt*islda,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(E_st0,mst*nkpt*islda,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(Box,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(InitTemp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(DesiredTemp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(TotalEn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(Ekin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(Ealphat,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(TS0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL mpi_bcast(E_dDrho0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    IF(itherm.eq.1) THEN
      CALL mpi_bcast(xi,7,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL mpi_bcast(vxi,7,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL mpi_bcast(axi,7,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    ENDIF

    CALL unfmtIO_r(18,focc_st,occ0,mst*nkpt*islda,0,1)

    DO iislda=1,islda
      filename=trim(adjustl(frho_out1(iislda)))//"."//trim(adjustl(fileindex))
      CALL rhoIO(AL,rho_storage_nL(1,iislda,1),mr_nL,2,n1L,n2L,n3L,filename)
      filename=trim(adjustl(frho_out2(iislda)))//"."//trim(adjustl(fileindex))
      CALL rhoIO(AL,rho_storage_nL(1,iislda,2),mr_nL,2,n1L,n2L,n3L,filename)
      filename=trim(adjustl(frho_out3(iislda)))//"."//trim(adjustl(fileindex))
      CALL rhoIO(AL,rho_storage_nL(1,iislda,3),mr_nL,2,n1L,n2L,n3L,filename)
    END DO

    cc_pp = cc_pp0

    IF(islda.eq.1.and.nkpt.eq.1) THEN
      ug_n_bp_0(:,:,1,1)=ug_n_bp(:,:)
    ELSE
!        ug_n_bp_0=ug_all
    END IF

    string = "load data"
    CALL timing_mpi(string,t_0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END IF

  itime = ntime_init
  CALL calc_psi()

  IF(inode.eq.1) THEN
    WRITE(17,*) "****************************************"
    CALL system_flush(17)
  END IF

  RETURN

END SUBROUTINE initial_TDDFT
!cccccccccccccccccccccccccccccccccccccccccccc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE output_psi()

  ! This subroutine will output:
  ! 1. td state |psi(jhole)|**2 in real space
  ! 2. the whole wavefunction file ug_psi
  ! 3. the project wavefunction file beta_psi (for PDOS plot)
  ! then call mpi_abort
  ! see README for details

  IF(inode.eq.1) THEN
    WRITE(6,*) "iCGmth(1) = -1, only do output."
    WRITE(6,*) "  time step is", jelec
    WRITE(6,*) "  rho_psi is for state", jhole
  END IF

  WRITE(fileindex,'(i)') jelec
  filename = trim(adjustl(fcc_st))//"."//trim(adjustl(fileindex))
  CALL unfmtIO_c(18, filename, cc_pp0, mst*mmn*nkpt*islda, 0, 1)
  IF(inode.eq.1) WRITE(6,*) "  read "//trim(adjustl(filename))

  DO iislda = 1, islda

    rho_tmp = 0.d0

    DO kpt = 1, nkpt

      IF(nkpt.ne.1) THEN
        CALL gen_G_comp(kpt, 0)
        CALL fftprep_comp(n1, n2, n3)
        ng_n = ngtotnod(inode, kpt)
      END IF

      CALL ugIOBP(ug_n_bp, kpt, 2, 0, iislda, -1, nkpt, islda)

      cpsi_td(:, :, kpt, iislda) = zero
      CALL zgemm('N', 'N', mg_nx, mmn, mst, one, ug_n_bp, mg_nx, &
                 cc_pp0(1,1,kpt,iislda), mst, one, &
                 cpsi_td(1,1,kpt,iislda), mg_nx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ONLY FOR nkpt = 1 AND islda = 1!!!!!!!!!
      IF(nkpt.eq.1 .and. islda.eq.1) THEN

        ug_n_bp_0(:,:,kpt,iislda) = ug_n_bp(:,:)
        ug_n_bp(:,1:mmn) = cpsi_td(:,1:mmn,kpt,iislda)

        WRITE(fileindex, '(i)') kpt
        filename = trim(adjustl(fpsi_ug))//"."//trim(adjustl(fileindex))//"."
        WRITE(fileindex, '(i)') iislda
        filename = trim(adjustl(filename))//trim(adjustl(fileindex))//"."
        WRITE(fileindex, '(i)') jelec
        filename = trim(adjustl(filename))//trim(adjustl(fileindex))

        CALL write_wg_BP(filename,AL,islda,nkpt)
        CALL MPI_barrier(MPI_COMM_WORLD,ierr)

        ! when k=1 and islda=1, write_wg_BP will destory ug_n_bp.
        ! value it again
        ug_n_bp(:,:) = ug_n_bp_0(:,:,kpt,iislda)
        ug_n_bp(:,1:mmn) = cpsi_td(:,1:mmn,kpt,iislda)

      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL d3fft_comp_block(cpsi_td(1,jhole,kpt,iislda), &
                            workr_n, -1, kpt, 1)

      rho_tmp(1:nr_n) = rho_tmp(1:nr_n) + abs(workr_n(1:nr_n))**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filename=trim(adjustl(focc))//"."
      WRITE(fileindex,'(i)') kpt
      filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
      WRITE(fileindex,'(i)') iislda
      filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
      WRITE(fileindex,'(i)') jelec
      filename=trim(adjustl(filename))//trim(adjustl(fileindex))

      IF(inode.eq.1) THEN
        OPEN(18,FILE=filename,STATUS='old',ACTION='read')
        REWIND(18)
        DO i=1,mst
          IF(i.le.mmn) THEN
            READ(18,'(3(E15.8,1x))') rtmp,rtmp,E_st(i,kpt,iislda)
          ELSE
            READ(18,'(3(E15.8,1x))') E_st(i,kpt,iislda),rtmp
          END IF
        END DO
        CLOSE(18)
        WRITE(6,*) "  read "//trim(adjustl(filename))
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END DO ! DO kpt

    CALL convert_SLvr(rho_tmp, rho_tmp_nL, 1)

    filename = trim(adjustl(fpsi))//"."
    WRITE(fileindex, '(i)') jhole
    filename = trim(adjustl(filename))//trim(adjustl(fileindex))//"."
    WRITE(fileindex, '(i)') iislda
    filename = trim(adjustl(filename))//trim(adjustl(fileindex))//"."
    WRITE(fileindex, '(i)') jelec
    filename = trim(adjustl(filename))//trim(adjustl(fileindex))

    CALL rhoIO(AL,rho_tmp_nL,mr_nL,1,n1L,n2L,n3L,filename)

    CALL MPI_barrier(MPI_COMM_WORLD,ierr)

   END DO ! DO iislda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(inode.eq.1) THEN
     OPEN(18, FILE="eigen_all.store", FORM="unformatted")
     REWIND(18)
     WRITE(18) islda,nkpt,mst,nref_tot,natom,nnodes
     DO iislda=1,islda
     DO kpt=1,nkpt
       WRITE(18) iislda,kpt,weighkpt(kpt),akx(kpt),aky(kpt),akz(kpt)
       WRITE(18) (E_st(i,kpt,iislda),i=1,mst)
     END DO
     END DO
     CLOSE(18)
   END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! call Etotcalc to do the projection
   fatom = 0.d0
   iforce_cal = 0
   ireplace = 0
   CALL Etotcalc(xatom, fatom, workr_n, Etot, iforce_cal, ido_rho, &
                 ido_vr, tolug, tolE, niter, nline, iCGmth, &
                 iscfmth, FermidE, itypeFermi, mCGbad, E_st, &
                 err_st, AL, nkpt, ntype, convergE, islda, igga, &
                 iwg_out, fwg_out, ivr_out, amx_mth, xgga, occ, &
                 ireplace, TS, E_dDrho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL MPI_barrier(MPI_COMM_WORLD,ierr)
   CALL MPI_abort(MPI_COMM_WORLD,ierr)

   STOP

END SUBROUTINE output_psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!cccccccccccccccccccccccccccccccccccccccccccc
!     Set up the time-evolved external electric field
!cccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE calc_vext()

  REAL(8) :: amp_ext, ome_ext, wid_ext, t0_ext, t1_ext, gau_ext, fav_ext, &
             d, E_ext_dt

  IF(ivext_dt.le.0) THEN

    vext_dt_n=0.d0
    vext_nL=0.d0
    ivext_in=0

  ELSEIF(ivext_dt.gt.0) THEN

    ivext_in=1
    IF(itime.eq.0) itime=1

    ! Change this part to set up your own potential

    amp_ext=1.d-3  ! Hart/Bohr
    ome_ext=2.d0*3.14159265d0/0.46625d0   ! rad/fs 
    wid_ext=2.d0   ! fs
    t0_ext=0.d0      ! fs
    t1_ext=2.d0    ! fs

!      gau_ext=amp_ext*exp(-(dble(itime-1)*dt_step-t0_ext)**2/
!     &                     2/wid_ext/wid_ext)
!     &        /dsqrt(2*3.14159265d0)/wid_ext

    IF((dble(itime-1)*dt_step).le.t0_ext) THEN
      gau_ext=0.d0
    ELSEIF((dble(itime-1)*dt_step).ge.t1_ext) THEN
      gau_ext=amp_ext
    ELSE
      gau_ext=dsin(3.14159265d0*(dble(itime-1)*dt_step-t0_ext)/ &
                   2.d0/(t1_ext-t0_ext))*amp_ext
    END IF

    fav_ext=dsin(ome_ext*dble(itime-1)*dt_step)*gau_ext

    IF(inode.eq.1) THEN
      WRITE(6,FMT1703) "  Real time        =",(itime-1)*dt_step
      WRITE(6,FMT1704) "  Gaussian envlp   =",gau_ext
      WRITE(6,FMT1704) "  E_ext            =",fav_ext
    END IF

    DO ii=1,nr_n

      jj=ii+(inode-1)*nr_n
      ix=(jj-1)/(n2*n3)+1
      iy=(jj-1-(ix-1)*n2*n3)/n3+1
      iz=jj-(ix-1)*n2*n3-(iy-1)*n3

      d=sqrt(dble(ix-1.d0-n1/2.d0)**2+dble(iy-1.d0-n2/2.d0)**2)
      E_ext_dt=1.d0-d*24.d0/(12.d0-3.14159265d0)/dble(n1)

      IF(E_ext_dt.lt.(-3.14159265d0/(12.d0-3.14159265d0))) THEN
        E_ext_dt=-3.14159265d0/(12.d0-3.14159265d0)
      ELSEIF(E_ext_dt.gt.1.d0) THEN
        E_ext_dt=1.d0
      END IF

      vext_dt_n(ii)=E_ext_dt*fav_ext

    END DO
!      vext_dt_n=fav_ext
    CALL convert_SLvr(vext_dt_n,vext_nL,1)

  END IF

  RETURN

END SUBROUTINE calc_vext

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!     At t=0, in real space, we could add an
!     exp[i * Kappa * X] term on every wavefunction,
!     which will act as an all-spectrum excitation.
!     It provides one possible way to calculate the
!     absorption spectrum.
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE calc_kappa()

  IF(ikappa.le.0) THEN

    DO iislda=1,islda
    DO kpt=1,nkpt
    DO imn=1,mmn
    DO m=1,mst2
      IF(m.eq.imn) cc_pp(m,imn,kpt,iislda)=one
      IF(m.ne.imn) cc_pp(m,imn,kpt,iislda)=zero
    END DO
    END DO
    END DO
    END DO

  ELSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO iislda=1,islda
    DO kpt=1,nkpt

      IF(nkpt.ne.1) THEN
        CALL gen_G_comp(kpt,0)
        CALL fftprep_comp(n1,n2,n3)
        ng_n=ngtotnod(inode,kpt)
      END IF
      CALL ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)

      workr=zero
      CALL d3fft_comp_block(ug_n_bp,workr,-1,1,nblock_band_mx)

      DO m=1,mmn

        workr_n(:)=workr(:,m)

        DO ii=1,nr_n
          jj=ii+(inode-1)*nr_n
          ix=(jj-1)/(n2*n3)+1
          iy=(jj-1-(ix-1)*n2*n3)/n3+1
          iz=jj-(ix-1)*n2*n3-(iy-1)*n3

          workr_n(ii) = workr_n(ii) * exp(cai*rkappa(1)*ix) &
                                    * exp(cai*rkappa(2)*iy) &
                                    * exp(cai*rkappa(3)*iz)

        END DO

        sum0=0.d0

        DO l=1,mst2

          cc=zero
          DO i=1,mr_n
            cc=cc+workr_n(i)*conjg(workr(i,l))
          END DO

          CALL global_sumc(cc)
          cc=cc*vol/(nr_nL*nnodes)
          cc_pp(l,m,kpt,iislda)=cc

          sum0=sum0+abs(cc_pp(l,m,kpt,iislda))**2

        END DO

        sum0=1.d0/dsqrt(sum0)

        IF(inode.eq.1) WRITE(6,*) "sum0=",sum0
        IF(inormal.eq.1) cc_pp(:,m,kpt,iislda)=cc_pp(:,m,kpt,iislda)*sum0

      END DO

    END DO
    END DO

    !CALL MPI_bcast(cc_pp,mst*mmn*nkpt*islda,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

  END IF

  RETURN

END SUBROUTINE calc_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!ccccccccccccccccccccccccccccccccccccccccccccccccc
!     When the integral inside [t1,t2] ends, we
!     will have c(i,j)s on phi(t1) basis set.
!     Then we need to change these coefficients
!     to phi(t2) basis set.
!ccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE cc_change_basis()

  DO iislda=1,islda
  DO kpt=1,nkpt

    IF(nkpt.ne.1) THEN
      CALL gen_G_comp(kpt,0)
      CALL fftprep_comp(n1,n2,n3)
      ng_n=ngtotnod(inode,kpt)
    END IF
    CALL ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)

    cdot_tmp=zero
    CALL dot_product_BP(ug_n_bp_0(1,1,kpt,iislda),ug_n_bp(1,1), &
                        ng_n,cdot_tmp)
    cdot_tmp=cdot_tmp*vol

    string = "<phi_t1|phi_t2>"
    CALL timing_mpi(string,t_0)

    call ortho_Gram(cdot_tmp, mst, mst)
    IF(mst2.lt.mst) cdot_tmp(mst2+1:mst,mst2+1:mst) = zero

    DO m2=1,mst2

      !cmax=zero
      !DO m1=1,mst2
      !  IF(abs(cdot_tmp(m1,m2)).gt.abs(cmax)) THEN
      !    cmax=cdot_tmp(m1,m2)
      !  END IF
      !END DO

      cmax = cdot_tmp( maxloc( abs(cdot_tmp(1:mst2, m2)), 1 ) , m2)
      cmax=cmax/abs(cmax)
      cdot_tmp(:,m2)=conjg(cmax)*cdot_tmp(:,m2)
      ug_n_bp(:,m2)=conjg(cmax)*ug_n_bp(:,m2)

    END DO

    cdot(:,:,kpt,iislda)=cdot_tmp(:,:)

    DO i=1,mst
      DO j=1,mst
        cdot_tmp(j,i)=conjg(cdot(i,j,kpt,iislda))
      END DO
    END DO

    cc_pp(:,:,kpt,iislda)=zero
    CALL zgemm('N','N',mst,mmn,mst,one, &
               cdot_tmp,mst,cc_pp0(1,1,kpt,iislda),mst,one, &
               cc_pp(1,1,kpt,iislda),mst)


    IF(inormal.eq.1) THEN

      cc_pp_tmp(:,1:mmn)=cc_pp(:,1:mmn,kpt,iislda)

      IF(mmn0.gt.1) THEN
        cc_pp_tmp(:,1:mmn0-1) = zero
        DO i = 1,mmn0-1
          cc_pp_tmp(i,i) = one
        ENDDO
      END IF

      call ortho_Gram(cc_pp_tmp, mst, mmn)
      IF(mst2.lt.mst) cc_pp_tmp(mst2+1:mst,1:mmn) = zero

      cc_pp(:,:,kpt,iislda)=cc_pp_tmp(:,:)

    END IF

    string = "orthonormal"
    CALL timing_mpi(string,t_0)

  END DO
  END DO

  RETURN

END SUBROUTINE cc_change_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!ccccccccccccccccccccccccccccccccccccccccccccccc
!     Calculate psi(t2) and rho(t2) from
!     c(i,j) and phi(t2)
!ccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE calc_psi()

  DO iislda=1,islda

    rho_tmp=0.d0
    rho_tmp2=0.d0
    rho_tmp3=0.d0

    DO kpt=1,nkpt

      IF(nkpt.ne.1) THEN
        CALL gen_G_comp(kpt,0)
        CALL fftprep_comp(n1,n2,n3)
        ng_n=ngtotnod(inode,kpt)
      END IF
      CALL ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)

      dos(:,kpt,iislda) = 0.d0
      DO imn=1,mmn
        EE=0.d0
        DO m=1,mst2
          EE=EE+abs(cc_pp(m,imn,kpt,iislda))**2*E_st(m,kpt,iislda)
          dos(m,kpt,iislda)=dos(m,kpt,iislda)+ &
                            abs(cc_pp(m,imn,kpt,iislda))**2*occ0(imn,kpt,iislda)
        END DO
        E_td(imn,kpt,iislda)=EE
      END DO

      cpsi_td(:,:,kpt,iislda)=zero
      CALL zgemm('N','N',mg_nx,mmn,mst,one, &
                 ug_n_bp,mg_nx,cc_pp(1,1,kpt,iislda),mst,one, &
                 cpsi_td(1,1,kpt,iislda),mg_nx)

      workr=zero
      CALL d3fft_comp_block(cpsi_td(1,1,kpt,iislda),workr,-1,kpt,mmn)

      DO m=1,mmn
        rho_tmp(1:nr_n)=rho_tmp(1:nr_n)+occ0(m,kpt,iislda)*abs(workr(1:nr_n,m))**2
      END DO

      ! Calculate specific psi charge density
      IF( iocc.gt.0 .and. (itime.eq.0 .or. mod((itime-1),10).eq.0) ) THEN
        rho_tmp2(1:nr_n)=rho_tmp2(1:nr_n)+abs(workr(1:nr_n,jelec))**2
        rho_tmp3(1:nr_n)=rho_tmp3(1:nr_n)+abs(workr(1:nr_n,jhole))**2
      ENDIF

    END DO ! nkpt

    CALL convert_SLvr(rho_tmp,rho_tmp_nL,1)
    IF( iocc.gt.0 .and. (itime.eq.0 .or. mod((itime-1),10).eq.0) ) THEN
      CALL convert_SLvr(rho_tmp2,rho_tmp2_nL,1)
      CALL convert_SLvr(rho_tmp3,rho_tmp3_nL,1)
    ENDIF

    sum0=0.d0
    DO i=1,nr_nL
      sum0=sum0+rho_tmp_nL(i)
    END DO

    CALL global_sumr(sum0)

    sum0=sum0*vol/(n1*n2*n3)
    sum0=totNel/sum0

    IF(inode.eq.1) THEN
      WRITE(6,*) "rho factor:",sum0
    END IF

    rho_tmp_nL=rho_tmp_nL*sum0
    rho_fac = sum0

    rho1(:,iislda)=rho_tmp_nL(:)

  END DO ! islda

  string = "psi_td calc"
  CALL timing_mpi(string,t_0)

  RETURN

END SUBROUTINE calc_psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Calculate the difference between input
!     rho1 and rho0, via
!     \sum_i abs(rho1(i) - rho0(i)) / \sum_i rho1(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE diff_rho()

  IF(inode.eq.1) THEN
    WRITE(6,*) "****************************************"
    WRITE(6,*) "CHECK (rho_out-rho_in)**2"
  END IF

  sum4=0.d0

  DO iislda=1,islda

    sum0=0.d0
    sum2=0.d0
    sum3=0.d0

    DO i=1,nr_nL
      sum0=sum0+abs(rho1(i,iislda)-rho0(i,iislda))
      sum2=sum2+rho1(i,iislda)
      sum3=sum3+rho0(i,iislda)
    END DO
    CALL global_sumr(sum0)
    CALL global_sumr(sum2)
    CALL global_sumr(sum3)

    IF(inode.eq.1) THEN
      WRITE(6,*) sum0,sum0/sum2,sum2,sum3
      WRITE(6,*) sum0*vol/(n1*n2*n3),sum2*vol/(n1*n2*n3),sum3*vol/(n1*n2*n3)
      WRITE(6,*)
    END IF

    sum4(iislda)=sum0*vol/(n1*n2*n3)

  END DO

  RETURN

END SUBROUTINE diff_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION check_converge()

        itmp = 0
        Do iislda = 1, islda
          If(sum4(iislda).lt.tolrho) itmp = itmp + 1
        Enddo

        If(itmp.eq.islda) Then
        ! Converged

          ! Timing
          CALL MPI_barrier(MPI_COMM_WORLD,ierr)
          IF(inode.eq.1) THEN
            t_loop1=mpi_wtime()
            WRITE(6,*) "Converged!"
            WRITE(17,'(a,i5,a,f17.10,a,e17.10,a,f10.3)') &
                       "  ", iiscf, &
                       ",  rho_fac =", rho_fac, &
                       ",  rho_diff =", (sum4(j), j=1,islda), &
                       ",  t = ", t_loop1-t_loop0
            WRITE(17,*) "------------------------------"
            WRITE(17,*) " Converged"
            t_loop0 = t_loop1
            CALL system_flush(17)
          END IF

          rho0 = rho1
          CALL recalc()

          check_converge = 1

        ELSEIF(iiscf.ge.mscf) THEN
!        ELSEIF(iiscf.ge.mscf .or. itime.eq.235) THEN
        ! iiscf reaches mscf

          IF(inode.eq.1) THEN
            t_loop1=mpi_wtime()
            WRITE(6,*) "LOOP NUM >",mscf
            WRITE(17,'(a,i5,a,f17.10,a,e17.10,a,f10.3)') &
                       "  ", iiscf, &
                       ",  rho_fac =", rho_fac, &
                       ",  rho_diff =", (sum4(j), j=1,islda), &
                       ",  t = ", t_loop1-t_loop0
            WRITE(17,*) "------------------------------"
            WRITE(17,*) " LOOP NUM >",mscf
            CALL system_flush(17)
          END IF

          check_converge = -1

        ELSE
        ! Not converged, do charge mixing then ontinue from do Loop_iiscf

          CALL mch_pulay_wrap(rho0, rho1, iiscf, islda)

          DO iislda = 1, islda
            CALL mch_kerk(rho0(1,iislda), rho1(1,iislda), amx_mth(5), &
                          0.8d0, 0.25d0, 0.1d0)
          END DO

          CALL MPI_barrier(MPI_COMM_WORLD,ierr)
          IF(inode.eq.1) THEN
            t_loop1 = mpi_wtime()
            WRITE(17,'(a,i5,a,f17.10,a,e17.10,a,f10.3)') &
                       "  ", iiscf, &
                       ",  rho_fac =", rho_fac, &
                       ",  rho_diff =", (sum4(j), j=1,islda), &
                       ",  t = ", t_loop1-t_loop0
            CALL system_flush(17)
          END IF

          check_converge = 0

        END IF ! If converge

END FUNCTION check_converge


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     For electron excitations, intrinsic Ehrenfest dynamics 
!     always under-estimates cooling down speeds.
!     Add Boltzmann factor to correct the speed.
!     Still in testing
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE calc_boltz()

  IF(inode.eq.1) THEN
    WRITE(6,*)
    WRITE(6,*) "------------------------------"
    WRITE(6,*) "| Add Boltzmann factor, Temp(K)=",temperature
    WRITE(6,*) "------------------------------"
    WRITE(6,*)
    t_0=mpi_wtime()
  END IF

  kT = temperature * 8.6173324d-5 / Hart

  cdot_tmp=zero
  S_boltz=one
  S_boltz_f=one

  DO k=1,mst2
  DO l=1,mst2

    IF(l.ne.k) THEN

      x_ltok=E1(l,kpt,iislda)-E1(k,kpt,iislda)
      rtmp=x_ltok*(exp((abs(H_T(l,k,kpt,iislda))-abs(x_ltok))/kT)+1.d0)
      S_boltz(l,k)=H_T(l,k,kpt,iislda)/rtmp

!      IF(iiscf.eq.1.or.itime.le.(ntime_init+3)) THEN

        DO i=1,mmn
          cdot_tmp(l,k)=cdot_tmp(l,k)+ &
                        conjg(cc_pp0(l,i,kpt,iislda))*cc_pp0(k,i,kpt,iislda)* &
                        occ0(i,kpt,iislda)*0.5d0
        END DO

        R_ltok_1(l,k)=cdot_tmp(l,k)*H_T(l,k,kpt,iislda)

!      END IF

      IF(real(R_ltok_1(l,k)).gt.0.d0) THEN
        r_factor(l,k)=exp(-abs(x_ltok/kT))
      ELSE
        r_factor(l,k)=1.d0
      END IF

      S_boltz_f(l,k)=S_boltz(l,k)*r_factor(l,k)

      IF((imag(R_ltok_1(l,k))*x_ltok).gt.0.d0) THEN
        H_T(l,k,kpt,iislda)=exp(-abs(x_ltok/kT))*H_T(l,k,kpt,iislda)
      END IF

    END IF

    S_boltz_st(l,k)=S_boltz(l,k)
    U_boltz_st(l,k)=conjg(cdot(k,l,kpt,iislda))

  END DO
  END DO

  CALL zgetrf(mst,mst,S_boltz,mst,ipiv,info)
  IF(info.ne.0) THEN
    WRITE(6,*) "Matrix inverse ERROR1, node:",inode-1
    CALL MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  END IF

  CALL zgetri(mst,S_boltz,mst,ipiv,work,10*mst,info)
  IF(info.ne.0) THEN
    WRITE(6,*) "Matrix inverse ERROR2, node:",inode-1
    CALL MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  END IF

  cdot_tmp=zero
  CALL zgemm('N','N',mst,mst,mst,one,S_boltz,mst, &
             S_boltz_f,mst,one,cdot_tmp,mst)
  S_boltz=cdot_tmp

  cdot_tmp=zero
  CALL zgemm('C','N',mst,mst,mst,one,S_boltz,mst, &
             cdot(1,1,kpt,iislda),mst,one,cdot_tmp,mst)

!  IF(inode.eq.1) THEN
!    t_1=mpi_wtime()
!    WRITE(6,*) "| Calc Factor,",t_1-t_0
!    t_0=mpi_wtime()
!  END IF

!cccccccccccccccccccccccccccccccccccccccccccccc
!      force orthonormal
  call ortho_Gram(cdot_tmp, mst, mst)
  IF(mst2.lt.mst) cdot_tmp(mst2+1:mst,mst2+1:mst) = zero
!cccccccccccccccccccccccccccccccccccccccccccccccc

  cdot(:,:,kpt,iislda)=cdot_tmp(:,:)

!  IF(inode.eq.1) THEN
!    t_1=mpi_wtime()
!    WRITE(6,*) "| Orthonormal 2,",t_1-t_0
!    WRITE(6,*) "------------------------------"
!    t_0=mpi_wtime()
!  END IF

END SUBROUTINE calc_boltz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Because of the finite accuracy of self-consistent
!     leapfrog, the rho(t2)->phi(t2)->psi(t2)->rho(t2)
!     loop could not be absolutely correct. We must decide
!     which one should be more precise, the rho or
!     the psi. Here we choose psi.
!     Here we use psi to generate rho, THEN use rho to
!     generate phi', THEN expand psi onto phi' to get
!     c(i,j)'. Finally we keep rho, phi' and c'.
!     Also in this subroutine the E_tot and force are
!     calculated.
!!!!! ^^^I G N O R E  T H O S E  C O N T E N T^^^ !!!!!
!     MODIFIED Nov. 2015, now we choose rho rather than psi,
!     so no further operations need to do.
!     Now this subroutine is only to calculate Etot and atom force.
!     Maybe we need to rename it.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE recalc()

!  implicit double precision (a-h,o-z)
  IMPLICIT NONE

  REAL(8) :: ewald, E_Hxc, E_coul, E_ion, E_extV, E_IVext, E_psiV, &
             E_rhoVext, E_cor, E_NSC, E_tmp

!      real*8 rho_ttmp(mr_nL,islda)
!      complex*16 ugbak(mg_nx,mst)

  IF(inode.eq.1) THEN
    WRITE(6,*) 
    WRITE(6,*) "------------------------------"
    WRITE(6,*) "| Recalc Adiabatic States, Etot & Force"
    WRITE(6,*) "------------------------------"
    WRITE(6,*)
  END IF

!ccccccccc
! MODIFIED Nov. 2015
  GOTO 233
!ccccccccc

  cc_pp0=cc_pp
  E_st0=E_st

  IF(islda.eq.1.and.nkpt.eq.1) THEN
    ug_n_bp_0(:,:,1,1)=ug_n_bp(:,:)
  ELSE
!       ug_n_bp_0=ug_all
  END IF

  string = "before Etotcalc"
  CALL timing_mpi(string,t_0)
  t_loop0=t_0

  iforce_cal=0
  ireplace=0
  rho_nL=rho0
  nnn=1
  IF(itime.eq.1.and.iocc.gt.0) nnn=2
  CALL Etotcalc(xatom, fatom, workr_n, Etot, iforce_cal, 0, 1, tolug, &
                tolE, nnn, nline, iCGmth_td, iscfmth_td, FermidE_td, &
                itypeFermi_td, mCGbad, E_st, err_st, AL, nkpt, ntype, &
                convergE, islda, igga, iwg_out, fwg_out, ivr_out, &
                amx_mth, xgga, occ, ireplace, TS, E_dDrho)
  string = "Etotcalc"
  CALL timing_mpi(string,t_0)

  CALL cc_change_basis()
  CALL calc_psi()
  CALL diff_rho()

!      CALL MPI_barrier(MPI_COMM_WORLD,ierr)
!      CALL MPI_abort(MPI_COMM_WORLD,ierr)
!      stop

233   CONTINUE

  IF(islda.eq.1.and.nkpt.eq.1) THEN
    ug_n_bp_0(:,:,1,1)=ug_n_bp(:,:)
  ELSE
!       ug_n_bp_0=ug_all
  END IF

  E_st0=E_st
  cc_pp0=cc_pp
  TS0=TS
  E_dDrho0=E_dDrho

  rho_storage_nL(:,:,1)=rho_storage_nL(:,:,2)
  rho_storage_nL(:,:,2)=rho_storage_nL(:,:,3)
  rho_storage_nL(:,:,3)=rho0(:,:)

!      IF(itime.ge.3) THEN
!       rho_storage_nL(:,:,4)=
!     &      3*rho_storage_nL(:,:,3)-
!     &      3*rho_storage_nL(:,:,2)+
!     &        rho_storage_nL(:,:,1)
!      ELSE
!       rho_storage_nL(:,:,4)=
!     &       rho_storage_nL(:,:,3)
!      END IF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      rho_ttmp=rho_nL
!      ugbak=ug_n_bp

  fatom=0.d0
  ewald=0.d0
  vvion=0.d0
  vvionT=0.d0
  rho_n=0.d0
  rrhocr=0.d0
  vvr_nL=0.d0
  vvr_n=0.d0
  vv0=0.d0
  E_Hxc=0.d0
  E_coul=0.d0
  E_ion=0.d0
  E_extV=0.d0
  E_IVext=0.d0

  CALL convert_SLvr(rho_n,rho_nL,-1)

  IF(icoul.eq.0) THEN
    CALL getewald(fatom,xatom,AL,ityatom,ewald)
  END IF

  IF(ivext_in.eq.1.or.ivext_dt.gt.0) THEN
    CALL getEextV(fatom,xatom,AL,ityatom,E_IVext)
  END IF

  CALL getVrhoL(AL,vvion,vvionT,xatom,ntype,iatom,rrhocr,totNel,0,islda)

  IF(ivext_in.eq.1.or.ivext_dt.gt.0) THEN
    vvion=vvion+vext_nL
  END IF

  IF(islda.eq.1.and.igga.eq.0) THEN
    CALL getpot2L(rho_nL,vvion,rrhocr,vvr_nL,vv0,E_Hxc,E_coul,E_ion)
  ELSEIF(islda.eq.2.and.igga.eq.0) THEN
    CALL getpot3L(rho_nL,vvion,rrhocr,vvr_nL,vv0,E_Hxc,E_coul,E_ion)
  ELSEIF(islda.eq.1.and.igga.eq.1) THEN
    CALL getpot4L(rho_nL,vvion,rrhocr,vvr_nL,vv0,E_Hxc,E_coul,E_ion,xgga)
  ELSEIF(islda.eq.2.and.igga.eq.1) THEN
    CALL getpot5L(rho_nL,vvion,rrhocr,vvr_nL,vv0,E_Hxc,E_coul,E_ion,xgga)
  END IF

  DO iislda=1,islda
    CALL convert_SLvr(vvr_n(1,iislda),vvr_nL(1,iislda),-1)
  END DO

  DO iislda=1,islda
    CALL mpi_bcast(vvr_n(1,iislda),nr_n,MPI_REAL8,0,MPI_COMM_N,ierr)
  END DO

  E_psiV=0.d0
  DO iislda=1,islda
    DO i=1,nr_n
      E_psiV=E_psiV+rho_n(i,iislda)*vvr_n(i,iislda)   ! rho_n is not symmetrize, but okay
    END DO
  END DO
  CALL global_sumr(E_psiV)

  E_rhoVext=0.d0
  DO iislda=1,islda
    DO i=1,nr_nL
     E_rhoVext=E_rhoVext+rho_nL(i,iislda)*vext_nL(i)
    END DO
  END DO
  CALL global_sumr(E_rhoVext)

  E_psiV=E_psiV*vol/(nr-0.d0)
  E_rhoVext=E_rhoVext*vol/(nrL-0.d0)

  E_cor=E_ion-E_rhoVext-E_psiV-E_dDrho0      ! E_ion=Vion*rho, includes E_rhoVext 

  E_extV=E_IVext+E_rhoVext    ! sum0 over ionic part and electronic part

!      rho_nL=rho_ttmp

  vvr_nL=0.d0
  vvionT=0.d0

  IF(islda.eq.1.and.igga.eq.0) THEN

    DO i=1,nr_nL
      vvr_nL(i,1)=rho_nL(i,1)
      vvionT(i)=UxcCA(rho_nL(i,1)+rrhocr(i),uxc2)
    END DO

  ELSEIF(islda.eq.2.and.igga.eq.0) THEN

    DO i=1,nr_nL
      vvr_nL(i,1)=rho_nL(i,1)+rho_nL(i,2)
      CALL UxcCA2(rho_nL(i,1)+rrhocr(i)*0.5d0, &
                  rho_nL(i,2)+rrhocr(i)*0.5d0, &
                  vxc1,vxc2,uxc1,uxc2)
      vvionT(i)=(vxc1+vxc2)/2.d0
    END DO

  ELSEIF(islda.eq.1.and.igga.eq.1) THEN

    rtmp = 0.d0
    DO i = 1, nr_nL
      vvr_nL(i,1) = rho_nL(i,1)
      rtmp = rtmp + dabs(rrhocr(i))
    ENDDO
    CALL global_sumr(rtmp)
    rtmp = rtmp * vol / nrL
    IF(rtmp.gt.1.D-5) THEN
      call getpot4_force(rho_nL,vvionT,rrhocr,xgga)
    ELSE
      vvionT=0.d0
    ENDIF

  ELSEIF(islda.eq.2.and.igga.eq.1) THEN

    rtmp = 0.d0
    DO i = 1, nr_nL
      vvr_nL(i,1) = rho_nL(i,1) + rho_nl(i,2)
      rtmp = rtmp + dabs(rrhocr(i))
    ENDDO
    CALL global_sumr(rtmp)
    rtmp = rtmp * vol / nrL
    IF(rtmp.gt.1.D-5) THEN
      call getpot5_force(rho_nL,vvionT,rrhocr,xgga)
    ELSE
      vvionT=0.d0
    ENDIF

  END IF

  IF(ifatom.eq.1) THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(icoul.eq.0) THEN
      CALL forcLC(AL,vvr_nL(1,1),vvionT,xatom,ntype,iatom,fatom)
    END IF

    IF(ilocal.eq.3) THEN

      ALLOCATE(fatom_tmp(3,natom))
      ALLOCATE(fatom_tmp2(3,natom))

      fatom_tmp=0.d0
      fatom_tmp2=0.d0

      DO iislda=1,islda
      DO kpt=1,nkpt

        IF(nkpt.ne.1) THEN
          CALL gen_G_comp(kpt,0)
          CALL fftprep_comp(n1,n2,n3)
          ng_n=ngtotnod(inode,kpt)
        END IF

        IF((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and. &
           (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) THEN

          IF(inode_tot.eq.1) WRITE(6,*) "Read TD ug"
          DO i=1,nblock_band_mx
            IF(i.gt.mmn) THEN
              ug_n_bp(:,i)=zero
            ELSE
              ug_n_bp(:,i)=cpsi_td(:,i,kpt,iislda)
            END IF
          END DO
          CALL forcNLq(fatom_tmp,occ0,kpt,nkpt,iislda,islda,E_td)

        END IF

      END DO
      END DO

      CALL mpi_allreduce(fatom_tmp,fatom_tmp2,3*natom, &
                         MPI_REAL8,MPI_SUM,MPI_COMM_K2,ierr)
      DO i=1,natom
        fatom(1,i)=fatom(1,i)+fatom_tmp2(1,i)
        fatom(2,i)=fatom(2,i)+fatom_tmp2(2,i)
        fatom(3,i)=fatom(3,i)+fatom_tmp2(3,i)
      END DO

      DEALLOCATE(fatom_tmp)
      DEALLOCATE(fatom_tmp2)

    END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DEBUG: Calc psi_td projection on atomic orbitals
!  IF( iocc.gt.0 .and. mod((itime-1),50).eq.0 ) THEN
!    call Proj1_comp( ilocal, nline, tolug, E_st(1,1,1), &
!                     err_st(1,1,1), 0.d0, vvr_n(1,1), &
!                     workr_n, 0, 1.d0, 1, 1)
!    filename='beta_psi.'
!    WRITE(fileindex,'(i)') itime-1
!    filename=trim(adjustl(filename))//trim(adjustl(fileindex))//"."
!    call beta_psiIO_TD( beta_psi, 1, 1, 0, 1, filename)
!  ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!  ONLY FOR K=1 !!!!!
  IF(islda.eq.1.and.nkpt.eq.1) THEN
    ug_n_bp(:,:)=ug_n_bp_0(:,:,1,1)
  END IF

  rho_nL=rho0

  E_NSC=0.d0

  DO iislda=1,islda
    DO kpt=1,nkpt
      DO i=1,mmn

        E_tmp=0.d0
        DO j=1,mst2
          E_tmp=E_tmp+abs(cc_pp0(j,i,kpt,iislda))**2* &
                E_st0(j,kpt,iislda)
        END DO

        E_NSC=E_NSC+E_tmp*occ0(i,1,iislda)

      END DO
    END DO
  END DO

  Etot=E_NSC+E_cor+E_Hxc+ewald+Ealphat-TS0
  Etot=Etot+E_extV

  string = "Etot and force calc"
  CALL timing_mpi(string,t_0)

  IF(inode.eq.1) THEN

    WRITE(6,*)
    WRITE(6,*) "---------------------------------------------------"
    WRITE(6,'(a,e21.14)') " Ewald        = ", ewald
    WRITE(6,'(a,e21.14)') " Alpha        = ", Ealphat
    WRITE(6,'(a,e21.14)') " E_NSC        = ", E_NSC
    WRITE(6,'(a,e21.14)') " E[-rho*V_Hxc]= ", E_cor
    WRITE(6,'(a,e21.14)') " E_Hxc        = ", E_Hxc
    WRITE(6,'(a,e21.14)') " E_extV       = ", E_extV
    WRITE(6,'(a,2(e17.10,1x))') " E_rhoVext,E_IVext     =", E_rhoVext,E_IVext
    WRITE(6,'(a,e21.14)') " -TS          = ", -TS
    WRITE(6,'(a,e21.14)') " E_tot        = ", Etot
    WRITE(17,FMT1704) "  E_tot            =",Etot
    CALL system_flush(17)
    WRITE(6,*) "---------------------------------------------------"
    WRITE(6,*) "E_ion,E_coul,E_xc",E_ion,E_coul,E_Hxc-E_coul
    WRITE(6,*) "ave(vtot)=ave_(V_xc)=: v0=",vv0
    WRITE(6,*) "-------------------------------------------"

  END IF

  RETURN

END SUBROUTINE recalc




END SUBROUTINE TDDFT
