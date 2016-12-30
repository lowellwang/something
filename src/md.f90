subroutine VVMDinit(Iseed, myid, iscale, InitTemp, Etot, TotalEn, Ekin, &
                    ivlct, ntemp, Vi, V_output)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: Iseed, myid, iscale, ivlct, ntemp
  real(8), dimension(3,matom) :: Vi, V_output
  real(8)                     :: InitTemp, Etot, TotalEn, Ekin

  integer               :: i, j, k
  real(8), dimension(3) :: vc, vd
  real(8)               :: mv2, Temp, total_mass, scaling, v_width
  real(8)               :: Ekin_imp
  real(8), external     :: RANG
  real(8), parameter    :: Hart = 27.211396D0
  real(8), parameter    :: Boltz = 0.023538D0 / 273.15D0 / 27.211396D0

  if(ivlct.le.0) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize velocities with maxwell-boltzmann distribution

    Vi = 0.d0
    mv2 = 0.d0
    CALL RANSET(Iseed)

    do i = 1, ntemp
      v_width = sqrt(InitTemp * Boltz / MDatom(i))
      do j = 1, 3
        Vi(j,i) = RANG(Iseed, 0.d0, v_width)
      enddo
      mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo

    ! scale velocities according to temperature, for the first time
    if(InitTemp.le.0.d0) then
      scaling = 0.d0
    else
      scaling = dabs(InitTemp * Boltz * 3.d0 * DBLE(ntemp) / mv2)
    endif
    if(myid.eq.1) then
      write(6,*) '*********************************'
      write(6,*) "Initialize velocities: first scaling is"
      write(6,*) scaling
    endif
    Vi(:,:) = Vi(:,:) * sqrt(scaling)
    mv2 = 0.d0
    do i = 1, ntemp
      mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endif

  if(ntemp.eq.1) goto 100

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! remove any spurious drift
  vd = 0.d0
  do i = 1, ntemp
    do j = 1, 3
      vd(j) = vd(j) + Vi(j,i)
    enddo
  enddo
  vd(:) = -vd(:) / ntemp
  if(myid.eq.1) then
    write(6,*) "total velocity drift:"
    write(6,*) vd(1), vd(2), vd(3)
  endif
  do i = 1, ntemp
    Vi(:,i) = Vi(:,i) + vd(:)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! make sure total momentum is zero
  vc = 0.D0
  total_mass = 0.d0
  do i = 1, ntemp
    do j = 1, 3
      vc(j) = vc(j) + Vi(j,i) * MDatom(i)
    enddo
    total_mass = total_mass + MDatom(i)
  enddo
  vc(:) = -vc(:) / total_mass
  if(myid.eq.1) then
    write(6,*) "total momentum drift:"
    write(6,*) vc(1), vc(2), vc(3)
  endif

  mv2 = 0.D0
  do i = 1, ntemp
    Vi(:,i) = Vi(:,i) + vc(:)
    mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! scale velocities for the second time
  if(ivlct.le.0) then
    if(InitTemp.le.0.d0) then
      scaling = 0.d0
    else
      scaling = dabs(InitTemp * Boltz * 3.d0 * DBLE(ntemp) / mv2)
    endif
    if(myid.eq.1) then
      write(6,*) "Initialize velocities: second scaling is"
      write(6,*) scaling
    endif
    Vi(:,:) = Vi(:,:) * sqrt(scaling)
    mv2 = 0.D0
    do i = 1, ntemp
      mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo
  endif

100 continue

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc current temperature
  Ekin = 0.5d0 * mv2
  Temp = 2.d0 * Ekin / (3.d0 * DBLE(ntemp) * Boltz)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc incidence ion kinetic energy, if they exist
  Ekin_imp = 0.d0
  if(ntemp.lt.natom) then
    do i = ntemp + 1, natom
      Ekin_imp = Ekin_imp + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo
    Ekin_imp = 0.5d0 * Ekin_imp
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! output
  V_output(:,:) = Vi(:,:)
  TotalEn = Etot + Ekin + Ekin_imp
  if(myid.eq.1) then
    write(6,*) '*********************************'
    write(6,*) "IniTotEn,Etot,Ekin,Ekin_imp"
    write(6,*) TotalEn*Hart, Etot*Hart, Ekin*Hart, Ekin_imp*Hart
    write(6,*) "Now Temp=", Temp, "Scaling=", scaling

    open(11,file='plot_MD.txt',access='append')
    write(11,'(A)') "  Time,             E_tot,            E_elec,             E_ion,         E_ion_imp,              Temp,           Scaling"
    write(11,666)   0, TotalEn*Hart, Etot*Hart, Ekin*Hart, Ekin_imp*Hart, Temp, scaling
    close(11)
666 format(i6,1x,12(f18.9,1x))
  endif

  return

end

subroutine VVMD(istep, myid, iscale, Etot, TotalEn, Ekin, &
                ntemp, Vi, V_output)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: istep, myid, iscale, ntemp
  real(8), dimension(3,matom) :: Vi, V_output
  real(8)                     :: Etot, TotalEn, Ekin

  integer               :: i, j, k
  real(8), dimension(3) :: vc, vd, vt
  real(8)               :: mv2, Temp, total_mass, scaling
  real(8)               :: Ekin_imp, Ekin_old
  real(8), parameter    :: Hart = 27.211396D0
  real(8), parameter    :: Boltz = 0.023538D0 / 273.15D0 / 27.211396D0

  if(ntemp.eq.1) goto 110

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! make sure total momentum is zero
  vc = 0.D0
  total_mass = 0.d0
  do i = 1,natom
    vc(:) = vc(:) + Vi(:,i) * MDatom(i)
    total_mass = total_mass + MDatom(i)
  enddo
  vc(:) = -vc(:) / total_mass

  mv2 = 0.D0
  do i = 1, ntemp
    Vi(:,i) = Vi(:,i) + vc(:)
    mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
  enddo

110 continue

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc incidence ion kinetic energy, if they exist
  Ekin_imp = 0.d0
  if(ntemp.lt.natom) then
    do i = ntemp + 1, natom
      Ekin_imp = Ekin_imp + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo
    Ekin_imp = 0.5d0 * Ekin_imp
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc energies, temp and scaling factor
  Ekin_old = 0.5d0 * mv2
  temperature = 2.d0 * Ekin_old / (3.d0 * DBLE(ntemp) * Boltz)
  Ekin = TotalEn - Etot - Ekin_imp
  scaling = dabs(Ekin / Ekin_old)

  if(myid.eq.1) then
    write(6,*) '*********************************'
    write(6,*) 'Time step', istep
    write(6,*) 'Before scaling'
    write(6,*) '*********************************'
    write(6,*) 'Ekin,Should', Ekin_old*Hart, Ekin*Hart
    write(6,*) 'Now Temp=', temperature, 'Scaling=', scaling
    write(6,*) "total momentum drift:"
    write(6,*) vc(1), vc(2), vc(3)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! scale vel if needed
  vt = 0.D0
  mv2 = 0.D0
  if(iscale.gt.0) Vi(:,:) = Vi(:,:) * sqrt(scaling)
  do i = 1, ntemp
    vt(:) = vt(:) + Vi(:,i) * MDatom(i)
    mv2 = mv2 + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
  enddo
  vt(:) = -vt(:) / total_mass

  Ekin = 0.5d0 * mv2
  temperature = 2.d0 * Ekin / (3.d0 * DBLE(ntemp) * Boltz)
  TotalEn = Etot + Ekin + Ekin_imp
  V_output(:,:) = Vi(:,:)

  if(myid.eq.1) then
    write(6,*) '*********************************'
    if(iscale.gt.0) write(6,*) 'After scaling'
    if(iscale.le.0) write(6,*) 'No scaling'
    write(6,*) '*********************************'
    write(6,*) 'Ekin', Ekin*Hart
    write(6,*) 'Now Temp=', temperature
    write(6,*) "total momentum drift:"
    write(6,*) vt(1), vt(2), vt(3)

    open(11,file='plot_MD.txt',access='append')
    write(11,666) istep, TotalEn*Hart, Etot*Hart, Ekin*Hart, Ekin_imp*Hart, temperature, scaling
    close(11)
666 format(i6,1x,12(f18.9,1x))
  endif

  return

end

subroutine nose_hoover_chain(DesiredTemp, Q, M, dt, ntemp, xi, vxi, axi, Ekin, Vi)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: M, ntemp
  real(8)                     :: dt, Ekin, DesiredTemp
  real(8), dimension(7)       :: Q, xi, vxi, axi
  real(8), dimension(3,matom) :: Vi
  real(8), parameter          :: Hart = 27.211396D0
  real(8), parameter          :: Boltz = 0.023538D0 / 273.15D0 / 27.211396D0

  integer                     :: i, j, k, nomega, ntick, nf
  real(8)                     :: dts, dts2, dts4, dts8, aa, kT
  real(8), dimension(7)       :: omega

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Yoshida Suzuki integration scheme
  ! nomega=3 && ntick=4 will be good enough in most cases
  ! Notice that omega changes with nomega
  ! If more stable integration is required, try nomega=5 or larger ntick
  nomega = 3
  omega(1) = 1.d0 / (2.d0 - 2.d0**(1.d0/3.d0))
  omega(2) = 1.d0 - 2.d0 * omega(1)
  omega(3) = omega(1)
  ntick = 4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  kT = DesiredTemp * Boltz
  nf = 3 * ntemp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(inode_tot.eq.1) write(6,*) "nf, nomega, ntick, tmpend", nf, nomega, ntick, DesiredTemp

  do k = 1, ntick
  do j = 1, nomega
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dts = omega(j) * dt / dble(ntick)
    dts2 = dts / 2.d0
    dts4 = dts / 4.d0
    dts8 = dts / 8.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = M, 2, -1
      axi(i) = (Q(i-1)*vxi(i-1)*vxi(i-1) - kT) / Q(i)
      vxi(i) = vxi(i) + dts4*axi(i)
      vxi(i-1) = vxi(i-1) * exp(-dts8*vxi(i))
    enddo
    axi(1) = (2.d0*Ekin - nf*kT) / Q(1)
    vxi(1) = vxi(1) + dts4*axi(1)
    vxi(1) = vxi(1) * exp(-dts8*vxi(2))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    aa = exp(-dts2 * vxi(1))
    Vi(:,:) = Vi(:,:) * aa
    Ekin = Ekin * aa * aa
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, M
      xi(i) = xi(i) + dts2*vxi(i)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vxi(1) = vxi(1) * exp(-dts8*vxi(2))
    axi(1) = (2.d0*Ekin - nf*kT) / Q(1)
    vxi(1) = vxi(1) + dts4*axi(1)
    do i = 2, M
      vxi(i-1) = vxi(i-1) * exp(-dts8*vxi(i))
      axi(i) = (Q(i-1)*vxi(i-1)*vxi(i-1) - kT) / Q(i)
      vxi(i) = vxi(i) + dts4*axi(i)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo
  enddo

  return

end


    
