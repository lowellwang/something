subroutine VVMDinit(Iseed, myid, iscale, InitTemp, Etot, TotalEn, Enki, &
                    ivlct, ntemp, Vi, V_output)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: Iseed, myid, iscale, ivlct, ntemp
  real(8), dimension(3,matom) :: Vi, V_output
  real(8)                     :: InitTemp, Etot, TotalEn, Enki

  integer               :: i, j, k
  real(8), dimension(3) :: vc, vd
  real(8)               :: mv2, Temp, total_mass, scaling, v_width
  real(8)               :: Enki_imp
  real(8), external     :: RANG
  real(8), parameter    :: HFs2vB2M = 27.211396D0
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
  Enki = 0.5d0 * mv2 * HFs2vB2M
  Temp = 2.d0 * Enki / (3.d0 * DBLE(ntemp) * Boltz * HFs2vB2M)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc incidence ion kinetic energy, if they exist
  Enki_imp = 0.d0
  if(ntemp.lt.natom) then
    do i = ntemp + 1, natom
      Enki_imp = Enki_imp + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo
    Enki_imp = 0.5d0 * Enki_imp * HFs2vB2M
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! output
  V_output(:,:) = Vi(:,:)
  TotalEn = Etot * HFs2vB2M + Enki + Enki_imp
  if(myid.eq.1) then
    write(6,*) '*********************************'
    write(6,*) "IniTotEn,Etot,Enki,Enki_imp"
    write(6,*) TotalEn, Etot*HFs2vB2M, Enki, Enki_imp
    write(6,*) "Now Temp=", Temp, "Scaling=", scaling

    open(11,file='plot_MD.txt',access='append')
    write(11,'(A)') " Time,             E_tot,            E_elec,             E_ion,         E_ion_imp,              Temp,           Scaling"
    write(11,666)   0, TotalEn, Etot*HFs2vB2M, Enki, Enki_imp, Temp, scaling
    close(11)
666 format(i5,1x,12(f18.9,1x))
  endif

  return

end

subroutine VVMD(istep, myid, iscale, Etot, TotalEn, Enki, &
                ntemp, Vi, V_output)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: istep, myid, iscale, ntemp
  real(8), dimension(3,matom) :: Vi, V_output
  real(8)                     :: Etot, TotalEn, Enki

  integer               :: i, j, k
  real(8), dimension(3) :: vc, vd, vt
  real(8)               :: mv2, Temp, total_mass, scaling
  real(8)               :: Enki_imp, Enki_old
  real(8), parameter    :: HFs2vB2M = 27.211396D0
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
  Enki_imp = 0.d0
  if(ntemp.lt.natom) then
    do i = ntemp + 1, natom
      Enki_imp = Enki_imp + MDatom(i) * (Vi(1,i)**2 + Vi(2,i)**2 + Vi(3,i)**2)
    enddo
    Enki_imp = 0.5d0 * Enki_imp * HFs2vB2M
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc energies, temp and scaling factor
  Enki_old = 0.5d0 * mv2 * HFs2vB2M
  temperature = 2.d0 * Enki_old / (3.d0 * DBLE(ntemp) * Boltz * HFs2vB2M)
  Enki = TotalEn - Etot * HFs2vB2M - Enki_imp
  scaling = dabs(Enki / Enki_old)

  if(myid.eq.1) then
    write(6,*) '*********************************'
    write(6,*) 'Time step', istep
    write(6,*) 'Before scaling'
    write(6,*) '*********************************'
    write(6,*) 'Enki,Should', Enki_old, Enki
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

  Enki = 0.5d0 * mv2 * HFs2vB2M
  temperature = 2.d0 * Enki / (3.d0 * DBLE(ntemp) * Boltz * HFs2vB2M)
  TotalEn = Etot * HFs2vB2M + Enki + Enki_imp
  V_output(:,:) = Vi(:,:)

  if(myid.eq.1) then
    write(6,*) '*********************************'
    if(iscale.gt.0) write(6,*) 'After scaling'
    if(iscale.le.0) write(6,*) 'No scaling'
    write(6,*) '*********************************'
    write(6,*) 'Enki', Enki
    write(6,*) 'Now Temp=', temperature
    write(6,*) "total momentum drift:"
    write(6,*) vt(1), vt(2), vt(3)

    open(11,file='plot_MD.txt',access='append')
    write(11,666) istep,TotalEn,Etot*HFs2vB2M,Enki,Enki_imp,temperature,scaling
    close(11)
666 format(i5,1x,12(f18.9,1x))
  endif

  return

end

subroutine nose_hoover_chain(DesiredTemp, Q, M, dt, ntemp, xi, vxi, axi, Enki, Vi)

  use data
  use data_TDDFT
  implicit none
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: M, ntemp
  real(8)                     :: dt, Enki, DesiredTemp
  real(8), dimension(7)       :: Q, xi, vxi, axi
  real(8), dimension(3,matom) :: Vi
  real(8), parameter          :: HFs2vB2M = 27.211396D0
  real(8), parameter          :: Boltz = 0.023538D0 / 273.15D0 / 27.211396D0

  integer                     :: i, j, k, nomega, nf
  real(8)                     :: dts, dts2, dts4, dts8, aa, kT
  real(8), dimension(7)       :: omega

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Yoshida Suzuki integration scheme
  nomega = 3
  omega(1) = 1.d0 / (2.d0 - 2.d0**(1.d0/3.d0))
  omega(2) = 1.d0 - 2.d0 * omega(1)
  omega(3) = omega(1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  kT = DesiredTemp * Boltz
  nf = 3 * ntemp

  do j = 1, nomega
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dts = omega(j) * dt
    dts2 = dts / 2.d0
    dts4 = dts / 4.d0
    dts8 = dts / 8.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calc thermostat acceleration
    ! axi = d^2(xi)/dt^2
    axi(1) = (2.d0*Enki - nf*kT) / Q(1)
    do i = 2, M
      axi(i) = (Q(i-1)*vxi(i-1)*vxi(i-1) - kT) / Q(i)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update thermostat velocity
    ! vxi = d(xi)/dt
    vxi(M) = vxi(M) + dts4*axi(M)
    do i = 1, M-1
      aa = exp(-dts8 * vxi(M+1-i))
      vxi(M-i) = (vxi(M-i) + dts4*axi(M-i)) * aa
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update particle velocity
    aa = exp(-dts2 * vxi(1))
    Vi(:,:) = Vi(:,:) * aa
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update thermostat position xi
    xi(1:M) = xi(1:M) + dts2 * vxi(1:M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update ion kinetic energy
    Enki = Enki * aa * aa
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update axi and vxi
    axi(1) = (2.d0*Enki - nf*kT) / Q(1)
    do i = 2, M
      axi(i) = (Q(i-1)*vxi(i-1)*vxi(i-1) - kT) / Q(i)
    enddo
    do i = 1, M-1
      aa = exp(-dts8 * vxi(i+1))
      vxi(i) = (vxi(i) + dts4*axi(i)) * aa
      axi(i+1) = (Q(i)*vxi(i)*vxi(i) - kT) / Q(i+1)
    enddo
    vxi(M) = vxi(M) + dts4 * axi(M)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo

  return

end


    
