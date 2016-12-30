subroutine MVATOMS(itime, iscale, Delt, Box, xatom, fatom0, fatom1, &
                   AL, Etot, InitTemp, DesiredTemp, TotalEn, Ekin, &
                   ivlct, ntemp, Vi, V_output, &
                   itherm, nhchain, qmass, xi, vxi, axi)

  use fft_data
  use load_data
  use data
  use data_TDDFT
  implicit double precision (a-h,o-z)
  include 'mpif.h'
  include 'param.escan_real_f90'

  integer                     :: itime, iscale, ivlct, ntemp, itherm, nhchain
  real(8), dimension(3,matom) :: xatom, fatom0, fatom1, Vi, V_output
  real(8), dimension(3,3)     :: AL
  real(8), dimension(7)       :: qmass, xi, vxi, axi
  real(8)                     :: InitTemp, DesiredTemp, TotalEn, Etot, Ekin, Box, Delt

  integer                     :: Iseed, istep
  real(8)                     :: delth, temp0, temp1
  real(8), dimension(3)       :: xx
  real(8), parameter          :: Hart = 27.211396D0
  real(8), parameter          :: Boltz = 0.023538D0 / 273.15D0 / 27.211396D0

  if(inode_tot.eq.1) then
    write(6,*) "START MD"
  endif

  if(itime.eq.1) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! init MD
    Iseed=12345
    call VVMDinit(Iseed, inode_tot, iscale, InitTemp, &
                  Etot, TotalEn, Ekin, ivlct, ntemp, Vi, V_output)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif(itime.gt.1) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! continue MD
    istep = itime - 1
    delth = Delt * 0.5d0
    do i = 1,natom
      Vi(:,i) = Vi(:,i) - delth * fatom0(:,i) / MDatom(i)
    enddo
    call VVMD(istep, inode_tot, iscale, &
              Etot, TotalEn, Ekin, ntemp, Vi, V_output)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calc vel and xatom
  delth = Delt * 0.5d0

  do i = 1,natom
    Vi(:,i) = Vi(:,i) - delth * fatom1(:,i) / MDatom(i)
  enddo

  if(itherm.eq.1) then
    temp0 = 2.d0 * Ekin / (3.d0 * DBLE(ntemp) * Boltz)
    call nose_hoover_chain(DesiredTemp, qmass, nhchain, Delt, ntemp, &
                           xi, vxi, axi, Ekin, Vi)
    temp1 = 2.d0 * Ekin / (3.d0 * DBLE(ntemp) * Boltz)
    if(inode_tot.eq.1) then
      write(6,*) '*********************************'
      write(6,*) 'Thermostat: Nose-Hoover chain'
      write(6,*) qmass
      write(6,*) xi
      write(6,*) vxi
      write(6,*) axi
      write(6,'(A,I7,A,E11.4,A,E11.4)') &
      '  itime = ', itime-1, ", T0 = ", temp0, ", T1 = ", temp1
    endif
  endif

  do i = 1, natom
    xx(:) = Delt * Vi(:,i)
    xatom(1,i) = xatom(1,i) + xx(1)*ALI(1,1) + xx(2)*ALI(2,1) + xx(3)*ALI(3,1)
    xatom(2,i) = xatom(2,i) + xx(1)*ALI(1,2) + xx(2)*ALI(2,2) + xx(3)*ALI(3,2)
    xatom(3,i) = xatom(3,i) + xx(1)*ALI(1,3) + xx(2)*ALI(2,3) + xx(3)*ALI(3,3)
    xatom(:,i) = xatom(:,i) - BOX * int(xatom(:,i) / BOX)
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  return

end
