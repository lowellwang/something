subroutine bandshift0(xatom, iatom, AL)

  use fft_data
  use load_data
  use data
  use data_TDDFT
  implicit none
  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

  real(8), dimension(3, matom) :: xatom
  integer, dimension(matom) :: iatom
  real(8), dimension(3, 3) :: AL

  integer i, j, isurf1, isurf2
  real(8) :: dis_qd_CtoS(2), dis_qd_ml, diameter, e_SC
  real(8), dimension(3) :: rtmp, rtmp1, rtmp2, rtmp3
  real(8), dimension(3, 2) :: mcenter
  real(8), dimension(2) :: totmass
  real(8), parameter :: e_SC_bulk = 6.2d0
  real(8), parameter :: e_sol = 2.17d0
  real(8), parameter :: Bohr = 0.529177d0
  real(8), parameter :: Hart = 27.211396d0

  isurf1 = 62
  isurf2 = 63
  !isurf1 = 1
  !isurf2 = 2

  e_SC = 4.71d0
  !e_SC = 4.35d0

  mcenter = 0.d0
  totmass = 0.d0
  gz_1 =-10000.d0
  gz_2 = 10000.d0

  do i = 1, natom
    if(iatom(i) .eq. 1 .or. &
       iatom(i) .eq. 6 .or. &
       iatom(i) .eq. 7) then
      totmass(1) = totmass(1) + MDatom(i)
      do j = 1, 3
        mcenter(j, 1) = mcenter(j, 1) + MDatom(i) * xatom(j, i)
      enddo
      if(xatom(3, i) .gt. gz_1) gz_1 = xatom(3, i)
    elseif(iatom(i) .eq. 34 .or. &
           iatom(i) .eq. 48) then
      totmass(2) = totmass(2) + MDatom(i)
      do j = 1, 3
        mcenter(j, 2) = mcenter(j, 2) + MDatom(i) * xatom(j, i)
      enddo
      if(xatom(3, i) .lt. gz_2) gz_2 = xatom(3, i)
    endif
  enddo

  mcenter(:,1) = mcenter(:,1) / totmass(1)
  mcenter(:,2) = mcenter(:,2) / totmass(2)

  gz_1 = gz_1 * AL(3,3)
  gz_2 = gz_2 * AL(3,3)
  if(gz_1 .lt. gz_2) then
    if(gz_1.lt.0.d0) then
      write(6,*) "only 2"
      gz_1 = 0.d0
      gz_2 = 0.d0
    elseif(gz_2.gt.AL(3,3)) then
      write(6,*) "only 1"
      gz_1 = AL(3,3)
      gz_2 = AL(3,3)
    elseif((gz_1 + 1.889d0) .lt. (gz_2 - 1.889d0)) then
      write(6,*) "good"
      ! gz_1 = gz_1 + 1.889d0
      ! gz_2 = gz_2 - 1.889d0
    else 
      write(6,*) "soso"
    endif
  else
    write(6,*) "bad"
    rtmp(1) = (gz_1 + gz_2) / 2.d0
    gz_1 = rtmp(1) - 0.945d0
    gz_2 = rtmp(1) + 0.945d0
  endif
  gz_1 = gz_1 / AL(3,3) * n3
  gz_2 = gz_2 / AL(3,3) * n3

  do i = 1, 3
    rtmp1(i) = abs(mcenter(i, 1) - mcenter(i, 2))
    rtmp2(i) = abs(xatom(i, isurf1) - mcenter(i, 2))
    rtmp3(i) = abs(xatom(i, isurf2) - mcenter(i, 2))
  enddo
  dis_qd_ml = (AL(1, 1) * rtmp1(1) + AL(1, 2) * rtmp1(2) + AL(1, 3) * rtmp1(3))**2 + &
              (AL(2, 1) * rtmp1(1) + AL(2, 2) * rtmp1(2) + AL(2, 3) * rtmp1(3))**2 + &
              (AL(3, 1) * rtmp1(1) + AL(3, 2) * rtmp1(2) + AL(3, 3) * rtmp1(3))**2
  dis_qd_ml = sqrt(dis_qd_ml)
  dis_qd_CtoS(1) = (AL(1, 1) * rtmp2(1) + AL(1, 2) * rtmp2(2) + AL(1, 3) * rtmp2(3))**2 + &
                   (AL(2, 1) * rtmp2(1) + AL(2, 2) * rtmp2(2) + AL(2, 3) * rtmp2(3))**2 + &
                   (AL(3, 1) * rtmp2(1) + AL(3, 2) * rtmp2(2) + AL(3, 3) * rtmp2(3))**2
  dis_qd_CtoS(1) = sqrt(dis_qd_CtoS(1))
  dis_qd_CtoS(2) = (AL(1, 1) * rtmp3(1) + AL(1, 2) * rtmp3(2) + AL(1, 3) * rtmp3(3))**2 + &
                   (AL(2, 1) * rtmp3(1) + AL(2, 2) * rtmp3(2) + AL(2, 3) * rtmp3(3))**2 + &
                   (AL(3, 1) * rtmp3(1) + AL(3, 2) * rtmp3(2) + AL(3, 3) * rtmp3(3))**2
  dis_qd_CtoS(2) = sqrt(dis_qd_CtoS(2))

  diameter = dis_qd_CtoS(1) + dis_qd_CtoS(2) + 1.d0 * 1.84d0 / Bohr
  !e_SC = 1.d0 + (e_SC_bulk - 1.d0) / (1.d0 + (7.5 / diameter / Bohr)**1.2d0)

  dP = 0.d0
  do i = 0, 100
    dP = dP - (e_SC - e_sol) * i / (e_SC * i + e_sol * (i+1)) * &
              (diameter/2.d0/dis_qd_ml)**(2*i+1)
  enddo
  dP = dP / 2.d0 / e_sol / dis_qd_ml

  write(6,*) "Diameter = ", diameter * Bohr
  write(6,*) "dielectric constant = ", e_SC
  write(6,*) "Dis = ", dis_qd_ml * Bohr
  write(6,*) "dP = ", dP * Hart
  write(6,*) "z_1 = ",gz_1
  write(6,*) "z_2 = ",gz_2

end subroutine bandshift0
