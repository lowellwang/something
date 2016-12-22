subroutine calc_dipole(flag, AL, nkpt, islda, frac, dipole, workr_phi, tot, ist, ied)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use fft_data
  use load_data
  use data
  use data_TDDFT

  implicit none

  include 'mpif.h'
  include 'param.escan_real_f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: flag, nkpt, islda
  integer :: ist, ied
  real(8) :: tot

  real(8), dimension(2, mst) :: frac, rtmp
  real(8), dimension(3, 3) :: AL
  real(8), dimension(mst, mst, nkpt, islda) :: dipole

  complex(8), dimension(mr_n, mst) :: workr_phi

  integer :: kpt, iislda, ii, jj, m, n, ierr, ilumo
  integer, dimension(3) :: igrid, ngrid

  real(8) :: omega
  real(8), dimension(3) :: rn
  !real(8), dimension(3, mst) :: rc, tmp

  complex(8) :: wf2, wf2_2, wf2_psi, wf2_psi_2
  complex(8), dimension(3, mst, mst) :: P, ctmp
  !complex(8), dimension(mr_n, mst) :: workr_psi

  real(8), parameter :: PI = 3.14159265d0
  complex(8), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0), ci = (0.d0, 1.d0)

  ngrid(1) = n1
  ngrid(2) = n2
  ngrid(3) = n3
  ilumo=floor(tot/2.d0+0.1)+1

  dipole = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc states distribution, for blend system
  ! Calc dipole matrix <psi_i|r|psi_j>
  do iislda = 1, islda
    do kpt = 1, nkpt

      if(nkpt.ne.1) then
        call gen_G_comp(kpt, 0)
        call fftprep_comp(n1, n2, n3)
        ng_n = ngtotnod(inode, kpt)
      endif

      workr_phi = zero
      !workr_psi = zero
      frac = 0.d0
      P = zero

      call d3fft_comp_block(ug_n_bp(1, ist), workr_phi, -1, kpt, ied - ist + 1)

      do ii = 1, nr_n

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        jj = ii + (inode - 1) * nr_n
        igrid(1) = (jj - 1) / (n2 * n3) + 1
        igrid(2) = (jj - 1 - (igrid(1) - 1) * n2 * n3) / n3 + 1
        igrid(3) = jj - (igrid(1) - 1) * n2 * n3 - (igrid(2) - 1) * n3

        if((igrid(3) - 1).le.gz_1) then
          omega = 1.d0
        elseif((igrid(3) - 1).ge.gz_2) then
          omega = 0.d0
        else
          omega = cos(PI * (igrid(3) - 1 - gz_1) / 2.D0 / (gz_2 - gz_1))**2
        endif

        if(flag.eq.1) then
          rn(:) = AL(:, 1) * (dble(igrid(1)) - dble(ngrid(1))/2.d0 - 1.d0) / dble(ngrid(1)) + &
                  AL(:, 2) * (dble(igrid(2)) - dble(ngrid(2))/2.d0 - 1.d0) / dble(ngrid(2)) + &
                  AL(:, 3) * (dble(igrid(3)) - dble(ngrid(3))/2.d0 - 1.d0) / dble(ngrid(3))
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do m = ist, ied
          wf2 = workr_phi(ii, m - ist + 1)
          frac(1, m) = frac(1, m) + omega * abs(wf2)**2
          frac(2, m) = frac(2, m) + (1.d0 - omega) * abs(wf2)**2
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(flag.eq.1) then
          do m = ist+1, ied
            do n = ist, m-1
              wf2   = workr_phi(ii, m - ist + 1)
              wf2_2 = workr_phi(ii, n - ist + 1)
              P(:, n, m) = P(:, n, m) + rn(:) * conjg(wf2) * wf2_2
            enddo
          enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      call mpi_allreduce(frac, rtmp, 2*mst, MPI_REAL8, MPI_SUM, MPI_COMM_K, ierr)
      do m = 1, mst
        if((rtmp(1, m) + rtmp(2, m)).eq.0.d0) then
          frac(1, m) = 0.5d0
          frac(2, m) = 0.5d0
        else
          frac(1, m) = rtmp(1, m) / (rtmp(1, m) + rtmp(2, m))
          frac(2, m) = rtmp(2, m) / (rtmp(1, m) + rtmp(2, m))
        endif
      enddo

      if(flag.eq.1) then
        call mpi_allreduce(P, ctmp, 3*mst*mst, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_K, ierr)
        P = ctmp
        do m = ist+1, ied
          do n = ist, m-1
            dipole(n, m, kpt, iislda) = abs(P(1,n,m))**2 + abs(P(2,n,m))**2 + abs(P(3,n,m))**2
          enddo
        enddo
      endif

    enddo
  enddo

end subroutine calc_dipole
