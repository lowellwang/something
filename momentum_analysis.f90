SUBROUTINE momentum_analysis(AL, rc, nkpt, islda, lmax, SPD, workr, ist, ied)

  USE fft_data
  USE load_data
  USE data
  USE data_TDDFT

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CHARACTER(LEN=25) :: subname = 'momentum_analysis'

  INTEGER :: nkpt, islda, lmax, ist, ied
  REAL(8), DIMENSION(3) :: rc
  REAL(8), DIMENSION(3, 3) :: AL
  REAL(8), DIMENSION(-lmax:lmax, 0:lmax, ist:ied) :: SPD
  COMPLEX(8), DIMENSION(mr_n, mst) :: workr

  INTEGER :: kpt, iislda
  INTEGER :: ii, jj, nn, i, j, k, l, m, ierr
  INTEGER, DIMENSION(3) :: igrid, ngrid
  REAL(8), DIMENSION(3) :: rn
  REAL(8), ALLOCATABLE, DIMENSION(:) :: rad, theta, phi
  COMPLEX(8) :: wf2
  COMPLEX(8), ALLOCATABLE, DIMENSION(:, :, :) :: ylm_a, Blm, Blm_tmp

  REAL(8), PARAMETER :: pi = 4.d0 * ATAN(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ngrid(1) = n1; ngrid(2) = n2; ngrid(3) = n3
  nn = n1 * n2 * n3
  ierr = 0

  ALLOCATE(rad(nn), theta(nn), phi(nn))
  ALLOCATE(ylm_a(-lmax:lmax, 0:lmax, 1:nn))
  ALLOCATE(Blm(-lmax:lmax, 0:lmax, ist:ied))

  ylm_a = (0.d0, 0.d0)

  IF(inode.EQ.1) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(lmax.gt.2) THEN
      PRINT *, "!!ERROR: lmax > 2 in ",subname
      CALL MPI_abort(MPI_COMM_WORLD,ierr)
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ii = 1
    DO i = 0, ngrid(1) - 1
      DO j = 0, ngrid(2) - 1
        DO k = 0, ngrid(3) - 1

          rn(:) = AL(:, 1) * dble(i) / dble(ngrid(1)) + &
                  AL(:, 2) * dble(j) / dble(ngrid(2)) + &
                  AL(:, 3) * dble(k) / dble(ngrid(3))
          rn(:) = rn(:) - rc(:)
          IF(ABS(rn(1)).LT.1.D-10) rn(1) = 0.d0
          IF(ABS(rn(2)).LT.1.D-10) rn(2) = 0.d0
          IF(ABS(rn(3)).LT.1.D-10) rn(3) = 0.d0

          rad(ii) = SQRT(rn(1)**2 + rn(2)**2 + rn(3)**2)
          IF(rad(ii).lt.1.D-10) THEN
            rad(ii) = 0.d0; theta(ii) = 0.d0; phi(ii) = 0.d0
          ELSE
            theta(ii) = ACOS(rn(3)/rad(ii))
            IF(theta(ii).lt.1.D-10) THEN
              theta(ii) = 0.d0; phi(ii) = 0.d0
            ELSEIF(ABS(pi - theta(ii)).lt.1.D-10) THEN
              theta(ii) = pi; phi(ii) = 0.d0
            ELSE
              IF(zn(1).eq.0.d0 .AND. zn(2).eq.0.d0) THEN
                phi(ii) = 0.d0
              ELSE
                phi(ii) = ATAN2(rn(2), rn(1))
              ENDIF
            ENDIF
          ENDIF

          DO l = 0, lmax
            DO m = -l, l
              ylm_a(m, l, ii) = YLM(l, m, theta(ii), phi(ii))
            ENDDO
          ENDDO

          ii = ii + 1

        ENDDO
      ENDDO
    ENDDO

    IF(ii.NE.(nn+1)) THEN
      PRINT *, "!!ERROR: Grid num error in ",subname
      CALL MPI_abort(MPI_COMM_WORLD,ierr)
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENDIF

  CALL MPI_BCAST(rad, nn, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(theta, nn, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(phi, nn, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(ylm_a, (2*lmax+1)*(lmax+1)*nn, &
                 MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
        
  DO iislda = 1, islda
    DO kpt = 1, nkpt

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(nkpt.ne.1) THEN
        call gen_G_comp(kpt, 0)
        call fftprep_comp(n1, n2, n3)
        ng_n = ngtotnod(inode, kpt)
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      workr = (0.d0, 0.d0)
      CALL d3fft_comp_block(ug_n_bp(1, ist), workr, -1, kpt, ied - ist + 1)

      DO i = ist, ied
        DO l = 0, lmax
          DO m = -l, l
            Blm(m, l, i) = (0.d0, 0.d0)
            DO ii = 1, nr_n
              jj = ii + (inode - 1) * nr_n
              wf2 = workr(ii, i - ist + 1)
              Blm(m, l, i) = Blm(m, l, i) + &
                 CONJG(ylm_a(m, l, jj)) * wf2 * SIN(theta(jj))
            ENDDO

          ENDDO
        ENDDO
      ENDDO

      ALLOCATE(Blm_tmp(-lmax:lmax, 0:lmax, ist:ied))
      CALL MPI_ALLREDUCE(Blm, Blm_tmp, (2*lmax+1)*(lmax+1)*nn, &
                 MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_K, ierr)
      Blm = Blm_tmp
      DEALLOCATE(Blm_tmp)

          
        


END SUBROUTINE
