SUBROUTINE bandshift(kpt, iislda, E_st, AL, tot, iatom, xatom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Correct band energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Assume there are two components A & B in the system
!     1. Generate the mask function mask(r),
!            phi(iA) = mask(r) * phi(i)
!            phi(iB) = (1 - mask(r)) * phi(i)
!     2. Calculate the dot_product:
!            dot_i_iN(i, iN) = <phi(i)|phi(iN)>, N = A,B
!     3. Generate the corrected eigen energy:
!            E_corrected(iN, iM) = 
!                      E_st(i) + Delta1, if N = M = A;
!                      E_st(i) + Delta2, if N = M = B;
!                      E_st(i) + (Delta1 + Delta2)/2, if N != M;
!     4. The corrected Hamiltonian should be (phi(i) basis):
!            H_corrected(i, j) = 
!                      \sum_{k,N,M} dot_i_iN(i, kN) * &
!                                   E_corrected(kN, kM) * &
!                                   conjg(dot_i_iN(j, kM))
!     5. Solve H_corrected to get ug_corrected
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE fft_data
  USE load_data
  USE data
  USE data_TDDFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Parameters passed in
  integer :: iislda, kpt, ierr
  real(8), dimension(mst) :: E_st
  real(8), dimension(3,3) :: AL
  real(8) :: tot   ! Total electron charge
  integer, dimension(matom) :: iatom
  real(8), dimension(3,matom) :: xatom

  ! mask function
  real(8), dimension(mr_n) :: mask

  ! corrected eigen energy
  complex(8), dimension(2*mst, 2*mst) :: E_corrected
  !real(8), dimension(2, 2, mst) :: E_corrected

  ! store phi(i) and phi(iN) in real space
  complex(8), dimension(mr_n, mst) :: workr_phi
  complex(8), dimension(mr_n, mst) :: workr_phi_iN

  ! dot product
  complex(8), dimension(mst, mst) :: dot_i_iN_half, dot_i_iN_tmp
  complex(8), dimension(mst, 2*mst) :: dot_i_iN

  ! corrected Hamiltonian
  complex(8), dimension(mst, mst) :: H_corrected

  ! corrected wft
  complex(8), dimension(mg_nx, mst) :: ug_corrected

  ! energy correction
  real(8) :: dV_1_CB, dV_1_VB, dV_2_CB, dV_2_VB, dV1, dV2

  ! the border of A and B, only available in Z axis now
  real(8) :: z1, z2

  ! normalization factor for workr_phi
  real(8) :: rfac

  ! timing
  real(8) :: time_begin, time_end

  ! tmp variables
  integer :: imn, nlumo
  integer :: i, j, k, info
  integer :: ix, iy, iz, ii, jj, iii, jjj
  real(8) :: rsum
  complex(8) :: csum

  real(8), dimension(mst, mst) :: ortho
  real(8), dimension(mst) :: W
  real(8), dimension(10*mst) :: rwork
  complex(8), dimension(10*mst) :: work
  complex(8), dimension(2*mst, mst) :: MTMP

  ! constant param, better to write them in a module
  REAL(8),    PARAMETER :: BOHR   = 0.529177d0
  REAL(8),    PARAMETER :: HART   = 27.211396d0
  REAL(8),    PARAMETER :: PI     = 4.d0 * ATAN(1.d0)
  COMPLEX(8), PARAMETER :: CI     = (0.d0, 1.d0)
  COMPLEX(8), PARAMETER :: ONE    = (1.d0, 0.d0)
  COMPLEX(8), PARAMETER :: ZERO   = (0.d0, 0.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get values from bandshift0.f

!  if(inode.eq.1) write(6,*) gz_1,gz_2,dP
  z1 = gz_1
  z2 = gz_2

  ! energy corrections
  ! better to read them from external files
  !dV_1_VB =-1.0295D0/HART - dP
  dV_1_VB =-0.6865D0/HART
  dV_2_VB =-0.6865D0/HART
  dV_1_CB = 1.0323D0/HART + dP
  dV_2_CB = 0.7920D0/HART

  ! find the LUMO state
  nlumo = floor(tot / 2.d0 + 0.1) + 1
  if(inode.eq.1) write(6,*) "LUMO=",nlumo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(inode.eq.1) then
    write(6,*)
    write(6,*) "------------------------------"
    write(6,*) "| Band Shift Test Start"
!    write(6,*) dV_1_VB, dV_2_VB, dV_1_CB, dV_2_CB
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! mask function

  mask = 0.d0
  do ii = 1, nr_n
    jj = ii + (inode-1)*nr_n
    ix = (jj-1)/(n2*n3) + 1
    iy = (jj-1- (ix-1)*n2*n3) / n3 + 1
    iz = jj - (ix-1)*n2*n3 - (iy-1)*n3
    if((iz-1).le.z1) then
      mask(ii) = 1.D0
    elseif((iz-1).gt.z2) then
      mask(ii) = 0.D0
    elseif((iz-1).gt.z1.and.(iz-1).le.z2) then
      mask(ii) = cos(PI * (iz-1-z1) / 2.D0 / (z2-z1))**2 ! use cos function if z1<z<z2
    endif
  enddo

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| Mask Calc Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! FFT phi(i) to real space

  ng_n = ngtotnod(inode, kpt)

  workr_phi = ZERO
  call d3fft_comp_block(ug_n_bp, workr_phi, -1, kpt, mst)

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| FFT Calc Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get phi(iN) in real space

  rfac = vol / (nr_nL * nnodes)
  do i=1,mst
    workr_phi_iN(:, i) = mask(:) * workr_phi(:, i) * rfac
  enddo

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| Phi(iN) Calc Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! dot_product <phi(i)|phi(i, iN)>

  ! first only do half to save time
  dot_i_iN_half = ZERO
  call zgemm('C', 'N', mst, mst, mr_n, ONE, &
              workr_phi, mr_n, workr_phi_iN, mr_n, &
              ONE, dot_i_iN_half, mst)

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| Dot Product Calc Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

  ! allreduce from mpi_comm_world? need to check~~~~~!
  dot_i_iN_tmp = ZERO
  call MPI_allreduce(dot_i_iN_half, dot_i_iN_tmp, mst*mst, &
                     MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
  dot_i_iN_half = dot_i_iN_tmp

  ! <phi(i)|phi(j)> = delta(i,j),
  ! so <phi(iB)|phi(i)> = delta(i,j) - <phi(iA)|phi(i)>
  dot_i_iN(1:mst, 1:mst) = dot_i_iN_half(1:mst, 1:mst)
  do i = 1, mst
    do j = 1, mst
      if(j.eq.i) then
        dot_i_iN(j, i+mst) = ONE - dot_i_iN(j, i)
      else
        dot_i_iN(j, i+mst) = ZERO - dot_i_iN(j, i)
      endif
    enddo
  enddo

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| Allreduce finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  if(inode.eq.1) then
!    open(11,file='cdot_check')
!    rewind(11)
!    do i=140,150
!      write(11,*) "i=",i
!      do j=1,320
!        write(11,*) dot_i_iN(j,i),dot_i_iN(j,i+mst)
!      enddo
!    enddo
!    close(11)
!  endif

!   call MPI_abort(MPI_COMM_WORLD,ierr)
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! correct eigen values

  E_corrected = ZERO
  do i = 1,mst
    if(i.ge.nlumo) then ! correct CB states
      dV1 = dV_1_CB
      dV2 = dV_2_CB
!    elseif(i.eq.(nlumo+1)) then ! additional correction for higher states
!      dV1 = dV_1_CB + 0.234d0/HART
!      dV2 = dV_2_CB
!    elseif(i.eq.(nlumo+2)) then ! additional correction for higher states
!      dV1 = dV_1_CB + 0.234d0/HART
!      dV2 = dV_2_CB
!    elseif(i.ge.(nlumo+3)) then ! additional correction for higher states
!      dV1 = dV_1_CB + 0.234d0/HART
!      dV2 = dV_2_CB
    else ! VB states
      dV1 = dV_1_VB
      dV2 = dV_2_VB
    endif
    E_corrected(i,     i)     = E_st(i) + dV1
    E_corrected(i+mst, i)     = E_st(i) + (dV1+dV2)/2.d0
    E_corrected(i,     i+mst) = E_st(i) + (dV1+dV2)/2.d0
    E_corrected(i+mst, i+mst) = E_st(i) + dV2
  enddo

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| E_corrected Calc Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calculate H_corrected

  ! MTMP = E_correct * conjg(dot_i_iN)
  MTMP = ZERO
  call zgemm('N', 'C', 2*mst, mst, 2*mst, ONE, &
              E_corrected, 2*mst, dot_i_iN, mst, &
              ONE, MTMP, 2*mst)

  ! H_corrected = dot_i_iN * MTMP
  H_corrected = ZERO
  call zgemm('N', 'N', mst, mst, 2*mst, ONE, &
              dot_i_iN, mst, MTMP, 2*mst, &
              ONE, H_corrected, mst)

!  if(inode.eq.1) then
!    write(6,*) "H_corrected_1:"
!    write(6,*) (H_corrected(j, 1), j = 1, 5)
!  endif

!  H_corrected = ZERO
!  do i = 1, mst
!    do j = 1, mst
!      csum = ZERO
!      do k = 1, mst
!        do ii = 1, 2
!          do jj = 1, 2
!            iii = k + (ii-1)*mst
!            jjj = k + (jj-1)*mst
!            csum = csum + dot_i_iN(i,iii) * E_corrected(iii,jjj) * &
!                          conjg(dot_i_iN(j,jjj))
!!            if(inode.eq.1) write(6,*) dot_i_iN(i,k+(ii-1)*mst) * E_corrected(ii,jj,k) * &
!!                          conjg(dot_i_iN(j,k+(jj-1)*mst))
!          enddo
!        enddo
!      enddo
!      H_corrected(i,j) = csum
!    enddo
!  enddo
!
!  if(inode.eq.1) then
!    write(6,*) "H_corrected_1:"
!    write(6,*) (H_corrected(j, 1), j = 1, 5)
!  endif
!
!  call MPI_abort(MPI_COMM_WORLD,ierr)
!  stop

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "| H' Calc finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! solve H_corrected

  call zheev('V', 'U', mst, H_corrected, mst, W, work, 5*mst, rwork, info)

  if(info.ne.0) then
    write(6,*) "error",info,inode
    stop
  endif

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "|  Diagnose Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! force orthonormal
  call ortho_Gram(H_corrected, mst, mst)
!  ! check orthonormal
!
!  do j = 1, mst
!    do i = 1, mst
!      csum = ZERO
!      do k = 1, mst
!        csum = csum + conjg(H_corrected(k,i)) * H_corrected(k,j)
!      enddo
!      ortho(i,j) = abs(csum) !*vol
!    enddo
!  enddo
!
!  do j = 1, mst
!    rsum = 0.d0
!    do i = 1, mst
!      if(i.ne.j.and.ortho(i,j).gt.rsum) rsum = ortho(i,j)
!    enddo
!    if(inode.eq.1) write(6,*) "i,ii,ij", j, ortho(j,j), rsum
!  enddo
!
!  if(inode.eq.1) then
!    do j=1,mst
!      write(6,*) 27.211396D0*W(j),27.211396D0*E_st(j)
!    enddo
!  endif

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "|  Orthonormal Finished,",time_end-time_begin
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calc new wft

  ug_corrected = ZERO
  call zgemm('N', 'N', mg_nx, mst, mst, ONE, &
             ug_n_bp, mg_nx, H_corrected, mst, ONE, ug_corrected, mg_nx)

  ug_n_bp = ug_corrected
  E_st(:) = W(:)

  if(inode.eq.1) then
    time_end=mpi_wtime()
    write(6,*) "|  New ug Finished,",time_end-time_begin
    write(6,*) "-------------------------------"
    time_begin=mpi_wtime()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return

END SUBROUTINE bandshift
