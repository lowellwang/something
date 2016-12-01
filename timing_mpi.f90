SUBROUTINE timing_mpi(string,t_0)
!-----------------------------------------

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

  INTEGER :: ierr
  REAL(8) :: t_0, t_1
  CHARACTER(50) :: string

  CALL MPI_barrier(MPI_COMM_WORLD,ierr)

  IF(inode.eq.1) THEN

    t_1 = mpi_wtime()

    WRITE(6,*)
    WRITE(6,*) "<-------------------------------->"
    WRITE(6,'(a,f12.4)') " Timing: "//trim(adjustl(string)), t_1-t_0
    WRITE(6,*) "<-------------------------------->"
    WRITE(6,*)

    t_0 = mpi_wtime()

  END IF

  RETURN

END SUBROUTINE timing_mpi
