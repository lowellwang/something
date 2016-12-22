      module data_pulay_6

      real*8,allocatable,dimension(:,:,:)      :: dw_6,dR_6
      real*8,allocatable,dimension(:,:)        :: R0_6,w_in0_6
      real*8,allocatable,dimension(:,:)        :: AA_pulay_6
      integer npulay_max_6,nint_pulay_6,num_dwdR_6,nreset_pulay_6
      
      contains
      
      subroutine data_allocate_6(mr_nL,islda, npulay_max_6)

      implicit none

      integer mr_nL,islda, npulay_max_6
      allocate(dw_6(mr_nL,npulay_max_6,islda))    ! store them on disk later
      allocate(dR_6(mr_nL,npulay_max_6,islda))    ! store them on disk later
      allocate(R0_6(mr_nL,islda))
      allocate(w_in0_6(mr_nL,islda))
      allocate(AA_pulay_6(npulay_max_6,npulay_max_6))
      nint_pulay_6=0
      
      end subroutine data_allocate_6


      subroutine data_deallocate_6()

      implicit none
      
      deallocate(dw_6)
      deallocate(dR_6)
      deallocate(R0_6)
      deallocate(w_in0_6)
      deallocate(AA_pulay_6)
      
      end subroutine data_deallocate_6

      end module data_pulay_6
