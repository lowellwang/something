      module data_TDDFT

      complex*16,allocatable,dimension(:,:,:)  :: cpsi_td
      real*8,allocatable,dimension(:,:,:) :: occ0
      real*8,allocatable,dimension(:,:,:) :: E_td

      integer mst2,mmn
      integer ntime,ntime_init,n_dt,mscf
      real*8 tolrho,tolintgl,init_time,imp_k

      contains
      
      subroutine data_allocate_TDDFT(mg_nx,mst,mmn,nkpt,islda)

      implicit none

      integer mg_nx,mst,mmn,nkpt,islda
      allocate(cpsi_td(mg_nx,mmn,islda))
      allocate(occ0(mst,nkpt,islda))
      allocate(E_td(mmn,nkpt,islda))

      end subroutine data_allocate_TDDFT


      subroutine data_deallocate_TDDFT()

      implicit none
      
      deallocate(cpsi_td)
      deallocate(occ0)
      deallocate(E_td)
      
      end subroutine data_deallocate_TDDFT

      end module data_TDDFT
