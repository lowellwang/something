      module data_TDDFT

      complex*16,allocatable,dimension(:,:,:,:)  :: cpsi_td
      real*8,allocatable,dimension(:,:,:) :: occ0
      real*8,allocatable,dimension(:,:,:) :: E_td
      real*8,allocatable,dimension(:)     :: vext_dt_n

      integer mst2,mmn
      integer ntime,ntime_init,mscf
      integer iocc,jhole,jelec,jxi
      integer noutput,inormal,ifatom,ikappa,ivext_dt
      integer iboltz,ibshift,i_scale,ibo_md
      real*8 tolrho,tolintgl,init_time,imp_k,rxi
      real*8 dP, gz_1, gz_2
      
      integer iMD
      real*8 temperature
      real*8 dtMD,MDatom(1800),MDtype(10)

      character(LEN=12),parameter :: FMT1701 = "(A,1X,I10)"
      character(LEN=14),parameter :: FMT1702 = "(A,1X,E10.3)"
      character(LEN=14),parameter :: FMT1703 = "(A,1X,F10.3)"
      character(LEN=15),parameter :: FMT1704 = "(A,1X,E23.16)"
      character(LEN=15),parameter :: FMT1705 = "(A,1X,F23.16)"

      contains
      
      subroutine data_allocate_TDDFT(mr_n,mg_nx,mst,mmn,nkpt,islda)

      implicit none

      integer mr_n,mg_nx,mst,mmn,nkpt,islda
      allocate(cpsi_td(mg_nx,mmn,nkpt,islda))
      allocate(occ0(mst,nkpt,islda))
      allocate(E_td(mmn,nkpt,islda))
      allocate(vext_dt_n(mr_n))

      end subroutine data_allocate_TDDFT


      subroutine data_deallocate_TDDFT()

      implicit none
      
      deallocate(cpsi_td)
      deallocate(occ0)
      deallocate(E_td)
      deallocate(vext_dt_n)
      
      end subroutine data_deallocate_TDDFT

      end module data_TDDFT
