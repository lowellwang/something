      subroutine mch_pulay_wrap(w_in,w_out,nint0,islda)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
ccccc  Now, this is in n1L,n2L,n3L
ccccc   This pulay wrap around the dw_6, dR_6, so it can use those
ccccc   values from the previous MD or relaxation time steps
cccccc  to use this subroutine, the nint_pulay_6=0 must be declaired
cccccc  before the MD simulation, this is done in data_allocate_6


      
c      use data_fft_L2_14
c      use data_fft_L_13
c      use data_fft_12
      
c      use data_Gcolumn_wave_9
c      use data_Gcolumn_L_10
c      use data_Gcolumn_L2_11
c      use data
c      use param_escan
c      use data_nonlocal_8
c      use data_wave_7
      use data_pulay_6
c      use data_vrho_5
c      use data_vrho_G_3

c      use data_commun_0
      
      use fft_data
      use load_data
      use data
      use data_TDDFT

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include 'mpif.h'

      real*8 w_in(mr_nL,islda),w_out(mr_nL,islda)

      real*8 tmp(npulay_max_6)
      real*8, allocatable, dimension(:,:)  :: AA1
      real*8, allocatable, dimension(:)  :: B
      integer nreset,ierr

      allocate(AA1(npulay_max_6,npulay_max_6))
      allocate(B(npulay_max_6))

      alpha2=1.d0   
      ! We could set alpha2=1.1, 1.2 etc, the code should work
cccc alpha2 controls how many recent charge densities
cccc to be used. If alpha2 > 1, then, the very old
cccc charge density is not used effectivly. 
cccc We find that alpha2=1 is O.K.
******************************************************
      nint_pulay_6=nint_pulay_6+1   
      ! nint_pulay_6, the total number of time pulay has been called.
      ! nint_pulay_6=0 will reset everything
      if(nint_pulay_6.eq.1) num_dwdR_6=0 
      ! initialize num_dwdR_6, which is the number of pairs of dR_6,
      ! dw_6 
      ! nint0 is the index in SCF iteration of each move ({R}). 
c      if(nint0.eq.1) nreset=0
c      nint=nint0-nreset
c      if(nint.gt.npulay_max_6) then
c      write(6,*) "restart pulay, nint0,npulay_max_6",
c     &    nint0,npulay_max_6
c      nreset=nint0-1
c      nint=1
c      endif

      if(nint0.eq.1) then

      do iislda=1,islda
      do i=1,nr_nL
      R0_6(i,iislda)=w_out(i,iislda)-w_in(i,iislda)
      w_in0_6(i,iislda)=w_in(i,iislda)
      enddo
      enddo

      endif
******************************************************
      if(nint0.gt.1) then   ! generates a new dw_6,dR_6 pair
!   num_dwdR_6 is the number of pairs of dw_6,dR_6 ever produced
      num_dwdR_6=num_dwdR_6+1  ! another ever increasing ind
      nposit=mod(num_dwdR_6-1,npulay_max_6)+1
      do iislda=1,islda
      do i=1,nr_nL
      dw_6(i,nposit,iislda)=w_in(i,iislda)-w_in0_6(i,iislda)
      dR_6(i,nposit,iislda)=w_out(i,iislda)-w_in(i,iislda)-
     &                 R0_6(i,iislda)
      R0_6(i,iislda)=w_out(i,iislda)-w_in(i,iislda)
      w_in0_6(i,iislda)=w_in(i,iislda)
      enddo
      enddo

cccccccccccccccoccccccccccccccccccccccc
      if(num_dwdR_6.le.npulay_max_6) then
      num=num_dwdR_6
      else
      num=npulay_max_6
      endif

cccccccccccccccccccccccccccccccccccccccc
c      do m=1,nint-1
c      s=0.d0
c      do iislda=1,islda
c      do i=1,nr_nL
c      s=s+dR_6(i,m,iislda)*R0_6(i,iislda)
c      enddo
c      enddo
c      call global_sumr(s)

c      s=s*vol/nrL
c      B(m)=-s
c      enddo

******************************************************
c      do m=1,nint-1
c      s1=0.d0

c      do iislda=1,islda
c      do i=1,nr_nL
c      s1=s1+dR_6(i,m,iislda)*dR_6(i,nint-1,iislda)
c      enddo
c      enddo

c      call global_sumr(s1)

c      s1=s1*vol/nrL
c      AA(m,nint-1)=s1
c      AA(nint-1,m)=s1
c      enddo

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       tmp(1:num)=0.d0
       fact=vol/nrL
       do iislda=1,islda
       call dgemv('T',nr_nL,num,fact,dR_6(1,1,iislda),
     & mr_nL,dR_6(1,nposit,iislda),1,1.d0,tmp,1)
       enddo

       call mpi_allreduce(tmp,AA_pulay_6(1,nposit),num,
     & MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)   

       do m=1,num
       AA_pulay_6(nposit,m)=AA_pulay_6(m,nposit)
       enddo

       endif    ! end of updating dw_6, dR_6 and AA_pulay_6, nint0.gt.1  


cccccccccc pulay optimization
**********************************************************
      if(num_dwdR_6.gt.0) then    ! using Pulay

      if(num_dwdR_6.le.npulay_max_6) then
      num=num_dwdR_6
      else
      num=npulay_max_6
      endif

      s_0=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      s_0=s_0+R0_6(i,iislda)**2
      enddo
      enddo

      call global_sumr(s_0)

      s_0=s_0*vol/nrL

      tmp(1:num)=0.d0
      fact=-vol/nrL
      do iislda=1,islda
      call dgemV('T',nr_nL,num,fact,dR_6(1,1,iislda),
     & mr_nL,R0_6(1,iislda),1,1.d0,tmp,1)
      enddo  

       call mpi_allreduce(tmp,B,num,
     & MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)   

ccccccccccccccccccccccccccccccccccccccc

      do m1=1,num
      do m2=1,num
      AA1(m1,m2)=AA_pulay_6(m1,m2)
      enddo
      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nposit=mod(num_dwdR_6-1,npulay_max_6)+1
      w=1.0d0
      do m=1,num
      npositx=nposit-m+1
      npositx=mod(npositx-1+2*npulay_max_6,npulay_max_6)+1
      AA1(npositx,npositx)=AA1(npositx,npositx)*w
      w=w*alpha2
      enddo


**********************************************************
      call gaussj(AA1,num,npulay_max_6,B,1,1)

      do iislda=1,islda
      do i=1,nr_nL
      w_out(i,iislda)=R0_6(i,iislda)
      enddo
      enddo

      
c      do m=1,nint-1
c      do iislda=1,islda
c      do i=1,nr_nL
c      w_in(i,iislda)=w_in(i,iislda)+B(m)*dw_6(i,m,iislda)
c      w_out(i,iislda)=w_out(i,iislda)+B(m)*dR_6(i,m,iislda)
c      enddo
c      enddo
c      enddo

       do iislda=1,islda
       call dgemv('N',nr_nL,num,1.d0,dw_6(1,1,iislda),
     & mr_nL,B,1,1.d0,w_in(1,iislda),1)
       call dgemv('N',nr_nL,num,1.d0,dR_6(1,1,iislda),
     & mr_nL,B,1,1.d0,w_out(1,iislda),1)
       enddo    

      s_1=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      w_out(i,iislda)=w_out(i,iislda)+w_in(i,iislda)
      s_1=s_1+(w_out(i,iislda)-w_in(i,iislda))**2
      enddo
      enddo

      call global_sumr(s_1)
      s_1=s_1*vol/nrL

      if(inode_tot.eq.1) then
      write(6,200) s_0,s_1
200   format("mch_pulay,drho_in,drho_out ", 2(E10.4,1x))
      endif

      endif     ! num_dwdR_6.gt.0

      deallocate(AA1)
      deallocate(B)

      return
      end
      

