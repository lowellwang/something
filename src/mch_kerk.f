      subroutine mch_kerk(w_in,w_out,amx_mth,alpha,beta,kerk_min)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
ccc   now, in n1L,n2L,n3L 


      use fft_data
      use load_data
      use data

      ! implicit double precision (a-h,o-z)
      implicit none
      include 'param.escan_real'

      real*8 w_in(mr_nL),w_out(mr_nL)

*************************************************************
***   zhiwang2016, replace the kerker mixing factor 
***   by passed in parameters
***
***   kerk(q) = alpha * max[ q^2 / (q^2 + beta), kerk_min ]
***
***   the original value here is alpha=0.8, beta=0.5
*************************************************************
      real*8 amx_mth,alpha,beta,kerk_min

      integer i,ng2_nL
      real*8 kerk_fac

      real*8, allocatable,dimension (:)  :: workr_n,vi

      allocate(workr_n(mr_nL))
      allocate(vi(mr_nL))

      kerk_fac=0.d0
      ng2_nL=ngtotnod2L(inode)
********************************************************
*** Kerker mixing, this is like the precoditioning,
*** and it adds the new component w_out-w_in into w_in
********************************************************
      do i=1,nr_nL
      w_out(i)=w_out(i)-w_in(i)
      vi(i)=0.d0
      enddo

      workr_n = w_out

      call d3fft_real2L(vi,workr_n,1,0)

      do i=1,ng2_nL

********** in SLDA, v(q=0) is not zero, but has the information of
********** relative highs between up and down, so this point should be
********** treated carefully.

       if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then
         kerk_fac=3.d0/8.d0
         write(6,*) "kerk, q=0",inode,i,vi(i*2),vi(i*2-1)
       else
         kerk_fac=max(gkk2_nL(i)/(gkk2_nL(i)+beta),kerk_min)
         !write(6,*) "kerk, q!=0",inode,i,kerk_fac
       endif
       vi(i*2)  =alpha*(kerk_fac-1.d0)*vi(i*2)
       vi(i*2-1)=alpha*(kerk_fac-1.d0)*vi(i*2-1)

cccccccc   the -0.8d0 is used to cancel the 0.8*w_out in real space, but left with 
cccccccc   the gkk components outside Ecut2L
      enddo

      call d3fft_real2L(vi,workr_n,-1,0)

cccccc  The 0.8*w_out cancels the 0.8*vi(i) except for gkk output Ecut2L

      do i=1,nr_nL
      w_in(i)=w_in(i)+workr_n(i)+alpha*w_out(i)
      enddo
**************************************************
      deallocate(workr_n)
      deallocate(vi)
      return
      end
      

