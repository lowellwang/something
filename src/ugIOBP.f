      subroutine ugIOBP(ug_n_tmp,kpt,iflag,istep,iislda,
     &   iflag2,nkpt,islda)
******************************************
c     c     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
*     *  copyright (c) 2003, The Regents of the University of California,
*     *  through Lawrence Berkeley National Laboratory (subject to receipt of any
*     *  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
cccccc     is all the icolor groups are writing or reading, 
cccccc     or only the group with kpt inside the [kpt_dis(1),kpt_dis(2)] are writing or reading ?

******************************************

****************************************
ccccc iflag=1, write, iflag=2, read
****************************************
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)

      complex*16 ug_n_tmp(mg_nx,nblock_band_mx)
      character*8 fname
*************************************************
      if(nkpt.eq.1.and.islda.eq.1.and.istep.eq.0) return
ccccc don't do anything in this case
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iflag2.lt.0) then         ! read ug according to icolor_k and kpt relationship
       if((iislda-1)*nkpt+kpt.lt.kpt_slda_dis(1).or.
     &    (iislda-1)*nkpt+kpt.gt.kpt_slda_dis(2)) return
       endif 

       if(iflag2.ge.0) then      ! read ug only for icolor_k.eq.iflag
       if(icolor_k.ne.iflag2) return
       endif 
      
*************************************************
      if(iflag.eq.2) then       ! read the wavefunction

        ug_n_tmp(:,:)=ug_all(:,:,kpt,iislda)

      endif                     ! for the iflag, read
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(iflag.eq.1) then       ! write the wavefunction

         ug_all(:,:,kpt,iislda)=ug_n_tmp(:,:)

      endif                     ! for the iflag, write
      
      return
      end

