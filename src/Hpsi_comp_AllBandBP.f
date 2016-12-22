      subroutine Hpsi_comp_AllBandBP(wg,wgh,nblock,
     &  ilocal,vr,workr_n,kpt,isij,swg,sumdum_m,iislda)

******************************************
c     c     Written by Byounghak Lee, Jan 28, 2008.
*     *  copyright (c) 2003, The Regents of the University of California,
*     *  through Lawrence Berkeley National Laboratory (subject to receipt of any
*     *  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

      use fft_data
      use load_data
      use data
!!!!!!!!!!!!!!!!!!!
      use data_TDDFT
!!!!!!!!!!!!!!!!!!!

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      parameter (mem_size=500000000)     ! 500 MB extra memory for fft

      real*8 vr(mr_n)
      complex*16 workr_n(mr_n)
      complex*16 wg(mg_nx,nblock), wgh(mg_nx,nblock)
      complex*16 swg(mg_nx,nblock)
      complex*16 wg_dum(mg_nx)
      complex*16 workr2_n(mr_n,nblock)
      complex*16, allocatable, dimension(:,:) :: workr3_n

      integer nmap(matom)

      complex*16 sumdum_m(mref,natom,nblock)
      complex*16 sumdumtmp_m(mref,natom,nblock)

      complex*16 sumdum(mref,natom),sumdum2(mref),sumdum3(mref)
      complex*16 sumdumtmp(mref,natom)

      complex*16 s,cc,s2,cc2

      integer isNLa(9,matom),ityatom(matom),ipsp_type(mtype)
      complex*16 cy

      real*8 occ_t(mtype)
      integer iiatom(mtype),icore(mtype),numref(matom)
      real*8 Dij0(32,32,mtype),Qij(32,32,mtype)

*************************************************
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
      common /comnmap/nmap
      common /comNL2/occ_t,iiatom,icore,numref,ityatom

      mblock=mem_size/(mr_n*16*3)              !  3 for three extra buffers inside fft
      if(mblock.gt.nblock) mblock=nblock
      mloop=nblock*1.d0/mblock+0.9999999999d0

      ng_n=ngtotnod(inode,kpt)

cccccccccccccccccccccccccccccccccccccccccccccccc
      do iloop=1,mloop
      mblock_tmp=mblock
      if(iloop.eq.mloop)
     &     mblock_tmp=nblock-(iloop-1)*mblock
      ip=(iloop-1)*mblock+1
      
      call d3fft_comp_block(wg(1,ip),workr2_n(1,ip),
     &  -1,kpt,mblock_tmp)
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccc


      if (ilocal==2) then
         call nonlocal_realsp_m()
      else             
         do m=1, nblock
            do i=1,nr_n
c               workr2_n(i,m)=(vr(i)+vext_dt_n(i))*workr2_n(i,m)
               workr2_n(i,m)=vr(i)*workr2_n(i,m)
            enddo      
         enddo   
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccc
      do iloop=1,mloop
      mblock_tmp=mblock
      if(iloop.eq.mloop)
     &     mblock_tmp=nblock-(iloop-1)*mblock
      ip=(iloop-1)*mblock+1
      call d3fft_comp_block(wgh(1,ip),workr2_n(1,ip),
     &     1,kpt,mblock_tmp)
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc

      
      do m=1, nblock
         do i=1,ng_n
            wgh(i,m)=wgh(i,m)+gkk_n(i,kpt)*wg(i,m)
         enddo
      enddo         

      if(ilocal.eq.3) then     
         call nonlocal_qsp_m()
      endif    

      return

      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nonlocal_qsp_m()
      implicit double precision (a-h,o-z)
      complex*16 cc1,cc0,sumy_m(nref_tot,nblock)
      complex*16, allocatable, dimension(:,:) :: sumy3_m


      if(isij.eq.1.and.ipsp_all.eq.2) then 
cccccc could need large memory 
      allocate(sumy3_m(nref_tot,nblock))
      endif

      cc1=dcmplx(1.d0,0.d0)
      cc0=dcmplx(0.d0,0.d0)

      sumdum_m = cc0

      call zgemm('c','n',nref_tot,nblock,ng_n,cc1,wqmask,
     &     mg_nx,wg,mg_nx,cc0,sumy_m,nref_tot)

      do ia=1,natom
         nref=numref(ia)
         iref_start2=iref_start(ia)
         
         if(ng_n.gt.0) then
            do m1=1, nblock         
               do jj=1,nref          
                  sumdum_m(jj,ia,m1)=dconjg(sumy_m(jj+iref_start2,m1))
               enddo
            enddo
         endif
      enddo

      call mpi_allreduce(sumdum_m,sumdumtmp_m,natom*mref*nblock,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
      sumdum_m = sumdumtmp_m*vol
      

      do ia=1,natom
         nref=numref(ia)
         iref_start2=iref_start(ia)
         iitype=ityatom(ia)

         do m1=1, nblock         
            do jj1=1,nref               

            if(ipsp_type(iitype).eq.2) then
            cc=dcmplx(0.d0,0.d0)
            do jj2=1,nref
            cc=cc+Dij(jj2,jj1,ia,iislda)*sumdum_m(jj2,ia,m1)
            enddo
            sumy_m(jj1+iref_start2,m1)=dconjg(cc)
            else
            sumy_m(jj1+iref_start2,m1)=Dij(jj1,jj1,ia,iislda)*
     &              dconjg(sumdum_m(jj1,ia,m1))
            endif

            if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
            cc=0.d0
            do jj2=1,nref
            cc=cc+Qij(jj2,jj1,iitype)*sumdum_m(jj2,ia,m1)
            enddo
            sumy3_m(jj1+iref_start2,m1)=dconjg(cc)
            endif

            enddo               
         enddo
      enddo   

cccccccccccccccccccccccccccccccccccccccccccccccccc

      
      call zgemm('n','n',ng_n,nblock,nref_tot,cc1,
     &     wqmask,mg_nx,sumy_m,nref_tot,cc1,wgh,mg_nx)                 

      
      if(isij.eq.1.and.ipsp_all.eq.2) then 

      swg=wg
      call zgemm('n','n',ng_n,nblock,nref_tot,cc1,
     &     wqmask,mg_nx,sumy3_m,nref_tot,cc1,swg,mg_nx)                 


      deallocate(sumy3_m)
      endif

      return
      end subroutine nonlocal_qsp_m
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine nonlocal_realsp_m()
      implicit double precision (a-h,o-z)
      complex*16 y,cc1,cc0
      parameter (nmap_max=13000)
      real*8 y_tmp1_m(nmap_max,nblock),y_tmp2_m(nmap_max,nblock)
      real*8 y_tmp31_m(nmap_max,nblock),y_tmp32_m(nmap_max,nblock)
      real*8 sumy1_m(mref,nblock),sumy2_m(mref,nblock)
      real*8 sumy31_m(mref,nblock),sumy32_m(mref,nblock)
      integer ico_m(nblock)

      cc1=dcmplx(1.d0,0.d0)
      cc0=dcmplx(0.d0,0.d0)

      sumdum = dcmplx(0.d0,0.d0)
      sumdum_m = dcmplx(0.d0,0.d0)
      ico0=0
      ico1=0
      ico2=0



      do ia=1,natom
       nref=numref(ia)
       if(nref.eq.0.or.nmap(ia).eq.0) goto 999
         
         if(nmap(ia).gt.nmap_max) then
      write(6,*) "nmap.gt.nmap_max, stop", nmap(ia),nmap_max,ia
            call  mpi_abort(MPI_COMM_WORLD,ierr)
         endif
         
         ico_m = ico0
         do m=1, nblock   
            do i=1,nmap(ia)
               ico_m(m)=ico_m(m)+1               
               y=dconjg(workr2_n(indm(ico_m(m)),m)*cphase(ico_m(m)))
               y_tmp1_m(i,m)=dreal(y)
               y_tmp2_m(i,m)=dimag(y)
            enddo            
         enddo
         ico0 = ico_m(1)


         call dgemm('N','N',nref,nblock,nmap(ia),1.d0,
     &        wmask(ico2+1),nref,y_tmp1_m,nmap_max,0.d0,sumy1_m,mref)

         
         call dgemm('N','N',nref,nblock,nmap(ia),1.d0,
     &        wmask(ico2+1),nref,y_tmp2_m,nmap_max,0.d0,sumy2_m,mref)


         if(nmap(ia).gt.0) then 
            do m=1, nblock
               do jj=1,nref
                  sumdum_m(jj,ia,m)=dcmplx(sumy1_m(jj,m),sumy2_m(jj,m))
               enddo
            enddo
            ico2=ico2+nmap(ia)*nref
         endif

999   continue
      enddo    ! ia=1,natom


      call mpi_allreduce(sumdum_m,sumdumtmp_m,natom*mref*nblock,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_B1,ierr)
      sumdum_m = sumdumtmp_m*vol/(n1*n2*n3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(isij.eq.1.and.ipsp_all.eq.2) then
      allocate(workr3_n(mr_n,nblock))
      workr3_n=workr2_n
      endif

      do m=1, nblock
         do i=1,nr_n
            workr2_n(i,m)=dble(vr(i))*workr2_n(i,m)
         enddo
      enddo
      

      ico1=0
      ico2=0
      do ia=1,natom
         nref=numref(ia)

         if(nref.eq.0.or.nmap(ia).eq.0) goto 998

         iitype=ityatom(ia)

        if(ipsp_type(iitype).eq.2) then
         do m=1,nblock
         do jj1=1,nref
         cc=dcmplx(0.d0,0.d0) 
         do jj2=1,nref
         cc=cc+Dij(jj2,jj1,ia,iislda)*sumdum_m(jj2,ia,m)
         enddo
         sumy1_m(jj1,m)=dreal(cc)
         sumy2_m(jj1,m)=-dimag(cc)
         enddo
         enddo
         else        ! ipsp_type(iitype).eq.2
         do m=1, nblock
            do jj=1,nref
               y=Dij(jj,jj,ia,iislda)*dconjg(sumdum_m(jj,ia,m))
               sumy1_m(jj,m)=dreal(y)
               sumy2_m(jj,m)=dimag(y)              
            enddo
         enddo         
         endif
ccccccccccccccccccccccccccccccccccc

        if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
         do m=1,nblock
         do jj1=1,nref
         cc2=dcmplx(0.d0,0.d0)
         do jj2=1,nref
         cc2=cc2+Qij(jj2,jj1,iitype)*sumdum_m(jj2,ia,m)
         enddo
         sumy31_m(jj1,m)=dreal(cc2)
         sumy32_m(jj1,m)=-dimag(cc2)
         enddo
         enddo
         endif
cccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccc
         
         call dgemm('T','N',nmap(ia),nblock,nref,1.d0,
     &    wmask(ico2+1),nref,sumy1_m,mref,0.d0,y_tmp1_m,nmap_max)

         call dgemm('T','N',nmap(ia),nblock,nref,1.d0,
     &    wmask(ico2+1),nref,sumy2_m,mref,0.d0,y_tmp2_m,nmap_max)

        if(isij.eq.1.and.ipsp_type(iitype).eq.2) then

         call dgemm('T','N',nmap(ia),nblock,nref,1.d0,
     &    wmask(ico2+1),nref,sumy31_m,mref,0.d0,y_tmp31_m,nmap_max)

         call dgemm('T','N',nmap(ia),nblock,nref,1.d0,
     &   wmask(ico2+1),nref,sumy32_m,mref,0.d0,y_tmp32_m,nmap_max)
         endif


        if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
         ico2=ico2+nmap(ia)*nref
         do i=1,nmap(ia)
            ico1=ico1+1
            do m=1, nblock
               workr2_n(indm(ico1),m)=workr2_n(indm(ico1),m)+
     &              dcmplx(y_tmp1_m(i,m),y_tmp2_m(i,m))
     &              *dconjg(cphase(ico1))
               workr3_n(indm(ico1),m)=workr3_n(indm(ico1),m)+
     &              dcmplx(y_tmp31_m(i,m),y_tmp32_m(i,m))
     &              *dconjg(cphase(ico1))
            enddo   
         enddo

        else
         ico2=ico2+nmap(ia)*nref
         do i=1,nmap(ia)
            ico1=ico1+1
            do m=1, nblock
               workr2_n(indm(ico1),m)=workr2_n(indm(ico1),m)+
     &              dcmplx(y_tmp1_m(i,m),y_tmp2_m(i,m))
     &              *dconjg(cphase(ico1))
            enddo   
         enddo
         endif

998   continue
      enddo   ! ia=1,natom  


      if(isij.eq.1.and.ipsp_all.eq.2) then
      call d3fft_comp_block(swg,workr3_n,1,kpt,nblock)
      deallocate(workr3_n)
      endif

     
      return
      end subroutine nonlocal_realsp_m



      end
