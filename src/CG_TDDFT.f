      subroutine CG_TDDFT(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt,iislda)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001. 
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
****************************************
****   The main conjugate gradient routine
****   mbad=0, normal iCG=1 run, ave_err not used
****   mbad>0, called from CG_new, fix the states: mx-mbad,mx to ave_err
******************************************
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
***********************************************
       complex*16 pg(mg_nx),ugh(mg_nx),pgh(mg_nx)
       complex*16 pg_old(mg_nx),ughh_old(mg_nx)
       complex*16,allocatable,dimension(:)  ::  spg

       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       integer isNLa(9,matom),ipsp_type(mtype)

       real*8 vr(mr_n),vr0(mr_n),vr_pert(mr_n)
       real*8 prec(mg_nx)
       complex*16 ug_n_1(mg_nx),ug_n_2(mg_nx)
       complex*16 workr_n(mr_n)
       complex*16 cc,cc_store(mx)
**********************************************
**** if iopt=0 is used, pghh_old can be deleted
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst),err2_st(mst)
	real*8 Ef
        complex*16 sumdum(mref,matom),sumdum2(mref,matom),cai
        real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)

       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comNL2/occ_t,iiatom,icore,numref,ityatom


       cai=dcmplx(0.d0,1.d0)
       pi=4*datan(1.d0)

       ng_n=ngtotnod(inode,kpt)

       if(ipsp_all.eq.1) then
       allocate(spg(1))
       else
       allocate(spg(mg_nx))
       endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc first, prepare the initial states
 

      vr0=vr     ! need to define vr0
      do i=1,nr_n
      jj=i+(inode-1)*nr_n
      i1=(jj-1)/(n2*n3)
      j1=(jj-1-(i1-1)*n2*n3)/n3+1
      k1=jj-(i1-1)*n2*n3-(j1-1)*n3
      x=(i1-1)*1.d0/n1
      y=(j1-1)*1.d0/n2
      z=(k1-1)*1.d0/n3
      vr_pert(i)=3*sin((x-0.2)*2*pi)/27.211396         ! 3 eV amplitude oscillations

      workr_n(i)=dexp(-((x-0.5)**2+(y-0.5)**2+(x-0.5)**2)/0.2**2)

      enddo

      call d3fft_comp(ug_n_2,workr_n,1,kpt)


      ug_n_1=dcmplx(0.d0,0.d0)
 
      m_init=16
      m_end=26
      cc_store=dcmplx(0.d0,0.d0)

      do m=m_init,m_end
      cc=dcmplx(0.d0,0.d0)
      do i=1,ng_n
      cc=cc+ug_n_2(i)*dconjg(ug_n_bp(i,m))
      enddo
      
      call global_sumc(cc)

      cc=cc*vol

      cc_store(m)=cc
      do i=1,ng_n
      ug_n_1(i)=ug_n_1(i)+cc*ug_n_bp(i,m)
      enddo
      enddo
cccccccccccccccccccccccc

      sum=0.d0
      do i=1,ng_n
      sum=sum+abs(ug_n_1(i))**2
      enddo

      call global_sumr(sum)

      fac=1.d0/dsqrt(sum*vol)
      ug_n_1=ug_n_1*fac
      cc_store=cc_store*fac
ccccccccccccccccccccccccccccccc


      hbar=0.65822d0         ! fs*eV

      hbar2=hbar/27.211396   ! fs*hartree
c      ntime=100000
c      dt=1.d0/100000           ! unit fs
      ntime=200000
      dt=1.d0/200000           ! unit fs
      Ephoton=1.d0           ! eV
      omega=Ephoton/hbar     ! 1/fs
      dt_fac=dt/hbar2

      if(inode.eq.1) then
      open(99,file="TD_DFT.out")
      rewind(99)
      open(98,file="TD_DFT.out2")
      rewind(98)
      endif 

      do 3000 itime=1,ntime

cccccccccccccccccccccccccccccccccccccccccccccccccc
      time=(itime-1)*dt
      time1=(itime+0.5-1)*dt

      if(mod(itime-1,1000).eq.0) then
      sum=0.d0
      do m=1,mx
      cc=dcmplx(0.d0,0.d0)
      do i=1,ng_n
      cc=cc+ug_n_1(i)*dconjg(ug_n_bp(i,m))
      enddo
      
      call global_sumc(cc)
      cc=cc*vol
      cc_store(m)=cc
      sum=sum+cdabs(cc_store(m))**2
      enddo

      if(inode.eq.1) then
      write(98,*) itime,time,sum
      write(98,202) (abs(cc_store(m))**2,m=1,mx)
      write(98,*) "*************"
      write(98,202) (real(cc_store(m))/abs(cc_store(m)),m=1,mx)
      call flush(98)
      endif
202   format(5(E11.4,1x))
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccc


      vr=vr0+vr_pert*sin(omega*time1) ! something like that

        call Hpsi_comp_tmp(ug_n_1,ugh,ilocal,vr,workr_n,kpt,
     &   1,sug_n(1,m),sumdum,iislda,ave_X,dX)

      sum=0.d0
      do i=1,ng_n
      sum=sum+real(ugh(i)*dconjg(ug_n_1(i)))
      enddo

      call global_sumr(sum)

      E=sum*vol
      
      ug_n_1=ug_n_1-cai*ugh*dt_fac

      sum=0.d0
      do i=1,ng_n
      sum=sum+abs(ug_n_1(i))**2
      enddo

      call global_sumr(sum)

      fact=1/dsqrt(sum*vol)
      ug_n_1=ug_n_1*fact

      if(inode.eq.1.and.mod(itime-1,100).eq.0) then
      write(99,200) itime,time,E,ave_X,dX,fact
200   format(i8,5(E13.6,1x))
      call flush(99)
      endif

3000  continue

      if(inode.eq.1) then
      close(99)
      endif

      call MPI_barrier(MPI_COMM_WORLD,ierr)
      call MPI_abort(MPI_COMM_WORLD,ierr)

***********************************************

      return
      end

