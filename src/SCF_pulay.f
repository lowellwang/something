      subroutine SCF_pulay(w_r_in,w_r_out,nrint0,RR,nrreset,islda,
     &                     w_r_in0,Rr0,dwr,dRr,npulay_r_max)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
ccccc  Now, this is in n1L,n2L,n3L


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 w_r_in(mr_nL,islda),w_r_out(mr_nL,islda)
      real*8 w_r_in0(mr_nL,islda),Rr0(mr_nL,islda)
      real*8 dwr(mr_nL,npulay_r_max,islda)
      real*8 dRr(mr_nL,npulay_r_max,islda)
      real*8 RR(npulay_r_max,npulay_r_max)
      real*8, allocatable, dimension(:,:)  :: RR1
      real*8, allocatable, dimension(:)  :: B
      integer nrreset,nrint,nrint0,npulay_r_max

      allocate(RR1(npulay_r_max,npulay_r_max))
      allocate(B(npulay_r_max))

      alpha2=1.d0
cccc alpha2 controls how many recent charge densities
cccc to be used. If alpha2 > 1, then, the very old
cccc charge density is not used effectivly. 
cccc We find that alpha2=1 is O.K.
******************************************************
      if(nrint0.eq.1) nrreset=0
      nrint=nrint0-nrreset

      if(nrint.gt.npulay_r_max) then
      write(6,*) "restart pulay, nrint0,npulay_r_max",
     &    nrint0,npulay_r_max
      nrreset=nrint0-1
      nrint=1
      endif

      if(nrint.eq.1) then

      do iislda=1,islda
      do i=1,nr_nL
      Rr0(i,iislda)=w_r_out(i,iislda)-w_r_in(i,iislda)
      w_r_in0(i,iislda)=w_r_in(i,iislda)
      enddo
      enddo

      endif
******************************************************
      if(nrint.gt.1) then

      s_0=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      dwr(i,nrint-1,iislda)=w_r_in(i,iislda)-w_r_in0(i,iislda)
      dRr(i,nrint-1,iislda)=w_r_out(i,iislda)-w_r_in(i,iislda)-
     &                 Rr0(i,iislda)
      Rr0(i,iislda)=w_r_out(i,iislda)-w_r_in(i,iislda)
      w_r_in0(i,iislda)=w_r_in(i,iislda)
      s_0=s_0+Rr0(i,iislda)**2
      enddo
      enddo

      if(inode.eq.1) write(6,*) "DEBUG00"

      call global_sumr(s_0)

      if(inode.eq.1) write(6,*) "DEBUG01"

      s_0=s_0*vol/nrL


      do m=1,nrint-1
      s=0.d0

      do iislda=1,islda
      do i=1,nr_nL
      s=s+dRr(i,m,iislda)*Rr0(i,iislda)
      enddo
      enddo

      call global_sumr(s)

      if(inode.eq.1) write(6,*) "DEBUG02"

      s=s*vol/nrL
      B(m)=-s
      enddo

******************************************************
      do m=1,nrint-1
      s1=0.d0

      do iislda=1,islda
      do i=1,nr_nL
      s1=s1+dRr(i,m,iislda)*dRr(i,nrint-1,iislda)
      enddo
      enddo

      call global_sumr(s1)

      if(inode.eq.1) write(6,*) "DEBUG03"

      s1=s1*vol/nrL
      RR(m,nrint-1)=s1
      RR(nrint-1,m)=s1
      enddo

cccccccccc pulay optimization
**********************************************************
      do m1=1,nrint-1
      do m2=1,nrint-1
      RR1(m1,m2)=RR(m1,m2)
      enddo
      enddo

      w=1.0d0
      do m=nrint-1,1,-1
      RR1(m,m)=RR1(m,m)*w
      w=w*alpha2
      enddo

**********************************************************

      call gaussj(RR1,nrint-1,npulay_r_max,B,1,1)

      if(inode.eq.1) write(6,*) "DEBUG04"

      do iislda=1,islda
      do i=1,nr_nL
      w_r_out(i,iislda)=Rr0(i,iislda)
      enddo
      enddo

      
      do m=1,nrint-1
      do iislda=1,islda
      do i=1,nr_nL
      w_r_in(i,iislda)=w_r_in(i,iislda)+B(m)*dwr(i,m,iislda)
      w_r_out(i,iislda)=w_r_out(i,iislda)+B(m)*dRr(i,m,iislda)
      enddo
      enddo
      enddo

      s_1=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      w_r_out(i,iislda)=w_r_out(i,iislda)+w_r_in(i,iislda)
      s_1=s_1+(w_r_out(i,iislda)-w_r_in(i,iislda))**2
      enddo
      enddo

      call global_sumr(s_1)
      s_1=s_1*vol/nrL

      if(inode_tot.eq.1) then
      write(6,*) "mch_pulay,nrint,dv_in,dv_out", nrint,s_0,s_1
      endif

      endif

      deallocate(RR1)
      deallocate(B)

      return
      end

