      subroutine MVATOMS(itime,iscale,Delt,Box,xatom,
     & fatom0,fatom1,AL,Etot,DesiredTemp,TotalEn,ivlct,ntemp,
     & VX,VY,VZ,V_output)

ccccccccccccccccccccccccccccccccccccccccccc
      use fft_data
      use load_data
      use data
      use data_TDDFT

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
ccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccc
      integer itime,iscale,ivlct,ntemp
      real*8 xatom(3,matom),fatom0(3,matom),fatom1(3,matom)
      real*8 AL(3,3)
      real*8 DesiredTemp,TotalEn,Etot,Box,Delt
      real*8 VX(matom),VY(matom),VZ(matom),V_output(3,matom)
ccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccc
      real*8 delth
      real*8 xx,yy,zz
      integer Iseed,istep
      parameter(HFs2vB2M=0.930895D0)
      parameter(Hdt=0.023538D0/273.15D0/27.211396D0)
ccccccccccccccccccccccccccccccccccccccccccc

      if(inode_tot.eq.1) then
        write(6,*) "START MD"
      endif

      if(itime.eq.1) then

cccccccccccccccccccccccccccc
c     INIT MD
cccccccccccccccccccccccccccc
        Iseed=12345
        irf = 1
        call VVMDinit(inode_tot,iscale,xatom,Iseed,DesiredTemp,
     &   Etot,TotalEn,ivlct,ntemp,VX,VY,VZ,V_output)

        delth = Delt*0.5d0
        do i=1,natom
         Vx(i)=Vx(i)-delth*fatom1(1,i)/MDatom(i)
         Vy(i)=Vy(i)-delth*fatom1(2,i)/MDatom(i)
         Vz(i)=Vz(i)-delth*fatom1(3,i)/MDatom(i)

         xx=Delt*Vx(i)
         yy=Delt*Vy(i)
         zz=Delt*Vz(i)
         xatom(1,i)=xatom(1,i)+xx*ALI(1,1)+yy*ALI(2,1)+zz*ALI(3,1)
         xatom(2,i)=xatom(2,i)+xx*ALI(1,2)+yy*ALI(2,2)+zz*ALI(3,2)
         xatom(3,i)=xatom(3,i)+xx*ALI(1,3)+yy*ALI(2,3)+zz*ALI(3,3)

c == Periodic boundary conditions,seems can be remove by test
         xatom(1,i)=xatom(1,i)-BOX*int(xatom(1,i)/BOX)
         xatom(2,i)=xatom(2,i)-BOX*int(xatom(2,i)/BOX)
         xatom(3,i)=xatom(3,i)-BOX*int(xatom(3,i)/BOX)
        enddo

        call mpi_barrier(MPI_COMM_WORLD,ierr)
!     INITIAL DONE
ccccccccccccccccccccccccccccccccccccccc

      elseif(itime.gt.1) then

ccccccccccccccccccccccccccccccccccccccc
        istep=itime-1
        call VVMD(iscale,Delt,BOX,istep,inode_tot,
     &   fatom0,fatom1,AL,xatom,
     &   DesiredTemp,Etot,TotalEn,ntemp,VX,VY,VZ,V_output)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
ccccccccccccccccccccccccccccccccccccccc

      endif

      call MPI_barrier(MPI_COMM_WORLD,ierr)

      return
      end
