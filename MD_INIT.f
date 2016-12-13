      subroutine VVMDinit(myid,iscale,xatom,Iseed,DesiredTemp,
     &          Etot,TotalEn,ivlct,ntemp,VX,VY,VZ,V_output)

      use data
      use data_TDDFT

c      implicit double precision (a-h,o-z)
      implicit none

      include 'mpif.h'
      include 'param.escan_real'

      real*8 xatom(3,matom)
      real*8 VX(matom),VY(matom),VZ(matom),V_output(3,matom)
      integer Iseed,i,myid,iscale,ivlct,ntemp
      real*8 mv2,vcx0,vcy0,vcz0,vx0t,vy0t,vz0t
      real*8 RANF,scaling,DesiredTemp,Etot
      real*8 TotalEn,Temp,Enki_Desired,total_mass
      real*8 Enki,Enki_imp
      real*8 HFs2vB2M,Bohr,Hdt
      parameter(HFs2vB2M=27.211396d0)
      parameter(Bohr=0.529177d0)
      parameter(Hdt=0.023538D0/273.15D0/27.211396D0)

      if(ivlct.le.0) then
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ivlct <= 0, initialize velocities with random numbers
       call RANTEST(Iseed)
       do i = 1, ntemp
        VX(i)=RANF(Iseed)-0.5D0
        VY(i)=RANF(Iseed)-0.5D0
        VZ(i)=RANF(Iseed)-0.5D0
       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif

      vcx0 = 0.D0
      vcy0 = 0.D0
      vcz0 = 0.D0
      total_mass=0.d0
      do i = 1, ntemp
       vcx0=vcx0+VX(i)*MDatom(i)
       vcy0=vcy0+VY(i)*MDatom(i)
       vcz0=vcz0+VZ(i)*MDatom(i)
       total_mass=total_mass+MDatom(i)
      enddo
      vcx0 = vcx0/total_mass
      vcy0 = vcy0/total_mass
      vcz0 = vcz0/total_mass
      vx0t = 0.D0
      vy0t = 0.D0
      vz0t = 0.D0
ccccccccccccccccccccccccc
c     modify the velocity center
      mv2 = 0.D0
      do i = 1, ntemp
        VX(i)=VX(i)-vcx0
        VY(i)=VY(i)-vcy0
        VZ(i)=VZ(i)-vcz0
        mv2=mv2+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      scaling = 1.d0
      if(ivlct.le.0) then
cccccccccccccccccccccccccccccccccc
c     use temperature from input file 
c     to modify current temperature
        scaling = dabs(DesiredTemp*Hdt*3.d0*DBLE(ntemp)/mv2)
        mv2 = 0.D0
        do i = 1, ntemp
         VX(i)=VX(i)*sqrt(scaling)
         VY(i)=VY(i)*sqrt(scaling)
         VZ(i)=VZ(i)*sqrt(scaling)
         vx0t=vx0t+VX(i)*MDatom(i)
         vy0t=vy0t+VY(i)*MDatom(i)
         vz0t=vz0t+VZ(i)*MDatom(i)
         mv2=mv2+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
        enddo
cccccccccccccccccccccccccccccccccc
      endif

cccccccccccccccccccccccccccccccccc
c     calc current temperature, then make the output
      Enki=0.5d0*mv2*HFs2vB2M
      Temp=2.d0*Enki/(3.d0*DBLE(ntemp)*Hdt*HFs2vB2M)

      vx0t = vx0t/total_mass
      vy0t = vy0t/total_mass
      vz0t = vz0t/total_mass

      Enki_imp=0.d0
      if(ntemp.lt.natom) then
       do i = ntemp+1, natom
        Enki_imp=Enki_imp+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
       enddo
       Enki_imp=0.5d0*Enki_imp*HFs2vB2M
      endif

      V_output(1,:)=VX(:)
      V_output(2,:)=VY(:)
      V_output(3,:)=VZ(:)

      TotalEn=Etot*HFs2vB2M+Enki+Enki_imp

      if(myid.eq.1) then
        write(6,*) "IniTotEn,Etot,Enki,Enki_imp"
        write(6,*) TotalEn,Etot*HFs2vB2M,Enki,Enki_imp
        write(6,*) "Now Temp=",Temp,"Scaling=",scaling
        write(6,*) "Velocity Center:",vx0t,vy0t,vz0t

        open(11,file='plot_MD.txt',access='append')
        write(11,'(A)') " Time,             E_tot,            E_elec,
     &             E_ion,         E_ion_imp,              Temp,
     &           Scaling"
        write(11,666) 0,TotalEn,Etot*HFs2vB2M,Enki,Enki_imp,
     &   Temp,scaling
        close(11)
666     format(i5,1x,12(f18.9,1x))
      endif

      return
      end



      subroutine VVMD(iscale,Delt,Box,istep,myid,
     &  fatom0,fatom1,AL,xatom,
     &  DesiredTemp,Etot,TotalEn,ntemp,VX,VY,VZ,V_output)

      use data
      use data_TDDFT
c      implicit double precision (a-h,o-z)
      implicit none

      include 'mpif.h'
      include 'param.escan_real'
******************************************
      real*8 xatom(3,matom)
      real*8 fatom0(3,matom),fatom1(3,matom)
      real*8 VX(matom),VY(matom),VZ(matom),V_output(3,matom)
      integer i,myid,istep,iscale,ntemp
      real*8 total_mass
      real*8 mv2,vcx0,vcy0,vcz0,vx0t,vy0t,vz0t,xx,yy,zz
      real*8 scaling,DesiredTemp,Delt,delth,newtemp
      real*8 TotalEn,Etot,Enki,Enki_imp,Enki_old,BOX,AL(3,3)
      real*8 HFs2vB2M,Bohr,Hdt
      parameter(HFs2vB2M=27.211396d0)
      parameter(Bohr=0.529177d0)
      parameter(Hdt=0.023538D0/273.15D0/27.211396D0)

      delth = Delt*0.5d0

      mv2=0.d0
      vcx0=0.d0
      vcy0=0.d0
      vcz0=0.d0
      total_mass=0.d0
      do i=1,natom
       Vx(i)=Vx(i)-delth*fatom0(1,i)/MDatom(i)
       Vy(i)=Vy(i)-delth*fatom0(2,i)/MDatom(i)
       Vz(i)=Vz(i)-delth*fatom0(3,i)/MDatom(i)
      enddo
      do i=1,ntemp
       mv2=mv2+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
       vcx0=vcx0+Vx(i)*MDatom(i)
       vcy0=vcy0+Vy(i)*MDatom(i)
       vcz0=vcz0+Vz(i)*MDatom(i)
       total_mass=total_mass+MDatom(i)
      enddo
      vcx0=vcx0/total_mass
      vcy0=vcy0/total_mass
      vcz0=vcz0/total_mass

      mv2 = 0.D0
      do i = 1, ntemp
       VX(i)=VX(i)-vcx0
       VY(i)=VY(i)-vcy0
       VZ(i)=VZ(i)-vcz0
       mv2=mv2+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
      enddo

      Enki_imp=0.d0
      if(ntemp.lt.natom) then
       do i = ntemp+1, natom
        Enki_imp=Enki_imp+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
       enddo
       Enki_imp=0.5d0*Enki_imp*HFs2vB2M
      endif


      Enki_old=0.5d0*mv2*HFs2vB2M
      temperature=2.d0*Enki_old/(3.d0*DBLE(ntemp)*Hdt*HFs2vB2M)
      Enki=TotalEn-Etot*HFs2vB2M-Enki_imp
      scaling=dabs(Enki/Enki_old)

      if(myid.eq.1) then
       write(6,*) 'Before scaling',istep
       write(6,*) '*********************************'
       write(6,*) 'Enki,Should',Enki_old,Enki
       write(6,*) 'Now Temp=',temperature,'Scaling=',scaling
       write(6,*) 'Velocity Center:',vcx0,vcy0,vcz0
       write(6,*) '*********************************'
       write(6,*)
      endif

      vx0t = 0.D0
      vy0t = 0.D0
      vz0t = 0.D0
      mv2 = 0.D0
      do i = 1, ntemp
       if(iscale.gt.0) then
        VX(i)=VX(i)*sqrt(scaling)
        VY(i)=VY(i)*sqrt(scaling)
        VZ(i)=VZ(i)*sqrt(scaling)
       endif
       vx0t=vx0t+VX(i)*MDatom(i)
       vy0t=vy0t+VY(i)*MDatom(i)
       vz0t=vz0t+VZ(i)*MDatom(i)
       mv2=mv2+MDatom(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
      enddo

      Enki=0.5d0*mv2*HFs2vB2M
      temperature=2.d0*Enki/(3.d0*DBLE(ntemp)*Hdt*HFs2vB2M)

      vx0t = vx0t/total_mass
      vy0t = vy0t/total_mass
      vz0t = vz0t/total_mass

      TotalEn=Etot*HFs2vB2M+Enki+Enki_imp

      V_output(1,:)=VX(:)
      V_output(2,:)=VY(:)
      V_output(3,:)=VZ(:)

      if(myid.eq.1) then
       if(iscale.gt.0) write(6,*) 'After scaling',istep
       if(iscale.le.0) write(6,*) 'No scaling',istep
       write(6,*) '*********************************'
       write(6,*) 'Enki',Enki
       write(6,*) 'Now Temp=',temperature
       write(6,*) 'Velocity Center:',vx0t,vy0t,vz0t
       write(6,*) '*********************************'
       write(6,*)

       open(11,file='plot_MD.txt',access='append')
       write(11,666) istep,TotalEn,Etot*HFs2vB2M,Enki,Enki_imp,
     &   temperature,scaling
        close(11)
666     format(i5,1x,12(f18.9,1x))
      endif

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

      return
      end
