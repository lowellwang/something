SUBROUTINE input_TDDFT(tot, m, mmn0, iwg_in, imax, rkappa, ivlct, ntemp, VX, VY, VZ, &
                   finp_td, ftmp_st, fdx_td, t_timewall)

  ! Read TDDFT parameter from file finp_td 
  ! If needed, find correct ftmp_st then read it too

  USE data_TDDFT
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'

  INTEGER :: m,mmn0,ierr,ierr2,iline,imax,i,j,k,l,l2,iwg_in,ivlct,ntemp,itmp
  REAL(8) :: rtmp(10),tot,rkappa(3),t_timewall
  REAL(8), DIMENSION(natom) :: VX, VY, VZ
  LOGICAL :: bmass,bkappa
  CHARACTER(20) :: finp_td,ftmp_st,fdx_td
  CHARACTER(50) :: filename,fileindex,ctmp
  CHARACTER(255) :: line

  iline = 0
  bmass = .false.
  bkappa = .true.

  mst2 = -1
  mmn = -1
  mmn0 = 1
  dtMD = -1.d0
  iMD = -1
  ibo_md = -1
  temperature = -1.d0
  ntemp = natom
  mscf = -1
  tolrho = -1.d0
  tolintgl = -1.d0
  ntime_init = -1
  init_time = -1.d0
  imp_k = 0.d0
  ntime = -1
  noutput = -1
  ifatom = -1
  ivext_dt = -1
  ikappa = 0
  rkappa = 0.d0
  iboltz = -1
  ibshift = -1
  i_scale = -1
  iocc = -2
  jhole = -1
  jelec = -1
  jxi = -1
  rxi = -1.d0
  ivlct = 0
  t_timewall = -1.d0

  open(14,file=finp_td,status='old',action='read',iostat=ierr)
  if(inode_tot.eq.1.and.ierr.ne.0) then
   write(6,*) "--ERROR: Cannot open "//trim(adjustl(finp_td))
   call MPI_abort(MPI_COMM_WORLD,ierr)
   stop
  endif

10    continue

  iline=iline+1
  read(14,'(a)',iostat=ierr,end=666,err=666) line

  if(line(1:1).eq."!") goto 10

  l=len_trim(line)
  l2=index(line,'!')
  i=index(line,'=')

  if(i.ne.0) then
    if(l2.ne.0) then
      if(l2.le.i) then
        goto 10
      else
        read(line(i+1:l2-1),'(a)') ctmp
      endif
    else
      read(line(i+1:l),'(a)') ctmp
    endif
  else
    goto 10
  endif

  if(index(line,'nbasis').ne.0) then
    read(ctmp,*) mst2

  elseif(index(line,'nbands').ne.0) then
    read(ctmp,*) mmn

  elseif(index(line,'nband0').ne.0) then
    read(ctmp,*) mmn0

  elseif(index(line,'nstep').ne.0) then
    read(ctmp,*) ntime

  elseif(index(line,'istep').ne.0) then
    read(ctmp,*) ntime_init

  elseif(index(line,'nelm').ne.0) then
    read(ctmp,*) mscf

  elseif(index(line,'iexci').ne.0) then
    read(ctmp,*) iocc

  elseif(index(line,'hole').ne.0) then
    read(ctmp,*) jhole

  elseif(index(line,'elect').ne.0) then
    read(ctmp,*) jelec

  elseif(index(line,'nexci').ne.0) then
    read(ctmp,*) jxi

  elseif(index(line,'rexci').ne.0) then
    read(ctmp,*) rxi

  elseif(index(line,'nwrite').ne.0) then
    read(ctmp,*) noutput

  elseif(index(line,'iforce').ne.0) then
    read(ctmp,*) ifatom

  elseif(index(line,'rkappa').ne.0) then
    iline=iline+1
    bkappa=.false.
    read(14,'(a)',iostat=ierr2,end=666,err=666) line
    read(line,*,iostat=ierr2,end=666,err=666) (rkappa(j),j=1,3)
    bkappa=.true.

  elseif(index(line,'ivdt').ne.0) then
    read(ctmp,*) ivext_dt

  elseif(index(line,'iboltz').ne.0) then
    read(ctmp,*) iboltz

  elseif(index(line,'ibandshift').ne.0) then
    read(ctmp,*) ibshift

  elseif(index(line,'iscaling').ne.0) then
    read(ctmp,*) i_scale

  elseif(index(line,'rhodiff').ne.0) then
    read(ctmp,*) tolrho

  elseif(index(line,'intglerr').ne.0) then
    read(ctmp,*) tolintgl

  elseif(index(line,'starttime').ne.0) then
    read(ctmp,*) init_time

  elseif(index(line,'kinetic').ne.0) then
    read(ctmp,*) imp_k

  elseif(index(line,'imd').ne.0) then
    read(ctmp,*) iMD

  elseif(index(line,'ibo_md').ne.0) then
    read(ctmp,*) ibo_md

  elseif(index(line,'ntmprtr').ne.0) then
    read(ctmp,*) ntemp

  elseif(index(line,'tmprtr').ne.0) then
    read(ctmp,*) temperature

  elseif(index(line,'mdmass').ne.0) then
    iline=iline+1
    read(14,'(a)',iostat=ierr2,end=666,err=666) line
    ierr2=1
    imax=11
    do while(ierr2.ne.0.and.imax.gt.1)
      imax=imax-1
      read(line,*,iostat=ierr2) (rtmp(j),j=1,imax)
    enddo
    if(ierr2.ne.0) goto 666
    if(.not.bmass) then
     bmass=.true.
    else
     goto 666
    endif
    MDtype(1:imax)=rtmp(1:imax)

  elseif(index(line,'dt').ne.0) then
    read(ctmp,*) dtMD

  elseif(index(line,'velocity').ne.0) then
    ivlct = 1
    if(inode_tot.eq.1) &
      write(6,*) "Read velocity from "//trim(adjustl(finp_td))
    do i = 1, natom
      iline=iline+1
      read(14,*,iostat=ierr2,end=666,err=666) itmp, VX(i), VY(i), VZ(i)
      if(inode_tot.eq.1) &
        write(6,'(i4,2x,3(E15.8,1x))') itmp, VX(i), VY(i), VZ(i)
    end do

  elseif(index(line,'dxatom_file').ne.0) then
    ifatom = 2
    read(ctmp,*) fdx_td
    if(inode_tot.eq.1) &
      write(6,*) "Read dxatom from "//trim(adjustl(fdx_td))

  elseif(index(line,'timewall').ne.0) then
    read(ctmp,*) t_timewall

  elseif(index(line,'=').ne.0.and.inode_tot.eq.1) then
    write(6,'(a,i5,a)') " --WARNING: Line ",iline," in "// &
                            trim(adjustl(finp_td))// &
        " may contain some values, but cannot be recognized"

  endif

  goto 10

  iline=0
666   continue
  close(14)

  if(ierr2.ne.0.and.inode_tot.eq.1.and.iline.ne.0) then
   write(6,*) "--ERROR: line",iline
   call MPI_abort(MPI_COMM_WORLD,ierr)
   stop

  elseif(ierr.ge.0.and.inode_tot.eq.1.and.iline.ne.0) then
   write(6,*) "--ERROR: line",iline
   call MPI_abort(MPI_COMM_WORLD,ierr)
   stop

  elseif(inode_tot.eq.1.and.(.not.bmass)) then
   write(6,*) "--ERROR: Cannot read nuclei masses"
   call MPI_abort(MPI_COMM_WORLD,ierr)
   stop

  elseif(inode_tot.eq.1.and.(.not.bkappa)) then
   write(6,*) "--ERROR: Cannot read rkappa in line",iline
   call MPI_abort(MPI_COMM_WORLD,ierr)
   stop

  endif

  if(mst2.lt.0) mst2=m
  if(mmn.lt.0.and.mst2.lt.0) mmn=m
  if(mmn.lt.0.and.mst2.gt.0) mmn=mst2

  if(mmn0.le.0.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: nband0 <= 0"
    call MPI_abort(MPI_COMM_WORLD,ierr)
  elseif(mmn0.gt.mmn.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: nband0 > nbands"
    call MPI_abort(MPI_COMM_WORLD,ierr)
  endif

  if(dtMD.le.0.d0.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: Cannot read dt, or dt = 0.0"
    call MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  endif

  if(ntime_init.lt.0.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: Cannot read istep"
    call MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  elseif(ntime_init.gt.0.and.inode_tot.eq.1.and.iwg_in.eq.0) then
    write(6,*) "--ERROR: istep.ne.0 while iwg_in.eq.0, no &
                wavefunctions have been read"
    call MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  endif

  if(ibo_md.gt.0.and.iocc.gt.0.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: ibo_md > 0 while iexci > 0"
    call MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  endif

  if(ntime.le.0.and.inode_tot.eq.1) then
    write(6,*) "--ERROR: Cannot read nstep, or nstep = 0"
    call MPI_abort(MPI_COMM_WORLD,ierr)
    stop
  endif

  if(iMD.lt.0) iMD=1
  if(ibo_md.lt.0) ibo_md=0
  if(temperature.lt.0.d0) temperature=0.d0
  if(mscf.lt.0) mscf=100
  if(tolrho.lt.0.d0) tolrho=1.d-3
  if(tolintgl.lt.0.d0) tolintgl=1.d-4
  if(init_time.lt.0.d0.or.ntime_init.eq.0) init_time=0.d0
  if(noutput.lt.0) noutput=10
  if(ifatom.lt.0) ifatom=1
  if(ivext_dt.lt.0) ivext_dt=0

  if(rkappa(1).ne.0.d0.or.rkappa(2).ne.0.d0.or.rkappa(3).ne.0.d0) &
      ikappa=1

  if(iboltz.lt.0) iboltz=0
  if(ibshift.lt.0) ibshift=0
  if(i_scale.lt.0) i_scale=0
  if(iocc.lt.-1) iocc=0
  if(jhole.lt.0) jhole=floor(tot/2.d0+0.1)+1
  if(jelec.lt.0) jelec=floor(tot/2.d0+0.1)+2
  if(jxi.lt.0) jxi=1
  if(rxi.lt.0.d0) rxi=1.d0
  if(t_timewall.lt.0.d0) t_timewall = 0.d0


  IF(ntime_init.gt.0) THEN
  IF(inode.eq.1) THEN

    filename=trim(adjustl(ftmp_st))//"."
    WRITE(fileindex,'(i)') ntime_init-1
    filename=trim(adjustl(filename))//trim(adjustl(fileindex))

    OPEN(29,FILE=filename,FORM='unformatted',STATUS='old',ACTION='READ',IOSTAT=ierr)

    IF(ierr.ne.0) THEN

      CLOSE(29)

      CALL system('ls '//trim(adjustl(ftmp_st))//'.* > filelog')

      OPEN(31,FILE='filelog',ACTION='READ')

      ierr=0
      DO WHILE(ierr.eq.0)
        READ(31,'(a)',IOSTAT=ierr) filename
      END DO
      CLOSE(31)

      filename=trim(adjustl(filename))

      i=len_trim(filename)
      j=index(filename,'.')

      IF(j.eq.0) THEN
        WRITE(6,*) "--ERROR: Cannot read "// &
                   trim(adjustl(ftmp_st))
        CALL MPI_abort(MPI_COMM_WORLD,ierr)
        STOP
      END IF

      READ(filename(j+1:i),*) ntime_init
      ntime_init=ntime_init+1

      OPEN(29,FILE=filename,FORM='unformatted',STATUS='old',ACTION='READ',IOSTAT=ierr)
      REWIND(29)

    END IF

    READ(29) itmp,init_time

    IF(itmp.ne.ntime_init.or.ierr.ne.0) THEN
      WRITE(6,*) "--ERROR: Wrong filehead in "//filename
      CALL MPI_abort(MPI_COMM_WORLD,ierr)
      STOP
    END IF

    init_time=init_time+dtMD

    CLOSE(29)

  END IF
  END IF

  CALL mpi_bcast(ntime_init,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL mpi_bcast(init_time,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  return
END SUBROUTINE input_TDDFT
