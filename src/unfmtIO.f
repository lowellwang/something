      subroutine unfmtIO_r(iIO,fname,rwork,isize,iflag,ibcast)

      implicit none

      include 'mpif.h'
      include 'param.escan_real'

      integer iIO,iflag,ibcast,isize,ierr
      real*8 rwork(isize)
      character*20 fname
cccccccccccccccccccccccccccccccccccccccccccc
c     iflag.eq.0, read; iflag.eq.1, write
c     ibcast.eq.1: do mpi_bcast
cccccccccccccccccccccccccccccccccccccccccccc
      if(iflag.eq.1.and.inode_tot.ne.1) return
      if(iflag.eq.0.and.ibcast.eq.0.and.inode_tot.ne.1) return
cccccccccccccccccccccccccccccccccccccccccccc
      if(iflag.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccc
       if(inode_tot.eq.1) then
cccccccccccccccccccccccccccccccccccccccccccc
        open(iIO,file=fname,form="unformatted",
     &   status='old',action='read',iostat=ierr)
        if(ierr.ne.0) then
         write(6,*) "file ",fname," does not exist, stop"
         write(17,*) "ERROR: file ",fname," does not exist, stop"
         close(17)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         stop
        endif
        rewind(iIO)
        read(iIO) rwork
        close(iIO)
ccccccccccccccccccccccccccccccccccccccccccccc
       endif
cccccccccccccccccccccccccccccccccccccccccccccc
       if(ibcast.eq.1) then
        call mpi_bcast(rwork,isize,MPI_REAL8,
     & 0,MPI_COMM_WORLD,ierr)
       endif
ccccccccccccccccccccccccccccccccccccccccccccc
      elseif(iflag.eq.1) then
ccccccccccccccccccccccccccccccccccccccccccccc
       if(inode_tot.eq.1) then
        open(iIO,file=fname,form="unformatted")
        rewind(iIO)
        write(iIO) rwork
        close(iIO)
       endif
ccccccccccccccccccccccccccccccccccccccccccccc
      endif

      return

      end



      subroutine unfmtIO_c(iIO,fname,cwork,isize,iflag,ibcast)

      implicit none

      include 'mpif.h'
      include 'param.escan_real'

      integer iIO,iflag,ibcast,isize,ierr
      complex*16 cwork(isize)
      character*20 fname
cccccccccccccccccccccccccccccccccccccccccccc
c     iflag.eq.0, read; iflag.eq.1, write
c     ibcast: do mpi_bcast or not
cccccccccccccccccccccccccccccccccccccccccccc
      if(iflag.eq.1.and.inode_tot.ne.1) return
      if(iflag.eq.0.and.ibcast.eq.0.and.inode_tot.ne.1) return
cccccccccccccccccccccccccccccccccccccccccccc
      if(iflag.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccc
       if(inode_tot.eq.1) then
cccccccccccccccccccccccccccccccccccccccccccc
        open(iIO,file=fname,form="unformatted",
     &   status='old',action='read',iostat=ierr)
        if(ierr.ne.0) then
         write(6,*) "file ",fname," does not exist, stop"
         write(17,*) "ERROR: file ",fname," does not exist, stop"
         close(17)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         stop
        endif
        rewind(iIO)
        read(iIO) cwork
        close(iIO)
ccccccccccccccccccccccccccccccccccccccccccccc
       endif
cccccccccccccccccccccccccccccccccccccccccccccc
       if(ibcast.eq.1) then
        call mpi_bcast(cwork,isize,MPI_DOUBLE_COMPLEX,
     & 0,MPI_COMM_WORLD,ierr)
       endif
ccccccccccccccccccccccccccccccccccccccccccccc
      elseif(iflag.eq.1) then
ccccccccccccccccccccccccccccccccccccccccccccc
       if(inode_tot.eq.1) then
        open(iIO,file=fname,form="unformatted")
        rewind(iIO)
        write(iIO) cwork
        close(iIO)
       endif
ccccccccccccccccccccccccccccccccccccccccccccc
      endif

      return

      end
