      subroutine system_flush(ifile)
      integer ifile
ccccccc flush the buffer of output file: ifile

      call flush(ifile)

      return
      end
