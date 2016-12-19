subroutine ortho_Gram(A, m_line, n_col)

    implicit none

    integer, intent(in) :: m_line, n_col
    complex(8), dimension(m_line, n_col) :: A

    integer :: i, j, k
    real(8) :: t, tt
    complex(8) :: s
    complex(8), dimension(m_line, n_col) :: R
    logical :: orig

    real(8), parameter :: mc = 3.D-16
    complex(8), parameter :: ONE = (1.d0, 0.d0), ZERO = (0.d0, 0.d0)

    R = ZERO
    do k = 1, n_col
        t = 0.d0
        do j = 1, m_line
            t = t + abs(A(j,k))**2
        enddo
        t = sqrt(t)
        if(t.le.mc) t = 0.d0
        tt = t
        orig = .true.
        if(k.gt.1) then
            do while(orig)
                do i = 1, k-1
                    s = ZERO
                    do j = 1, m_line
                        s = s + conjg(A(j,i)) * A(j,k)
                    enddo
                    R(i,k) = R(i,k) + s
                    A(1:m_line,k) = A(1:m_line,k) - s * A(1:m_line,i)
                enddo
                tt = 0.d0
                do j = 1, m_line
                    tt = tt + abs(A(j,k))**2
                enddo
                tt = sqrt(tt)
                if(tt.gt.mc) then
                    if(tt.lt.(t/1.d1)) then
                        print *, "WARNING: ortho_Gram: reortho column:",k
                        orig = .true.
                        t = tt
                    else
                        orig = .false.
                    endif
                else
                    print *, "WARNING: ortho_Gram: vector norm is 0, column:",k
                    orig = .false.
                    tt = 0
                endif
            enddo
        endif
        R(k,k) = tt
        if(tt.gt.mc) tt = 1.d0 / tt
        A(1:m_line,k) = A(1:m_line,k) * tt
    enddo

    return

end subroutine ortho_Gram
