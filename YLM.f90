COMPLEX(8) FUNCTION YLM(l, m, theta, phi)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(8), INTENT(IN) :: theta, phi

    REAL(8), PARAMETER ::    pi = 4.d0 * ATAN(1.d0)
    COMPLEX(8), PARAMETER :: ci =   (0.d0, 1.d0), &
                             one =  (1.d0, 0.d0), &
                             zero = (0.d0, 0.d0)

    YLM = zero

    SELECT CASE (l)
        CASE (0)
            SELECT CASE (m)
                CASE (0)
                    YLM = 1.d0 / 2.d0 / SQRT(pi)
                CASE DEFAULT
                    GOTO 9999
            END SELECT
        CASE (1)
            SELECT CASE (m)
                CASE (0)
                    YLM = SQRT(3.d0) * COS(theta) / 2.d0 / SQRT(pi)
                CASE (-1)
                    YLM = SQRT(3.d0) * SIN(theta) * EXP(-ci * phi) / 2.d0 / SQRT(2.d0 * pi)
                CASE (1)
                    YLM = 0.d0 - SQRT(3.d0) * SIN(theta) * EXP(ci * phi) / 2.d0 / SQRT(2.d0 * pi)
                CASE DEFAULT
                    GOTO 9999
            END SELECT
        CASE (2)
            SELECT CASE (m)
                CASE (0)
                    YLM = SQRT(5.d0) * (3.d0 * (COS(theta))**2 - 1.d0) / 4.d0 / SQRT(pi)
                CASE (-1)
                    YLM = SQRT(15.d0) * SIN(theta) * COS(theta) * EXP(-ci * phi) / 2.d0 / SQRT(2.d0 * pi)
                CASE (1)
                    YLM = 0.d0 - SQRT(15.d0) * SIN(theta) * COS(theta) * EXP(ci * phi) / 2.d0 / SQRT(2.d0 * pi)
                CASE (-2)
                    YLM = SQRT(15.d0) * (SIN(theta))**2 * EXP(-ci * 2.d0 * phi) / 4.d0 / SQRT(2.d0 * pi)
                CASE (2)
                    YLM = SQRT(15.d0) * (SIN(theta))**2 * EXP(ci * 2.d0 * phi) / 4.d0 / SQRT(2.d0 * pi)
                CASE DEFAULT
                    GOTO 9999
            END SELECT
        CASE DEFAULT
            GOTO 9999
    END SELECT

9999 CONTINUE

END FUNCTION YLM
