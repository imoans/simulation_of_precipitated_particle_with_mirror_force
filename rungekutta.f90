
!!!*
! ルンゲクッタ法でサイクロトロン運動を求める
! @module CyclotronWithRungeKutta
!!!
module CyclotronWithRungeKutta

    ! h: 刻み幅
    double precision, parameter :: h = 1.0d-3

    ! q: 電荷
    double precision, parameter :: q = 1.0d0

    ! m: 質量
    double precision, parameter :: m = 1.0d0

    ! B: 磁場
    double precision, parameter :: B(3) = (/ 0, 0, 1 /)


    contains

        !!!*
        ! ルンゲクッタ法を繰り返し、各ステップにおける位置を計算
        !
        ! @subroutine run
        ! @public
        !!!
        subroutine run(init_vel)
            double precision init_vel

            ! t: 時間
            double precision :: t = 0

            ! r: 位置
            double precision :: r(3) = (/ 0, 1, 0 /)

            ! v: 速度
            double precision v(3)
            v(1) = init_vel
            v(2) = 0
            v(3) = 0

            do i = 1, 10000
                write(*,*) r(1), r(2)
                v = v + rungekutta(v, t, acceleration)
                r = r + v * h
            end do
        end


        !!!*
        ! 加速度
        ! @function acceleration
        !!!
        function acceleration(v, time, dr)
            double precision acceleration(3)
            double precision time
            double precision v(3)
            double precision dr
            acceleration = q / m * (E + cross(v, B(r + dr)))
        end


        !!!*
        ! 外積
        ! @function cross
        !!!
        function cross(x, y)
            double precision cross(3), x(3), y(3)

            cross(1) = x(2) * y(3) - x(3) * y(2)
            cross(2) = x(3) * y(1) - x(1) * y(3)
            cross(3) = x(1) * y(2) - x(2) * y(1)
        end


        !!!*
        ! ルンゲクッタ法
        ! @function rungekutta
        !!!
        function rungekutta(vec3, time, f)

            double precision rungekutta(3), vec3(3), k1(3), k2(3), k3(3), k4(3)
            double precision time

            ! 一次導関数微分方程式の右辺
            ! 第３引数にはΔ vecΔ h (deltaVecTime) をいれる
            interface
                function f(vec3, time, deltaVecTime)
                    double precision f(3)
                    double precision vec3(3)
                    double precision deltaVecTime(3)
                    double precision time
                end
            end interface

            k1 = h * f(vec3, time, 0)
            k2 = h * f(vec3 + k1 / 2, time + h / 2, k1 * h / 4)
            k3 = h * f(vec3 + k2 / 2, time + h / 2, k2 * h / 4)
            k4 = h * f(vec3 + k3, time + h, k3 * h)

            rungekutta = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        end
    !end contains
end

module ReadFileModule
    implicit none

    contains

        function readFile(filename, datanum)
            double precision readFile(datanum), data(datanum)
            character(len=20) filename
            integer datanum, i
            real value

            open(17, file=filename, status='old')

                do i = 1, datanum

                    read(17, *)value
                    data(i) = value

                end do

            close(17)
            readFile = data
        end

end

program main
    use CyclotronWithRungeKutta
    call run(1.7d-2)

end
