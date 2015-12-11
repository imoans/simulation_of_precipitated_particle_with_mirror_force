
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


    contains


        !!!*
        ! ルンゲクッタ法を繰り返し、各ステップにおける位置を計算
        !
        ! @subroutine run
        ! @public
        ! @param {double} init_vel x,y方向の初速
        !!!
        subroutine run(init_vel)
            double precision init_vel

            ! t: 時間
            double precision :: t = 0

            ! r: 位置
            double precision :: r(3) = (/ 5.12d0, 0d0, 0d0 /) ! x = 磁場の配列の中央の値

            ! v: 速度
            double precision v(3)



            v(1) = init_vel
            v(2) = init_vel
            v(3) = 0d0

            do i = 1, 10000
                write(*,*) r(1), r(2), r(3)
                v = v + rungekutta(v, t, acceleration)
                r = r + v * h
            end do
        end


        !!!*
        ! ある座標における磁場を取得
        ! @function B
        !!!
        function B(r)

            double precision B(3), r(3)
            double precision x100, Bx
            integer i, j
            integer, dimension(1024) :: BArr = (/ (i**2,i=1,1024) /)

            x100 = r(1) * 100 + 1 ! 粒子のx座標の100倍 (BArrの単位に合わせるため. 1を足すことで配列に対応)
            i    = int(x100) ! 粒子のいるBArrの左端
            j    = i + 1     ! 粒子のいるBArrの右端

            ! 粒子が磁場の右側の範囲外にいるとき
            if (i > 1024) then
                Bx = 0d0

            ! 粒子が磁場の左側の範囲外にいるとき
            elseif (i < 1) then
                Bx = BArr(1)

            ! 粒子のいるBArrの右端が未定義のときはそれを0とみなす
            elseif (i == 1024) then
                Bx = (j - x100) * BArr(i)

            else
                Bx = (j - x100) * BArr(i) + (x100 - i) * BArr(j)

            endif

            B(1) = Bx
            B(2) = 0d0
            B(3) = 0d0

        end


        !!!*
        ! ある座標における電場を取得
        ! @function E
        !!!
        function E(r)

            double precision E(3), r(3)
            double precision x100, Ez
            integer i, j
            integer, dimension(1024) :: EArr = 1

            x100 = r(3) * 100 + 1 ! 粒子のx座標の100倍 (BArrの単位に合わせるため. 1を足すことで配列に対応)
            i    = int(x100) ! 粒子のいるEArrの左端
            j    = i + 1     ! 粒子のいるEArrの右端

            ! 粒子が磁場の右側の範囲外にいるとき
            if (i > 1024) then
                Ez = 0d0

            ! 粒子が磁場の左側の範囲外にいるとき
            elseif (i < 1) then
                Ez = EArr(1)

            ! 粒子のいるEArrの右端が未定義のときはそれを0とみなす
            elseif (i == 1024) then
                Ez = (j - x100) * EArr(i)

            else
                Ez = (j - x100) * EArr(i) + (x100 - i) * EArr(j)

            endif

            E(1) = 0d0
            E(2) = 0d0
            E(3) = Ez

        end




        !!!*
        ! 加速度
        ! @function acceleration
        !!!
        function acceleration(v, time, dr)
            double precision acceleration(3), v(3), dr(3)
            double precision time

            acceleration = q / m * ( E(r + dr) + cross( v, B(r + dr) ) )
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
            ! d1〜d2 微小時間で移動した距離
            double precision d1(3), d2(3), d3(3), d4(3)
            double precision time

            ! 一次導関数微分方程式の右辺
            ! 第３引数には 微小時間でのvec3の積分値 (台形近似) を与える
            interface
                function f(vec3, time, deltaVecTime)
                    double precision f(3)
                    double precision vec3(3)
                    double precision time
                    double precision deltaVecTime(3)

                end
            end interface

            d1 = (/0d0, 0d0, 0d0/)
            k1 = f(vec3, time, d1)

            d2 = (vec3 + vec3 + k1 * h / 2) * h / 2 / 2 ! 台形
            k2 = f(vec3 + k1 * h / 2, time + h / 2, d2)

            d3 = (vec3 + vec3 + k2 * h / 2) * h / 2 / 2 ! 台形
            k3 = f(vec3 + k2 * h / 2, time + h / 2, d3)

            d4 = (vec3 + vec3 + k3 * h) * h / 2 ! 台形
            k4 = f(vec3 + k3 * h, time + h, d4)

            rungekutta = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        end
    !end contains
end


program main
    use CyclotronWithRungeKutta
    call run(1.7d-2)

end
