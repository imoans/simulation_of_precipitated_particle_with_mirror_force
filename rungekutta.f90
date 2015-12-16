!!!*
! ルンゲクッタ法でサイクロトロン運動を求める
! @module CyclotronWithRungeKutta
!!!
module CyclotronWithRungeKutta

    use mtmod

    ! h: 刻み幅
    double precision, parameter :: h = 1.0d-3

    type Particle
        ! r: 位置, v: 速度
        double precision :: r(3), v(3)

        ! q: 電荷
        double precision :: q = 1.0d0

        ! m: 質量
        double precision :: m = 1.0d0

    end type Particle



    contains


        !!!*
        ! ルンゲクッタ法を繰り返し、各ステップにおける位置を計算
        !
        ! @subroutine run
        ! @public
        !!!
        subroutine run

            ! t: 時間
            double precision :: t = 0

            integer i

            type(Particle) electron

            electron = createParticle((/5.12d0, 0d0, 0d0/), (/1.7d-2, 1.7d-2, 0d0/))

            do i = 1, 10000
                call update(electron, t)
            end do
        end


        !!!*
        ! 与えられた粒子の位置と磁場を更新
        ! @subroutine update
        ! @param {Particle} pt 粒子
        ! @param {double} time 時間
        !!!
        subroutine update(pt, time)

            type(Particle) pt
            double precision time

            pt%v = pt%v + rungekuttaParticle(pt, time, acceleration)
            pt%r = pt%r + pt%v * h

        end


        !!!*
        ! 与えられた位置と速度から粒子を生成
        ! @function createParticle
        ! @param {double(3)} 粒子の位置
        ! @param {double(3)} 粒子の速度
        ! @return {Particle} 粒子
        !!!
        function createParticle(r, v)
            type(Particle) createParticle
            double precision r(3), v(3)

            createParticle%r = r
            createParticle%v = v
        end


        !!!*
        ! 与えられた長さと角度からベクトルをランダムに生成する
        ! @function createVectorByRandom
        ! @param {double} absLen 長さ
        !!!
        function createVectorByRandom(absLen)
            double precision createVectorByRandom(3)
            double precision absLen, angle, pi, randomNum

            call sgrnd(4357) ! 乱数生成開始

            randomNum = grnd()
            pi = atan(1.0) * 4.0
            angle = randomNum * 2 * pi

            createVectorByRandom(1) = cos(angle) * absLen
            createVectorByRandom(2) = sin(angle) * absLen
            createVectorByRandom(3) = 0d0
        end


        !!!*
        ! 与えられた粒子のラーマー半径を取得
        ! @function larmorRadius
        ! @param {Particle} particle 粒子
        ! @return {double} ラーマー半径
        !!!
        function larmorRadius(pt)

            type(Particle) pt
            double precision larmorRadius, currentB(3)
            integer :: lort = 1 ! ローレンツファクター

            currentB = B(pt%r)

            v_perp_len = sqrt( pt%v(2)**2 + pt%v(3)**2 )
            larmorRadius = v_perp_len / currentB(1) * lort

        end


        !!!*
        ! ある座標において粒子の作る動径方向の磁場の大きさを取得
        ! @function radialB
        ! @param {double(3)} r 粒子の位置
        ! @return {double} 磁場の大きさ
        !!
        function radialB(r)

            double precision radialB, r(3)
            integer, dimension(1024) :: radialBArr = 1.0d0

            radialB = getLinearly(radialBArr, r(1) * 100)

        end


        !!!*
        ! ある座標における磁場を取得
        ! @function B
        ! @param {double(3)} r 粒子の位置
        ! @return {double(3)} 磁場の大きさ
        !!!
        function B(r)

            double precision B(3), r(3)
            integer, dimension(1024) :: BArr = (/ (i**2,i=1,1024) /)

            B(1) = getLinearly(BArr, r(1) * 100)
            B(2) = 0d0
            B(3) = 0d0
        end


        !!!*
        ! ある座標における電場を取得
        ! @function E
        ! @param {double(3)} r 粒子の位置
        ! @return {double(3)} 電場の大きさ
        !!!
        function E(r)

            double precision E(3), r(3)
            integer, dimension(1024) :: EArr = 0

            E(1) = 0d0
            E(2) = 0d0
            E(3) = getLinearly(EArr, r(3) * 100)

        end




        !!!*
        ! 加速度
        ! @function acceleration
        ! @param {Particle} pt 粒子
        ! @param {double} time 時間
        ! @return {double(3)} 加速度
        !!!
        function acceleration(pt, time)

            type(Particle) pt
            double precision acceleration(3), correctedB(3)
            double precision time

            correctedB = B(pt%r) + correctB(pt)

            acceleration = pt%q / pt%m * ( E(pt%r) + cross( pt%v, correctedB ) )
        end



        !!!*
        ! 動径方向の磁場を考慮した磁場の補正分を計算
        ! @function correctB
        ! @param {Particle} pt 粒子
        ! @return {double(3)} 磁場の補正分の大きさ
        !!!
        function correctB(pt)
            type(Particle) pt
            double precision correctB(3), v_perp
            double precision currentBry, currentBrz, currentBr
            integer :: lort = 1 ! ローレンツファクター

            v_perp = sqrt( pt%v(2)**2 + pt%v(3)**2 )

            currentBry = radialB(pt%r) * pt%v(3) / v_perp
            currentBrz = radialB(pt%r) * pt%v(2) / v_perp
            currentBr  = radialB(pt%r) * larmorRadius(pt) / lort

            correctB(1) = - v_perp  * currentBr
            correctB(2) = - pt%v(1) * currentBry
            correctB(3) =   pt%v(1) * currentBrz

        end




        !!!*
        ! 外積
        ! @function cross
        ! @param {double(3)} x
        ! @param {double(3)} y
        ! @return {double(3)} 結果
        !!!
        function cross(x, y)
            double precision cross(3), x(3), y(3)

            cross(1) = x(2) * y(3) - x(3) * y(2)
            cross(2) = x(3) * y(1) - x(1) * y(3)
            cross(3) = x(1) * y(2) - x(2) * y(1)
        end




        !!!*
        ! 与えられたスカラー配列を連続値とみなし、インデックスから値を得る
        ! 内分した値を利用
        ! @function getLinearly
        ! @param {dimension} arr スカラーの配列
        ! @param {double} val インデックス (0スタート)
        !!!
        function getLinearly(arr, val)

            double precision getLinearly, val
            integer left, right
            integer, dimension(1024) :: arr

            val   = val + 1   ! インデックスを 1スタートに変更
            left  = int(val)  ! 左端
            right = left + 1  ! 右端

            ! 値が配列の右端より右にあるときは0を返す
            if (right > 1024) then
                getLinearly = 0d0

            ! 値が配列の左端より左にあるときはarr(1)を返す
            elseif (right < 1) then
                getLinearly = arr(1)

            ! valを挟む整数の右端が未定義のときはそれを0とみなす
            elseif (left == 1024) then
                getLinearly = (right - val) * arr(left)

            else
                getLinearly = (right - val) * arr(left) + (val - left) * arr(right)

            endif
        end



        !!!*
        ! 粒子の速度に対してルンゲクッタ法を適用する
        ! 粒子の加速度の計算に粒子の位置が必要なため粒子の構造体を渡す設計にしている
        !
        ! @function rungekuttaParticle
        ! @param {Particle} pt 粒子
        ! @param {double} time 時間
        ! @param {function} accel 加速度の式 (一次導関数微分方程式の右辺)
        ! @return {double(3)} 速度の変分
        !!!
        function rungekuttaParticle(pt, time, accel)
            type(Particle) pt, p1, p2, p3
            double precision rungekuttaParticle(3), time
            double precision k1(3), k2(3), k3(3), k4(3)
            double precision d2(3), d3(3), d4(3) ! d2〜d4 微小時間で移動した距離

            ! 一次導関数微分方程式の右辺
            interface
                function accel(p, t)
                    import :: Particle
                    type(Particle) p
                    double precision accel(3)
                    double precision t
                end
            end interface

            k1 = accel(pt, time)
            p1 = createParticle(pt%r, pt%v + k1 * h / 2)

            d2 = (pt%v + pt%v + k1 * h / 2) * h / 2 / 2 ! 台形
            k2 = accel(p1, time + h / 2)
            p2 = createParticle(p1%r + d2, p1%v + k2 * h / 2)

            d3 = (pt%v + pt%v + k2 * h / 2) * h / 2 / 2 ! 台形
            k3 = accel(p2, time + h / 2)
            p3 = createParticle(p2%r + d3, p2%v + k3 * h)

            d4 = (pt%v + pt%v + k3 * h) * h / 2 ! 台形
            k4 = accel(p3, time + h)

            rungekuttaParticle = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        end
    !end contains
end


program main
    use CyclotronWithRungeKutta
    call run

end
