!!!*
! ルンゲクッタ法でサイクロトロン運動を求める
! @module CyclotronWithRungeKutta

! xyz座標系について 原点： 地表の北極点 (地球の中心ではない！)
! x: 赤道平面と平行な軸 [m]
! y: (xとzに垂直な方向) [m]
! z: 天頂方向 [m]
!!!
module CyclotronWithRungeKutta

    ! 円周率
    double precision, parameter:: pi = 3.14159265358979323d0

    ! 自然対数の底
    double precision, parameter:: e  = 2.71828182845904524d0

    ! 光速 [m/s]
    double precision, parameter:: c = 2.99792458d8

    ! spaceRange: グリッドの範囲
    integer, parameter :: spaceRange = 1024

    ! h: 時間の刻み幅 [s]
    double precision, parameter :: h = 0.5d-2 * 7.15d-7 ! 電離圏でのサイクロトロン周期(仮)

    ! initial_r: 粒子の初期位置
    double precision, parameter :: initial_r(3) = (/0d0, 0d0, 300d3/) ! 北極方向に地表から300km地点

    ! earthRadius: 地球の半径 [m]
    double precision, parameter :: earthRadius = 6371d3

    ! V_LEN: 粒子の速度の絶対値 [m/s]
    double precision, parameter :: V_LEN = 6.24d-2 * c

    ! minZ: 最低高度 [m]
    double precision, parameter :: minZ = 50d3

    ! maxZ : 最高高度 [m]
    double precision, parameter :: maxZ = 300d3

    ! 粒子
    type Particle
        ! r: 位置([m], [m], [m]), v: 速度 ([m/s], [m/s], [m/s])
        double precision :: r(3), v(3)

        ! 初期位置でのピッチ角 [°]
        double precision :: initialPA

        ! Mirror Forceを無視するかどうか
        logical:: ignoreMF = .true.

        ! q: 電荷 [C]
        double precision :: q = -1.602176d-19

        ! m: 質量 [kg]
        double precision :: m = 9.10938356d-31

    end type Particle


    ! 球面座標
    type Spherical

        ! 地球の中心からの距離 [m]
        double precision :: r

        ! 赤道平面とのなす角 [rad]
        double precision :: lambda

    end type Spherical


    contains


        !!!*
        ! ルンゲクッタ法を繰り返し、各ステップにおける位置を計算
        !
        ! @subroutine run
        ! @public
        !!!
        subroutine run

            ! t: 時間 [s]
            double precision :: t = 0

            integer i

            type(Particle) electron

            electron = createInitialParticle(45d0)

            do i = 1, 10000

                !if (collisionOccurred(electron, h)) then
                !    exit
                !endif
                !write(*,*) electron%r
                write(*,*) i, vlen(electron%v)
                !write(*,*) electron%v

                t = t + h
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
            double precision time, v(3), r(3), runge(3)

            runge = rungekuttaParticle(pt, time, acceleration)

            pt%v = pt%v + runge
            pt%r = pt%r + pt%v * h

        end


        !!!*
        ! 粒子の速度から運動エネルギーを計算
        ! @function kineticEnergy
        ! @param {Particle} 粒子
        ! @return {double} 運動エネルギー
        !!!
        function kineticEnergy(pt)
            double precision kineticEnergy, v_len
            type(Particle) pt

            kineticEnergy = 0.5d0 * pt%m * vlen(pt%v) ** 2
        end


        !!!*
        ! 与えられた数だけ粒子を生成する
        ! @function createManyParticles
        ! @param {integer} num 粒子の数
        ! @return {Particle(num)} 粒子を要素に持つ配列
        !!!
        function createManyParticles(num)
            integer num
            type(Particle) createManyParticles(num)

            do i = 1, num
                createManyParticles(i) = createInitialParticle(89d0)
            end do

        end


        !!!*
        ! 中性大気との衝突が起きたかどうか、高度(x)に応じてランダムに結果を生成
        ! @function collisionOccurred
        ! @param {Particle} pt 粒子
        ! @param {double} dt 時間の微小変化
        ! @return {logical} 中性大気との衝突したかどうか
        !!!
        function collisionOccurred(pt, dt)
            type(Particle) pt
            double precision randomNum, dt
            logical collisionOccurred

            call random_number(randomNum)
            collisionOccurred = (randomNum < collisionProbabilityByZ(pt%r(3), dt))

        end


        !!!*
        ! 与えられたz座標における中性大気との衝突確率
        ! @function collisionProbabilityByZ
        ! @param {double} z
        ! @param {double} dt 時間の微小変化
        ! @return {double} 確率
        !!!
        function collisionProbabilityByZ(z, dt)
            double precision z, dt, collisionProbabilityByZ

            collisionProbabilityByZ = e**(-z * 10000 * dt)

        end


        !!!*
        ! ピッチ角を与え、初期の粒子を得る
        ! @function createInitialParticle
        ! @return {Particle} 粒子
        !!!
        function createInitialParticle(angle)
            type(Particle) createInitialParticle

            double precision initialPA, initialBx, angle, initialVel(3)


            initialVel(1) = 0d0
            initialVel(2) = sin(angle * pi / 180) * V_LEN
            initialVel(3) = - cos(angle * pi / 180) * V_LEN

            createInitialParticle = Particle(initial_r, initialVel, initialPA)

        end


        !!!*
        ! 粒子のピッチ角を求める
        ! @function pitchAngle
        ! @param {Particle} pt
        ! @return {double} ピッチ角 [°]
        !!!
        function pitchAngle(pt)
            type(Particle) pt
            double precision absVelocity, pitchAngle

            absVelocity = vlen(pt%v)
            pitchAngle = 180 - acos(pt%v(3) / absVelocity) * 180 / pi

        end



        !!!*
        ! ある座標におけるdipole磁場を、divB = 0 となるように補正したもの
        ! @function B
        ! @param {double(3)} r 粒子の位置
        ! @return {double(3)} 磁場
        !!!
        function B(r)

            double precision B(3), r(3), Bxy, k, sinL, cosL, spR3, xyAngle
            double precision, parameter :: permeability = 4d-7 * pi ! 真空での透磁率 [H/m]
            double precision, parameter :: dipoleMoment = 8.05d22   ! Earth's dipole moment [m]
            type(Spherical) sp

            !sp = toSpherical(r) ! xyz座標を球面座標に変換

            !sinL = sin(sp%lambda)
            !cosL = cos(sp%lambda)
            !spR3 = sp%r**3

            !k = permeability * dipoleMoment / 4 / pi ! 係数

            !B(3) = k * (-3*sinL**2 + 1) / spR3

            !! divB = 0となるように調整された、平行な磁場ベクトルのxyに平行な成分
            !Bxy = k * cosL * sinL * ( 9 - 6 * sinL**2) / (2 * spR3)

            !xyAngle = atan2(r(2), r(1))

            !B(1) = Bxy * cos(xyAngle)
            !B(2) = Bxy * sin(xyAngle)
            B(1) = 0d0
            B(2) = 0d0
            B(3) = 50000d-9

        end


        !!!*
        ! ある座標における電場を取得
        ! @function El
        ! @param {double(3)} r 粒子の位置
        ! @return {double(3)} 電場の大きさ
        !!!
        function El(r)

            double precision El(3), r(3)

            El(1) = 0d0
            El(2) = 0d0
            El(3) = 0d0

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
            double precision acceleration(3), time

            acceleration = pt%q / pt%m * ( El(pt%r) + cross( pt%v, B(pt%r) ))

        end function



        !!!*
        ! 磁場の大きさ、粒子の質量、電荷からサイクロトロン周期を求める
        ! @function period
        ! @param {Particle} pt 粒子
        ! @return {double} サイクロトロン周期 [s]
        !!!
        function cyclotronPeriod(pt)
            type(Particle) pt
            double precision cyclotronPeriod

            cyclotronPeriod = 2 * pi * pt%m / abs(pt%q) / vlen(B(pt%r))

        end function


        !!!*
        ! 与えられた3次元ベクトルの長さを取得
        ! @function vlen
        ! @param {double(3)} vec
        ! @return {double} length
        !!!
        function vlen(vec)
            double precision vlen, vec(3)
            vlen = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
        end function


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
            integer, dimension(spaceRange) :: arr

            val   = val + 1   ! インデックスを 1スタートに変更
            left  = int(val)  ! 左端
            right = left + 1  ! 右端

            ! 値が配列の右端より右にあるときは0を返す
            if (right > spaceRange) then
                getLinearly = 0d0

            ! 値が配列の左端より左にあるときはarr(1)を返す
            elseif (right < 1) then
                getLinearly = arr(1)

            ! valを挟む整数の右端が未定義のときはそれを0とみなす
            elseif (left == spaceRange) then
                getLinearly = (right - val) * arr(left)

            else
                getLinearly = (right - val) * arr(left) + (val - left) * arr(right)

            endif
        end



        !!!*
        ! 与えられたxyz座標を球面座標に変換
        ! 球面座標の中心は地球の中心 (xyz座標の原点とは異なる)
        ! @function toSpherical
        ! @param {double(3)} xyz
        ! @param {Spherical}
        !!!
        function toSpherical(xyz)

            type(Spherical) toSpherical

            double precision xyz(3), earthCenteredXYZ(3), r, lambda, xylen

            earthCenteredXYZ(1) = xyz(1)
            earthCenteredXYZ(2) = xyz(2)
            earthCenteredXYZ(3) = xyz(3) + earthRadius

            r = vlen(earthCenteredXYZ)

            xylen = sqrt(xyz(1)**2 + xyz(2)**2)

            lambda = acos(xylen / r)

            toSpherical = Spherical(r, lambda)

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
            double precision time
            double precision rungekuttaParticle(3), k1(3), k2(3), k3(3), k4(3)
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
            p1 = Particle(pt%r, pt%v + k1 * h / 2, pt%initialPA, pt%ignoreMF)

            d2 = (pt%v + pt%v + k1 * h / 2) * h / 2 / 2 ! 台形
            k2 = accel(p1, time + h / 2)
            p2 = Particle(p1%r + d2, p1%v + k2 * h / 2, pt%initialPA, pt%ignoreMF)

            d3 = (pt%v + pt%v + k2 * h / 2) * h / 2 / 2 ! 台形
            k3 = accel(p2, time + h / 2)
            p3 = Particle(p2%r + d3, p2%v + k3 * h, pt%initialPA, pt%ignoreMF)

            d4 = (pt%v + pt%v + k3 * h) * h / 2 ! 台形
            k4 = accel(p3, time + h)

           rungekuttaParticle = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

           return
        end
    !end contains
end


program main
    use CyclotronWithRungeKutta

    call run

end
