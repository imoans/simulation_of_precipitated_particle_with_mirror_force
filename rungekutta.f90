!!!*
! ルンゲクッタ法でサイクロトロン運動を求める
! @module CyclotronWithRungeKutta
!!!
module CyclotronWithRungeKutta

    ! spaceRange: グリッドの範囲
    integer, parameter :: spaceRange = 1024

    ! h: 時間の刻み幅
    double precision, parameter :: h = 1.0d-3

    ! initial_r: 粒子の初期位置
    double precision, parameter :: initial_r(3) = (/5.12d0, 0d0, 0d0/) ! 磁場の配列の長さの中央の値
    ! TODO: 配列の長さが可変になった時、xの値も調節する

    ! V_LEN: 粒子の速度の絶対値
    double precision, parameter :: V_LEN = 1.7d-2

    type Particle
        ! r: 位置, v: 速度
        double precision :: r(3), v(3)

        ! 初期位置でのピッチ角 [°]
        double precision :: initialPA

        ! Mirror Forceを無視する
        logical:: ignoreMF = .true.

        ! q: 電荷
        double precision :: q = -1.0d0

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
            double precision energy

            electron = createInitialParticleByRandom()

            do i = 1, 1000

                !write(*,*) i, electron%r(1)
                !write(*,*) electron%v
                !write(*,*) i, pitchAngle(electron)
                write(*,*) t, kineticEnergy(electron)

                if (collisionOccurred(electron, h)) then
                    exit
                endif

                t = t + h
                call update(electron, t)
                energy = kineticEnergy(electron)

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
        ! @param {double(3)} 速度
        ! @return {double} 運動エネルギー
        !!!
        function kineticEnergy(pt)
            double precision kineticEnergy, v_len
            type(Particle) pt

            v_len = sqrt(pt%v(1) ** 2 + pt%v(2) **2 + pt%v(3) **2)
            kineticEnergy = 0.5d0 * pt%m * v_len ** 2
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
                createManyParticles(i) = createInitialParticleByRandom()
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
            collisionOccurred = (randomNum < collisionProbabilityByX(pt%r(1), dt))

        end


        !!!*
        ! 与えられたx座標における中性大気との衝突確率
        ! x > 0を想定
        ! @function collisionProbabilityByX
        ! @param {double} x
        ! @param {double} dt 時間の微小変化
        ! @return {double} 確率
        !!!
        function collisionProbabilityByX(x, dt)
            double precision x, dt, collisionProbabilityByX
            double precision, parameter:: e = 2.71828182845904524d0 ! 自然対数の底

            collisionProbabilityByX = e**(-x * 10000 * dt)

        end


        !!!*
        ! 初期位置、初期速度粒子の初速をランダムに取得
        ! ピッチ角(0°~90°)をランダムな一様分布として与える
        ! @function createInitialParticleByRandom
        ! @return {Particle} 粒子
        !!!
        function createInitialParticleByRandom()
            type(Particle) createInitialParticleByRandom

            double precision initialPA, initialBx, angle, pi, randomNum, initialVel(3)

            pi = atan(1.0) * 4.0

            call random_number(randomNum) ! 乱数生成開始

            initialPA = randomNum * 90 ! 初期位置でのPitch Angle

            initialVel(1) = - cos(initialPA * pi / 180) * V_LEN
            initialVel(2) = sin(initialPA * pi / 180) * V_LEN
            initialVel(3) = 0d0

            createInitialParticleByRandom = Particle(initial_r, initialVel, initialPA)


        end


        !!!*
        ! 粒子のピッチ角を求める
        ! @function pitchAngle
        ! @param {Particle} pt
        ! @return {double} ピッチ角 [°]
        !!!
        function pitchAngle(pt)
            type(Particle) pt
            double precision absVelocity, pitchAngle, pi

            pi = atan(1.0) * 4.0
            absVelocity = sqrt( pt%v(1)**2 + pt%v(2)**2 + pt%v(3)**2 )
            pitchAngle = 180 - acos(pt%v(1) / absVelocity) * 180 / pi

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
            integer, dimension(spaceRange) :: radialBArr = (/ (0.5d0*(spaceRange+1-i),i=1,spaceRange) /)

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
            integer, dimension(spaceRange) :: BArr = (/ ((spaceRange+1-i),i=1,spaceRange) /)

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
            integer, dimension(spaceRange) :: EArr = 0

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
            double precision acceleration(3), time, accel(3), mf(3)

            if (pt%ignoreMF) then
                mf = (/0d0, 0d0, 0d0/) ! 磁場の配列の長さの中央の値
            else
                mf = correctLortForce(pt)
            endif

            acceleration = pt%q / pt%m * ( E(pt%r) + cross( pt%v, B(pt%r) ) + mf)
        end function acceleration



        !!!*
        ! 動径方向の磁場を考慮したローレンツ力の補正分を計算
        ! @function correctLortForce
        ! @param {Particle} pt 粒子
        ! @return {double(3)} 補正分の大きさ
        !!!
        function correctLortForce(pt)
            type(Particle) pt
            double precision correctLortForce(3), v_perp
            double precision currentBry, currentBrz, currentBr
            integer :: lort = 1 ! ローレンツファクター

            v_perp = sqrt( pt%v(2)**2 + pt%v(3)**2 )

            currentBry = radialB(pt%r) * pt%v(3) / v_perp
            currentBrz = - radialB(pt%r) * pt%v(2) / v_perp
            currentBr  = radialB(pt%r) * larmorRadius(pt) / lort

            correctLortForce(1) = - v_perp  * currentBr
            correctLortForce(2) = - pt%v(1) * currentBrz
            correctLortForce(3) =   pt%v(1) * currentBry

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
