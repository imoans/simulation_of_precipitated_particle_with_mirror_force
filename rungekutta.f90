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

    ! h: 時間の刻み幅
    double precision, parameter :: h = 5d-2

    ! initial_r: 粒子の初期位置 (184mを1として、z = 300[km])
    double precision, parameter :: initial_r(3) = (/0d0, 0d0, 1.63d3/)

    ! V_LEN: 粒子の速度の絶対値
    double precision, parameter :: V_LEN = 6.24d-2

    ! 粒子
    type Particle
        ! r: 位置, v: 速度
        double precision :: r(3), v(3)

        ! Mirror Forceを無視するかどうか
        logical:: ignoreMF = .false.

        ! q: 電荷
        double precision :: q = -1d0

        ! m: 質量
        double precision :: m = 1d0

    end type Particle


    ! 単一粒子運動シミュレーションに関わる変数(結果含む)を格納
    type SimulationScope

        ! 初期位置でのピッチ角 [°]
        double precision :: initialPA

        ! 終了するまでのイテレーション回数
        integer iterations

        ! 何回イテレーションするか
        integer:: ITERATION_TIMES = 20000000

        ! 衝突したのかどうか
        logical:: collided = .false.

        ! 衝突した場合はその位置
        double precision collidedR(3)

        ! 衝突した場合そのときのピッチ角[°]
        double precision collidedPA

        ! ミラーポイントに達したのかどうか
        logical:: bounded = .false.

        ! ミラーポイント (あれば)
        double precision mirrorPoint(3)

    end type SimulationScope


    contains


        !!!*
        ! ルンゲクッタ法を繰り返し、各ステップにおける位置を計算
        ! シミュレーション結果を返す
        !
        ! @subroutine run
        ! @param {double} initialPA 初期のピッチアングル
        ! @param {logical} ignoreMF ミラー力を無視するかどうか
        ! @return {SimulationScope} simscope
        ! @public
        !!!
        function run(initialPA, ignoreMF)

            ! t: 時間
            double precision :: t = 0

            double precision initialPA, prevVz
            logical ignoreMF

            type(SimulationScope) sim, run
            type(Particle) electron

            sim%initialPA = initialPA

            electron = createInitialParticle(sim%initialPA, ignoreMF)

            do i = 1, sim%ITERATION_TIMES

                sim%iterations = i

                if (collisionOccurred(electron, h)) then
                    sim%collided = .true.
                    sim%collidedR(1) = electron%r(1)
                    sim%collidedR(2) = electron%r(2)
                    sim%collidedR(3) = electron%r(3)
                    sim%collidedPA = pitchAngle(electron)
                    exit
                endif

                t = t + h

                prevVz = electron%v(3)
                call update(electron, t)

                ! ミラーポイントの判定
                if (prevVz <= 0 .AND. electron%v(3) >= 0) then
                    sim%bounded = .true.
                    sim%mirrorPoint(1) = electron%r(1)
                    sim%mirrorPoint(2) = electron%r(2)
                    sim%mirrorPoint(3) = electron%r(3)
                endif
            end do

            run = sim

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
            double precision z, dt, collisionProbabilityByZ, lambda


            ! 75d3/184のとき0 , 125d3/184のとき-3.5
            ! lambda: サイクロトロン1周期のうちに何回衝突するか
            lambda = 10 ** (-1.288d-2 * z + 1.288d-2 * 75d3 / 184)

            collisionProbabilityByZ = 1 - e ** (-lambda * dt)

            write(*,*) 184d-3 * z, collisionProbabilityByZ

        end


        !!!*
        ! ピッチ角を与え、初期の粒子を得る
        ! @function createInitialParticle
        ! @param {double} angle ピッチ角[°]
        ! @param {logical} ignoreMF ミラー力を無視するかどうか
        ! @return {Particle} 粒子
        !!!
        function createInitialParticle(angle, ignoreMF)
            type(Particle) createInitialParticle

            double precision initialBx, angle, initialVel(3)
            logical ignoreMF


            initialVel(1) = 0d0
            initialVel(2) = sin(angle * pi / 180) * V_LEN
            initialVel(3) = - cos(angle * pi / 180) * V_LEN

            createInitialParticle = Particle(initial_r, initialVel, ignoreMF)

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
        ! @param {logical} ignoreMF ミラー力を無視するかどうか
        ! @return {double(3)} 磁場
        !!!
        function B(r, ignoreMF)

            double precision B(3), r(3), B0, r0, z, xylen, Br
            logical ignoreMF

            z = r(3)
            r0 = 3.46d4! 地球表面の座標
            B0 = 1d0 ! 地球表面の磁場
            B(3) = - B0 * (1 + z / r0)**(-3)

            if (ignoreMF) then
                B(1) = 0d0
                B(2) = 0d0

            else
                xylen = sqrt(r(1)**2 + r(2)**2)
                xyangle = atan2(r(2), r(1))

                Br = B(3) * 3 * xylen / (2 * (r0 + z))

                B(1) = Br * cos(xyangle)
                B(2) = Br * sin(xyangle)
            endif

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

            acceleration = pt%q / pt%m * ( El(pt%r) + cross( pt%v, B(pt%r, pt%ignoreMF) ))

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

            cyclotronPeriod = 2 * pi * pt%m / abs(pt%q) / vlen(B(pt%r, pt%ignoreMF))

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

        ! 配列の値の総和 を、 denomiで割った値を返す
        function arrMean(arr, length, denomi)

            integer length, denomi

            double precision, dimension(length):: arr
            double precision arrMean, total

            if (denomi == 0) then
                arrMean = 0d0

            else

                total = 0d0

                do i = 1, length
                    total = total + arr(i)
                end do

                arrMean = total / denomi

            endif
        end function




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
            type(Particle) pt, p1, p2, p3, p4
            double precision time
            double precision rungekuttaParticle(3), k1(3), k2(3), k3(3), k4(3)
            double precision d1(3), d2(3), d3(3), d4(3) ! d1〜d3 微小時間で移動した距離

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
            d1 = (pt%v + pt%v + k1 * h / 2) * h / 2 / 2 ! 台形
            p1 = Particle(pt%r + d1, pt%v + k1 * h / 2, pt%ignoreMF)

            k2 = accel(p1, time + h / 2)
            d2 = (p1%v + p1%v + k2 * h / 2) * h / 2 / 2 ! 台形
            p2 = Particle(p1%r + d2, p1%v + k2 * h / 2, pt%ignoreMF)

            d3 = (pt%v + pt%v + k2 * h / 2) * h / 2 / 2! 台形
            p3 = Particle(pt%r + d3, pt%v + k2 * h / 2, pt%ignoreMF)
            k3 = accel(p3, time + h / 2)

            d4 = (pt%v + pt%v + k3 * h) * h / 2 ! 台形
            p4 = Particle(pt%r + d5, pt%v + k3 * h, pt%ignoreMF)

            k4 = accel(p4, time + h)

           rungekuttaParticle = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

           return
        end
    !end contains
end




























program main

    use CyclotronWithRungeKutta

    implicit none

    integer, parameter:: TRIAL_NUM = 3

    ! 100回のシミュレーション結果
    type MultiSimulationScope

        ! 最初のピッチ角
        double precision initialPA

        ! 終了するまでのイテレーション回数
        integer, dimension(TRIAL_NUM):: iterationsArr = 0d0

        ! 衝突した回数
        integer:: collidedNum = 0

        ! 衝突した位置(z)の配列
        double precision, dimension(TRIAL_NUM):: collidedZArr = 0d0

        ! 衝突したピッチ角の配列
        double precision, dimension(TRIAL_NUM):: collidedPAArr = 0d0

        ! ミラーポイントに達した回数
        integer:: boundedNum = 0

        ! ミラーポイント(z)の配列
        double precision, dimension(TRIAL_NUM):: mirrorZArr = 0d0

    end type MultiSimulationScope

    type(MultiSimulationScope) msim
    type(SimulationScope) sim

    ! ミラー力を無視するかどうか
    logical:: ignoreMF = .true.

    double precision:: angle = 0

    integer i,j

    do i = 0, 0

        !angle = 0.02 * i + 72
        angle = i

        msim = MultiSimulationScope(angle)

        do j = 1, TRIAL_NUM
            sim = run(angle, ignoreMF)

            msim%iterationsArr(j) = sim%iterations
            if (sim%collided) then
                msim%collidedNum = msim%collidedNum + 1
                msim%collidedZArr(j) = sim%collidedR(3)
                msim%collidedPAArr(j) = sim%collidedPA
            endif

            if (sim%bounded) then
                msim%boundedNum = msim%boundedNum + 1
                msim%mirrorZArr(j) = sim%mirrorPoint(3)
            endif
        end do

        !write(*, *) angle, msim%boundedNum, msim%collidedNum, &
        !    & 184d-3 * arrMean(msim%mirrorZArr, TRIAL_NUM, msim%boundedNum), &
        !    & 184d-3 * arrMean(msim%collidedZArr, TRIAL_NUM, msim%collidedNum), &
        !    & arrMean(msim%collidedPAArr, TRIAL_NUM, msim%collidedNum)

    end do


end
