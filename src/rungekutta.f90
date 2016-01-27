!!!*
! solving cyclotron motion with Runge-Kutta method
! @module CyclotronWithRungeKutta
!
! xyz coordinate system
! oringin: center of Earth
! x: parallel to the equatorial plane [m]
! z: zenith direction [m]
!!!
module CyclotronWithRungeKutta

    ! Pi
    double precision, parameter:: pi = 3.14159265358979323d0

    ! base of natural logarithm
    double precision, parameter:: e  = 2.71828182845904524d0

    ! h: step size of time
    double precision, parameter :: h = 5d-2

    ! initial_r: initial position of a particle (184m as one unit, z = 300[km])
    double precision, parameter :: initial_r(3) = (/0d0, 0d0, 1.63d3/)

    ! KM_PER_UNIT: distance unit [km]
    ! distance unit is speed of light multiplied by cyclotron period at surface of Earth
    double precision, parameter :: KM_PER_UNIT = 184d-3

    ! V_LEN: absolute value of particle velocity
    double precision, parameter :: V_LEN = 6.24d-2

    ! particle
    type Particle
        ! r: position, v: velocity
        double precision :: r(3), v(3)

        ! ignoreMF: whether Mirror Force is ignored
        logical ignoreMF

        ! q: charge (unit is charge of electron [Q])
        double precision :: q = -1d0

        ! m: mass (unit is mass of electron [kg])
        double precision :: m = 1d0

    end type Particle


    ! SimulationScope: variables and results involved single particle simulation
    type SimulationScope

        ! pitch angle at initial position [°]
        double precision :: initialPA

        ! number of iterations until exit this program
        integer iterations

        ! number of iterations
        integer:: ITERATION_TIMES = 20000000

        ! whether collision has occured
        logical:: collided = .false.

        ! collision position when collision occured
        double precision collidedR(3)

        ! pitch angle when collision occured [°]
        double precision collidedPA

        ! whether a particle reaches mirror point
        logical:: bounded = .false.

        ! mirror point
        double precision mirrorPoint(3)

    end type SimulationScope


    contains


        !!!*
        ! calcurate position of particle in each step with Runge-Kutta method
        ! return result of simulation
        !
        ! @function run
        ! @param {double} initialPA pitch angle of initial position
        ! @param {logical} ignoreMF whether Mirror Force is ignored
        ! @return {SimulationScope} simscope
        ! @public
        !!!
        function run(initialPA, ignoreMF) result(sim)

            ! t: time
            double precision :: t = 0

            double precision initialPA, prevVz
            logical ignoreMF

            type(SimulationScope) sim
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

                ! exit when a particle reaches initial position
                elseif (electron%r(3) > initial_r(3)) then
                    exit

                ! exit when a particle reaches surface of Earth
                elseif (electron%r(3) < 0) then
                    exit
                endif

                t = t + h

                prevVz = electron%v(3)
                call update(electron, t)

                ! save mirror point
                if (prevVz <= 0 .AND. electron%v(3) >= 0) then
                    sim%bounded = .true.
                    sim%mirrorPoint(1) = electron%r(1)
                    sim%mirrorPoint(2) = electron%r(2)
                    sim%mirrorPoint(3) = electron%r(3)
                endif
            end do

        end function


        !!!*
        ! update position of a particle and magnitude of magnetic field
        ! @subroutine update
        ! @param {Particle} pt particle
        ! @param {double} time time
        !!!
        subroutine update(pt, time)

            type(Particle) pt
            double precision time, v(3), r(3), runge(3)

            runge = rungekuttaParticle(pt, time, acceleration)

            pt%v = pt%v + runge
            pt%r = pt%r + pt%v * h

        end subroutine

        !!!*
        ! calculate kinetic energy by velocity of a particle
        ! @function kineticEnergy
        ! @param {Particle} particle
        ! @return {double} kinetic energy
        !!!
        function kineticEnergy(pt)
            double precision kineticEnergy, v_len
            type(Particle) pt

            kineticEnergy = 0.5d0 * pt%m * vlen(pt%v) ** 2
        end function


        !!!*
        ! determine whether collision with particles in a neutral atmosphere has occurred by height(z value) of particle
        ! @function collisionOccurred
        ! @param {Particle} pt particle
        ! @param {double} dt minor change of time
        ! @return {logical} whether collision has occurred
        !!!
        function collisionOccurred(pt, dt)
            type(Particle) pt
            double precision randomNum, dt
            logical collisionOccurred

            call random_number(randomNum)
            collisionOccurred = (randomNum < collisionProbabilityByHeight(pt%r(3), dt))

        end function


        !!!*
        ! collision probability at certain height(z value)
        ! @function collisionProbabilityByHeight
        ! @param {double} z
        ! @param {double} dt minor change of time
        ! @return {double} probability
        !!!
        function collisionProbabilityByHeight(z, dt)
            double precision z, dt, collisionProbabilityByHeight, lambda

            lambda = 10 ** (-1.288d-2 * z + 1.288d-2 * 75d3 / 184)

            collisionProbabilityByHeight = 1 - e ** (-lambda * dt)

        end function


        !!!*
        ! generate a initial particle by pitch angle
        ! @function createInitialParticle
        ! @param {double} angle pitch angle[°]
        ! @param {logical} ignoreMF whether Mirror Force is ignored
        ! @return {Particle} particle
        !!!
        function createInitialParticle(angle, ignoreMF) result(pt)
            type(Particle) pt

            double precision initialBx, angle, initialVel(3)
            logical ignoreMF


            initialVel(1) = 0d0
            initialVel(2) = sin(angle * pi / 180) * V_LEN
            initialVel(3) = - cos(angle * pi / 180) * V_LEN

            pt = Particle(initial_r, initialVel, ignoreMF)

        end function


        !!!*
        ! get pitch angle of a particle
        ! @function pitchAngle
        ! @param {Particle} pt
        ! @return {double} pitchAngle [°]
        !!!
        function pitchAngle(pt)
            type(Particle) pt
            double precision absVelocity, pitchAngle

            absVelocity = vlen(pt%v)
            pitchAngle = 180 - acos(pt%v(3) / absVelocity) * 180 / pi

        end function



        !!!*
        ! get dipole magnetic field (divB = 0) by position of a particle
        ! @function B
        ! @param {double(3)} r position of a particle
        ! @param {logical} ignoreMF whether Mirror Force is ignored
        ! @return {double(3)} magnitude of magnetic field
        !!!
        function B(r, ignoreMF)

            double precision B(3), r(3), B0, r0, z, xylen, Br
            logical ignoreMF

            z = r(3)
            r0 = 3.46d4! position at surface of Earth
            B0 = 1d0 ! magnitude of magnetic field at surface of Earth
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

        end function


        !!!*
        ! get electric field at certian position
        ! @function El
        ! @param {double(3)} r position of a particle
        ! @return {double(3)} magnitude of electric field
        !!!
        function El(r)

            double precision El(3), r(3)

            El(1) = 0d0
            El(2) = 0d0
            El(3) = 0d0

        end function



        !!!*
        ! get acceleration
        ! @function acceleration
        ! @param {Particle} pt particle
        ! @param {double} time time
        ! @return {double(3)} acceleration
        !!!
        function acceleration(pt, time)

            type(Particle) pt
            double precision acceleration(3), time

            acceleration = pt%q / pt%m * ( El(pt%r) + cross( pt%v, B(pt%r, pt%ignoreMF) ))

        end function



        !!!*
        ! get cyclotron period with magnitude of magnetic field, mass and charge of a particle
        ! @function period
        ! @param {Particle} pt particle
        ! @return {double} cyclotron period
        !!!
        function cyclotronPeriod(pt)
            type(Particle) pt
            double precision cyclotronPeriod

            cyclotronPeriod = 2 * pi * pt%m / abs(pt%q) / vlen(B(pt%r, pt%ignoreMF))

        end function


        !!!*
        ! get length of three-dimensional vector
        ! @function vlen
        ! @param {double(3)} vec
        ! @return {double} length
        !!!
        function vlen(vec)
            double precision vlen, vec(3)
            vlen = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
        end function


        !!!*
        ! cross product
        ! @function cross
        ! @param {double(3)} x
        ! @param {double(3)} y
        ! @return {double(3)} result
        !!!
        function cross(x, y)
            double precision cross(3), x(3), y(3)

            cross(1) = x(2) * y(3) - x(3) * y(2)
            cross(2) = x(3) * y(1) - x(1) * y(3)
            cross(3) = x(1) * y(2) - x(2) * y(1)
        end function


        !!!*
        ! get value divided by sum of value of array in denomi
        ! @function arrMean
        !!!
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
        ! get particle velocity with Runge-Kutta method
        !
        ! @function rungekuttaParticle
        ! @param {Particle} pt particle
        ! @param {double} time time
        ! @param {function} accel acceleration (right-hand side of the first derivative differential equation)
        ! @return {double(3)} variation of speed
        !!!
        function rungekuttaParticle(pt, time, accel)
            type(Particle) pt, p1, p2, p3, p4
            double precision time
            double precision rungekuttaParticle(3), k1(3), k2(3), k3(3), k4(3)
            double precision d1(3), d2(3), d3(3), d4(3) ! d1〜d3 distance traveled by minor time

            ! The right-hand side of the first derivative differential equation
            interface
                function accel(p, t)
                    import :: Particle
                    type(Particle) p
                    double precision accel(3)
                    double precision t
                end function
            end interface

            k1 = accel(pt, time)
            d1 = (pt%v + pt%v + k1 * h / 2) * h / 2 / 2 ! area of trapezoid
            p1 = Particle(pt%r + d1, pt%v + k1 * h / 2, pt%ignoreMF)

            k2 = accel(p1, time + h / 2)
            d2 = (p1%v + p1%v + k2 * h / 2) * h / 2 / 2
            p2 = Particle(p1%r + d2, p1%v + k2 * h / 2, pt%ignoreMF)

            d3 = (pt%v + pt%v + k2 * h / 2) * h / 2 / 2
            p3 = Particle(pt%r + d3, pt%v + k2 * h / 2, pt%ignoreMF)
            k3 = accel(p3, time + h / 2)

            d4 = (pt%v + pt%v + k3 * h) * h / 2
            p4 = Particle(pt%r + d5, pt%v + k3 * h, pt%ignoreMF)

            k4 = accel(p4, time + h)

           rungekuttaParticle = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

           return
        end function
    !end contains
end module





!!!*
! various simulations of cyclotron motion
! @module CyclotronSimulation
!!!
module CyclotronSimulation

    use CyclotronWithRungeKutta

    implicit none

    integer, parameter:: TRIAL_NUM = 100

    ! 100 times of simulation results
    type MultiSimulationScope

        ! pitch angle at initial position
        double precision initialPA

        ! number of iterations until the end
        integer, dimension(TRIAL_NUM):: iterationsArr = 0d0

        ! number of collision
        integer:: collidedNum = 0

        ! array of collision position
        double precision, dimension(TRIAL_NUM):: collidedZArr = 0d0

        ! array of pitch angle when collision occurred
        double precision, dimension(TRIAL_NUM):: collidedPAArr = 0d0

        ! number of times paritcle reaches mirror point
        integer:: boundedNum = 0

        ! array of mirror point
        double precision, dimension(TRIAL_NUM):: mirrorZArr = 0d0

    end type MultiSimulationScope


    contains

        !!!*
        ! display result of 100 times simulation increasing pitch angle by 1°
        ! @param {logical} ignoreMF whether Mirror Force is ignored
        !
        ! @result
        ! {pitch angle}\t{nuber of times particle reachegs mirror point}\t{number of collision}\t{mean of mirror point [km]}\t{mean
        ! of collision height [km]}\t{mean of pitch angle when collision occurred[°]}
        !!!
        subroutine everyAngle(ignoreMF)

            type(MultiSimulationScope) msim
            type(SimulationScope) sim
            integer i,j
            logical ignoreMF
            double precision angle

            do i = 0, 90

                angle = i

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

                write(*, *) angle, msim%boundedNum, msim%collidedNum, &
                    & KM_PER_UNIT * arrMean(msim%mirrorZArr, TRIAL_NUM, msim%boundedNum), &
                    & KM_PER_UNIT * arrMean(msim%collidedZArr, TRIAL_NUM, msim%collidedNum), &
                    & arrMean(msim%collidedPAArr, TRIAL_NUM, msim%collidedNum)
            end do
        end subroutine


        !!!*
        ! get collision position by random pitch angle obtained from certian distribution
        ! @param {Integer} iterations number of iterations
        ! @param {logical} ignoreMF whether Mirror Force is ignored
        ! @result {pitch angle}\t{collision altitude[km]}
        !!!
        subroutine collidedHeights(iterations, ignoreMF)

            integer i, iterations
            logical ignoreMF
            double precision angle
            type(SimulationScope) sim

            do i = 1, iterations

                angle = uniformPA(0, 90)

                sim = run(angle, ignoreMF)

                if (sim%collided) then
                    write(*, *) angle, sim%collidedR(3) * KM_PER_UNIT
                else
                    write(*, *) angle, 0
                end if

            end do
        end subroutine


        !!!*
        ! get pitch angle from normal distribution [minPA, maxPA) by random
        !!!
        function uniformPA(minPA, maxPA) result(angle)
            integer minPA, maxPA
            double precision angle
            call random_number(angle)
            angle = minPA + angle * (maxPA - minPA)
        end function

end module


program main

    use CyclotronSimulation

    implicit none
    integer :: iterations = 1000
    logical :: ignoreMF = .false.

    integer :: k, argc
    character :: argv*10
    character(10) meanStr, sdStr

    argc = iargc()

    do k = 1, argc
        call getarg(k, argv)

        if (argv == '--ignoreMF') then
            ignoreMF = .true.
        end if
    end do

    write (*, *) '# ignoreMF:', ignoreMF

    call collidedHeights(iterations, ignoreMF)


end program
