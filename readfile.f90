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
    use ReadFileModule

    integer,parameter :: initVelNum = 10
    double precision initVels(initVelNum), initV

    initVels = readFile('init_vel.tsv', initVelNum)

    print *, initVels

end
