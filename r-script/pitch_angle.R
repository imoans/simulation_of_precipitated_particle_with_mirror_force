pa <- function(pitchAngle_eq) {

    B_300 = 0.00005075586826
    B_500 = 0.000046451446
    B_eq  = 0.0000001255362951

    return (asin((B_500 / B_eq)^(1/2) * sin(pitchAngle_eq * pi / 180)) * 180 / pi )

}

x = 0

#i <- 0:28506
i <- 0:29799
transformed <- unlist(lapply(i * 0.0001, pa))

#plot(pa, 0, 3, xlab='pitch angle at equatorial plane', ylab='pitch angle at 300km')
