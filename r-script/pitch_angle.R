pa <- function(pitchAngle_eq) {

    B_300 = 0.00005075586826
    B_eq  = 0.0000001255362951

    return (asin((B_300 / B_eq)^(1/2) * sin(pitchAngle_eq * pi / 180)) * 180 / pi )

}
