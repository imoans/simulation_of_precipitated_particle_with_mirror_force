collisionProbability <- function(z) {

    x = z / 0.184
    dt = 0.05
    return (1 - exp(-10^(-0.01288 * x + 0.01288 * 75000 / 184) * dt))

}

collisionProbability(300)
