B <- function(z) {

    x = z / 0.184

    r0 = 34600
    B0 = 1

    return (B0 * (1 + x / r0)^(-3))

}

B(300)
