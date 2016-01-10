# lambda: サイクロトロン1周期のうちに何回衝突するか
lambda <- function(x) {
    z = x / 0.184
    return ( 10 ^ (-1.288e-2 * z + 1.288e-2 * 75e3 / 184) )
}

colliP <- function(z) {

    dt = 0.05
    return (1 - exp(-lambda(z) * dt))
}


