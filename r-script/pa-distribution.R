
paDensity <- function(x) {
    a <- 0.0375 / 2
    b <- -20.25
    alpha = 1 / 2.151235e-16
    return ( alpha * (10^(a * x^3 + b)) )
}

i <- 0:29799
den <- unlist(lapply(i * 0.0001, paDensity))

sumden <- cumsum(den)


# sumdenの逆関数を求める
randomDen <- function(uni) {

    return (which( abs(sumden - uni) == min(abs(sumden - uni)) )* 0.0001)

}
