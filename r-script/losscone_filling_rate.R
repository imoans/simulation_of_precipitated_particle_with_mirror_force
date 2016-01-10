lossconeDensity <- function(mean, sd) {

    return (norm(mean, sd))
}

lossconeDistribution <- function(x) {

    return (1/2*(1+erf(x-mean)/sqrt(2*dispersion^2)))

}
