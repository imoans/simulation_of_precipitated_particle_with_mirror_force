mf <- read.table('data/mf_trial100.d')
mf <- subset(mf, mf[,5] > 0)

nmf <- read.table('data/nmf_trial100.d')
nmf <- subset(nmf, nmf[,5] > 0)

mirrorPoint <- function(angle300) {

    Bmirror = 0.00005075586826 / (sin(angle300 * pi / 180))^2

    C = 8050000000000000

    dist = (C / Bmirror * sqrt(1 + 3* (sin(66 * pi / 180))^2 ))^(1/3)

    return (dist * 10^(-3) - 6371)


}

altiPlot <- function() {
    cols <- c(rgb(1,0,0), rgb(0,0,1), rgb(0.3, 0.3, 0.3))
    lwds <- c(3, 3, 1)
    labels <- c('mirror force ON', 'mirror force OFF', 'mirror point')

    plot(mf[,1], mf[,5], xlim = c(0, 90), ylim=c(0,300), type='l', col=rgb(1,0,0), xlab='pitch angle at 300km [°]', ylab='altitude [km]', lwd=3)
    par(new=T)
    plot(nmf[,1], nmf[,5], xlim = c(0, 90), ylim=c(0,300), type='l', col=rgb(0,0,1), xlab='', ylab='', lwd=3)
    par(new=T)
    plot(mirrorPoint, xlim = c(0, 90), ylim=c(0, 300), xlab='', ylab='', col=rgb(0.3, 0.3, 0.3))


    legend('topleft', legend = labels, col= cols, lty = c(1,1,1), lwd = lwds)
}
