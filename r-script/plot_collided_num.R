# 文字列連結
L <- function(..., f) {
    if (missing(f)) 
        f <- paste(rep("%s", length(c(...))), collapse = "")
    sprintf(fmt = f, ...)
}

showAltitudeHists <- function(dir, mf = T) {

    basedir = 'data/collidedHeights/'
    fileName = ifelse(mf == T, 'mf.d', 'nmf.d')

    normHs <- read.table(L(basedir, dir, '/', fileName))
    normHsEx0 = subset(normHs, normHs[,2] > 0)

    uniformHs <- read.table(L(basedir, 'uniform/', fileName))
    uniformHsEx0 = subset(uniformHs, uniformHs[,2] > 0)

    dataSize = min(dim(uniformHsEx0)[1], dim(normHsEx0)[1])

    uniformHsEx0 = uniformHsEx0[1:dataSize,]
    normHsEx0    = normHsEx0[1:dataSize,]

    print(dataSize)
    print(wilcox.test(uniformHsEx0[,2], normHsEx0[,2]))
    print(t.test(log(uniformHsEx0[,2]), log(normHsEx0[,2])))

    hist(normHsEx0[,2],
         breaks = 30,
         col    = "#ff00ff40",
         border = "#ff00ff",
         xlab   = 'altitude [km]',
         ylab   = 'number of electron',
         main   = L('histgram of collided altitude (N=', dataSize, ' ,', ifelse(mf, 'with', 'ignoring'), ' Mirror Force)'),
         ylim   = c(0,120),
         xlim   = c(90, 180)
    )

    par(new=T)


    hist(uniformHsEx0[,2],
         breaks = 30,
         col    = "#0000ff40",
         border = "#0000ff",
         xlab   = '',
         main   = '',
         ylab   = '',
         ylim   = c(0,120),
         xlim   = c(90, 180)
    )


    legend('topright',
        legend = c('normal distribution', 'uniform distribution'),
        col    = c('#ff00ff40', '#0000ff40'),
        lty    = c(1,1),
        lwd    = c(5, 5)
    )
}
