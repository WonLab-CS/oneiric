tinkerbell <- function(time = 6000,
    a = 0.9,
    b = -0.6013,
    c = 2,
    d = 0.5,
    x_0 = -0.72,
    y_0 = -0.64){
    x <- c(x_0, rep(0, time - 1))
    y <- c(y_0, rep(0, time - 1))

    for (i in seq(1, time - 1)) {
        x[i + 1] <- x[i]^2 - y[i]^2 + a * x[i] + b * y[i]
        y[i + 1] <- 2 * x[i] * y[i] + c * x[i] + d * y[i]
    }
    return(data.frame(x,y))
}