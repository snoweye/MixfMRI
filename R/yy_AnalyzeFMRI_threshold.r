### Copy some functions from "AnalyzeFMRI/R/threshold.R" since
### "AnalyzeFMRI" is archived by CRAN.

Threshold.Bonferroni <- function(p.val, n, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL) {

    ## calculate the Bonferroni threshold for n iid tests to give a p-value of p.val
    ## type specifies the univariate distribution of the test statistics under consideration
    
    if(type[1] == "Normal") return(qnorm(1 - p.val / n))
    if(type[1] == "t") return(qt(1 - p.val / n, df = df1))
    if(type[1] == "F") return(qf(1 - p.val / n, df1 = df1, df2 = df2))

}


Threshold.RF <- function(p.val, sigma, voxdim = c(1, 1, 1), num.vox, type = c("Normal", "t"), df = NULL) {

    ## calculates the Random Field theory threshold to give a p-value of p.val
    ## type specifies the marginal distribution of the field
    
    EC.func <- function(u , sigma, voxdim, num.vox, p.val) (EC.3D(u, sigma, voxdim, num.vox, type[1], df) - p.val)^2

    threshold <- optimize(f = EC.func, interval = c(2,10), sigma = sigma, voxdim = voxdim, num.vox = num.vox, p.val = p.val, tol = 1e-9)$minimum
    
    return(threshold)
}


EC.3D <- function(u, sigma, voxdim = c(1, 1, 1), num.vox, type = c("Normal", "t"), df = NULL) {

    ## The Expectation of the Euler Characteristic for a 3D Random Field above a threshold u
    ## type specifies the marginal distribution of the field
    
    V <- prod(voxdim) * num.vox
    S <- sqrt(det(solve(2 * sigma)))
    EC <- switch(type[1],
                 Normal = V * S * (u^2 - 1) * exp(-u^2 / 2) / (4 * pi^2),
                 t = V * S * (1 + (u^2 / df))^(-(df - 1) / 2) * (u^2 * (df - 1) / df - 1) / (4 * pi^2)
                 )
    return(EC)
}


Threshold.FDR <- function(x, q, cV.type = 2, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL) {

    ## calculates the FDR threshold for a vector of p-values in x
    ## q specifies the desired FDR
    ## cV specfies the type of FDR threshold used (See Genovese et al. (2002))
    
    if(type[1] == "Normal") p <- sort(1 - pnorm(x))
    if(type[1] == "t") p <- sort(1 - pt(x, df = df1))
    if(type[1] == "F") p <- sort(1 - pf(x, df1 = df1, df2 = df2))

    V <- length(p)
    cV <- switch(cV.type, 1, log(V) + 0.5772)
    
    i <- 1
    while (p[i] <= (i * q) / (V * cV)) i <- i + 1
    i <- max(i - 1, 1)

    if(type[1] == "Normal") thr <- qnorm(1 - p[i])
    if(type[1] == "t") thr <- qt(1 - p[i], df = df1)
    if(type[1] == "F") thr <- qf(1 - p[i], df1 = df1, df2 = df2)

    return(thr)

}


