
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                              FUNCTIONS
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
get.alpha <- function(n, beta, delta, sigma, interval=c(0,1)) {
    f.diff <- function(alpha) {
        n.alpha <- get.planned.size(delta, sigma, alpha, beta);
        n.alpha - n
    }

    rst <- tryCatch({
        uniroot(f.diff, interval=interval)$root;
    }, error = function (e) {
        print(e);
        NA;
    })

    rst
}

##pretty print
pp <- function(x) {
    format(x, digits = 4);
}

##get z.alpha, z.beta
get.z.ab <- function(alpha=0.05, beta=0.2) {
    z.alpha <- qnorm(1 - alpha/2);
    z.beta  <- qnorm(1 - beta);
    c(z.alpha, z.beta);
}


##sample size calculation per group
##rr: randomization ratio n1/n0
get.planned.size <- function(delta, sigma, alpha=0.05, beta=0.2, rr=1) {

    zab     <- get.z.ab(alpha=alpha, beta=beta);
    z.alpha <- zab[1];
    z.beta  <- zab[2];
    N       <- (z.alpha + z.beta)^2/delta^2;
    N       <- N * sigma^2 * (1+1/rr);

    ##return
    ceiling(N);
}


##use summary information to get zscore
get.zscore <- function(y1, y0) {
    eff.size <- mean(y1) - mean(y0);
    eff.sd   <- sqrt(var(y1)/length(y1) + var(y0)/length(y0));
    zscore   <- eff.size / eff.sd;
    pval     <- 2*(min(1 - pnorm(zscore), pnorm(zscore)));
    c(ey=eff.size, sd=eff.sd, z=zscore, pval=pval)
}


##-------------------------------------------------------------------
##
##              SCENARIO I
##
##-------------------------------------------------------------------
simu.trial <- function(delta1, delta2, n.perarm, pi1, sigma) {
    s1 <- rbinom(2, size = n.perarm, pi1);

    ##arm1
    y11 <- rnorm(s1[1],          delta1, sigma);
    y12 <- rnorm(n.perarm-s1[1], delta2, sigma);

    ##arm0
    y01 <- rnorm(s1[2],          0, sigma);
    y02 <- rnorm(n.perarm-s1[2], 0, sigma);


    ##summary
    zs.1 <- get.zscore(y11, y01);
    zs.2 <- get.zscore(y12, y02);
    zs   <- get.zscore(c(y11, y12), c(y01, y02));

    c(zs.1["ey"],   zs.2["ey"],   zs["ey"],
      zs.1["pval"], zs.2["pval"], zs["pval"],
      zs.1["z"],    zs.2["z"],    zs["z"]);
}


simu.single.setting <- function(pi1,
                                n.rep,
                                alpha,
                                delta1,
                                delta2,
                                n.perarm,
                                sigma,
                                method=NULL, n.cores=4,
                                rej.regions = NULL) {

    ##true
    te <- c(delta1, delta2, delta1*pi1 + delta2*(1-pi1));

    ##simulation
    rst <- parallel::mclapply(1:n.rep,
                              function(x) {
                         simu.trial(delta1, delta2, n.perarm, pi1, sigma);
                     }, mc.cores=n.cores);

    rst <- t(simplify2array(rst));

    ##summary
    srst <- NULL;
    for (i in 1:3) {
        cur.bias <- mean(rst[,i]) - te[i];
        cur.mse  <- mean((rst[,i]-te[i])^2);
        cur.sd   <- sd(rst[,i]);
        srst     <- c(srst, cur.bias, cur.mse, cur.sd);
    }

    ## three hypothesis
    if ("mb" == method) {
        s.rej <- apply(rst, 1, function(x) {
            cur.r <- NULL;
            for (j in c("1", "2", "0")) {
                cur.r <- c(cur.r,
                           IfInRejRegion(rej.regions[[j]], x[c(9,7:8)]));
            }
            cur.r
        })
        s.rej <- t(s.rej);
    } else {
        s.rej <- NULL;
        for (i in 1:3) {
            cur.rej <- rst[,3+i] < alpha[i];
            s.rej   <- cbind(s.rej, cur.rej);
        }
    }

    srst <- c(srst, apply(s.rej, 2, mean));

    ## any hypothesis
    srst <- c(srst, mean(apply(s.rej, 1, function(x) {any(x)})));
    srst;
}
