##-----------------------------------------------------------------------------------
##                     run Shiny Gui
##-----------------------------------------------------------------------------------
#' Run Web-Based StressTest application
#'
#' Call Shiny to run \code{stresstest} as a web-based application
#'
#'
#'
#' @export
#'
stShiny <- function() {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Shiny needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("shinythemes", quietly = TRUE)) {
        stop("shinythemes needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("DT", quietly = TRUE)) {
        stop("shinythemes needed for this function to work. Please install it.",
             call. = FALSE)
    }


    appDir <- system.file("shiny", package = "stresstest")
    if (appDir == "") {
        stop("Could not find Shiny directory. Try re-installing `stresstest`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal");
}


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


##optimized and typical bonferonni
r.getn.bf <- function(delta1, delta2, sigma, pi1,
                      method = c("bon", "bon.opt", "mb"),
                      alpha=0.05,
                      beta=0.2,
                      alpha.min   = 0.0001,
                      alpha.input = NULL) {

    method <- match.arg(method);

    ##get sample sizes for all populations
    f.bf <- function() {
        delta  <- delta1 * pi1 + delta2 * (1-pi1);
        size.1 <- get.planned.size(delta1, sigma, alpha=alpha/3, beta=beta);
        size.2 <- get.planned.size(delta2, sigma, alpha=alpha/3, beta=beta);
        size.c <- get.planned.size(delta,  sigma, alpha=alpha/3, beta=beta);

        c(max(ceiling(size.1/pi1),
              ceiling(size.2/(1-pi1)),
              size.c),
          size.1,
          size.2,
          size.c);
    }

    ##optimized bonferonni
    f.n.bo <- function(cur.n) {
        cur.n1 <- round(cur.n * pi1);
        cur.n2 <- cur.n - cur.n1;

        alpha1 <- get.alpha(cur.n1, beta, delta1, sigma);
        alpha2 <- get.alpha(cur.n2, beta, delta2, sigma);
        alphac <- get.alpha(cur.n, beta, delta, sigma);
        a3     <- c(alpha1, alpha2, alphac);

        if (any(is.na(a3)) |
            any(a3 < alpha.min) |
            sum(a3, na.rm=TRUE) > alpha)
            return(NULL);

        c(cur.n1, cur.n2, cur.n, a3);
    }

    delta  <- delta1 * pi1 + delta2 * (1-pi1);
    init.n <- f.bf();
    if ("bon" == method) {
        ##naive bonferonni adjustment
        rst <- list(N      = init.n[1],
                    alpha3 = rep(alpha/3,3),
                    Nall   = init.n[2:4]);
    } else if ("bon.opt" == method) {
        init.n.1  <- init.n[1];
        cur.range <- c(0, init.n.1);
        cur.n     <- ceiling(mean(cur.range));
        while (cur.n < cur.range[2]) {
            cur.rst <- f.n.bo(cur.n);
            if (is.null(cur.rst)) {
                cur.range[1] <- cur.n;
            } else {
                cur.range[2] <- cur.n;
                last.rst     <- cur.rst;
            }

            ##mid point
            cur.n <- ceiling(mean(cur.range));
        }

        rst <- list(N      = last.rst[3],
                    alpha3 = last.rst[4:6],
                    Nall   = c(last.rst[1:3]));

    } else if ("mb" == method) {
        cur.n <- SampleSizeMB.fix.alpha(pi.1    = pi1,
                                        delta.1 = delta1,
                                        delta.2 = delta2,
                                        sigma   = sigma,
                                        alpha   = alpha,
                                        beta    = beta,
                                        alpha.input = alpha.input)$sample.size;
        rej.regions <- NULL;
        for (j in c("0", "1", "2")) {
            rej.regions[[j]] <- RejRegion(alpha.input, j);
        }

        rst <- list(N      = cur.n,
                    alpha3 = alpha.input,
                    Nall   = c(round(cur.n * pi1),
                               cur.n - round(cur.n * pi1),
                               cur.n),
                    rej.regions = rej.regions);
    }
    rst$method <- method;
    rst$pars   <- c(beta   = beta,
                    alpha  = alpha,
                    delta1 = delta1,
                    delta2 = delta2,
                    sigma  = sigma,
                    pi1    = pi1);

    rst
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

simu.single.setting <- function(n.rep, alpha, delta1, delta2,
                                n.perarm, pi1, sigma,
                                method=NULL, n.cores=4, rej.regions = NULL) {

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


##-------------------------------------------------------------------
##
##              PLOT
##
##-------------------------------------------------------------------

plot.rst <- function(simu.rst, y.var, x.var = "Pi", n = 10) {

    f.con <- function(v, simu.rst) {
        if (is.factor(simu.rst[[v]])) {
            simu.rst[[v]] <- as.numeric(levels(simu.rst[[v]]))[simu.rst[[v]]];
        }
        simu.rst
    }

    simu.rst <- f.con(x.var, simu.rst);
    simu.rst <- f.con(y.var, simu.rst);


    p <- ggplot2::ggplot(simu.rst, ggplot2::aes_string(x.var, y.var, color="Design", group = "Design")) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = n)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = n));
    p
}
