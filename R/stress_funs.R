#' Get design parameters
#' @export
stDesign <- function(delta1, delta2, sigma, pi1,
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
                                        beta    = 1-beta,
                                        alpha.input = alpha.input)$sample.size;
        ##per arm
        cur.n <- ceiling(cur.n/2);
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


#' Simulate trials
#' @export
stSimu <- function(sizen,
                   all.pi,
                   progress=NULL,
                   cnames = c("Design",
                              "Pi",
                              "Bias1", "MSE1", "SD1",
                              "Bias2", "MSE2", "SD2",
                              "BiasC", "MSEC", "SDC",
                              "Rej1", "Rej2", "RejC", "RejAny"),
                   n.rep  = 1000,
                   delta1 = 0,
                   delta2 = 0,
                   n.cores = 4) {
    rst <- NULL;
    mtd <- NULL;
    for (j in 1:length(sizen)) {
        for (i in 1:length(all.pi)) {
            if (!is.null(progress))
                progress$set(value  = (length(all.pi)*(j-1)+i)/length(sizen)/length(all.pi),
                             detail = paste(names(sizen)[[j]],
                                            ", ",
                                            "pi=", all.pi[i], sep=""));

            cur.rst <- simu.single.setting(all.pi[i],
                                           n.rep       = n.rep,
                                           delta1      = delta1,
                                           delta2      = delta2,
                                           n.cores     = n.cores,
                                           alpha       = sizen[[j]]$alpha3,
                                           n.perarm    = sizen[[j]]$N,
                                           sigma       = sizen[[j]]$pars['sigma'],
                                           method      = sizen[[j]]$method,
                                           rej.regions = sizen[[j]]$rej.regions);

            rst  <- rbind(rst, c(all.pi[i], cur.rst));
        }
        mtd  <- c(mtd, rep(names(sizen)[j], length(all.pi)));
    }

    colnames(rst)    <- cnames[-1];
    rst              <- data.frame(rst);
    rst[[cnames[1]]] <- mtd;
    rst
}



##-------------------------------------------------------------------
##
##              PLOT
##
##-------------------------------------------------------------------

#' Plot Results
#'
#' @export
stPlot <- function(simu.rst, y.var, x.var = "Pi", n = 10) {
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


#' Plot Results
#'
#' @export
stPlotStress <- function(rst.d, power=0.8, ..., cols=c("blue", "red", "brown")) {
    ss <- c("C" = "Power: Combined Pop. Null Hyp.",
            "1" = "Power: Subpop. 1 Null Hyp.",
            "2" = "Power: Subpop. 2 Null Hyp.");

    vec.pi     <- unique(rst.d$Pi);
    vec.method <- unique(rst.d$Design);
    min.power  <- min(c(rst.d$Rej1, rst.d$Rej2, rst.d$RejC));

    par(mfrow=c(1,3), ...);
    for (i in c("C", "1", "2")) {
        plot(NULL, xlim=range(vec.pi), ylim=c(min.power,1),
             xlab = "Prop. of Subpop. 1",
             ylab="Power",
             main = ss[i]);

        cur.lst <- NULL;
        for (j in 1:length(vec.method)) {
            cur.sub <- subset(rst.d, Design == vec.method[j]);
            cur.y   <- as.numeric(cur.sub[, paste("Rej", i, sep="")]);
            cur.p   <- as.numeric(cur.sub[, "Pi"]);
            cur.lo  <- loess(cur.y ~ cur.p);

            cur.lst[[j]] <- list(cur.y=cur.y, cur.p=cur.p, cur.lo=cur.lo);
        }

        ##lines(cur.p1, cur.y1, type="l", col="blue");
        ##lines(cur.p2, cur.y2, type="l", col="black");
        lines(c(0,1), c(power, power), lty=3, lwd=2, col="gray");

        ##polygon
        for (j in 1:length(vec.method)) {
            cur.p  <- cur.lst[[j]]$cur.p;
            cur.lo <- cur.lst[[j]]$cur.lo;
            p.xx   <- seq(min(vec.pi), max(vec.pi), 0.001);
            p.yy   <- predict(cur.lo, p.xx);

            inx    <- which(p.yy < power);
            if (0 == length(inx))
                next;

            ss.xx  <- p.xx[inx];
            ss.yy  <- p.yy[inx];
            p.xx   <- c(min(ss.xx), ss.xx, max(ss.xx));
            p.yy   <- c(0, ss.yy, 0);
            polygon(p.xx, p.yy, col=cols[j], density=30, border=NULL);

            ##text(mean(ss.xx), 0.6, "Stress-Test", cex=1);
            ##text(mean(ss.xx), 0.56, "Failure Zone", cex=1);
            ##text(mean(ss.xx), 0.52, "(<80% Power)", cex=1);
        }

        ##curve
        for (j in 1:length(vec.method)) {
            cur.p  <- cur.lst[[j]]$cur.p;
            cur.lo <- cur.lst[[j]]$cur.lo;
            lines(cur.p, predict(cur.lo, cur.p), lty=1, lwd=2, col=cols[j]);
        }

        if ("C" == i)
            legend("bottomleft", legend=vec.method,
                   lty = 1,
                   col = cols, lwd=2, cex=1.1, bty="n");
    }
}


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

