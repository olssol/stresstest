
#' Naive estimator
#'
#' @export
#'
stEstNaive <- function(data, ...) {
    data   <- stData(data, ...);
    rst.lm <- lm("Y~A", data = data$data);

    coefficients(rst.lm)["A"]
}

#' ANCOVA estimator without interaction
#'
#'
#' @export
#'
stEstANCOVA <- function(data, ...) {
    data   <- stData(data, ...);
    fml.A  <- update(data$fml.W, "Y ~ A + .");
    rst.lm <- lm(fml.A, data = data$data);

    coefficients(rst.lm)["A"]
}


#' ANCOVA estimator with interaction
#'
#'
#' @export
#'
stEstANCOVAII <- function(data, ...) {
    data   <- stData(data, ...);
    fml.A  <- update(data$fml.W, "Y ~ A + . + A:.");
    rst.lm <- lm(fml.A, data = data$data);

    ys <- NULL;
    for (a in 0:1) {
        d0    <- data$data;
        d0$A  <- a;
        cur.y <- predict(rst.lm, newdata = d0);
        ys    <- c(ys, mean(cur.y));
    }

    ys[2] - ys[1];
}

#' Koch estimator
#'
#'
#' @export
#'
stEstKoch <- function(data, ...) {
    data      <- stData(data, ...);
    fml.noint <- update(data$fml.W, "~ . -1");
    wmat      <- model.matrix(fml.noint, data$data);

    ybar <- NULL;
    wbar <- NULL;
    vyw  <- 0;
    vww  <- 0;

    for (a in 0:1) {
        cur.inx <- which(a == data$data$A);
        cur.n   <- length(cur.inx);
        cur.y   <- data$data$Y[cur.inx];
        cur.w   <- wmat[cur.inx,,drop = FALSE];

        ybar    <- c(ybar, mean(cur.y));
        wbar    <- rbind(wbar, apply(cur.w, 2, mean));
        vyw     <- vyw + cov(cur.y, cur.w)/cur.n;
        vww     <- vww + cov(cur.w)/cur.n;
    }

    (ybar[2] - ybar[1]) - vyw %*% solve(vww) %*% (wbar[2,] - wbar[1,]);
}

#' Leon estimator
#'
#'
#' @export
#'
stEstLeon <- function(data, ...) {
    data   <- stData(data, ...);
    wmat   <- model.matrix(data$fml.W, data$data);

    A      <- data$data$A;
    Y      <- data$data$Y;
    inx.0  <- which(0 == A);
    inx.1  <- which(1 == A);
    Y0bar  <- mean(Y[inx.0]);
    Y1bar  <- mean(Y[inx.1]);
    gn1    <- mean(A);

    S3  <- apply(wmat, 1, function(x) {as.vector(x %*% t(x))});
    S3  <- apply(S3, 1, sum);
    S3  <- matrix(S3, nrow = ncol(wmat));

    S4  <- apply((A-gn1) * wmat, 2, sum);

    S1  <- apply((Y[inx.1] - Y1bar) * wmat[inx.1,], 2, sum);
    S2  <- apply((Y[inx.0] - Y0bar) * wmat[inx.0,], 2, sum);

    (Y1bar - Y0bar) - length(A) * (S1/sum(A)^2 + S2/sum(1-A)^2) %*% solve(S3) %*% S4;
}

#' TMLE estimator
#'
#'
#' @export
#'
stEstTMLE<- function(data, ...) {
    data    <- stData(data, ...);
    fml.1   <- update(data$fml.W, "Y ~ .");
    A       <- data$data$A;
    gn1     <- mean(A);

    data$data$weights <- A /gn1^2 + (1-A) / (1-gn1)^2;

    ##first lm
    lm.1 <- lm(fml.1, data = data$data, weights = weights);

    ## offset Y
    Y2 <- data$data$Y;
    Y2 <- Y2 - predict(lm.1)
    Y2 <- Y2 + coefficients(lm.1)["(Intercept)"];

    ##2nd lm
    lm.2 <- lm("Y2 ~ A", data = data.frame(Y2 = Y2, A = A));

    coefficients(lm.2)["A"]
}


#' IPW estimator
#'
#'
#' @export
#'
stEstIPW <- function(data, ...) {
    data  <- stData(data, ...);

    A     <- data$data$A;
    fml.A <- update(data$fml.W, "A ~ .");
    glm.A <- glm(fml.A, data = data$data, family = "binomial");
    gw    <- predict(glm.A, type="response");

    inx.1 <- 1 == A;
    ybar1 <- sum(data$data$Y[inx.1]/gw[inx.1]);
    ybar1 <- ybar1 / sum(1/gw[inx.1]);

    inx.0 <- 0 == A;
    ybar0 <- sum(data$data$Y[inx.0]/(1 - gw[inx.0]));
    ybar0 <- ybar0 / sum(1/(1 - gw[inx.0]));

    ybar1 - ybar0
}

#' Model standardization estimator
#'
#'
#' @export
#'
stEstMdlStd <- function(data, ...) {
    data  <- stData(data, ...);

    A     <- data$data$A;
    fml.y <- update(data$fml.W, "Y ~ .");

    ybar <- NULL;
    for (a in 0:1) {
        cur.d <- data$data[a == A,];
        if (data$isbinary) {
            glm.y  <- glm(fml.y, data = cur.d, family = "binomial");
            pred.y <- predict(glm.y, newdata = data$data, type="response");
        } else {
            lm.y   <- lm(fml.y, data = cur.d);
            pred.y <- predict(lm.y, newdata = data$data);
        }
        ybar  <- c(ybar, mean(pred.y));
    }

    ybar[2] - ybar[1]
}

#'  doubly-robust weighted least squares  estimator
#'
#'
#' @export
#'
stEstDRWLS <- function(data, ...) {
    data  <- stData(data, ...);

    fml.A <- update(data$fml.W, "A ~ .");
    glm.A <- glm(fml.A, data = data$data, family = "binomial");
    gw    <- predict(glm.A, type="response");

    pred.y <- get.dr.pred(data, gw);
    ybar   <- apply(pred.y, 2, mean);

    ybar[2] - ybar[1]
}

#'  PLEASE estimator
#'
#'
#' @export
#'
stEstPLEASE <- function(data, ...) {
    data   <- stData(data, ...);

    stopifnot(data$isbinary);

    fml.A  <- update(data$fml.W, "A ~ .");
    glm.A  <- glm(fml.A, data = data$data, family = "binomial");
    gw     <- predict(glm.A, type="response");

    ##1st double robust
    pred.y <- get.dr.pred(data, gw);

    ##augmentation
    data$data$u0 <- pred.y[,1] - mean(pred.y[,1]);
    data$data$u1 <- pred.y[,2] - mean(pred.y[,2]);
    fml.A.2      <- update(data$fml.W, "A ~ . + u0 + u1");
    glm.A.2      <- glm(fml.A.2, data = data$data, family = "binomial");
    gw.2         <- predict(glm.A.2, type="response");

    ##2nd double robust
    pred.y.2     <- get.dr.pred(data, gw.2);
    ybar         <- apply(pred.y.2, 2, mean);

    browser();

    ybar[2] - ybar[1]
}
