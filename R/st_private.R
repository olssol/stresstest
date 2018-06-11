#' Double Robust function
#'
#'
#'
get.dr.pred <- function(data, gw) {
    fml.y <- update(data$fml.W, "Y ~ .");
    A     <- data$data$A;
    rst   <- NULL;
    for (a in 0:1) {
        inx.a   <- which(a == A);
        weights <- gw[inx.a];
        if (0 == a)
            weights <- 1 - weights;
        weights <- 1/weights;

        cur.d         <- data$data[inx.a,];
        cur.d$weights <- weights;
        if (data$isbinary) {
            glm.y  <- glm(fml.y, data = cur.d, weights = weights, family = "binomial");
            pred.y <- predict(glm.y, newdata = data$data, type="response");
        } else {
            lm.y   <- lm(fml.y, data = cur.d, weights = weights);
            pred.y <- predict(lm.y, newdata = data$data);
        }
        rst <- cbind(rst, pred.y);
    }

    rst
}
