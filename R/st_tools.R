#' Check data format for stresstest
#'
#' @export
#'
stData <- function(data, Y = "Y", A = "A", fml.W = NULL) {
    if (inherits(data, "STDATA"))
        return(data);

    stopifnot(inherits(data, "data.frame"));

    cnames <- colnames(data);
    wnames <- cnames[!(cnames %in% c(Y,A))];

    ## Y an A
    stopifnot(all(c(Y, A) %in% cnames));
    data$Y <- data[[Y]];
    data$A <- data[[A]];
    stopifnot(is.numeric(data$Y));
    isbinary <- all(Y %in% c(0,1));


    trts <- unique(data$A);
    stopifnot(2 == length(trts));
    data$A <- as.numeric(trts[1] == data$A);

    if (!is.null(fml.W)) {
        fml.W  <- formula(fml.W);
        w.vars <- all.vars(fml.W);
        stopifnot(all(w.vars %in% wnames));
    } else {
        fml.W <- paste("~", paste(wnames, collapse = "+"));
        fml.W <- as.formula(fml.W);
    }

    rst <- list(data = data, fml.W = fml.W, isbinary = isbinary);
    class(rst) <- "STDATA";

    return(rst)
}


