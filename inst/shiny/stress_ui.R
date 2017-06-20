require(stresstest);

##----------------------------------------------------------------------------
##                           CONSTANTS
##----------------------------------------------------------------------------
get.consts <- reactive({
    list(rstnames=c("Design",
                    "Pi",
                    "Bias1", "MSE1", "SD1",
                    "Bias2", "MSE2", "SD2",
                    "BiasC", "MSEC", "SDC",
                    "Rej1", "Rej2", "RejC", "RejAny"),
         designs=c("Bonferroni"           = "bon",
                   "Optimized Bonferroni" = "bon.opt",
                   "Optimized MB"         = "mb"
                   ))
})

##-------------------------------------------------------------
##           UI DEFINITIONS
##-------------------------------------------------------------

##show different type of messages
msg.box <- function(contents, type="info") {
    switch(type,
           info    = cls <- "cinfo",
           warning = cls <- "cwarning",
           success = cls <- "csuccess",
           error   = cls <- "cerror");
    rst <- '<div class="';
    rst <- paste(rst, cls, '">');
    rst <- paste(rst, contents);
    rst <- paste(rst, "</div>");
    HTML(rst);
}

##----------------------------------------------------------------------------
##                           PANELS
##----------------------------------------------------------------------------
tab.utility <- reactive({
    tabPanel("Note",
             withMathJax(includeHTML("www/paper_stress.html")))
})

tab.advanced <- reactive({
    tabPanel("Advanced",
             wellPanel(
                 fluidRow(
                     column(4,
                            h4("Uncertainty"),
                            sliderInput(inputId = "inDelP",
                                        label=withMathJax("$$\\Delta (\\pi_1)$$"),
                                        value = 0.1, min = 0, max = 0.5, step = 0.01),
                            sliderInput(inputId = "inNP",
                                        label="Number of grids",
                                        value = 10, min = 5, max = 50, step = 1)
                            ),
                     column(4,
                            h4("Other"),
                            sliderInput(inputId = "inNCore", label="Number of cores",
                                        value = 1, min = 1,
                                        max = (parallel::detectCores()-1), step = 1),
                            sliderInput(inputId = "inNRep", label="Number of replications",
                                        value = 1000, min = 100, max = 10000, step = 100)
                            ),
                     column(4,
                            h4("Simulation Scenario"),
                            radioButtons(inputId = "inH0", label="",
                                         choices = c("Under global null hypothesis" = 1,
                                                     "Under alternative hypothesis" = 2),
                                         selected = 2)
                            )
                 )),
             wellPanel(
                 fluidRow(
                     column(3,
                            h4("Sample Size"),
                            radioButtons(inputId = "inMulti", label="Multiplicity adjustment",
                                         choices = get.consts()$designs)
                            ##div(actionButton("btnSize", "Compute"),
                            ##    style="margin-bottom:10px"),
                            ),
                     column(6, htmlOutput("txtSmpSize"))
                 )),
             wellPanel(
                 h4("Overall Results"),
                 fluidRow(
                     column(3,
                            selectInput(inputId = "inVarName", label="Select variable",
                                        choices = get.consts()$rstnames[-(1:2)])),
                     column(9,
                            plotOutput("pltrst"))
                 )
             ))
})

tab.rst <- reactive({
    tabPanel("Stress Test",
             wellPanel(
                 fluidRow(
                     column(4,
                            h4("Data Generation"),
                            numericInput(inputId = "inDelta1",
                                         label=withMathJax("$$E(Y|A=1,T=1) - E(Y|A=1,T=0)$$"),
                                         value=3),
                            numericInput(inputId = "inDelta2",
                                         label=withMathJax("$$E(Y|A=2,T=1) - E(Y|A=2,T=0)$$"),
                                         value=3),
                            numericInput(inputId = "inSigma",
                                         label="Standard Deviation of Y",
                                         value=5)
                            ),
                     column(4,
                            h4("Data Generation"),
                            sliderInput(inputId = "inPi",
                                        label=withMathJax("$$\\pi_1=P(A=1)$$"),
                                        value = 0.4, min = 0, max = 1, step = 0.05),
                            h4("Constraints"),
                            selectInput(inputId = "inAlpha", label="Two-Sided alpha level",
                                        choices=c(0.05, 0.1, 0.2)),
                            sliderInput(inputId = "inPower", label="Power",
                                        value = 0.8, min = 0, max = 1, step = 0.05)
                            ),
                     column(4,
                            div(msg.box("After settings of the study are specified, click the button below to
                                     conduct the stress-test."),
                                div(actionButton("btnRst", "STRESS-TEST"),
                                    style="margin: 0 auto; width: 30%"),
                                style="padding-top: 200px;")
                            )
                 )),
             uiOutput("uiRst"))
})

tab.main <- reactive({
    panels <- list(type = "pills",
                   id="mainpanel",
                   selected="Stress Test",
                   tab.utility(),
                   ##tab.setting(),
                   tab.rst(),
                   tab.advanced(),
                   tabPanel("Download",
                            wellPanel(h4("Download Simulation Results"),
                                      downloadButton('btnRstDload')),
                            DT::dataTableOutput("tblrst")
                            )
                   );
    ##generate
    do.call(tabsetPanel, panels);
})




##----------------------------------------------------------------------------
##                           FUNCTIONS
##----------------------------------------------------------------------------

get.par <- reactive({
    userLog$lpar;
})

get.N <- reactive({
    userLog$N;
})

get.allpi <- reactive({
    pi1 <- get.par()$pi1;
    dpi <- input$inDelP;
    npi <- input$inNP;

    seq(max(0, pi1 - dpi),
        min(1, pi1 + dpi),
        length.out = npi+1);
})

get.lpar <- function() {
    lpar <- list(delta1 = input$inDelta1,
                 delta2 = input$inDelta2,
                 sigma  = input$inSigma,
                 alpha  = as.numeric(input$inAlpha),
                 beta   = 1 - input$inPower,
                 pi1    = input$inPi
                 );

}

get.sample.sizes <- function(lpar, methods){

    ##Create a Progress object
    progress <- shiny::Progress$new(session, min=0, max=1);
    progress$set(message = "Calculating sample sizes...", value=0);
    ##Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close());

    rst <- NULL;
    for (multi in 1:length(methods)) {
        progress$set(value  = (multi-1)/length(methods),
                     detail = names(methods)[multi]);

        if ("mb" == methods[multi]) {
            alpha.input <- rst[["bon.opt"]]$alpha3[c(3,1:2)]/2;
        } else {
            alpha.input <- NULL;
        }

        rst[[methods[multi]]]  <- stDesign(lpar$delta1,
                                           lpar$delta2,
                                           lpar$sigma,
                                           lpar$pi1,
                                           alpha  = lpar$alpha,
                                           beta   = lpar$beta,
                                           method = methods[multi],
                                           alpha.input = alpha.input
                                           );
    }
    rst
}


##get utility without bcratio
get.simu.rst <- reactive({

    if (is.null(input$btnRst))
        return(NULL);

    if (0 == input$btnRst)
        return(NULL);

    isolate({
        lpar  <- userLog$lpar <- get.lpar();
        sizen <- userLog$N <- get.sample.sizes(lpar,
                                               get.consts()$designs);
        all.pi <- get.allpi();
        if (1 == input$inH0) {
            delta1 <- 0;
            delta2 <- 0;
        } else {
            delta1 <- sizen[[1]]$pars["delta1"];
            delta2 <- sizen[[1]]$pars["delta2"];
        }

        ##Create a Progress object
        progress <- shiny::Progress$new(session, min=0, max=1);
        progress$set(message = "Stress-test in progress...", value=0);
        ##Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close());

        rst <- stSimu(sizen,
                      all.pi,
                      progress,
                      cnames  = get.consts()$rstnames,
                      n.rep   = input$inNRep,
                      delta1  = delta1,
                      delta2  = delta2,
                      n.cores = input$inNCore);
    })
})

