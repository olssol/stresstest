

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

##----------------------------------------------------------------------------
##                           PANELS
##----------------------------------------------------------------------------
tab.utility <- reactive({
    tabPanel("Note",
             withMathJax(includeHTML("www/paper_stress.html")))
})

tab.setting <- reactive({
    tabPanel("Setting",
             wellPanel(
                 fluidRow(
                     column(4,
                            h4("Data Generation"),
                            sliderInput(inputId = "inPi",
                                        label=withMathJax("$$\\pi_1=P(A=1)$$"),
                                        value = 0.4, min = 0, max = 1, step = 0.05),
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
                            h4("Constraints"),
                            selectInput(inputId = "inAlpha", label="Two-Sided alpha level",
                                        choices=c(0.05, 0.1, 0.2)),
                            sliderInput(inputId = "inPower", label="Power",
                                        value = 0.8, min = 0, max = 1, step = 0.05),
                            radioButtons(inputId = "inH0", label="Simulation Scenario",
                                         choices = c("Under global null hypothesis" = 1,
                                                     "Under alternative hypothesis" = 2),
                                         selected = 2)
                            ),
                     column(4,
                            h4("Sample Size"),
                            radioButtons(inputId = "inMulti", label="Multiplicity adjustment",
                                         choices = get.consts()$designs),
                            div(actionButton("btnSize", "Compute"),
                                style="margin-bottom:10px"),
                            htmlOutput("txtSmpSize")
                            )
                 )),
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
                            )
                 ))
             )
})

tab.rst <- reactive({
    tabPanel("Stress Test",
             div(actionButton("btnRst", "Conduct Simulation"),
                 style="margin-bottom:10px"),
             wellPanel(
                 h4("Overall Results"),
                 fluidRow(
                     column(3,
                            selectInput(inputId = "inVarName", label="Select variable",
                                        choices = get.consts()$rstnames[-(1:2)])),
                     column(9,
                            plotOutput("pltrst"))
                 )
             ),
             wellPanel(
                 h4("Power Loss"),
                 plotOutput("pltpower")
             ))
})

tab.main <- reactive({
    panels <- list(type = "pills",
                   id="mainpanel",
                   selected="Setting",
                   tab.utility(),
                   tab.setting(),
                   tab.rst(),
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

##get utility without bcratio
get.simu.rst <- reactive({

    if (is.null(input$btnRst))
        return(NULL);

    if (0 == input$btnRst)
        return(NULL);

    isolate({
        sizen  <- get.N();
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
        progress$set(message = "Computation in progress...", value=0);
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

