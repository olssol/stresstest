
shinyServer(function(input, output, session) {

    source("stress_ui.R", local=TRUE);

    userLog          <- reactiveValues();
    userLog$uid      <- 1;
    userLog$model    <- -1;
    userLog$lpars    <- NULL;
    userLog$N        <- NULL;

    ##------------------------------------
    ##---------main page------------------
    ##------------------------------------
    output$mainpage <- renderUI({
        tab.main();
    })

    ##--------------------------------------
    ##---------exit-------------------------
    ##--------------------------------------
    observeEvent(input$close, {stopApp()});

    ##------------------------------------
    ##---------simulation results---------
    ##------------------------------------

    observeEvent(input$btnSize, {

        lpar <- list(delta1 = input$inDelta1,
                     delta2 = input$inDelta2,
                     sigma  = input$inSigma,
                     alpha  = as.numeric(input$inAlpha),
                     beta   = 1 - input$inPower,
                     pi1    = input$inPi
                     );

        ##Create a Progress object
        progress <- shiny::Progress$new(session, min=0, max=1);
        progress$set(message = "Computation in progress...", value=0);
        ##Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close());

        methods <- get.consts()$designs;
        rst     <- NULL;
        for (multi in 1:length(methods)) {
            progress$set(value  = (multi-1)/length(methods),
                         detail = names(get.consts()$designs)[multi]);

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
        userLog$lpar <- lpar;
        userLog$N    <- rst;
    }, ignoreInit=TRUE);

    ## design details
    output$txtSmpSize <- renderText({

        multi <- input$inMulti;
        drst  <- get.N();

        if (is.null(multi) | is.null(drst))
            return(NULL);

        sn <- get.N()[[multi]];
        sp <- sn$pars;

        if ("bon" == multi) {
            rst <- sprintf("<h5> <p>To achieve %5.0f power at alpha level %5.3f,
                  with Bonferroni multiplicity control, a total of %5.0f subjects per arm
                  are needed for subgroup 1, a total of %5.0f subjects per arm
                  are needed for subgroup 2, and a total of %5.0f subjects per arm
                  are needed for the overall population. </p><p> With the subgroup
                  1 proportion %5.2f, a total of <b>%5.0f</b> subjects per arm are needed to
                  achieve %5.0f%% power for all the hypothesis. </p> </h5>",
                  100*(1-sp['beta']), sp['alpha'],
                  sn$Nall[1], sn$Nall[2], sn$Nall[3], sp['pi1'],
                  sn$N, 100*(1-sp['beta']));
        } else if ("bon.opt" == multi | "mb" == multi) {
            rst <- sprintf("<h5> <p>With %s multiplicity control,
                   a total of %5.0f subjects per arm
                   are needed for subgroup 1 to achieve %5.0f%% power at alpha level %5.3f,
                   a total of %5.0f subjecs per arm
                   are needed for subgroup 2 to achieve %5.0f%% power at alpha level %5.3f,
                   and a total of %5.0f subjects per arm
                   are needed for the overall population to achieve %5.0f%% power at alpha level %5.3g. </p>
                   <p> With the subgroup 1 proportion %5.2f, a total of <b>%5.0f</b> subjects per arm are needed to
                   achieve %5.0f%% power for all the hypothesis under the specific type I error level
                   and control the family-wise
                   type I error at %5.3f level </p> </h5>",
                   names(get.consts()$designs)[which(multi == get.consts()$designs)],
                   sn$Nall[1], 100*(1-sp['beta']), sn$alpha3[1],
                   sn$Nall[2], 100*(1-sp['beta']), sn$alpha3[2],
                   sn$Nall[3], 100*(1-sp['beta']), sn$alpha3[3],
                   sp['pi1'], sn$N, 100*(1-sp['beta']), sp['alpha']);
        } else {
            rst <- NULL;
        }

        rst
    })


    output$pltrst <- renderPlot({
        cur.rst <- get.simu.rst();
        if (is.null(cur.rst))
            return(NULL);

        cur.var <- input$inVarName;
        stPlot(cur.rst, cur.var);
    }, bg="transparent")


    output$pltpower <- renderPlot({
        cur.rst <- get.simu.rst();
        if (is.null(cur.rst))
            return(NULL);
        if (1 == input$inH0)
            return(NULL);

        stPlotStress(cur.rst, power=1-userLog$lpar$beta, cex.main=1.3, cex.lab = 1.3, mar=c(4,4,2,1));
    }, bg="transparent")

    ##------------------------------------
    ##---------download-------------------
    ##------------------------------------

    ##results
    output$tblrst<- DT::renderDataTable({
                            get.simu.rst();
                        },
                        rownames=NULL,
                        options=list(pageLength=100,
                                     scrollX=TRUE))

    ##download result data
    output$btnRstDload <- downloadHandler(
        filename=function() {
        paste('rst_',
              format(Sys.time(), "%m%d%Y%H%M%S"),
              '.txt',sep="")
    },

    content=function(file) {
        tfile <- tempfile();
        write.table(get.simu.rst(), tfile, row.names=FALSE);
        bytes <- readBin(tfile, "raw", file.info(tfile)$size);
        writeBin(bytes, file);
    })

})
