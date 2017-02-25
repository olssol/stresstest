
shinyServer(function(input, output, session) {

    source("stress_ui.R", local=TRUE);

    userLog          <- reactiveValues();
    userLog$uid      <- 1;
    userLog$model    <- -1;


    ##------------------------------------
    ##---------main page------------------
    ##------------------------------------
    output$mainpage <- renderUI({
        tab.main();
    })

    ##------------------------------------
    ##---------simulation results---------
    ##------------------------------------

    output$txtSmpSize <- renderText({

        multi <- as.numeric(input$inMulti);
        lpar  <- get.par();
        sn    <- get.N()[[multi]];
        if (1 == multi) {
            rst <- sprintf("<h5> <p>To achieve %5.0f%% power at alpha level %5.3f,
                  with %s multiplicity control, a total of %5.0f subjects per arm
                  are needed for subgroup 1, a total of %5.0f subjects per arm
                  are needed for subgroup 2, and a total of %5.0f subjects per arm
                  are needed for the overall population. </p><p> With the subgroup
                  1 proportion %5.2f, a total of <b>%5.0f</b> subjects per arm are needed to
                  achieve %5.0f%% power for all the hypothesis. </p> </h5>",
                  100*(1-sn$beta), sn$alpha, sn$method,
                  sn$Nall[1], sn$Nall[2], sn$Nall[3], sn$pars['pi1'],
                  sn$N, 100*(1-sn$beta));
        } else if (2 == multi) {
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
                   sn$method,
                   sn$Nall[1], 100*(1-sn$beta), sn$alpha3[1],
                   sn$Nall[2], 100*(1-sn$beta), sn$alpha3[2],
                   sn$Nall[3], 100*(1-sn$beta), sn$alpha3[3],
                   sn$pars['pi1'], sn$N, 100*(1-sn$beta), sn$alpha);
        }
        rst
    })

    output$tblrst<- DT::renderDataTable({
                            get.simu.rst();
                        },
                        rownames=NULL,
                        options=list(pageLength=100,
                                     scrollX=TRUE))

    output$pltrst <- renderPlot({
        cur.rst <- get.simu.rst();
        if (is.null(cur.rst))
            return(NULL);

        cur.var <- input$inVarName;
        plot.rst(cur.rst, cur.var);
    }, bg="transparent")


    ##------------------------------------
    ##---------download-------------------
    ##------------------------------------

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
