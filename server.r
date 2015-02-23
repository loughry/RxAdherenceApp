library(shiny)
library(MASS)
library(ISLR)

##----------------------------------------------------------------
source("adherenceFunctions.r")


##--------------- Globals -------------------------------
## path to default data that loads at startup  !?!?
oncoData <- "/Users/loughry/RxREVU/Data/Oncology/SAMPLECXP_fileR_PHARMACY_CLAIMS.csv"
diabData <- "/Users/loughry/RxREVU/Data/Diabetes/SAMPLE10K_PHARMACY_CLAIMS.csv"
pharmClaimCols <- c("PT_ID","WRITTEN_DAY","FILLED_DAY","MED_SNM","QTY","NEW_REFILL","DAYS_SUPPLY")
naStrings <- c("NA","9X9X") ## 9X9X used for no value available -- what to do with SQL queries?
MAX.ROWS <- 1000 # max number of rows to display in data table

## ToDo
## Data
## Data paths for onco, diab, later setup SQL queries w UI mod's
## to select data for individuals cohort
## need to put in UI number in cohort to set sample size limits

## report statistic values for individual(s)
## report ... for aggregate, select number in sample
## report for repeated samplings

## What exactly do we want from this?  Evaluation tool for
## reliability of different adherence measures, should be stable
## over repeated sampling within cohort, all, etc.  Yes?

## Ultimately what we want is POC medication recommendations that take into
## account:
## individuals history with med's (adherence, efficacy, side effects
## condition being treated
## what measures of adherence are most relevant to efficacy for individuals
## cohort and relation of individuals adherence behavior so recommendation
## targets best measure of adherence that the individual can support effectively

## Need labelled data like from questionnaire:
##   is you condition currenctly being treated effectively?
##   are your side effects acceptable?

## Then can tie med's, adherence, cohorts to efficacy (and side effects)

## Possible adherence related measures (like risk measures):
## time between "write_day" and "fill_day"

shinyServer(function(input, output) {
    ## Full unconfigured data set
    fullData <- reactive({
        inFile <- input$selData
        if (is.null(inFile)){ # should happen
            print("Null input$selData")
        } else {
            ## chose data here
            fDat <- read.csv(oncoData,header=T,na.strings=naStrings) # whole pharm claims
            for(i in 1:dim(fDat)[2]){
                cDat <- fDat[,i]
                ## sel <- cDat == "9X9X"  # this still left a factor, but with no 9X9X.. using
                ## fDat[sel,i] <- NA      # na.strings instead but need a SQL solution
            }
            fDat[,pharmClaimCols]
        }
    })

    ## setSeed <- reactive({
    ##     if(input$randomSeed == 0)
    ##         set.seed(NULL)
    ##     else
    ##         set.seed(input$randomSeed)
    ## })

    fullDataDim <- reactive({
        dim(fullData())
    })

    fullColNms <- reactive({
        names(fullData())
    })

    fullDataName <- reactive({
        paste(input$selData,"pharmacy claims")
    })

    fullNumericColNms <- reactive({
        fullColNms()[as.numeric(input$numericCol)]
    })

    fullTextColNms <- reactive({
        fullColNms()[as.numeric(input$textCol)]
    })

    fullFactorColNms <- reactive({
        fullColNms()[as.numeric(input$factorCol)]
    })

    fullSelectColNms <- reactive({
        fullColNms()[as.numeric(input$selectCol)]
    })

    fullSampleIdx <- reactive({
        if(input$sampleSize >= fullDataDim()[1] || input$sampleSize == 0){
            1:fullDataDim()[1]
        } else {
            sample(fullDataDim()[1],input$sampleSize)
        }
    })

    fullColumnDataTypes <- reactive({
        types <- c(NA,fullDataDim()[2])
        types[as.numeric(input$numericCol)] <- "numeric"
        types[as.numeric(input$factorCol)] <- "factor"
        types[as.numeric(input$textCol)] <- "text"
        types
    })

    ## Reduce fullData to working data set.

    data <- reactive({
        fullData()[fullSampleIdx(),as.numeric(input$selectCol)]
    })

    colNms <- reactive({
        names(data())
    })

    dataName <- reactive({
        if (!input$getData){
            "Iris Data"  # the default data set !?!?
        } else {
            if(input$importData){
                paste("Imported Data from",input$file1)
            } else {
                input$dataSet
            }
        }
    })

    ## How to translate/compute col specificaions from raw to active data set?
    colTypes <- reactive({
        fullColumnDataTypes()[as.numeric(input$selectCol)]
    })

    numericColNms <- reactive({
        colNms()[as.numeric(input$numericCol)]
    })

    textColNms <- reactive({
        colNms()[as.numeric(input$textCol)]
    })

    factorColNms <- reactive({
        colNms()[as.numeric(input$factorCol)]
    })


    ##------------------- Evaluate Metrics ----------------------

    ptIDs <- reactive({
        unique(data()[,"PT_ID"])
    })

    medNms <- reactive({
        unique(as.character(data()[,"MED_SNM"]))
    })

    metricParms <- reactive({
        mets <- list(mpr=NA,lot=NA,gaps=NA,doc=NA,pers=NA)
        for(sm in names(mets)){
            switch(sm,
                   "lot" = mets[[sm]] <- NA,
                   "mpr" = mets[[sm]] <- list(type=input$mprType,termLen=input$mpr.TermLen),
                   "gaps"= mets[[sm]] <- list(evalType=input$gapType, evalCase=input$gapFun,
                                             byYear=input$byYear),
                   "doc"= mets[[sm]] <-  list(termLen=input$docLen),
                   "pers"= mets[[sm]] <- list(threshold=input$persThresh,
                                   percent=input$persPercent)
                   )
        }
        mets
    })

    compMets <- function(pmcD,adhrMets,metParms){
        mets <- list(mpr=NA,lot=NA,gaps=NA,doc=NA,pers=NA)
        ## Use isolate() to avoid dependencies
        for(sm in adhrMets){
            switch(sm,
                   "lot" = mets[[sm]] <- lenOfTherapy(pmcD),
                   "mpr" = mets[[sm]] <- medPosRat(pmcD,type=metParms$mpr$type,
                       termLen=metParms$mpr$termLen),
                   "gaps"= mets[[sm]] <- medGaps(pmcD, evalType=metParms$gaps$evalType,
                       evalCase=metParms$gaps$evalCase, byYear=metParms$gaps$byYear),
                   "doc"= mets[[sm]] <- daysOfCover(pmcD,termLen=metParms$doc$termLen),
                   "pers"= mets[[sm]] <- persistence(pmcD, threshold=metParms$pers$threshold,
                       percent=metParms$pers$percent),
                   mets[["ERROR"]] <- paste("metric not found: ",sm)
                   )
        }
        mets
    }

    getPmcDat <- function(dat,pt,med){
        ## dat <- data()
        ## ## ?!?! Need to select pt and med
        ## pt <- ptID
        ## med <- medNm
        sel <- dat[,"PT_ID"] == pt & dat[,"MED_SNM"] == med
        pDat <- dat[sel,c("QTY","DAYS_SUPPLY","FILLED_DAY")]
        pDat <- pDat[order(pDat[,"FILLED_DAY"],decreasing=F),]
        daysElapsed <- c(NA,diff(pDat[,"FILLED_DAY"]))
        pDat <- cbind(pDat,daysElapsed)
        names(pDat) <- c("qty","days","fillDay","diffDays")
        pDat
    }

    pmcDat <- reactive({
        dat <- data()
        ## ?!?! Need to select pt and med
        pt <- ptIDs()[1]
        med <- medNms()[1]
        getPmcDat(dat,pt,med)
    })

    aMets <- reactive({
        ## Take a dependency on input$goButton
        if(input$evaluate){ # do metric computations
            ## Use isolate() to avoid dependencies
            pmcD <- isolate(pmcDat())
            mtrcParms <- isolate(metricParms())
            selMets <- isolate(input$adhrMtrcs)
            compMets(pmcD,selMets,mtrcParms)
        }
    })

    ## For each sample (partial time series of pharm.claims data
    ## compute the specified metrics.  Returns matrix of samples
    ## (partial series) by metrics
    compMetSeries <- function(pmc,stepSz=356){
        pmc1 <- isolate(pmcDat()) # get patient data
        mets <- list(mpr=NA,lot=NA,gaps=NA,doc=NA,pers=NA) ## ?!?! c/b sparse
        nMets <- length(mets)
        ## compute interval masks
        fDays <- pmc1[,"fillDay"]
        nDays <- last(fDays) - first(fDays)
        nMsk <- floor(nDays/stepSz) + 1 # of samples, cieling includes partials
        msks <- list()
        for(m in nMsk:1){  ## get sample masks
            msks[[nMsk-m+1]] <- fDays <= last(fDays)-stepSz*(m-1)
        }
        ## Setup matrix and compute metrics for each sample
        metM <- array(NA,c(nMsk,nMets))  # matrix of <samples,metrics>
        dNames <- list(1:nMsk,names(mets))
        names(dNames) <- c("Sample","Metrics")
        dimnames(metM) <- dNames
        for(s in 1:nMsk){ # for each sample (mask) compute metrics)
            pmc1s <- pmc1[msks[[s]],]
            if(dim(pmc1s)[1] < 2 ) next  ## it it is 'empty' skip
            mtrcParms <- isolate(metricParms())
            selMets <- isolate(input$adhrMtrcs)
            mets <- compMets(pmc1s,selMets,mtrcParms)
#            metM[s,] <- as.numeric(compMets(pmc1s,selMets,mtrcParms))
            metM[s,] <- as.numeric(mets)

        }
        metM
    }

    aMetSeries <- reactive({
        ## Take a dependency on input$goButton
        if(input$evaluate){ # do metric computations
            ## Use isolate() to avoid dependencies
            pmcD <- isolate(pmcDat())
            mtrcParms <- isolate(metricParms())
            selMets <- isolate(input$adhrMtrcs)
            compMetSeries(pmcD)
        }
    })

    ## Now I need to aggregate metrics and stats and implement
    ## moving window.

    ## Compute stats for all patients in cohort
    ## cMets <- reactive({
    ##     mets <- list()
    ##     if(input$evaluate){ # do metric computations
    ##         pts <- unique(isolate(ptIDs()))  # specify cohort IDs here ?!?!
    ##         med <- isolate(medNms()[1]) # specify medNm here ?!?!
    ##         dat <- isolate(data())
    ##         ## Use isolate() to avoid dependencies
    ##         for(pt in pts[1]){
    ##             pmcD <- getPmcDat(dat,pt,med)
    ##             mtrcParms <- isolate(metricParms())
    ##             selMets <- isolate(input$adhrMtrcs)
    ##             mets[[pt]] <- compMets(pmcD,selMets,mtrcParms)
    ##         }

    ##     }
    ##     mets
    ## })

    ## I should create header text to be used below where
    ## there is redundant code
    ## Models Configuration
    ## Configure names to be used for analysis including formula's
    ## Should do some error checking here too, e.g. class==factor,
    ## and numeric == numeric, output not in input, ...
    ## Add * ^ : log(), etc.  ?!?!

    ## This is the generic form, not used anywhere I don't think
    modelNms <- reactive({
        inNms <- colNms()[as.numeric(input$inputCol)]
        outNms <- colNms()[as.numeric(input$outputCol)[1]]   #?!?!? forcing one only
        numericNms <- colNms()[as.numeric(input$numericCol)]
        classNms <- colNms()[as.numeric(input$factorCol)]
        if(!is.null(input$outputCol)){  # must have output
            if(!is.null(input$inputCol)){  # may have inputs
                form1 <- as.formula(paste(paste(outNms,"~", paste(inNms, collapse= "+"))))
            } else {
                form1 <- as.formula(paste(outNms,"~ 1"))
            }
        } else {
            form1=NA
        }
        list(In=inNms,Out=outNms,Numeric=numericNms,Class=classNms,Formula1=form1)
    })

    ## Linear Model
    modelNmsLM <- reactive({
        inNms <- colNms()[as.numeric(input$inputLM)]
        outNms <- colNms()[as.numeric(input$outputLM)[1]]   #?!?!? forcing one only
        if(!is.null(input$outputLM)){  # must have output
            if(!is.null(input$inputLM)){  # may have inputs
                form1 <- as.formula(paste(paste(outNms,"~", paste(inNms, collapse= "+"))))
            } else {
                form1 <- as.formula(paste(outNms,"~ 1"))
            }
        } else {
            form1=NA
        }
        list(In=inNms,Out=outNms,Formula1=form1)
    })

    ## Logistic Regression
    modelNmsLogReg <- reactive({
        inNms <- colNms()[as.numeric(input$inputLogReg)]
        outNms <- colNms()[as.numeric(input$outputLogReg)[1]]   #?!?!? forcing one only
        if(!is.null(input$outputLogReg)){  # must have output
            if(!is.null(input$inputLogReg)){  # may have inputs
                form1 <- as.formula(paste(paste(outNms,"~", paste(inNms, collapse= "+"))))
            } else {
                form1 <- as.formula(paste(outNms,"~ 1"))
            }
        } else {
            form1=NA
        }
        list(In=inNms,Out=outNms,Formula1=form1)
    })

    ## LDA
    modelNmsLDA <- reactive({
        inNms <- colNms()[as.numeric(input$inputLDA)]
        outNms <- colNms()[as.numeric(input$outputLDA)[1]]   #?!?!? forcing one only
        if(!is.null(input$outputLDA)){  # must have output
            if(!is.null(input$inputLDA)){  # may have inputs
                form1 <- as.formula(paste(paste(outNms,"~", paste(inNms, collapse= "+"))))
            } else {
                form1 <- as.formula(paste(outNms,"~ 1"))
            }
        } else {
            form1=NA
        }
        list(In=inNms,Out=outNms,Formula1=form1)
    })

    ## Clustering
    modelNmsClust <- reactive({
        inNms <- colNms()[as.numeric(input$inputClust)]
        outNms <- colNms()[as.numeric(input$outputClust)[1]]   #?!?!? forcing one only
        if(!is.null(input$outputClust)){  # must have output
            if(!is.null(input$inputClust)){  # may have inputs
                form1 <- as.formula(paste(paste(outNms,"~", paste(inNms, collapse= "+"))))
            } else {
                form1 <- as.formula(paste(outNms,"~ 1"))
            }
        } else {
            form1=NA
        }
        list(In=inNms,Out=outNms,Formula1=form1)
    })

    ## This should probly be split up into different models, dispatched elsewhere?
    ## Need these separated to get output display right
    ## modelResult <- reactive({
    ##     modelBuilds <- input$buildModel
    ##     if(modelBuilds == 0){  # User has selected to build model
    ##         mod <- NA
    ##     } else {
    ##         ## Use switch to select model type
    ##         ## using isolate to prevent rebuilds triggered by data changes
    ##         mod <- isolate(  ## disconnect the reactive reexecution link
    ##             switch(input$analysis,
    ##                    "linReg" = lm(modelNmsLM()$Formula1,data=data()),
    ##                    "logReg" = glm(modelNmsLogReg()$Formula1,data=data(),family=binomial),
    ##                    "lda" = lda(modelNmsLDA()$Formula1,data=data()),
    ##                    "cluster" = NA)
    ##         ) # close isolate
    ##     }
    ##     list(Model=mod,Builds=modelBuilds)
    ## })

    modResLM <- reactive({
        modelBuilds <- input$buildModel  ## Trigger rebuild
        if(modelBuilds == 0){  # User has selected to build model
            mod <- NA
        } else {
            ## using isolate to prevent rebuilds triggered by data changes
            mod <- isolate(  ## disconnect the reactive reexecution link
                lm(modelNmsLM()$Formula1,data=data())
            ) # close isolate
        }
        list(Model=mod,Builds=modelBuilds)
    })

    modResLogReg <- reactive({
        modelBuilds <- input$buildModel  ## Trigger rebuild
        if(modelBuilds == 0){  # User has selected to build model
            mod <- NA
        } else {
            ## using isolate to prevent rebuilds triggered by data changes
            mod <- isolate(  ## disconnect the reactive reexecution link
                glm(modelNmsLogReg()$Formula1,data=data(),family=binomial)
            ) # close isolate
        }
        list(Model=mod,Builds=modelBuilds)
    })

    modResLDA <- reactive({
        modelBuilds <- input$buildModel  ## Trigger rebuild
        if(modelBuilds == 0){  # User has selected to build model
            mod <- NA
        } else {
            ## using isolate to prevent rebuilds triggered by data changes
            mod <- isolate(  ## disconnect the reactive reexecution link
                lda(modelNmsLDA()$Formula1,data=data())
            ) # close isolate
        }
        list(Model=mod,Builds=modelBuilds)
    })

    modResClust <- reactive({
        modelBuilds <- input$buildModel  ## Trigger rebuild
        if(modelBuilds == 0){  # User has selected to build model
            mod <- NA
        } else {
            ## using isolate to prevent rebuilds triggered by data changes
            mod <- isolate(  ## disconnect the reactive reexecution link
                NA  ## insert clustering here ?!?!
            ) # close isolate
        }
        list(Model=mod,Builds=modelBuilds)
    })

    ## Get numeric plot columns -- should probable decouple this like I have
    ## for the analysis display, especially as it could be complexified
    plotCols <- reactive({
        switch(input$plotType,
                       "plotX" = as.numeric(input$plotXCol),
                       "hist" = as.numeric(input$histCol),
                       "plotXY" = as.numeric(input$plotXYCol),
                       "pairs" = as.numeric(input$pairsCol),
                       "contour" = as.numeric(input$contourCol),
                       "image" = as.numeric(input$imageCol),
                       "persp" = as.numeric(input$perspCol),
                       "clust" = as.numeric(input$clustCol),
                       1  # default column 1, should never happen, make NA ?!?!
                       )
    })

    ##=======================================================================
    ##=============== Begin Output Items ====================================
    ##=======================================================================


    ##----------------------------------------------------------------------
    ##----------------  Data Config Panels  --------------------------------
    ##----------------------------------------------------------------------

    ## Generate a summary of the data
    ## Fix the formatting and content here
    output$configInfo <- renderPrint({
        ## cat("------------------  Data Info ----------------------\n")
        cat("Full Data dimensions: ")
        cat(dim(fullData()))
        cat("\n----------------------------------------\n")
        cat("Selected Data dimensions: ")
        cat(dim(data()))
        cat("\n----------------------------------------\n")
        cat("Column Names:")
        cLen <- length(fullColNms())
        cat(paste("\n",1:cLen,fullColNms()))
        ## list column types
        cat("\n----------------------------------------\n")
        cat("Column Classes Found:\n")
        colCls <- NULL
        for(i in 1:7)    colCls <- c(colCls,class(fullData()[,i]))
        cat(colCls)
        cat("\n")
        cat(fullNumericColNms())
        cat("\n----------------------------------------\n")
        cat("Numeric Columns Specified:\n")
        cat(fullNumericColNms())
        cat("\n----------------------------------------\n")
        cat("Class/Factor Columns Specified:\n")
        cat(fullFactorColNms())
        cat("\n----------------------------------------\n")
        cat("Text Columns: Specfied\n")
        cat(fullTextColNms())
        cat("\n----------------------------------------\n")
        cat("Selected Columns:\n")
        cat(fullSelectColNms())
        cat("\n----------------------------------------\n")
        cat("Sample Index:\n")
        if(length(fullSampleIdx()) >  100){
            cat("  Greater than 100 samples so not displayed")
        } else {
            cat(fullSampleIdx())
        }
        cat("\n----------------------------------------\n")
    })

    output$configSummary <- renderPrint({
        summary(data())
    })

    output$configDebug <- renderPrint({
    })

    output$configDataTable <- renderTable({
        dat <- data()
        if(dim(dat)[1] > MAX.ROWS){
            dat[1:MAX.ROWS,] # crop to max row
        } else {
            dat
        }
    })

    ##-----------------------------------------------------------------
    ##   Adherence Eval Setup  Panels
    ##-----------------------------------------------------------------

    output$adhrInfo <- renderPrint({
        cat("----------------  Working Data Info ------------------\n")
        cat("Working Data dimensions: ")
        cat("\n----------------------------------------\n")
        cat("Selected Adherence Metrics:\n")
        cat(input$adhrMtrcs)
        cat("\n----------------------------------------\n")
        cat("Metric Configurations :\n")
        cat("MPR: ")
        cat(input$mprType,"   ",input$mrp.TermLen)
        cat("\n")
        cat("LoT: n/a")
        cat("\n")
        cat("Persistence: ")
        cat(input$persThresh,"   ",input$persPercent)
        cat("\n")
        cat("DoC: ")
        cat(input$docLen)
        cat("\n")
        cat("Gaps: ")
        cat(input$gapFun, "   ", input$gapType, "   ", input$byYear)
        cat("\n")
        cat("\n----------------------------------------\n")
        cat("parameters in metricParms()\n")
        print(metricParms())
        cat("\n")
        cat("\n----------------------------------------\n")
    })

    output$adhrSummary <- renderPrint({
        cat("Adherence Metrics Computed: \n")
        cat(names(aMets()), "\n")
        paste(aMets())
        ## cat("cMets()", "\n")
        ## paste(cMets())
    })

    output$adhrDebug <- renderPrint({
        cat("Evaluate Counter: ")
        cat(paste("evaluate:", input$evaluate))
        cat("\n----------------------------------------\n")
        cat("input$adhrMtrcs:\n")
        cat(input$adhrMtrcs)
        cat("\n----------------------------------------\n")
        cat("\nPatient IDs: \n")
        cat(ptIDs())
        cat("\n----------------------------------------\n")
        cat("Medication Names: \n")
        mLen <- length(medNms())
        cat(paste(medNms(),","))
        cat("\n----------------------------------------\n")
    })

    output$adhrMetricsTable <- renderTable({
        aMetSeries()
        ## if(dim(dat)[1] > MAX.ROWS){
        ##     dat[1:MAX.ROWS,] # crop to max row
        ## } else {
        ##     dat
        ## }
    })

    output$adhrDataTable <- renderTable({
        dat <- pmcDat()
        if(dim(dat)[1] > MAX.ROWS){
            dat[1:MAX.ROWS,] # crop to max row
        } else {
            dat
        }
    })

    ##---------------------------------------------------------------------
    ##-------------------- Data Eval Panels -----------------------------
    ##---------------------------------------------------------------------

    ## output$dataNames <- renderUI({
    ##     dataNms <<- names(data())
    ##     dataNms
    ## })

    ## Generate an HTML table view of the data
    output$table <- renderTable({
        data.frame(data())
    })

    ## Generate a summary of the data
    output$summary <- renderPrint({
        summary(data())
    })

    ## Fix the formatting and content here
    output$info <- renderPrint({
        print("------------------  Data Info ----------------------")
        print("Column Names:")
        print(colNms())
        print("Column Data Types:")
        print(colTypes())
        print("dim(data)=")
        print(dim(data()))
        print("--------------- Plot Info -------------------")
        print("Plot Type:")
        print(input$plotType)
        print("Ploting Columns Selected")
        print(colNms()[plotCols()])
        print("Contour Function")
        print(input$contourFun)
    })

    output$debug <- renderPrint({
    })

    ## Probly want to break this apart ?!?!
    output$dataPlot <- renderPlot({
        switch(input$plotType,
               ## Should deal with 0 & > 2 columns ?!?!
               "plotX" = plot(data()[,plotCols()[1]],main="Observations Plot",
                   sub=NULL,
                   xlab= "Observation",
                   ylab= colNms()[plotCols()[1]]),
               "hist" = hist(data()[,plotCols()[1]],
                   main="Histogram",
                   sub=NULL,
                   xlab=colNms()[plotCols()[1]],
                   ylab="Frequency" ),
               "plotXY" = plot(data()[,plotCols()[1:2]],main="X-Y Plot",
                   sub=NULL,
                   xlab= colNms()[plotCols()[1]],
                   ylab= colNms()[plotCols()[2]]),
               "pairs" = pairs(data()[,plotCols()],main="Pairs (scatter) Plot"),
               "contour" = contour(
                   outer(data()[,plotCols()[1]],
                         data()[,plotCols()[2]], FUN=input$contourFun),
                   main="Contour Plot"),
               ## "image" = pairs(data()[,plotCol],main="Pairs (scatter) Plot"),
               ## "persp" = pairs(data()[,plotCol],main="Pairs (scatter) Plot"),
               ## "clust" = pairs(data()[,plotCol],main="Pairs (scatter) Plot")
               plot(1:10)  # default plot should display message
               ) # end switch
    })

    ##===========================================
    #           Analysis Panels
    #============================================

    ##--------------- LinReg -------------------
    output$infoLinReg <- renderPrint({
        print("Linear Regression Details")
        print(input$analysis)
        print(colNms())
        print("Inputs:")
        ## print(colNms()[input$inputCol])
        print(modelNmsLM()$In)
        print("Outputs:")
        print(modelNmsLM()$Out)
        print("Formula:")
        print(modelNmsLM()$Formula1)
        print("Model Summary:")
        summary(modResLM()$Model)
    })

    output$plotLinReg <- renderPlot({
        if(modResLM()$Builds == 0){
            plot(1:10)  ## ?!?! Replace this with message
        } else {
            plot(modResLM()$Model,which=input$plotNumber)
        }
    })

    output$debugLinReg <- renderPrint({
        print("LinReg Debug info")
        print(input$analysis)
        print("Model Builds")
        print(as.numeric(modResLM()$Builds))
    })

    ##--------------- LogReg -------------------
    output$infoLogReg <- renderPrint({
        print("Logistic Regression Details")
        print(input$analysis)
        print(colNms())
        print("Inputs:")
        ## print(colNmsLogReg()[input$inputCol])
        print(modelNmsLogReg()$In)
        print("Outputs:")
        print(modelNmsLogReg()$Out)
        print("Formula:")
        print(modelNmsLogReg()$Formula1)
        print("Model Summary:")
        summary(modResLogReg()$Model)
    })

    output$plotLogReg <- renderPlot({
        if(modResLogReg()$Builds == 0){
            plot(1:10)  ## ?!?! Replace this with message
        } else {
            plot(modResLogReg()$Model)  # ?!?! which=...
        }
    })

    output$debugLogReg <- renderPrint({
        print("LogReg Debug Info")
        print(input$analysis)
        print("Model Builds:")
        print(as.numeric(modResLogReg()$Builds))
    })

    ##--------------- LDA -------------------
    output$infoLDA <- renderPrint({
        print("Linear Discriminant Analysis Details")
        print(input$analysis)
        print(colNms())
        print("Inputs:")
        ## print(colNms()[input$inputCol])
        print(modelNmsLDA()$In)
        print("Outputs:")
        print(modelNmsLDA()$Out)
        print("Formula:")
        print(modelNmsLDA()$Formula1)
        print("Model Summary:")
        summary(modResLDA()$Model)
    })

    output$plotLDA <- renderPlot({
        if(modResLDA()$Builds == 0){
            plot(1:10)  ## ?!?! Replace this with message
        } else {
            plot(modResLDA()$Model)
        }
    })

    output$debugLDA <- renderPrint({
        print("LDA Debug Info")
        print(input$analysis)
        print("Model Builds:")
        print(as.numeric(modResLDA()$Builds))
    })

    ##--------------- Cluster -------------------
    output$infoCluster <- renderPrint({
        print("Clustering Details")
        print(input$analysis)
        print(colNms())
        print("Inputs:")
        ## print(colNms()[input$inputCol])
        print(modelNmsClust()$In)
        print("Outputs:")
        print(modelNmsClust()$Out)
        print("Formula:")
        print(modelNmsClust()$Formula1)
        print("Model Summary:")
        summary(modResClust()$Model)
    })

    output$plotCluster <- renderPlot({
        if(modResClust()$Builds == 0){
            plot(1:10)  ## ?!?! Replace this with message
        } else {
            plot(modResClust()$Model)
        }
    })

    output$debugCluster <- renderPrint({
        print("Clustering Debug Info")
        print(input$analysis)
        print("Model Builds:")
        print(modResClust()$Builds)
    })


})

