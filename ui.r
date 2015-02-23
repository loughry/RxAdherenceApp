library(shiny)
library(MASS)
library(ISLR)

colRng <- 1:100  # hard coded for now, needs to be based on data size
##dataNms <- c("one","two")

shinyUI(fluidPage(
    titlePanel("RxREVU Adherence Analysis 1"),
    sidebarLayout(
        sidebarPanel(
            ##----------------------- Select Data(base) to use -------------------
            selectInput("selData", "Select Data Source",
                        c("Oncology100" = "onco100",
                          "Diabetes10k" = "diab10k"),
                        selected="onco100"),
            ##------------------ Config data -------------------------------------
            checkboxInput("dataConfig",
                          label = "Specify Data Types and Columns", value = TRUE),
            conditionalPanel(  ## Data speficification panel
                condition = "input.dataConfig == true",
                helpText("Specify the data types found in each column, etc."),
                ##                            dataNms <<- uiOutput("dataNames"),
                selectInput('numericCol', 'Numeric Columns', colRng,
                            multiple=TRUE, selectize=TRUE),
                selectInput('factorCol', 'Class Columns', colRng,
                            multiple=TRUE, selectize=TRUE),
                selectInput('textCol', 'Text Columns', colRng,
                            multiple=TRUE, selectize=TRUE),
                selectInput('selectCol', 'Select Columns', colRng,
                            multiple=TRUE, selectize=TRUE,selected=c(1:5,7)),
                ## numericInput('randomSeed',"Set Random Seed (0 for no seed)",1),
                numericInput("sampleSize",
                             "Sample Size", 0)
            ),
            ##--------------- Configure Adherence Metrics -------------------------
            checkboxInput("configMetrics",
                          label = "Configure Adherence Metrics", value = FALSE),
            conditionalPanel(  ## Data speficification panel
                condition = "input.configMetrics == true",
                selectInput("metric", "Adherence Metric",
                            c("Medication Possision Ratio" = "mpr",
                              "Length Of Therapy" = "lot",
                              "Persistence" = "pers",
                              "Medication Gaps" = "gap",
                              "Days of Coverage" = "doc"),
                            selected="mpr"),
                conditionalPanel(  ## MPR Config Panel
                    condition = "input.metric == 'mpr'",
                    helpText("Configure MPR"),
                    selectInput("mprType", "Select MPR Type",
                                c("Length Of Therapy" = "lot",
                                  "Recent Time Period" = "last"),
                                selected="lot"),
                    conditionalPanel(
                        condition = "input.mprType == 'last'",
                        sliderInput("mpr.TermLen", "Period Length in Years",
                                    min=.1, max=5.0, value=1.0) )
                ), # Close MPR config panel
                conditionalPanel(  ## LOT Config Panel
                    condition = "input.metric == 'lot'",
                    helpText("Nothing to configure for LoT")
                ), # Close LOT config panel
                conditionalPanel(  ## Persistence Config Panel
                    condition = "input.metric == 'pers'",
                    helpText("Configure Persistence"),
                    numericInput("persThresh", "Persistence Threshold in Days", 100),
                    checkboxInput("persPercent", label = "Percent",
                                  value = FALSE)
                ), # Close Persistence config panel
                conditionalPanel(  ## Gap Config Panel
                    condition = "input.metric == 'gap'",
                    helpText("Configure Gap"),
                    selectInput("gapFun", "Reporting Statistic",
                                c("Mean" = "mean",
                                  "Max" = "max",
                                  "Sum"="sum",
                                  "Median"="median",
                                  "Min"="min"),
                                selected="mean"),
                    radioButtons("gapType", "Gap Measure",
                                 c("Count" = "count",
                                   "Values" = "value"),
                                 selected="value"),
                    checkboxInput("byYear", label = "Stats By Year",
                                  value = FALSE)
                ), # Close GAP config panel
                conditionalPanel(  ## Days Of Coverage Config Panel
                    condition = "input.metric == 'doc'",
                    helpText("Configure DoC, period ends at last data point"),
                    numericInput("docLen", "Period Length in Years", 1.0,min=0.1,
                                 max=5.0,step=0.1)

                ) # Close DoC config panel
            ), # end config metrics panels

            ## ----------- Setup Adherence Eval ---------------------
            checkboxInput("setAdhr",
                          label = "Setup Adherence Metrics Evaluation", value = FALSE),
            conditionalPanel(  # If evaluating adherence metrics
                condition = "input.setAdhr == true",
                selectInput("adhrMtrcs", "Adherence Metrics",  ## Select metrics
                            c("Medication Possesion Ratio" = "mpr",
                              "Length of Therapy" = "lot",
                              "Persistence" = "pers",
                              "Medication Gaps" = "gaps",
                              "Days of Coverage" = "doc"),
                            multiple=TRUE
                            ),
                ## One or more ID's, one or more metrics, range of metric parameters,
                ## repeated
                ## sampling ID's, moving endpoint time series, Cohort samples,
                ## Single ID:
                ##   Single Metric
                ##   Multiple Metric
                ## Multiple ID:
                ##   group selection, sampling from
                ##   stats over ID collection or over resampling
                numericInput("smpStepSz", "Sampling Step Size (0 none)",0,min=0,
                                 max=720,step=30),
                helpText("Coming: cohort selection, ways to view and analyze"),
                actionButton("evaluate", "Evaluate")
            ), # Close adherence metric's eval panel
            hr(),

            ##----------------------- Explore Data -----------------------------
            checkboxInput("showDataPlot",
                          label = "Explore Data (show data & plot)", value = FALSE),
            conditionalPanel(  # If showing and plotting data, select columns and plot type
                condition = "input.showDataPlot == true",
                selectInput("plotType", "Plot Type",  ## Select analysis/model type
                            c("Observations" = "plotX",
                              "Histogram" = "hist",
                              "XY Plot" = "plotXY",
                              "Pairs/Scatter" = "pairs",
                              "Contour" = "contour",
                              "Image (3D)" = "image",
                              "Perspective" = "persp",
                              "Cluster" = "clust"),
                            selected="plotX"),
                conditionalPanel(  ## PlotX Config Panel
                    condition = "input.plotType == 'plotX'",
                    helpText("Configure Observations Plot"),
                    selectInput('plotXCol', 'Plot Column', colRng,selected=1,
                                multiple=FALSE)
                    ), # Close PlotX panel
                conditionalPanel(  ## Hist Config Panel
                    condition = "input.plotType == 'hist'",
                    helpText("Configure Histogram Plot"),
                    selectInput('histCol', 'Plot Column', colRng,selected=1,
                                multiple=FALSE)
                    ), # Close Hist panel
                conditionalPanel(  ## PlotXY Config Panel
                    condition = "input.plotType == 'plotXY'",
                    helpText("Configure XY Plot"),
                    selectInput('plotXYCol', 'Plot Column', colRng,selected=1:2,
                                multiple=TRUE,selectize=TRUE)
                    ), # Close PlotXY panel
                conditionalPanel(  ## Pairs Config Panel
                    condition = "input.plotType == 'pairs'",
                    helpText("Configure Pairs/Scatter Plot"),
                    selectInput('pairsCol', 'Plot Columns', colRng,selected=1:2,
                                multiple=TRUE,selectize=TRUE)
                    ), # Close Pairs panel
                conditionalPanel(  ## Contour Config Panel
                    condition = "input.plotType == 'contour'",
                    helpText("Configure Contour Plot"),
                    selectInput('contourCol', 'Plot Columns X & Y',
                                colRng,selected=1:3,
                                multiple=TRUE,selectize=TRUE),
                    selectInput("contourFun", "Select Function (X,Y)",
                                c("Multiply" = "*",
                                  "Divide" = "/",
                                  "Add" = "+",
                                  "Subtract" = "-"),
                                  selected="*")
                ), # Close Contour panel
                conditionalPanel(  ## Image Config Panel
                    condition = "input.plotType == 'image'",
                    helpText("Configure Image Plot"),
                    selectInput('imageCol', 'Plot Columns X,Y & Z',
                                colRng,selected=1:3,
                                multiple=TRUE,selectize=TRUE)
                    ), # Close Image panel
                conditionalPanel(  ## Perspective Config Panel
                    condition = "input.plotType == 'persp'",
                    helpText("Configure Perspective Plot"),
                    selectInput('perspCol', 'Plot Columns X,Y & Z',
                                colRng,selected=1:3,
                                multiple=TRUE,selectize=TRUE)
                    ), # Close Perspective panel
                conditionalPanel(  ## Cluster Config Panel
                    condition = "input.plotType == 'clust'",
                    helpText("Configure Cluster Plot"),
                    selectInput('clustCol', 'Plot Column',
                                colRng,selected=1:3,
                                multiple=TRUE,selectize=TRUE)
                    ), # Close Cluster panel
                ## selectInput('plotCol', 'Plot Columns', colRng,selected=1:2,
                ##             multiple=TRUE,selectize=TRUE),
                hr()
            ),
            ##------------------- Analysis ------------------------
            checkboxInput("analysisSetup",
                          label = "Analysis", value = FALSE),
            conditionalPanel(  ## Analysis Panels
                condition = "input.analysisSetup == true",
                selectInput("analysis", "Analysis",  ## Select analysis/model type
                            c("Linear Regression" = "linReg",
                              "Logistic Regression" = "logReg",
                              "Linear Discriminant Analysis" = "lda",
                              "Cluster Analysis" = "cluster"),
                            selected="linReg"),
                conditionalPanel(  ## Linear Regression Config Panel
                    condition = "input.analysis == 'linReg'",
                    helpText("Configure Linear Regression"),
                    selectInput('inputLM', 'Input Columns',colRng,
                                multiple=TRUE, selectize=TRUE),
                    selectInput('outputLM', 'Output Column', colRng,
                                multiple=TRUE, selectize=TRUE),
                    numericInput('plotNumber','Plot Number',1,min=1,max=6,step=1)
                    ), # Close Linear Regression Panel
                conditionalPanel(  ## LDA Config Panel
                    condition = "input.analysis == 'lda'",
                    helpText("Configure LDA"),
                    selectInput('inputLDA', 'Input Column',colRng,
                                multiple=TRUE, selectize=TRUE),
                    selectInput('outputLDA', 'Output Columns', colRng,
                                multiple=TRUE, selectize=TRUE)
                ), # Close LDA panel
                conditionalPanel(  ## Logistic Regression Config Panel
                    condition = "input.analysis == 'logReg'",
                    helpText("Configure Logistic"),
                    selectInput('inputLogReg', 'Input Columns',colRng,
                                multiple=TRUE, selectize=TRUE),
                    selectInput('outputLogReg', 'Output Column', colRng,
                                multiple=TRUE, selectize=TRUE),
                    numericInput('plotNumber','Plot Number',1,min=1,max=6,step=1)
                ), # Close Logistic panel
                conditionalPanel(  ## Cluster Config Panel
                    condition = "input.analysis == 'cluster'",
                    helpText("Configure Clustering -- NA"),
                    selectInput('inputClust', 'Columns',colRng,
                                multiple=TRUE, selectize=TRUE)
                ), # Close Cluster panel
                actionButton('buildModel','Build Model')
            ) # Close Analysis Panels
        ),
        mainPanel(
            conditionalPanel(  # Data Config Panels
                condition = "input.dataConfig == true",
                tabsetPanel(type = "tabs",
                            tabPanel("Data Configuration",
                                     verbatimTextOutput("configInfo")),
                            ## tabPanel("Data Specs",
                            ##          verbatimTextOutput("configSummary")),
                            tabPanel("Data Debug",
                                     verbatimTextOutput("configDebug")),
                            tabPanel("Data Table",
                                     tableOutput("configDataTable"))
                            )
            ),
            conditionalPanel(  # If setting up evaluating adherence metrics, show this paneld
                condition = "input.setAdhr == true",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("adhrInfo")),
                            tabPanel("Summary & Stats",
                                     verbatimTextOutput("adhrSummary")),
                            tabPanel("Metric Series Table",
                                     tableOutput("adhrMetricsTable")),
                            tabPanel("Data Table",
                                     tableOutput("adhrDataTable")),
                            tabPanel("Debug",
                                     verbatimTextOutput("adhrDebug"))
                            )
            ), # Close Adherence Setup Panel
            conditionalPanel(  # Analysis Panels
                condition = "input.showDataPlot == true",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("info")),
                            tabPanel("Summary",
                                     verbatimTextOutput("summary")),
                            tabPanel("Plot",plotOutput('dataPlot')),
                            tabPanel("Debug",
                                     verbatimTextOutput("debug"))
                            )
            ),
            ## conditionalPanel(  # Analysis Panels
            ##     condition = "input.showDataPlot == false",
##                condition = "input.analysisSetup == true",
            conditionalPanel(  # Lin Reg Output
                condition = "input.analysisSetup == true && input.analysis == 'linReg'",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("infoLinReg")),
                            tabPanel("Plot",plotOutput('plotLinReg')),
                            tabPanel("Debug",
                                     verbatimTextOutput("debugLinReg"))
                            )
            ),
            conditionalPanel(  # Log Reg Output
                condition = "input.analysisSetup == true && input.analysis == 'logReg'",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("infoLogReg")),
                            ## Using plotLinReg since it is glm()
                            tabPanel("Plot",plotOutput('plotLogReg')),
                            tabPanel("Debug",
                                     verbatimTextOutput("debugLogReg"))
                            )
            ),
            conditionalPanel(  # LDA Output
                condition = "input.analysisSetup == true && input.analysis == 'lda'",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("infoLDA")),
                            tabPanel("Plot",plotOutput('plotLDA')),
                            tabPanel("Debug",
                                     verbatimTextOutput("debugLDA"))
                            )
            ),
            conditionalPanel(  # Cluster Output
                condition = "input.analysisSetup == true && input.analysis == 'cluster'",
                tabsetPanel(type = "tabs",
                            tabPanel("Info",
                                     verbatimTextOutput("infoCluster")),
                            ## tabPanel("Plot",plotOutput('plotLinReg'))
                            ## )
                            tabPanel("Debug",
                                     verbatimTextOutput("debugCluster"))
                            )
            )
##        )
        )  # mainPanel
    ) # sideBarLayout
) )  # fluidPage, shinyUI

