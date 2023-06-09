import random
import os


def WGCNA_shiny_app_creator(wgcna_expression_df, wgcna_trait_df):

    output_app_dir_name = str(random.randint(int(1e7), int(1e8)))
    output_app_dir = '/data1/Shiny_Apps/WGCNA/' + output_app_dir_name
    os.makedirs(output_app_dir, exist_ok=True)

    wgcna_expression_df.to_csv(output_app_dir + '/expression_FPKM.csv')
    wgcna_trait_df.to_csv(output_app_dir + '/trait_data.csv')

    wgcna_mainapp_text_raw = """
    
suppressMessages(library(devtools))
suppressMessages(library(ShinyWGCNA))
suppressMessages(library(shinyjs))
suppressMessages(library(dashboardthemes))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(shiny))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library(ape))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(reshape2))
suppressMessages(library(edgeR))
suppressMessages(library(shinythemes))
suppressMessages(library(ggplotify))
suppressMessages(library(ggprism))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(shinyjqui))
suppressMessages(library(ggpubr))
suppressMessages(library(conflicted))
options(shiny.maxRequestSize = 300*1024^2)
options(scipen = 6)
conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("rename","dplyr")
conflict_prefer("desc","dplyr")
conflict_prefer("cor","stats")

# Logo ----
customLogo <- shinyDashboardLogoDIY(

  boldText = "HEMU WGCNA Interactive Analysis Module"
  ,mainText = ""
  ,textSize = 14
  ,badgeText = "Powered by R::Shiny"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "gray"
  ,badgeBorderRadius = 3
)

# UI Frontend ----
ui <- shinyUI(
  navbarPage(customLogo,
             
             # Data_import_panel ----
             tabPanel(
               useShinyjs(),
               title = "I. Data Filtering",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "Sidebar",
                     sidebarPanel(
                       width = 4,
                       selectInput(
                         inputId = "method1",
                         label = "Normalized method",
                         choices = c(FPKM = "rawFPKM",
                                     lgFPKM = "lgFPKM"),
                         selected = "rawFPKM"
                       ),
                       hr(),
                       HTML('<font color = #FF6347  size = 3.2><b>Gene Expression Filter</b></font>'),
                       textInput(
                         inputId = "SamPer",
                         label = "Sample percentage",
                         value = "0.9"
                       ),
                       textInput(
                         inputId = "RCcut",
                         label = "Expression Cutoff",
                         value = "1"
                       ),
                       hr(),
                       # p("Noise removal, for instance, removing all genes that FPKM<1 in 90%+ provided samples",
                       #   style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       # br(),
                       # HTML('<font color = #FF6347 size = 3.2><b>Second Time filter</b></font>'),
                       # radioButtons(
                       #   inputId = "CutMethod",
                       #   label = "Filter Method",
                       #   choices = c("MAD","Var"),
                       #   selected = "MAD"
                       # ),
                       textInput(
                         inputId = "remain",
                         label = "Reserved genes Num.",
                         value = "4000"
                       ),
                     )# sidebarPanel
                 ),# div
                 
                 # Main panel ----
                 mainPanel(
                   fluidPage(
                     tabsetPanel(
                       tabPanel(title = "Input file check",height = "500px",width = "100%",
                                br(),
                                actionButton("action1", "Perform Filtering & Sample Clustering"),
                                htmlOutput("Inputcheck"),
                                htmlOutput("filter1"),
                                htmlOutput("filter2"),
                                
                       ),
                       tabPanel(title = "Preview of Input",height = "500px",width = "100%",
                                DT::dataTableOutput("Inputbl"),
                       ),
                       tabPanel(title = "SampleCluster",height = "500px",width = "100%",
                                jqui_resizable(
                                  plotOutput("clustPlot")
                                ),
                                downloadButton("downfig1","Download")
                       ),# tabPanel
                     ),
                   )# fluidPage
                 )#mainPanel
               ) # sidebarLayout
             ),##tabPanel
             
             
             
             # SFT_power_panel ----
             tabPanel(
               useShinyjs(),
               title = "II. Selection of SFT & Power",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "Sidebar2",
                     sidebarPanel(
                       width = 4,
                       sliderInput(
                         inputId = "CutoffR",
                         label = HTML('R<sup>2</sup> cutoff'),
                         min = 0,
                         max = 1,
                         value = 0.8
                       ),
                       br(),
                       HTML('<font size = 2.5 color = #7a8788><i>WGCNA will generate a recommended power value.<font color = blue>If you find that the R <sup>2</sup> value corresponding to experience power given by the software lower than your setting Threshold </font>,<font color = purple><b> please select a customized power based on the SFT plot.</b></i></font></font>'),
                       radioButtons(
                         inputId = "PowerTorF",
                         label = "Power type",
                         choices = c("Recommended","Customized"),
                         selected = "Recommended"
                       ),
                       
                       sliderInput(
                         inputId = "PowerSelect",
                         label = "Final Power Selection",
                         min = 1,
                         max = 33,
                         value = 6
                       )
                     )
                 ),
                 # Main panel ----
                 mainPanel(
                   fluidPage(
                     #### output field
                     #actionButton("toggleSidebar2",
                                  #"Toggle sidebar"),
                     tabsetPanel(
                       tabPanel(title = "Select Power",height = "500px",width = "100%",
                                br(),
                                actionButton("Startsft","Start Power Estimation"),
                                h4("Please conduct both power estimation and scale-free estimation before pushing on to step 3."),
                                htmlOutput("powerout"),
                                jqui_resizable(
                                  plotOutput("sftplot")
                                ),
                                textInput(inputId = "width2",
                                          label = "figure width",
                                          value = 10),
                                textInput(inputId = "height2",
                                          label = "figure height",
                                          value = 10),
                                actionButton("adjust2","Set fig size"),
                                downloadButton("downfig2","Download")
                                
                                
                       ),
                       tabPanel(title = "Information of sft table",height = "500px",width = "100%",
                                DT::dataTableOutput("sfttbl")
                       ),
                       tabPanel(title = "Scale-free estimation",height = "500px",width = "100%",
                                br(),
                                actionButton("Startcheck","Start scale-free estimation"),
                                p("This step may take a relatively long time when expression dataframe is large."),
                                jqui_resizable(
                                  plotOutput("sfttest")
                                ),
                                textInput(inputId = "width3",
                                          label = "width",
                                          value = 10),
                                textInput(inputId = "height3",
                                          label = "height",
                                          value = 10),
                                actionButton("adjust3","Set fig size"),
                                downloadButton("downfig3","Download")
                                
                       )## tabPanel
                     )## tabsetPanel
                   )## fluidPage
                 )
               )
             ),##tabPanel
             
             
             # Module-net_panel ----
             tabPanel(
               useShinyjs(),
               title = "III. Network Construction",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "Sidebar3",
                     sidebarPanel(
                       width = 4,
                       sliderInput(
                         inputId = "minMsize",
                         label = "Minimum Module Size",min = 0,max = 200,value = 30
                       ),
                       sliderInput(
                         inputId = "mch",
                         label = "Module cuttree height",
                         min = 0, max = 1, value = 0.25
                       ),
                       textInput(inputId = "blocksize",
                                 label = "Maximum block size",
                                 value = 6000),
                       p("Maximum Block Size (Max: 6000)",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                     )
                 ),
                 # Main panel ----
                 mainPanel(
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "Cluster",height = "500px",width = "100%",
                         br(),
                         actionButton("Startnet","Start k-means Clustering"),
                         jqui_resizable(
                           plotOutput("cluster")
                         ),
                         textInput(inputId = "width4",
                                   label = "figure width",
                                   value = 10),
                         textInput(inputId = "height4",
                                   label = "figure height",
                                   value = 10),
                         actionButton("adjust4","Set fig size"),
                         downloadButton("downfig4","Download"),
                         
                         br(),
                         br(),
                         tableOutput("m2num")
                       ),
                       tabPanel(
                         title = "Eigengene adjacency heatmap",height = "500px",width = "100%",
                         jqui_resizable(
                           plotOutput("eah")
                         ),
                         textInput(inputId = "width5",
                                   label = "figure width",
                                   value = 10),
                         textInput(inputId = "height5",
                                   label = "figure height",
                                   value = 10),
                         actionButton("adjust5","Set fig size"),
                         downloadButton("downfig5","Download")
                       ),
                       tabPanel(
                         title = "Gene to module",height = "500px",width = "100%",
                         DT::dataTableOutput("g2m"),
                         downloadButton("downtbl2","download")
                       )
                     )
                     
                   )
                 )
               )
             ),##tabPanel
             
             # Module-trait_panel ----
             tabPanel(
               useShinyjs(),
               title = "IV. Module-Trait Relation",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "Sidebar4",
                     sidebarPanel(
                       width = 4,
                       colourpicker::colourInput(inputId = "colormin",
                                                 label = "Color: Minimum",
                                                 value = "gray"),
                       colourpicker::colourInput(inputId = "colormid",
                                                 label = "Color: Middle",
                                                 value = "white"),
                       colourpicker::colourInput(inputId = "colormax",
                                                 label = "Color: Maximum",
                                                 value = "red"),
                       textInput(
                         inputId = "xangle",
                         label = "x axis label angle",
                         value = 0
                       ),
                     )
                 ),
                 # Main panel ----
                 mainPanel(
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "Module to trait",height = "500px",width = "100%",
                         br(),
                         actionButton("starttrait","Start Module-Trait Heatmap Generation and Analysis"),
                         jqui_resizable(
                           plotOutput("mtplot")
                         ),
                         textInput(inputId = "width6",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height6",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust6","Set fig size"),
                         downloadButton("downfig6","Download")
                       ),
                       tabPanel(
                         title = "Eigengene-based connectivities,KME",height = "500px",width = "100%",
                         DT::dataTableOutput("KME"),
                         
                         downloadButton("downtbl3","download")
                       ),
                     )
                   )
                 )
               )
             ),##tabPanel
             
             # Interested_module_panel ----
             tabPanel(
               useShinyjs(),
               title = "V. Module Analysis",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "sidebar5",
                     sidebarPanel(
                       width = 4,
                       selectInput(
                         inputId = "strait",
                         label = "Select trait",
                         choices = c("select a trait","trait2","..."),
                         selected = "select a trait",
                         multiple = F
                       ),
                       selectInput(
                         inputId = "smodule",
                         label = "Select module",
                         choices = c("red","black","..."),
                         selected = "red",
                         multiple = F
                       ),
                       actionButton("InterMode","Start Analysis")
                     )
                 ),
                 # Main panel ----
                 mainPanel(
                   # actionButton("toggleSidebar5",
                   #              "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "GS-Connectivity",height = "500px",width = "100%",
                         jqui_resizable(
                           plotOutput("GSCon")
                         ),
                         textInput(inputId = "width7",
                                   label = "figure width",
                                   value = 10),
                         textInput(inputId = "height7",
                                   label = "figure height",
                                   value = 10),
                         actionButton("adjust7","Set fig size"),
                         downloadButton("downfig7","Download")
                       ),
                       tabPanel(
                         title = "Heatmap",height = "500px",width = "100%",
                         jqui_resizable(
                           plotOutput("heatmap")
                         ),
                         textInput(inputId = "width8",
                                   label = "figure width",
                                   value = 10),
                         textInput(inputId = "height8",
                                   label = "figure height",
                                   value = 10),
                         actionButton("adjust8","Set fig size"),
                         downloadButton("downfig8","Download")
                       ),
                     )
                   )
                 )
               ),
             ),##tabPanel
             
             # Hubgene panel ----
             tabPanel(
               useShinyjs(),
               title = "VI. Hub Gene Identification",
               sidebarLayout(
                 # Sidebar ----
                 div(id = "sidebar6",
                     sidebarPanel(
                       width = 4,
                       selectInput(
                         inputId = "hubtrait",
                         label = "Select trait",choices = c("select a trait","..."),multiple = F
                       ),
                       selectInput(
                         inputId = "hubmodule",
                         label = "Select module",choices = c("select a module","..."),multiple = F
                       ),
                       actionButton("starthub","Start Analysis")
                     )
                 ),
                 # Main panel ----
                 mainPanel(
                   # actionButton("toggleSidebar6",
                   #              "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "Identification by kME and Gene Significance (GS)",
                         sliderInput(
                           inputId = "kMEcut",
                           label = "kME cutoff",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         sliderInput(
                           inputId = "GScut",
                           label = "GS cutoff",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         DT::dataTableOutput("kMEhub"),
                         downloadButton("downtbl4","download")
                       ),
                       tabPanel(
                         title = "Cytoscape output",
                         textInput(
                           inputId = "threshold",
                           label = "weight threshold",
                           value = 0.02
                         ),
                         actionButton("threadd","choose the threshold"),
                         DT::dataTableOutput("edgeFile"),
                         DT::dataTableOutput("nodeFile"),
                         downloadButton("downtbl5","download edgefile"),
                         downloadButton("downtbl6","download nodefile")
                       )
                     )
                   )
                 )
               )
             )##tabPanel
             
  )## navbarPage
)## UI


# Server Backend ----
server <- function(input, output, session){
  
  # Section0 - Data load ----
  data = reactive({
    file1 = "./expression_FPKM.csv"
    if(is.null(file1)){return()}
    
    file_main = read.delim(file=file1,
               sep=",",
               header = T,
               stringsAsFactors = F,
              )
    
    convert_numeric <- function(x){
      return (as.numeric(x))}
    
    rownames(file_main) = file_main[,1]
    file_main_new = file_main[,-1]
    
    
    file_main_new = as.data.frame(lapply(file_main_new,convert_numeric))
    rownames(file_main_new) = file_main[,1]
    file_main_new = file_main_new[which(rowSums(file_main_new) > 0),]
    file_main_new$ID = rownames(file_main_new)
    file_main_new = file_main_new %>% dplyr::select('ID',colnames(file_main_new)[1:dim(file_main_new)[2]-1],everything())
    # print(file_main_new[1:5,1:5])
    # print(summary(file_main_new))
    
    file_main_new
  })
  
  # Section1 - Check Input Format, for blank values ----
  output$Inputcheck = renderUI({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) == 0) {
      HTML(paste0('<br><font color = blue><b>[SYSTEM] </b></font>Expression matrix has been successfully loaded onto the HEMU WGCNA interactive analysis module.',
                  '<br><font color = blue><b>[SYSTEM] </b></font>Loaded ', nrow(data()), ' genes and ', ncol(data())-1, ' samples.'
                  ))
    } else {
      HTML(
        '<br><hr><font color = blue><b>
        [SYSTEM] </b></font>Expression matrix contains blank (NA) values or rows, please inspect data structure.
        ')
    }
  })
  
  # Section1 - Count Number Normalization ----
  fmt = reactive({
    "FPKM"
  })
  observe({
      updateSelectInput(session, "method1",choices = c(FPKM = "rawFPKM", logFPKM = "lgFPKM"))
      updateTextInput(session,"RCcut",value = 1)
  })
  
  mtd = reactive({
    input$method1
  })
  
  sampP = reactive({
    as.numeric(input$SamPer)
  })
  rccutoff = reactive({
    as.numeric(input$RCcut)
  })
  GNC = reactive({
    as.numeric(input$remain)
  })
  cutmethod = reactive({
    input$CutMethod
  })
  ## set reactiveValues
  exp.ds<-reactiveValues(data=NULL)
  downloads <- reactiveValues(data = NULL)
  
  # Section1 - Expression Data Filtration ----
  observeEvent(
    input$action1, # Update Information
    {
      if(is.null(data())){return()}
      if(length(which(is.na(data()))) != 0) {return()}
      
      exp.ds$table = data.frame()
      exp.ds$table2 = data.frame()
      exp.ds$param = list()
      exp.ds$layout = as.character(input$treelayout)
      
      output$filter1 = renderUI({
        input$action1
        p_mass = c("Processing step1, remove very low expressed genes",
                   paste("Processing step2, pick out high variation genes via",cutmethod()))
        withProgress(
          message = "Raw data normlization",
          value = 0,{
            for (i in 1:2) {
              incProgress(1/2,detail = p_mass[i]) # Progress bar configuration
              if(i == 1) {
                exp.ds$table = getdatExpr(rawdata = data(),
                                          RcCutoff = rccutoff(),samplePerc = sampP(),
                                          datatype = fmt(),method = mtd())
                # print("Completed, normalization, step1")
                
              } else if (i == 2){
                # exp.ds$table2 = getdatExpr2(datExpr = exp.ds$table,
                #                             GeneNumCut = 1-GNC()/nrow(exp.ds$table),cutmethod = cutmethod())
                exp.ds$table2 = t(exp.ds$table)
                # print("Completed, normalization, step2")
                exp.ds$param = getsampleTree(exp.ds$table2,layout = exp.ds$layout)
                # print("Completed, plotting sample clustering dendrogram")
                
              }
              Sys.sleep(0.1)
            }
          }
        )
        isolate(HTML(paste0('<font color = blue><b>[SYSTEM] </b></font>Filter threshold: FPKM < <font color = blue><b>',rccutoff(),'</b></font> in > <font color = blue> <b>',100*sampP(),'% </b></font>samples','<br/>',
                            '<font color = blue><b>[SYSTEM] </b></font>
Genes after filtering: ',ncol(exp.ds$table2))))
      })
    }
  )
  
  
  
  
  ## summary num ----
  output$Inputbl = DT::renderDataTable({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) != 0) {return()}
    as.data.frame(t(exp.ds$table2))
  })
  ## sample tree
  output$clustPlot = renderPlot({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) != 0) {return()}
    if(is.null(exp.ds$table2)){return()}
    plot(exp.ds$param$sampleTree,main = "Sample clustering dendrogram", sub = "", xlab = "")
  })
  
  ##remove outlier
  s_outlier = reactive({
    colnames(data())
  })
  
  observe({
    updateSelectInput(session, "outlier",choices = s_outlier())
  })
  observeEvent(
    input$Startremove,
    {
      exp.ds$outliersamples = as.character(input$outlier)
      if(!exp.ds$outliersamples %in% s_outlier()) {return()} else {
        exp.ds$rmoutlier = 
          mv_outlier(x = data(),y = exp.ds$outliersamples)
      }
      output$new_mat_preview = DT::renderDataTable({
        if(is.null(exp.ds$table2)){return()}
        if(!exp.ds$outliersamples %in% s_outlier()) {return()}
        as.data.frame(exp.ds$rmoutlier)
      })
    }
  )
  
  
  
  
  ## download sample tree ----
  
  rscut = reactive({
    as.numeric(input$CutoffR)
  })
  
  observeEvent(
    input$Startsft,
    {
      if(is.null(exp.ds$table2)){return()}
      exp.ds$sft = list()
      output$powerout = renderUI({
        sft_mess = c("pick soft threshold in processing ...",
                     "Finish.")
        withProgress(message = 'SFT selection', value = 0,
                     expr = {
                       for (i in 1:2) {
                         incProgress(1/2, detail = sft_mess[i] )
                         if (i == 1) {
                           exp.ds$sft = getpower(datExpr = exp.ds$table2,rscut = rscut())
                         } else {
                           return()
                         }
                         
                       }
                     })
        isolate(HTML(paste0('<font color = red> <b>The power recommended by WGCNA is:</b> </font><font color = bule><b>',exp.ds$sft$power,'</b></font> ','<br/>',
                            '<font color = gray> <i>If all power values lower than the R square threshold which you set, it means that the power value is an empirical value. At this time, you need to infer a power value based on the results on this plot and check whether it can form a scale-free network. </i> </font>')))
      })
    }
  )
  
  
  # outsft ----
  output$sftplot = renderPlot({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    exp.ds$sft$plot
  })
  ## outtbl
  output$sfttbl = DT::renderDataTable({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    as.data.frame(exp.ds$sft$sft)
  })
  ## test sft
  pcus = reactive({
    as.numeric(input$PowerSelect)
  })
  PowerTorF = reactive({
    input$PowerTorF
  })
  observeEvent(
    input$Startcheck,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$sft)){return()}
      sftcheck_mess = c("Checking scale free network ...",
                        "Finish.")
      exp.ds$power = exp.ds$sft$power
      exp.ds$cksft = list()
      withProgress(message = 'SFT selection', value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2, detail = sftcheck_mess[i] )
                       if (i == 1) {
                         if(PowerTorF() == "Recommended"){
                           exp.ds$power = exp.ds$sft$power
                           exp.ds$cksft = powertest(power.test = exp.ds$sft$power,datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         } else if (PowerTorF() == "Customized"){
                           exp.ds$power = pcus()
                           exp.ds$cksft = powertest(power.test = pcus(),datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         }
                       } else {
                         return()
                       }
                       
                     }
                   })
    }
  )
  
  output$sfttest = renderPlot({
    if(is.null(exp.ds$sft)){return()}
    input$Startcheck
    if(length(exp.ds$cksft ) == 0){return()}
    exp.ds$cksft
  })
  mms = reactive({
    as.numeric(input$minMsize)
  })
  mch = reactive({
    as.numeric(input$mch)
  })
  blocksize = reactive({
    as.numeric(input$blocksize)
  })
  observeEvent(
    input$Startnet,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$power)){return()}
      exp.ds$netout = getnetwork(datExpr = exp.ds$table2,power = exp.ds$power,
                                 minModuleSize = mms(),mergeCutHeight = mch(),maxBlocksize = blocksize())
      exp.ds$nSamples = nrow(exp.ds$table2)
      exp.ds$net = exp.ds$netout$net
      exp.ds$moduleLabels = exp.ds$netout$moduleLabels
      exp.ds$moduleColors = exp.ds$netout$moduleColors
      exp.ds$MEs_col = exp.ds$netout$MEs_col
      exp.ds$MEs = exp.ds$netout$MEs
      exp.ds$Gene2module = exp.ds$netout$Gene2module
    }
  )
  output$cluster = renderPlot({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  })
  output$m2num = renderTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    table(exp.ds$moduleColors)
  })
  output$eah = renderPlot({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2), plotDendrograms = T,
                          xLabelsAngle = 90)
  })
  
  output$g2m = DT::renderDataTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    exp.ds$Gene2module
  })
  
  phen <- reactive({
    file2 = "./trait_data.csv"
    if(is.null(file2)){return()}
    file_main2 = read.delim(file=file2,
                           sep=",",
                           header = T,
                           stringsAsFactors = F,
    )
    file_main2
  })
  
  
  observeEvent(
    input$starttrait,
    {
      if(is.null(phen())){return()}
      if (ncol(phen()) == 2) {
        x <- phen()
        Tcol = as.character(unique(x[,2]))
        b <- list()
        for (i in 1:length(Tcol)) {
          b[[i]] = data.frame(row.names = x[,1],
                              levels = ifelse(x[,2] == Tcol[i],1,0))
        }
        c <- bind_cols(b)
        c <- data.frame(row.names = x$name,
                        c)
        colnames(c) = Tcol
        rownames(c) = phen()[,1]
        exp.ds$phen<- c
      } else {
        exp.ds$phen = data.frame(row.names = phen()[,1],
                                 phen()[,-1])
      }
      exp.ds$phen =  exp.ds$phen[match(rownames(exp.ds$table2),rownames(exp.ds$phen)),]
      exp.ds$traitout = getMt(phenotype = exp.ds$phen,
                              nSamples = exp.ds$nSamples,moduleColors = exp.ds$moduleColors,datExpr = exp.ds$table2)
      exp.ds$xangle = as.numeric(input$xangle)
      exp.ds$c_min = as.character(input$colormin)
      exp.ds$c_mid = as.character(input$colormid)
      exp.ds$c_max = as.character(input$colormax)
      exp.ds$modTraitCor = exp.ds$traitout$modTraitCor
      exp.ds$modTraitP = exp.ds$traitout$modTraitP
      exp.ds$textMatrix = exp.ds$traitout$textMatrix
      exp.ds$KME = getKME(datExpr = exp.ds$table2,moduleColors = exp.ds$moduleColors,MEs_col = exp.ds$MEs_col)
      exp.ds$mod_color = gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitCor))
      exp.ds$mod_color_anno = setNames(exp.ds$mod_color,rownames(exp.ds$modTraitCor))
      exp.ds$Left_anno = rowAnnotation(
        Module = rownames(exp.ds$modTraitCor),
        col = list(
          Module = exp.ds$mod_color_anno
        ),
        show_legend = F,
        show_annotation_name = F
      )
    }
  )
  
  output$mtplot = renderPlot({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    Heatmap(
      matrix = exp.ds$modTraitCor,
      cluster_rows = F, cluster_columns = F,
      left_annotation = exp.ds$Left_anno,
      cell_fun = function(j,i,x,y,width,height,fill) {
        grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
      },
      row_names_side = "left",
      column_names_rot = exp.ds$xangle,
      heatmap_legend_param = list(
        at = c(-1,-0.5,0,0.5, 1),
        labels = c("-1","-0.5", "0","0.5", "1"),
        title = "",
        legend_height = unit(9, "cm"),
        title_position = "lefttop-rot"
      ),
      rect_gp = gpar(col = "black", lwd = 1.2),
      column_title = "Module-trait relationships",
      column_title_gp = gpar(fontsize = 15, fontface = "bold"),
      col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
    )
  })
  
  output$KME = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$KME)
  })
  
  observeEvent(
    input$run_filter,
    {
      if(is.null(phen())){return()};
      exp.ds$kme_method = as.numeric(input$inter_method);
      exp.ds$kme_cutoff = as.numeric(input$kme_cutoff);
      exp.ds$kme_outlist = iterative_out(
        g2m = exp.ds$Gene2module,
        rawMat = data(),
        tbl = as.data.frame(exp.ds$KME),
        method = exp.ds$kme_method,
        KME_cutoff = exp.ds$kme_cutoff
      )
      message("kme analysis finish")
    }
  )
  
  output$Iter_retained = DT::renderDataTable({
    input$run_filter
    if(is.null(phen())){return()}
    as.data.frame(exp.ds$kme_outlist$retain)
  })
  
  output$Iter_removed = DT::renderDataTable({
    input$run_filter
    if(is.null(phen())){return()}
    as.data.frame(exp.ds$kme_outlist$remove)
  })
  
  s_mod = reactive({
    gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitP))
  })
  
  observe({
    updateSelectInput(session, "smodule",choices = s_mod())
  })
  
  s_trait = reactive({
    colnames(exp.ds$modTraitP)
  })
  observe({
    updateSelectInput(session, "strait",choices = s_trait())
  })
  
  observeEvent(
    input$InterMode,
    {
      if(is.null(phen())){return()}
      if(is.null(exp.ds$phen)){return()}
      exp.ds$GSout = getMM(datExpr = exp.ds$table2,MEs_col = exp.ds$MEs_col,nSamples = exp.ds$nSamples,corType = "pearson")
      exp.ds$MM = exp.ds$GSout$MM
      exp.ds$MMP = exp.ds$GSout$MMP
      exp.ds$sml = as.character(input$smodule)
      exp.ds$st = as.character(input$strait)
      exp.ds$Heatmap = moduleheatmap(datExpr = exp.ds$table2,MEs = exp.ds$MEs_col,which.module = exp.ds$sml,
                                     moduleColors = exp.ds$moduleColors)
    }
  )
  
  output$GSCon = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                   traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                   geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples)
  })
  
  output$heatmap = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    exp.ds$Heatmap
  })
  
  output$GSMM.all = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    MMvsGSall(which.trait = exp.ds$st,
              traitData = exp.ds$phen,
              datExpr = exp.ds$table2,
              moduleColors = exp.ds$moduleColors,
              geneModuleMembership = exp.ds$MM,
              MEs = exp.ds$MEs_col,
              nSamples = exp.ds$nSamples)
  })
  
  observe({
    updateSelectInput(session, "hubmodule",choices = s_mod())
  })
  
  observe({
    updateSelectInput(session, "hubtrait",choices = s_trait())
  })
  
  observeEvent(
    input$starthub,
    {
      exp.ds$hubml = as.character(input$hubmodule)
      exp.ds$hubt = as.character(input$hubtrait)
      exp.ds$kMEcut = as.numeric(input$kMEcut)
      exp.ds$GScut = as.numeric(input$GScut)
      print(exp.ds$hubml)
      exp.ds$hub.all = hubgenes(datExpr = exp.ds$table2,
                                mdl = exp.ds$hubml,
                                power = exp.ds$power,
                                trt = exp.ds$hubt,
                                KME = exp.ds$KME,
                                GS.cut = exp.ds$GScut,
                                kME.cut =exp.ds$kMEcut,
                                datTrait = exp.ds$phen,
                                g2m = exp.ds$Gene2module
      )
    }
  )
  
  observeEvent(
    input$threadd,
    {
      exp.ds$threshold = as.numeric(input$threshold)
      exp.ds$cyt = cytoscapeout(datExpr = exp.ds$table2,
                                power = exp.ds$power,module = exp.ds$hubml,
                                moduleColors = exp.ds$moduleColors,
                                threshold = exp.ds$threshold)
    }
  )
  output$kMEhub = DT::renderDataTable({
    input$starthub
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$hubt)){return()}
    if(is.null(exp.ds$kMEcut)){return()}
    if(is.null(exp.ds$GScut)){return()}
    exp.ds$hub.all$hub3
  })
  
  output$edgeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[1]]
  })
  
  output$nodeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[2]]
  })
  
  observeEvent(
    input$Integrate_data,
    {
      exp.ds$data_interagrate = list(
        expmat = exp.ds$table2,
        traitmat = exp.ds$phen,
        expmat_format = as.character(fmt()),
        method = as.character(mtd()),
        samplePercentage = as.numeric(sampP()),
        rccutoff = as.numeric(rccutoff()),
        GNC = as.numeric(GNC()),
        cutmethod = as.character(cutmethod()),
        power = exp.ds$power,
        min_module_size = as.numeric(mms()),
        module_cuttree_height = as.numeric(mch()),
        blocksize = as.numeric(blocksize()),
        net = exp.ds$net,
        moduleLabels = exp.ds$moduleLabels,
        moduleCoolors = exp.ds$moduleColors,
        MEs = exp.ds$MEs,
        MEs_col = exp.ds$MEs_col,
        exp.ds$Gene2module,
        modTraitCor = exp.ds$modTraitCor,
        modTraitP = exp.ds$modTraitP,
        textMatrix = exp.ds$textMatrix,
        KME = exp.ds$KME
      )
      outrsd <<- exp.ds$data_interagrate
    }
  )
  
  # Plot and data download ----
  observeEvent(
    input$adjust1,
    {
      downloads$width1 <- as.numeric(input$width1)
      downloads$height1 <-  as.numeric(input$height1)
    }
  )
  observeEvent(
    input$adjust2,
    {
      downloads$width2 <- as.numeric(input$width2)
      downloads$height2 <-  as.numeric(input$height2)
    }
  )
  observeEvent(
    input$adjust3,
    {
      downloads$width3 <- as.numeric(input$width3)
      downloads$height3 <-  as.numeric(input$height3)
    }
  )
  observeEvent(
    input$adjust4,
    {
      downloads$width4 <- as.numeric(input$width4)
      downloads$height4 <- as.numeric(input$height4)
    }
  )
  observeEvent(
    input$adjust5,
    {
      downloads$width5 <- as.numeric(input$width5)
      downloads$height5 <-  as.numeric(input$height5)
    }
  )
  observeEvent(
    input$adjust6,
    {
      downloads$width6 <- as.numeric(input$width6)
      downloads$height6 <-  as.numeric(input$height6)
    }
  )
  observeEvent(
    input$adjust7,
    {
      downloads$width7 <- as.numeric(input$width7)
      downloads$height7 <-  as.numeric(input$height7)
    }
  )
  observeEvent(
    input$adjust8,
    {
      downloads$width8 <- as.numeric(input$width8)
      downloads$height8 <-  as.numeric(input$height8)
    }
  )
  observeEvent(
    input$adjust10,
    {
      downloads$width10 <- as.numeric(input$width10)
      downloads$height10 <-  as.numeric(input$height10)
    }
  )
  library(ape)
  output$downfig1 = downloadHandler(
    filename = function() {
      "01.SampleCluster.nwk"
    },
    content = function(file) {
      write.tree(phy = exp.ds$param$tree,file = file)
    }
  )
  output$downfig2 = downloadHandler(
    filename = function() {
      "02.SftResult.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$sft$plot,filename = file,width = downloads$width2,height = downloads$height2)
    }
  )
  output$downfig3 = downloadHandler(
    filename = function() {
      "03.CheckSft.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$cksft,filename = file,width = downloads$width3,height = downloads$height3)
    }
  )
  output$downfig4 = downloadHandler(
    filename = function() {
      "04.ClusterDendrogram.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width4, height = downloads$height4 )
      plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
    }
  )
  output$downfig5 = downloadHandler(
    filename = function() {
      "05.EigengeneadJacencyHeatmap.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width5, height = downloads$height5)
      plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = T,
                            xLabelsAngle = 90)
      dev.off()
    }
  )
  output$downfig6 = downloadHandler(
    filename = function() {
      "06.Module2Trait.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width6, height = downloads$height6)
      
      print(Heatmap(
        matrix = exp.ds$modTraitCor,
        cluster_rows = F, cluster_columns = F,
        left_annotation = exp.ds$Left_anno,
        cell_fun = function(j,i,x,y,width,height,fill) {
          grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        column_names_rot = exp.ds$xangle,
        heatmap_legend_param = list(
          at = c(-1,-0.5,0,0.5, 1),
          labels = c("-1","-0.5", "0","0.5", "1"),
          title = "",
          legend_height = unit(9, "cm"),
          title_position = "lefttop-rot"
        ),
        rect_gp = gpar(col = "black", lwd = 1.2),
        column_title = "Module-trait relationships",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
      ))
      
      dev.off()
    }
  )
  output$downfig7 = downloadHandler(
    filename = function() {
      paste0("07.GS",exp.ds$sml,"-",exp.ds$st,"-Connectivity.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width7, height = downloads$height7)
      print(getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                           traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                           geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples))
      dev.off()
    }
  )
  output$downfig8 = downloadHandler(
    
    filename = function() {
      paste0("08.",exp.ds$sml,"-",exp.ds$st,"MEandGeneHeatmap.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(exp.ds$Heatmap)
      dev.off()
    }
  )
  output$downfig10 = downloadHandler(
    
    filename = function() {
      "09.GSvsMM.all.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(MMvsGSall(which.trait = exp.ds$st,
                      traitData = exp.ds$phen,nSamples = exp.ds$nSamples,
                      datExpr = exp.ds$table2,
                      moduleColors = exp.ds$moduleColors,
                      geneModuleMembership = exp.ds$MM,MEs = exp.ds$MEs_col))
      dev.off()
    }
  )
  output$downtbl2 = downloadHandler(
    
    filename = function() {
      if(is.null(exp.ds$net)){return()}
      "01.Gene2Module.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$Gene2module,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl3 = downloadHandler(
    filename = function() {
      "02.KMEofAllGenes.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$KME,file = file,sep = "\t",row.names = T,quote = F)
    }
  )
  output$downtbl4 = downloadHandler(
    filename = function() {
      paste0("03.",exp.ds$hubml,"-",exp.ds$hubt,"hubgene_by_GS_MM.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$hub.all$hub3,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl5 = downloadHandler(
    filename = function() {
      paste0("04.",exp.ds$hubml,".edge.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[1]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl6 = downloadHandler(
    filename = function() {
      paste0("04.cyt",exp.ds$hubml,".node.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[2]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtabout = downloadHandler(
    filename = function() {
      paste0("00.Remove_outlier_Table.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$rmoutlier,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_Iter_1 = downloadHandler(
    filename = function() {
      paste0("00.Retained_GeneSet_for_Next_Round.xls")
    },
    content = function(file) {
      write.table(x = as.data.frame(exp.ds$kme_outlist$retain),file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_Iter_2 = downloadHandler(
    filename = function() {
      paste0("00.Removed_GeneSet.xls")
    },
    content = function(file) {
      write.table(x = as.data.frame(exp.ds$kme_outlist$remove),file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_EnvImg = downloadHandler(
    filename = function() {
      paste0("EnvImg.rds")
    },
    content = function(file) {
      saveRDS(outrsd,file = file)
    }
  )
}


shinyApp(ui,server)
    
    
    """

    try:
        with open(output_app_dir + '/app.R', mode='w') as out_fh:
            out_fh.write(wgcna_mainapp_text_raw)  # Create shiny application
        return output_app_dir, output_app_dir_name

    except:
        return None
