# To run this app (as a developer), start R and then:
# source("app.R")
# runApp()

library(shiny); library(readr); library(dplyr)
library(limma); library(Glimma); library(edgeR)
library(Mus.musculus)
library(RColorBrewer); library(gplots)
library(msigdbr)
library(DT); #library(plotly) # will likely need this later

getProjects <- function() {
  read.csv("projectDataSets.csv", stringsAsFactors = FALSE)
}

projects<- getProjects()


ui <- fluidPage(
  
  # App title ----
  titlePanel("RNAeasy123"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 3,
      selectInput("projectChooser", "Choose a project dataset",
                  choices = projects$projectName, selected = ""),
      selectInput("groupChooser", "Choose a dataset group",
                  choices = "", selected = ""
      ),
      actionButton(inputId = "submitButton",
                   label = "Submit")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Counts", DTOutput("mytable")),
                  tabPanel("Unsup. Clustering",
                           fluidRow(
                             # plot titles are hardcoded - TODO FIXME
                             column(6, h2(align="center", "A. SL-TBI"), plotOutput("leftUSPlot")),
                             column(6, h2(align="center", "B. Batch"), plotOutput("rightUSPlot"))
                           )
                  ),
                  tabPanel("Mean-variance trends",
                    fluidRow(
                      column(6, h2(align="center", "voom: Mean-variance trend"), plotOutput("leftMVPlot")),
                      column(6, h2(align="center", "Final model: Mean-variance trend"), plotOutput("rightMVPlot"))
                    )
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  # This gets the metadata for the chosen project
  getMetaData <- eventReactive(input$projectChooser, {
    # print(sessionInfo())
    whichData <- filter(projects, projects$projectName == input$projectChooser)
    studyDesign <- read.csv(whichData$studyDesign, stringsAsFactors = FALSE)
    return(studyDesign)
  })
  
  # This observes what the user chooses for the project and updates the values possible in the group dropdown
  observe({
           values <- unique(getMetaData()$cellType)
    # Can also set the label and select items
    updateSelectInput(session, "groupChooser",
                      choices = values)
  })
  
  # This goes and gets the right dataset after the user clicks submit, subsets and reformats that data to show only the specific data group
  getData <- eventReactive(input$submitButton, {
      # print(sessionInfo())
      whichData <- filter(projects, projects$projectName == input$projectChooser)
      counts <- read.table(whichData$countsData, header = TRUE)
      rownames(counts) <- counts$GeneSymbol
      counts$GeneSymbol <- NULL
      whichSamples <- getMetaData() %>% filter(cellType == input$groupChooser)
      ifelse(input$projectChooser == "TGR", dataToUse <- whichSamples$molecular_id,
             dataToUse <- whichSamples$sampleName)
      colnames(counts)<- gsub("\\.", "-", colnames(counts)) # WTF????  Why is it reading in the file and then teh column names have .'s in them instead of -'s?
      #print(head(rownames(counts)))
      countsSub <- counts[,dataToUse]
      rownames(countsSub) <- rownames(counts) #WTFx2!!!! This refused to replace the rowname indexes with actual text strings on one computer.  
      #print(head(rownames(countsSub)))
      return(countsSub)
    })
    
  # This creates a rendered counts table with Gene Symbol as a column
  output$mytable = DT::renderDT({
    getData()
  }, rownames = TRUE, colnames = c('Gene Symbol' = 1))
  

  makeMVPlot <- eventReactive(input$submitButton, {
    cells <- getData()
    metadata <- getMetaData()
    counts.m <- as.matrix(cells)
    data <- DGEList(counts.m)
    l <- makeUSPlot()
    #Normalising gene expression distributions
    data <- calcNormFactors(data, method = "TMM")
    design <- model.matrix(~0+metadata$timepoint)
    str(design)
    colnames(design) <- gsub("timepoint", "", colnames(design))
    contr.matrix <- makeContrasts(
      DC.D0vsD4 = D0 - D4, 
      DC.D0vsD7 = D0 - D7,
      DC.D4vsD7 = D4 - D7,
      levels = colnames(design))
    # par(mfrow=c(1,2))
    v <- voom(data, design, plot=TRUE)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    tfit <- treat(vfit, lfc=log2(1.5))
    print(paste("class of tfit is", class(tfit)))
    # plotSA(tfit, main="Final model: Mean-variance trend")

    return(list(tfit=tfit))
  })

  makeUSPlot <- eventReactive(input$submitButton, {
    cells <- getData()
    metadata <- getMetaData() %>% filter(cellType == input$groupChooser)
    counts.m <- as.matrix(cells)
    rownames(counts.m) <- rownames(cells)
    data <- DGEList(counts.m)
    colnames(data) <- metadata$sampleName
    data$samples$batch <- metadata$batch
    data$samples$timepoint <- metadata$timepoint
    data$genes <- suppressMessages(select(Mus.musculus, keys=rownames(cells), 
                         columns=c("ENTREZID", "TXCHROM"), 
                         keytype="SYMBOL") %>% filter(!duplicated(SYMBOL)))
    lcpm <- cpm(data, log=TRUE)
    col.batch <- as.factor(metadata$batch)
    levels(col.batch) <-  brewer.pal(nlevels(col.batch)+1, "Dark2")
    col.batch <- as.character(col.batch)
    return (list(lcpm=lcpm, batch=metadata$batch, col.batch=col.batch, data=data, dim=c(3,4)))
  })
  
  output$leftUSPlot <- renderPlot({
    left <- makeUSPlot()
    plotMDS(left$lcpm, labels=left$timepoint, col=left$col.timepoint)
  })
  
  output$rightUSPlot <- renderPlot({
    right <- makeUSPlot()
    plotMDS(right$lcpm, labels=right$batch, col=right$col.batch, dim=right$dim)
  })

  output$leftMVPlot <- renderPlot({
    par(mfrow=c(1,2))
    plotSA(makeMVPlot()$tfit, main="Final model: Mean-variance trend")
  })
  
}

shinyApp(ui = ui, server = server)

