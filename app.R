# To run this app (as a developer), start R and then:
# source("app.R")
# runApp()

library(shiny)
library(readr)
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(Glimma)
library(gplots)
library(msigdbr)
library(DT)


getFileChoices <- function() {
  files <- list.files()
  files[grep("counts.txt", files, fixed=TRUE)]
}

filesAvailable <- getFileChoices() # Get file list to populate the drop down menu


ui <- fluidPage(
  
  # App title ----
  titlePanel("RNAeasy123"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("fileChooser", "Choose a file",
                  choices = filesAvailable, selected = filesAvailable),
      # for now we are hardcoding the cell types
      # in future we want to get the distinct cell type values 
      # for each file when the file is selected
      # TODO fix that....
      selectInput("cellTypeChooser", "Choose a cell type",
                  choices = c("DC", "DP"), selected = "DC"
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
                  tabPanel("# of DE Genes",
                    selectInput("testChooser", "Choose a statistical test",
                      choices=c("efit", "tfit")),
                    DTOutput("testTable")
                  )
      )
    )
  )
)

# add an action button on the left side



server <- function(input, output, session) {
  # TODO add a function to get the data
  
  getData <- eventReactive(input$submitButton, {
    # print(sessionInfo())
    counts <- data.frame(read_delim(input$fileChooser,
                                    "\t", escape_double = FALSE, trim_ws = TRUE))
    theseCells <- dplyr::select(counts, contains(input$cellTypeChooser))
    rownames(theseCells) <- counts$GeneSymbol
    theseCells
  })
  getMetaData <- eventReactive(input$submitButton, {
    metadataFile <- gsub("\\.txt$", ".csv", input$fileChooser)
    metadataFile <- gsub("\\.counts\\.", ".metadata.", metadataFile)
    metadata <- data.frame(read.csv(metadataFile), stringsAsFactors = TRUE)
    thisMeta <- dplyr::filter(metadata, grepl(input$cellTypeChooser, sampleName)==TRUE)
    thisMeta
  })
  
  output$mytable = DT::renderDT({
    getData()
  }, rownames = TRUE, colnames = c('Gene Symbol' = 1))
  

  makeMVPlot <- eventReactive(input$submitButton, {
    cells <- getData()
    metadata <- getMetaData()
    counts.m <- as.matrix(cells)
    data <- DGEList(counts.m)
    l <- makeUSPlot()
    #keep.exprs <- filterByExpr(l$data, group=metadata$group)
    #data <- l$data[keep.exprs,, keep.lib.sizes=FALSE]
    #Normalising gene expression distributions
    data <- calcNormFactors(data, method = "TMM")
    #keep.exprs <- filterByExpr(cells, group=metadata$group)
    #data <- data[keep.exprs,, keep.lib.sizes=FALSE]
    #data <- calcNormFactors(data, method = "TMM")
    print(metadata$group)

    design <- model.matrix(~0+metadata$group)
    str(design)
    colnames(design) <- gsub("group", "", colnames(design))
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
    metadata <- getMetaData()
    counts.m <- as.matrix(cells)
    data <- DGEList(counts.m)
    colnames(data) <- metadata$sampleName
    data$samples$batch <- metadata$batch
    data$samples$group <- metadata$group
    data$genes <- select(Mus.musculus, keys=rownames(cells), 
                         columns=c("ENTREZID", "TXCHROM"), 
                         keytype="SYMBOL") %>% filter(!duplicated(SYMBOL))
    lcpm <- cpm(data, log=TRUE)
    # col.batch <- batch
    col.batch <- metadata$batch
    levels(col.batch) <-  brewer.pal(nlevels(metadata$batch), "Set2")
    col.batch <- as.character(col.batch)
    return (list(lcpm=lcpm, batch=metadata$batch, col.batch=col.batch, data=data, dim=c(3,4)))
  })
  
  
  output$leftUSPlot <- renderPlot({
    left <- makeUSPlot()
    plotMDS(left$lcpm, labels=left$group, col=left$col.group)
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

