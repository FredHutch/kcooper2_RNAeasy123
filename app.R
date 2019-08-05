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
  files[grep(".txt", files, fixed=TRUE)]
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
          filesAvailable),
        # for now we are hardcoding the cell types
        # TODO fix that....
        selectInput("cellTypeChooser", "Choose a cell type",
          c("DC", "DP")
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
            column(6, plotOutput("leftUSPlot")),
            column(6, plotOutput("rightUSPlot"))
          )
        )
      )
    )
  )
)

# add an action button on the left side



server <- function(input, output) {
  # TODO add a function to get the data

  getData <- eventReactive(input$submitButton, {
      print(sessionInfo())
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

  makeUSPlot <- eventReactive(input$cellTypeChooser, {
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
    return (list(lcpm=lcpm, batch=metadata$batch, col.batch=col.batch, dim=c(3,4)))
  })
  

  output$leftUSPlot <- renderPlot({
    left <- makeUSPlot()
    plotMDS(left$lcpm, labels=left$group, col=left$col.group)
  })

  output$rightUSPlot <- renderPlot({
    right <- makeUSPlot()
    plotMDS(right$lcpm, labels=right$batch, col=right$col.batch, dim=right$dim)
  })

}

shinyApp(ui = ui, server = server)

