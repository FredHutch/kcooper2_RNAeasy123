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


getCounts <- function() {
  counts <- read_delim("2019.07.12.counts.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
  # # rownames(counts) <- counts$GeneSymbol

}

getFileChoices <- function() {
  files <- list.files()
  files[grep(".txt", files, fixed=TRUE)]
}



ui <- fluidPage(

  # App title ----
  titlePanel("RNAeasy123"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
        selectInput("fileChooser", "Choose a file",
          getFileChoices()),
        # for now we are hardcoding the cell types
        # TODO fix that....
        selectInput("cellTypeChooser", "Choose a cell type",
          c("DC", "DP"))

    ),

    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type="tabs",
        tabPanel("Counts", dataTableOutput("mytable")),
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

makeUSPlot <- function(input) {
    counts <- read_delim(input$fileChooser, 
        "\t", escape_double = FALSE, trim_ws = TRUE)
    rownames(counts) <- counts$GeneSymbol

    cells <- dplyr::select(counts, contains(input$cellTypeChooser))
    counts.m <- as.matrix(cells)
    data <- DGEList(counts.m)
    samplenames <- substring(colnames(data), 1, nchar(colnames(data)))
    colnames(data) <- samplenames
    # batch is hardcoded atm - TODO FIXME
    batch <- as.factor(c("09182018", "09182018", 
                     "09182018", "09182018", "09182018", 
                     "10312018", "10312018"))
    data$samples$batch <- batch
    group <- as.factor(rep(c("D0", "D4", "D7"), c(2,3,2)))
    data$samples$group <- group
    geneid <- rownames(cells)
    geneid_from_shiny <- geneid
    save("geneid_from_shiny", file="geneid_from_shiny.rda")
    genes <- select(Mus.musculus, keys=geneid, columns=c("ENTREZID", "TXCHROM"), 
                    keytype="SYMBOL")
    genes <- genes[!duplicated(genes$SYMBOL),]
    data$genes <- genes
    cpm <- cpm(data)
    lcpm <- cpm(data, log=TRUE)
    L <- mean(data$samples$lib.size) * 1e-6
    M <- median(data$samples$lib.size) * 1e-6
    keep <- filterByExpr(data)
    edgeR.data <- data[keep, , keep.lib.sizes=FALSE]
    keep.exprs <- filterByExpr(data, group=group)
    data <- data[keep.exprs,, keep.lib.sizes=FALSE]
    data <- calcNormFactors(data, method = "TMM")
    cpm <- cpm(data)
    lcpm <- cpm(data, log=TRUE)
    par(mfrow=c(1,2))
    col.group <- group
    levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
    col.group <- as.character(col.group)
    col.batch <- batch
    levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
    col.batch <- as.character(col.batch)
    return (list(lcpm=lcpm, group=group, batch=batch, col.group=col.group, col.batch=col.batch, dim=c(3,4)))
}



server <- function(input, output) {
  output$mytable = shiny::renderDataTable({
      read_delim(input$fileChooser, 
        "\t", escape_double = FALSE, trim_ws = TRUE)
  })


  output$leftUSPlot <- renderPlot({
    left <- makeUSPlot(input)
    plotMDS(left$lcpm, labels=left$group, col=left$col.group)
  })

  output$rightUSPlot <- renderPlot({
    right <- makeUSPlot(input)
    plotMDS(right$lcpm, labels=right$batch, col=right$col.batch, dim=right$dim)
  })

}

shinyApp(ui = ui, server = server)

