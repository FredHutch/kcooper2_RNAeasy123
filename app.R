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

getData <- function(input) {
  # redundant reading of counts file:
  counts <- read_delim(input$fileChooser, 
        "\t", escape_double = FALSE, trim_ws = TRUE)
  cells <- dplyr::select(counts, contains(input$cellTypeChooser))

  counts.m <- as.matrix(cells)

  data <- DGEList(counts.m)
  samplenames <- substring(colnames(data), 1, nchar(colnames(data)))
  batch <- as.factor(c("09182018", "09182018", 
                     "09182018", "09182018", "09182018", 
                     "10312018", "10312018"))
  group <- as.factor(rep(c("D0", "D4", "D7"), c(2,3,2)))
  data$samples$group <- group
  data
}

server <- function(input, output) {
  output$mytable = shiny::renderDataTable({
      read_delim(input$fileChooser, 
        "\t", escape_double = FALSE, trim_ws = TRUE)
  })

  # cells <- reactive({dplyr::select(counts, contains(input$cellTypeChooser)) })
  # counts.m <- reactive({as.matrix(cells)})
  # data <- getData(input)

  output$leftUSPlot <- renderPlot({
    # redundant read of counts
    counts <- read_delim(input$fileChooser, 
        "\t", escape_double = FALSE, trim_ws = TRUE)
    cells <- dplyr::select(counts, contains(input$cellTypeChooser))
    counts.m <- as.matrix(cells)
    data <- DGEList(counts.m)
    samplenames <- substring(colnames(data), 1, nchar(colnames(data)))
    colnames(data) <- samplenames
    # batch is hardcoded atm
    batch <- as.factor(c("09182018", "09182018", 
                     "09182018", "09182018", "09182018", 
                     "10312018", "10312018"))
    data$samples$batch <- batch
    group <- as.factor(rep(c("D0", "D4", "D7"), c(2,3,2)))
    data$samples$group <- group
    geneid <- rownames(cells)
    print(paste("geneid is", geneid))
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
    plotMDS(lcpm, labels=group, col=col.group)
  })

  sd <- 0.3*sqrt(4/rchisq(1000,df=4))
  x <- matrix(rnorm(1000*6,sd=sd),1000,6)
  rownames(x) <- paste("Gene",1:1000)
  x[1:50,4:6] <- x[1:50,4:6] + 2
  # without labels, indexes of samples are plotted.
  mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
  # or labels can be provided, here group indicators:
  # plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))

  # output$leftUSPlot <- renderPlot({  plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))})





  # output$leftUSPlot <- renderPlot(plot(1:10))
  output$rightUSPlot <- renderPlot(plot(10:1))
  # source("funcs.R", max=Inf)

  # library(limma)
  # system.time({plotData =  getOutput()})

  # pl <- plotData$plot1

  # output$plot1 = renderPlot({  barcodeplot(pl$x, index=pl$index, index2=pl$index2, main=pl$main) }) 

}

shinyApp(ui = ui, server = server)

