# This is just a dummy hello world shiny app at the moment.

getCounts <- function() {
  library(readr)
  counts <- read_delim("2019.07.12.counts.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
  # # rownames(counts) <- counts$GeneSymbol

}


# getPlot1 <- function() {
#   library(limma)
#   # return(barcodeplot(1:10, index=1))
#   return (list(x=1:10, index=1))

# }

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("RNAeasy123"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of bins ----
    #   sliderInput(inputId = "bins",
    #               label = "Number of bins:",
    #               min = 1,
    #               max = 50,
    #               value = 30)

    ),

    # Main panel for displaying outputs ----
    mainPanel(
      h2("Counts"),
      DT::dataTableOutput("mytable"),
      plotOutput("plot1")


    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$mytable = DT::renderDataTable({getCounts()})
  source("funcs.R", max=Inf)

  library(limma)
  system.time({plotData =  getOutput()})

  pl <- plotData$plot1

  output$plot1 = renderPlot({  barcodeplot(pl$x, index=pl$index, index2=pl$index2, main=pl$main) }) 
  # output$plot1 = renderPlot({barcodeplot(1:10, index=1) })
  # output$plot1 = renderPlot({ plot(1:10)})

}

shinyApp(ui = ui, server = server)

