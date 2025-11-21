suppressPackageStartupMessages({
    library(shiny)
    library(ga.data)
  #  library(bettertrace)
    library(shinyjs)
})

shinyUI(fluidPage(

    useShinyjs(),
    ## Application title
    ## headerPanel("Plate Analyzer"),

    tags$head(
        tags$style("#diPlot {height:60vh !important;}"),
        tags$link( rel="stylesheet", type="text/css", href="extra-styles.css" )
    ),

    tabsetPanel(

        tabPanel(

            "Discovery Plotter",

            wellPanel(
                fluidRow(

                    column(
                        4,
                        fileInput("bc_file", label=h4("Plate Analyzer"), multiple=TRUE,
                                  accept="text/csv" )
                    ),
                    shinyjs::hidden(column(
                      4,
                      h4("Normalized data (with HYC01 and QCC29)"),
                      downloadButton("downloadResults","Download"),
                      id="download-button-column"
                    ))
                   
                )

            ),
           fluidRow(
                  column(12,
                         plotOutput("Discoveryplot", height="1000px")
                  )
                )

        ),
         
        tabPanel(
          "QC Overview",
          
          wellPanel(
            fluidRow(
              column( 4, div(tableOutput("qcTableST"),class="qcTable" ))
            ),
            fluidRow(
              column( 4, tableOutput("qcTableQCC29") )
            ),
            fluidRow(
              column( 4, tableOutput("qcTableQCC23") )
            ),
            fluidRow(
              column( 4, tableOutput("qcTableQCC33") )
            ),
            fluidRow(
              column( 4, tableOutput("qcTableZymo") )
            )
          )
          
          
        ),
        tabPanel(
          
          "PC plot",
          
          fluidRow(
            column(12,
                   plotOutput("PCplot", height="500px")
            )
          )
          
          
        ),
        tabPanel(

            "Probe Annotations",

            fluidRow(
                column(
                    12,
                    dataTableOutput("ProbeAnnotations")
                )
            )

        )
    )

    ## tags$script(HTML("console.log(output)"))


))
