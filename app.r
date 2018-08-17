library(shiny)
library(ggplot2)

###################################################################################################
############# User interface
###################################################################################################
ui <- fluidPage(
  titlePanel("Custom cluster analysis"),

  sidebarLayout(

    sidebarPanel(
        fileInput("embedfile", label = h3("Embedding table")),
        br(),
        fileInput("expfile", label = h3("Expression table"))

        ),


    mainPanel(
      # tableOutput("value"),
      # br(),
      # verbatimTextOutput("value"),
      # br(),
      tabsetPanel(
          tabPanel("Embedding",
            plotOutput("embeddingPlot", brush = brushOpts(id = "brush", delayType = "debounce", delay = "1000"),
                       height = 650)
          ),
          tabPanel("Expression",
            plotOutput("exprsPlot", height = 650)

          )
        ),
      br(),
      p("Number of cells selected: "),
      textOutput("brushed"),
      br(),
      DT::dataTableOutput("test"),
      br(),
      p("Select a gene above. It's expression will be plotted below:"),
      br()



      )

  )

)

###################################################################################################
############# Server
###################################################################################################


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  options(shiny.maxRequestSize=500*1024^2)

  ####
  # Reading in the data
  ####

  embedding = reactive({
      infile = input$embedfile

      if(is.null(infile)){
        if(file.exists("dummy_data/umap.csv")){
          df = read.table(
            "dummy_data/umap.csv",
            header = F,
            row.names = 1,
            sep = ","
          )
          colnames(df) = c("Dim1", "Dim2")
          return(df)
        }else{
          return(NULL)
        }

      }else{
        df = read.table(
            infile$datapath,
            header = F,
            row.names = 1,
            sep = ","
          )

        colnames(df) = c("Dim1", "Dim2")

        return(df)
      }




    })



  exprs = reactive({
      infile = input$expfile

      if(is.null(infile)){
        if(file.exists("dummy_data/anscombe.rds")){
          df = readRDS("dummy_data/anscombe.rds")
          return(df)
        }else{
          return(NULL)
        }

      }else{

        df = read.table(
            infile$datapath,
            header = T,
            row.names = 1,
            sep = ","
            )

        return(as.matrix(df))

      }



    })

  output$value = renderPrint({
      exprs.selected = exprs()[,selected.cells.bool()]
      exprs.selected = exprs.selected[rowSums(exprs.selected > 0) > 0, ]

      return(head(exprs.selected))
    })

  ####
  # Testing whether the reading works
  ####




  ####
  # Making the embedding plot
  ####


  output$embeddingPlot = renderPlot({
      em = embedding()

      if(is.null(em)){
        return(NULL)
      }else{
        return(

            ggplot(em, aes(x = Dim1,
                 y = Dim2)) +
              geom_point() +
              ylab("Dim2") +
              xlab("Dim1")

          )
      }



    })

  ####
  # Finding out which cells are selected on the plot
  ####

  sel.cells = reactive({
      em = embedding()

      brushedPoints(em, input$brush, allRows = TRUE)

    })

  output$brushed = renderText({

      sum(selected.cells.bool())

    })

  selected.cells.bool = reactive({
      sel.cells()$selected_
    })


  # Returns the gene with the highest summed expression in the selection
  output$test = DT::renderDataTable({

    if(is.null(embedding())){
      return(NULL)
    }else{
      return(DT::datatable(wilcox.out(),
                         options = list(pageLength = 10),
                         selection = 'single'))
    }




  })

  wilcox.out = reactive({

      exprs.selected = exprs()[,selected.cells.bool()]
      exprs.selected = exprs.selected[rowSums(exprs.selected > 0) > 0, ]



      exprs.background = exprs()[,!selected.cells.bool()]


      tested.gene = NULL
      p.value = NULL

      for(gene in row.names(exprs.selected)){

        ex.selected = exprs.selected[gene, ]
        print(ex.selected)
        ex.background = exprs.background[gene, ]

        test.out = wilcox.test(x = ex.selected , y = ex.background, alternative = "greater", paired = FALSE)

        tested.gene = c(tested.gene, gene)
        p.value = c(p.value, test.out$p.value)

      }

      out.df = data.frame(gene = tested.gene,
                          p.value = p.value)

      out.df$p.adjusted = p.adjust(out.df$p.value, method = "BH")

      return(out.df)

    })


  output$exprsPlot = renderPlot({

      s = input$test_rows_selected

      if(length(s) == 1){

          exprs.gene = wilcox.out()[s,1]

          exprs.vec = exprs()[row.names(exprs()) == exprs.gene, ]
          em = embedding()
          em$expression = exprs.vec

          ggplot(em, aes(x = Dim1,
                                  y = Dim2,
                                  color = expression)) +
              geom_point() +
              ylab("Dim2") +
              xlab("Dim1") +
              scale_color_continuous(low = "gray", high = "red", name = exprs.gene)


        }else{
          return(NULL)
        }

    })

}



###################################################################################################
############# Running the app
###################################################################################################

shinyApp(ui = ui, server = server)
