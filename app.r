library(shiny)
library(ggplot2)
library(purrr)
library(pbmcapply)
library(umap)

###################################################################################################
############# User interface
###################################################################################################
ui <- fluidPage(
  titlePanel("Embedding explorer"),

  sidebarLayout(

    sidebarPanel(
        radioButtons("input",
          label = h4("Kind of input"),
          choices = list("Raw input" = "raw",
                         "Embedding + Expression" = "embedding")
        ),
        conditionalPanel(condition = "input.input == 'raw'",
          fileInput("rawfile", label = h3("Raw UMI counts")),
          p("CSV file: rows x columns = genes x cells"),
          p("First column = geneIDs, first row = cellIDs"),
          numericInput("min_cells", label = h3("Number of cells to accept gene"), value = 20),
          br(),
          selectInput("normMethod", "Select normalization method",
            choices = c("Anscombe", "Log norm"), selected = "Anscombe"
          )
        ),
        conditionalPanel(condition = "input.input == 'embedding'",
          fileInput("embedfile", label = h3("Embedding table")),
          p("The embedding should be a csv file with rows x colums: cells x dimensions."),
          p("Make sure that the cellIDs are the first column (row names) and no column names should be present."),
          br(),
          fileInput("expfile", label = h3("Expression table")),
          p("The expression table should be a CSV file: geneIDs x cellIDs"),
          p("First column: GeneIDs, first row: CellIDs")),
        br(),
        br(),
        h4("Server status"),
        conditionalPanel(
            condition="($('html').hasClass('shiny-busy'))",
            h5("R is busy"))
        ),


    mainPanel(
      # tableOutput("value"),
      # br(),
      # verbatimTextOutput("value"),
      # br(),
      tabsetPanel(
          tabPanel("Dimensional Reduction",
            p("If raw UMI counts are provided, use this tab to do the dimensional reduction"),
            uiOutput("QCtab")
          ),
          tabPanel("Embedding",
            plotOutput("embeddingPlot", brush = brushOpts(id = "brush", delayType = "debounce", delay = "1000"),
                       height = 650),
            br(),
            p("Draw a rectangular selection around a cluster of cells, you are interested in."),
            p("The app will calculate a Kruskal-Wallis-Wilcox test for all genes expressed within this cluster."),
            p("Click on the \"Gene enrichment analysis\" tab and wait until the calculation finished."),
            br(),
            h5("Number of cells selected: "),
            textOutput("brushed")
          ),
          tabPanel("Expressed genes",
            br(),
            p("Use this drop-down menu to select a gene. You can also click inside the menu and hit once delete, then you can enter a gene yourself."),
            uiOutput("selectGenes"),
            br(),
            plotOutput("exprsPlot", height = 650)
          ),
          tabPanel("Gene enrichment analysis",
            tabsetPanel(
              tabPanel("Results table",
                h4("Gene enrichment analysis results"),
                p("Select a gene below. It's expression will be plottet in the tab \"Gene enrichment plot\"."),
                h5("Attention"),
                p("If there is no table below this text, select cells in the \"Embedding\t tab and _wait_ here. The calculation can take a couple of seconds to minutes."),
                DT::dataTableOutput("test")
              ),
              tabPanel("Gene enrichment plot",
                plotOutput("gePlot", height = 650)
            )
              

            )

          )

        # )
        ),
      br(),
      h4("Debugger"),
      verbatimTextOutput("debugger")
      
      



      )

  )

)


###################################################################################################
############# Server
###################################################################################################


# Define server logic 
server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)

  #####################################
  # Functions
  #####################################


  # Function for doing the anscombe transformation
  anscombe_transform = function(counts, sf){
    apply(counts[, ], 2, function(umi) {
    ans <- sqrt(umi+3/8) - sqrt(3/8)
    ans/sqrt(sum(ans^2)) })
  }

  # Function for calculating the size factor per cell
  calculate_size_factor = function(counts){
    Matrix::colSums(counts) / mean(Matrix::colSums(counts))
  }

  disp <- function(k) (var(k) - mean(k)) / mean(k) / mean(k) # var = mu + disp * mu^2


  fitNB <- function(k, sf, initialVals = c(round(mean(k), 10), round(disp(k), 10))) { 
      k = as.numeric(k)
      o = optim( 
        initialVals,  # mu, alpha
        function(x) {
          -sum( lgamma( k + 1/x[2] ) - lgamma( 1/x[2] ) - lgamma( k+1 ) - ( k + 1/x[2] ) * log( 1 + sf * x[2] * x[1] ) + k * log( sf * x[2] * x[1] ) )
          },
        
        function(x) c(
          -sum( ( k - sf * x[1] ) / ( x[1] + sf * x[2] * x[1]^2 ) ),
          -sum( ( x[2] * ( k - sf * x[1] ) / ( 1 + sf * x[2] * x[1] ) + log( 1 + sf * x[2] * x[1] ) - digamma( k + 1/x[2] ) + digamma( 1/x[2] ) ) / x[2]^2 ) ),
        hessian = TRUE,
        lower = c( 1e-10, 1e-10 ),
        method = "L-BFGS-B" )
      c( mean = o$par[1], 
         disp = o$par[2], 
         SE_m = 1 / sqrt(o$hessian[1,1]), 
         SE_d = 1 / o$hessian[2,2] ) }


  library(purrr) # for possibly
  safe_fitNB <- possibly(fitNB, otherwise = c(mean=NA, disp=NA, SE_m = NA, SE_d = NA))


  stack_disp_table = function(dispList){
      dtable = NULL
      for(i in 1:length(dispList)){
          tmp = dispList[[i]]
          dtable = rbind(dtable, tmp)
      }
      
      return(dtable)
  }

  calculate_disp_table = function(counts, sf){
    require(pbmcapply)


    disp.table = pbmclapply(rownames(counts), function(g){
                            NBparams = safe_fitNB(counts[g, ], sf)
                            data.frame(gene = g, t(NBparams))
                            })

    return(stack_disp_table(disp.table))
  }

  inspect_disptable = function(disptable){
    signif_dispersion = (disptable$disp / disptable$SE_d) > 2

    noNA = colSums(apply(disptable, 1, is.na)) == 0

    disptable$noNA = noNA
    disptable$signif_dispersion = signif_dispersion


    disp_plot = ggplot(disptable, aes(x = mean, y = disp)) +
                geom_point(size = 0.25, alpha = 0.2, aes(color = noNA & signif_dispersion))+
                scale_x_log10() +
                 scale_y_log10() +
                annotation_logticks(base = 10, size = 0.2)

    return(list("signif_dispersion" = signif_dispersion,
                "noNA" = noNA,
                "plot" = disp_plot))

  }



  #####################################
  # Reading in the data
  #####################################

  embedding = reactive({
      if(input$input == "embedding"){
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
            )[,1:2]
        }
        
      }else if(input$input == "raw"){
        df = own.embedding()[,1:2]
      }

      colnames(df) = c("Dim1", "Dim2")
      return(as.data.frame(df))



    })



  exprs = reactive({

    if(input$input == "embedding"){
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
      }else if(input$input == "raw"){
        return(as.matrix(normed.exprs()))
      }
      



    })

  # Reading in the raw UMI input,
  # using the input$rawfile$datapath variable.
  raw.umi.input = reactive({
      infile = input$rawfile

      if(is.null(infile)){
        return(NULL)
      }else{
        df = read.table(infile$datapath,
          header = T,
          row.names = 1,
          sep = ",")
      }

      return(as.matrix(df))

    })


  #####################################
  # QC
  #####################################


  # Using only UMI of genes which are expresed in more than
  # min_cells cells.
  raw.umi = reactive({

      if(is.null(input$rawfile)){
        return(NULL)
      }

      df = raw.umi.input()

      return(df[as.vector(Matrix::rowSums(df > 0)) >= input$min_cells, ])

    })

  # Calculating the number of genes expressed per cell
  n_genes = reactive({

      if(is.null(raw.umi())){
        return(NULL)
      }
      umi = raw.umi()

      expressed.genes = umi > 0
      expressed.genes = Matrix::colSums(expressed.genes)

      return(as.vector(expressed.genes))
    })

  # Calculating the number of UMI per cell
  n_umi = reactive({

      if(is.null(raw.umi())){
        return(NULL)
      }
      umi = raw.umi()

      umi.counts.per.cell = Matrix::colSums(umi)

      return(as.vector(umi.counts.per.cell))
    })


  # Rendering the dimensional reduction tabsetPanels
  output$QCtab = renderUI({
      if(is.null(raw.umi())){
        return(NULL)
      }else{
        tabsetPanel(
          tabPanel("QC",
            fluidRow(
                   column(4, uiOutput("selectMinMaxUMI")),
                   column(4, uiOutput("selectMinMaxGenes"))
                  ),
            fluidRow(column(8, plotOutput("n_genes_vs_n_umi", height = 400))),
            br(),
            fluidRow( column(4, plotOutput("violin_nUMI")),
                      column(4, plotOutput("violin_nGenes"))
            )
          ),
          tabPanel("Gene correlation",

            fluidRow(
              column(4, uiOutput("selectCorrGene1")),
              column(4, uiOutput("selectCorrGene2"))
            ),
            fluidRow(
              column(8, plotOutput("geneCorrPlot"))
            )
          ),
          tabPanel("NB fitting",
            verbatimTextOutput("head_disp_table"),
            br(),
            fluidRow(
              column(8, plotOutput("fitting_plot"))
            )
          ),
          tabPanel("Dim. red. plots",
            fluidRow(column(4, selectInput("whichDimRed", "What dim. reduction to perform?", choices = c("PCA", "UMAP"), selected = "PCA"))),
            fluidRow(column(8, plotOutput("dim_red_plot")))
          )
        
        )
        
      }
    })

  # Slider input for min and max number of UMI
  output$selectMinMaxUMI = renderUI({
      if(is.null(n_umi())){
        return(NULL)
      }else{
        sliderInput("MinMaxUMI", label = h3("UMI range"), min = min(n_umi()), max = max(n_umi()), value = c(min(n_umi()), max(n_umi())))
      }
      
    })

    # Slider input for min and max number of genes
  output$selectMinMaxGenes = renderUI({
      if(is.null(n_umi())){
        return(NULL)
      }else{
        sliderInput("MinMaxGenes", label = h3("Gene range"), min = min(n_genes()), max = max(n_genes()), value = c(min(n_genes()), max(n_genes())))
      }
      
  })

    # Find cells which fall into the criteria
    # selectMinMaxUMI and selectMinMaxGenes
  selected_cells = reactive({
        return(n_umi() >= input$MinMaxUMI[1] & n_umi() <= input$MinMaxUMI[2] & 
                         n_genes() >= input$MinMaxGenes[1] & n_genes() <= input$MinMaxGenes[2])
  })

    # Generate a color vector for selected cells
    # black cells are inside the selection window,
    # grey cells are outside.
  selected_cells_color = reactive({
        ifelse(selected_cells(),
                         "black", "gray")
  })


  # Scatterplot, showing n_genes vs n_umi
  # Window is selectable using the selectMinMaxUMI and
  # selectMinMaxGenes sliders.
  output$n_genes_vs_n_umi = renderPlot({
    if(is.null(raw.umi)){
      return(NULL)
    }else{

      colBool = selected_cells_color()

      return(
      ggplot(NULL, aes(x = n_umi(), y = n_genes())) +
        geom_point(color = colBool) +
        xlab("nUMI") +
        ylab("nGenes") +
        geom_vline(aes(xintercept = input$MinMaxUMI[1])) +
        geom_vline(aes(xintercept = input$MinMaxUMI[2])) +
        geom_hline(aes(yintercept = input$MinMaxGenes[1])) +
        geom_hline(aes(yintercept = input$MinMaxGenes[2]))
      )
    }
    
  })


  # Violinplot for all cells
  # passing the filter criteria showing the number
  # of UMI per cell
  output$violin_nUMI = renderPlot({
      ggplot(NULL, aes(y = n_umi()[selected_cells()], x = "nUMI")) +
        geom_violin() +
        ylab("UMI counts")
    })

  # As above, but showing the number of genes in the window
  output$violin_nGenes = renderPlot({
      ggplot(NULL, aes(y = n_genes()[selected_cells()], x = "nGenes")) +
        geom_violin()
    })

  normed.exprs = reactive({
      norm.choice = input$normMethod

      if(norm.choice == "Anscombe"){
          umi.counts = raw.umi()[ , selected_cells()]

          sf = calculate_size_factor(umi.counts)

          normed = anscombe_transform(umi.counts, sf)

          return(normed)
        }else if(norm.choice == "Log norm"){
          return(NULL)  ## Here I will perform the normal Seurat approach later. TODO
        }
    })
    
  ######
  # Making the gene correlation plot
  ######
  output$selectCorrGene1 = renderUI({
      selectInput("corrGene1", "Select gene A", choices = row.names(normed.exprs()), selected = rownames(normed.exprs())[1])
    })

  output$selectCorrGene2 = renderUI({
      selectInput("corrGene2", "Select gene B", choices = row.names(normed.exprs()), selected = rownames(normed.exprs())[1])
    })

  output$geneCorrPlot = renderPlot({
      ggplot(NULL, aes(x = as.vector(normed.exprs()[input$corrGene1, ]),
                       y = as.vector(normed.exprs()[input$corrGene2, ]) ) 
             ) +
        geom_point() +
        xlab(input$corrGene1) +
        ylab(input$corrGene2)
    })

  ######
  # Negative binomial fitting
  ######

  initial.nb.fit = reactive({
      umi.counts = raw.umi()[ , selected_cells()]


      calculate_disp_table(umi.counts, calculate_size_factor(umi.counts))

    })


  output$head_disp_table = renderPrint({
      head(initial.nb.fit())
    })

  significant.dispersion = reactive({
      disp.table = initial.nb.fit()

      noNA = colSums(apply(disp.table, 1, is.na)) == 0

      return(((disp.table$disp / disp.table$SE_d) > 2) & noNA)
    })

  output$fitting_plot = renderPlot({
      disp.table = initial.nb.fit()
      disp.table$signif_disp = significant.dispersion()

      ggplot(disp.table, aes(x = mean, y = disp)) +
        geom_point(size = 0.25, alpha = 0.2, aes(color = signif_disp)) +
        scale_x_log10() +
        scale_y_log10()

    })

  normed.exprs.signif = reactive({
      n.exprs = normed.exprs()

      n.exprs = n.exprs[significant.dispersion(), ]

      return(n.exprs)
    })

  pca.embedding = reactive({
      n.exprs = as.matrix(normed.exprs.signif())

      pca.embedding = prcomp(t(n.exprs))$x

      return(pca.embedding)
    })

  umap.embedding = reactive({
      pca.em = pca.embedding()

      umap.em = umap(pca.em[,1:10])

      return(umap.em$layout)
    })

  own.embedding = reactive({
      if(input$whichDimRed == "PCA"){
        return(pca.embedding())
      }else if(input$whichDimRed == "UMAP"){
        return(umap.embedding())
      }
    })

  # Plot the dimensional reduced data output
  # Deceide what to plot using a drop-down menu (whichDimRed)
  output$dim_red_plot = renderPlot({

    em = own.embedding()

    ggplot(NULL, aes(x = em[,1],
                       y = em[,2])) +
        geom_point(size = 0.25, alpha = 0.2) +
        xlab("Dim1") +
        ylab("Dim2")
    })


  ##########
  # Debugger print window
  ##########
  output$debugger = renderPrint({
      # if(!is.null(normed.exprs())){
      #   print(
      #       head(
      #             as.vector(normed.exprs()[input$corrGene2, ])
      #         ))
      # }

      return(NULL)
    })


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

    if(is.null(embedding()) | is.null(exprs())){
      return(NULL)
    }else{
      return(DT::datatable(wilcox.out(),
                         options = list(pageLength = 10),
                         selection = 'single'))
    }




  })

  wilcox.out = reactive({

      exprs.selected = exprs()[,selected.cells.bool()]
      exprs.selected = exprs.selected[Matrix::rowSums(exprs.selected > 0) > 0, ]



      exprs.background = exprs()[,!selected.cells.bool()]


      tested.gene = NULL
      p.value = NULL

      for(gene in row.names(exprs.selected)){

        ex.selected = exprs.selected[gene, ]
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


  output$gePlot = renderPlot({

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

  ####
  # Expression plot for selected genes
  ####

  # Creating a drop-down menu for the genes:
  output$selectGenes = renderUI({

      if(is.null(exprs())){
        h4("Please upload an expression table!")
      }else{
        genes.in.dataset = row.names(exprs())

        selectInput("selectGenes", "Select a gene", choices = genes.in.dataset, selected = genes.in.dataset[1])
      }
      
    })

  output$exprsPlot = renderPlot({

    if(is.null(input$selectGenes)){
      return(NULL)
    }else{
      em = embedding()
      em$geneExprs = exprs()[input$selectGenes, ]

      ggplot(em, aes(x = Dim1, y = Dim2, color = geneExprs)) +
        geom_point() +
        scale_color_continuous(high = "red", low = "gray", name = input$selectGenes)
    }

  })



}






###################################################################################################
############# Running the app
###################################################################################################

shinyApp(ui = ui, server = server)
