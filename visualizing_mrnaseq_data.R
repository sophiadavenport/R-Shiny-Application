## Author: Sophia Davenport
## sophiada@bu.edu
## BU BF591
## Final Project 591
library(tidyr)
library(bslib)
library(ggplot2)
library(colourpicker) 
library(DT)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(ggbeeswarm)

ui <- fluidPage(
  titlePanel("Gene Expression in Brain Tissue for Huntington's Disease vs Normal Individuals"),
  theme = bslib::bs_theme(bootswatch = "sandstone", version = 5),
  div(
    style = "margin-top: 10px;",
    h6(HTML("Utilizing data from GEO DataSet: GSE64810"))
  ),
  
  sidebarLayout(
    sidebarPanel(
      style = "background-color: #f0f0f0; padding: 20px;",
      fileInput("fileInput", label = "Upload Sample File"),
      fileInput("CountsfileInput", label = "Upload Counts Matrix"),
      fileInput("DeseqfileInput", label = "Upload DeSeq2 Results"),
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
                 tabsetPanel(
                   tabPanel("Summary", HTML("<h5><b>Summary for Sample Characteristics</b></h5>"), tableOutput("summary_table")),
                   tabPanel("Table", HTML("<h5><b>Sample Information</b></h5>"), DTOutput("table")),
                   tabPanel("Plot", radioButtons("selected_column", "Choose the Sample Characteristic to Visualize", 
                                                 choices = c("PMI", "Age of Death", "RIN")), 
                            plotOutput("violin_plot"))
                 )
        ),
        tabPanel("Counts",
                 tabsetPanel(
                   tabPanel("Normalized Counts", value = "Normalized Counts", 
                            sliderInput("varianceSlider", "Variance Percentile", min = 0, max = 100, value = 50),
                            sliderInput("nonZeroSlider", "Non-Zero Samples", min = 0, max = 69, value = 10),
                            HTML("<h5><b>Summary of the effect of filtering by variance and non-zero samples:</b></h5>"),
                            tableOutput("summary_counts")
                   ),
                   tabPanel("Scatterplots", HTML("<h5><b>Diagnostic Scatterplots</b></h5>"), plotOutput("diagnosis_scatter_plot")),
                   tabPanel("Heatmap", value =  "Heatmap", HTML("<h5><b>Heatmap of Genes Remaining after Filtering</b></h5>"), 
                            plotOutput("heatmap")),
                   tabPanel("Principal Component Analysis", value =  "Principal Component Analysis", 
                            numericInput("x_axis", "Select X-axis component (PC1-69)", value = 1, min = 1, max = 69),
                            numericInput("y_axis", "Select Y-axis component (PC1-69)", value = 2, min = 1, max = 69), 
                            plotOutput("pca")
                   )
                 )
        ),
        tabPanel("Differential Gene Expression", tabsetPanel(
          tabPanel("Table of DeSeq2 Results", HTML("<h5><b>DeSeq2 Results</b></h5>"), DTOutput("deseq_table")),
          tabPanel("Volcano Plot", 
                   fluidRow(
                     column(4,
                            radioButtons("x", "Choose the column for the x-axis",
                                         choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))),
                     column(4,
                            radioButtons("y", "Choose the column for the y-axis",
                                         choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))),
                     column(4,
                            sliderInput("pAdjustMagnitude", "Select the magnitude of the p-adjusted coloring:",
                                        min = -50, max = 1, value = -2))),
                   plotOutput("deseq_volcano"))
        )),
        tabPanel("Visualize Gene Expression", 
                 HTML("<h5><b>Visualize Expression Data for Individual Genes</b></h5>"),
                 #plot of the selected type with the normalized gene counts for the selected gene split out by the categorical variable chosen
                 fluidRow(column(4, 
                                 selectizeInput("selectedGene", "Select Gene:", choices = NULL, multiple = FALSE)),
                 column(4, radioButtons("category", "Select Categorical Variable:", choices = c("Diagnosis"))),
                 column(4, radioButtons("plot_type", "Select Plot Type:", 
                              choices = c("Bar plot", "Boxplot", "Violin plot", "Beeswarm plot")))),
                 actionButton("submitButton", "Create Plot"),
                 #tableOutput("genetable"),
                 plotOutput("geneplot")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30*1024^2) #increased so that counts matrix can be uploaded
  
  load_data <- reactive({
    req(input$fileInput)
    filepath <- input$fileInput$datapath # getting path from input button
    # Attempt to read as CSV
    data <- tryCatch({
      read.csv(file = filepath)
    }, error = function(e_csv) {
      # If reading as CSV fails, attempt to read as TSV
      tryCatch({
        read.delim(file = filepath, sep = "\t")
      }, error = function(e_tsv) {
        return("Please input a CSV or TSV file.") # if cannot be read as csv or tsv, return an error message 
      })
    })
    # Check if data is not NULL
    if (!is.null(data)) {
      col_names <- c("Sample Name", "Sample Geo Accession", "Diagnosis", "PMI", "Age of Death", "RIN", "mRNA Seq Reads")
      names(data) <- col_names
      return(data)
    } else {
      return("Error: Unable to read the file.")
    }
  })
  #Raw sample information in a sortable table
  output$table <- renderDT({
    data <- load_data()
    datatable(data, options = list(ordering = TRUE))
  })
  #summary of sample information table
  summary_table <- reactive({
    data <- load_data()
    data <- data[, c("Diagnosis", "PMI", "Age of Death", "RIN", "mRNA Seq Reads")]
    data$Diagnosis <- factor(data$Diagnosis, levels = c("Neurologically normal", "Huntington's Disease"))
    if (!is.null(data)) {
      summary_data <- data.frame(
        Column_Name = names(data),
        Type = sapply(data, typeof),
        Mean_or_Distinct = sapply(data, function(x) {
          if (is.numeric(x)) {
            paste(round(mean(x, na.rm = TRUE), 2), " (", round(sd(x, na.rm = TRUE), 2), ")", sep = "")
          } else { 
            "Neurologically normal, Huntington's Disease"
          }
        })
      )
      col_names <- c("Column Name", "Type", "Mean (sd) or Distinct Values")
      names(summary_data) <- col_names
      return(summary_data)
    } else {
      return(NULL)
    }
  })
  #Output for summary of sample information table
  output$summary_table <- renderTable({
    summary_data <- summary_table()
    return(summary_data)
  })
  
  create_violin_plot <- function(data, selected_column) {
    if (is.null(selected_column)) {
      return(print("Please select an attribute to plot."))
    }
    
    plot <- ggplot(data, aes(x = data$Diagnosis, y = data[[selected_column]], fill = Diagnosis)) +
      geom_violin(alpha = 0.7) +
      labs(x = "Diagnosis", y = selected_column, title = paste("Violin Plot of", selected_column, "by Diagnosis")) +
      theme_minimal()
    
    return(plot)
  }
  #Output of violin plot for samples tab
  output$violin_plot <- renderPlot({
    data <- load_data()
    selected_column <- input$selected_column
    
    if (!is.null(data) && !is.null(selected_column)) { #ensuring data is input and column is selected
      create_violin_plot(data, selected_column)
    }
  })
  
  ###Start of Counts Tab Code:
  normalized_counts <- reactive({
    req(input$CountsfileInput)
    counts_data <- read.csv(input$CountsfileInput$datapath, header = TRUE) #Read the CSV file
    names(counts_data)[1] <- "Gene" #change name from X to gene
    number_samples <- ncol(counts_data)-1 #total number of samples in matrix
    number_genes <- nrow(counts_data) #total number of genes in matrix
    
    with_var <- counts_data %>%
      mutate(row_variance = apply(select(., -Gene), 1, var)) %>%
      mutate(variance_percentile = rank(row_variance) / length(row_variance) * 100)
    
    zero_counts <- rowSums(counts_data[, -1] == 0)
    with_var$zero_count <- zero_counts
    #include genes with at least X percentile of variance
    filtered_data <- with_var[with_var$variance_percentile >= input$varianceSlider, ]
    #include genes with at least X non-zero samples
    filtered_data <- filtered_data[number_samples-(filtered_data$zero_count) >= input$nonZeroSlider, ]
    
    number_passing <- nrow(filtered_data) #number of genes passing the current filter
    number_not_passing <- (number_genes - number_passing) #number of genes NOT passing the current filter
    percent_passing <- (number_passing/number_genes)*100 #percent of genes passing filter
    percent_not_passing <- (number_not_passing/number_genes)*100 #percent of genes NOT passing filter
    
    #list so items can be called
    result <- list(number_genes = number_genes, number_samples = number_samples, number_passing = number_passing,
                   number_not_passing = number_not_passing, percent_passing = percent_passing, 
                   percent_not_passing = percent_not_passing, filtered_data = filtered_data, counts_data = counts_data)
    
    return(result)
  })
  
  output$summary_counts <- renderTable({
    result <- normalized_counts()
    summary_data <- data.frame(
      "Attribute" = c("Total Number of Genes", "Total Number of Samples", "Number of Genes Passing Filter", 
                      "Number of Genes Not Passing Filter", "Percent of Genes Passing Filter", 
                      "Percent of Genes Not Passing Filter"),
      "Value" = c(result$number_genes, result$number_samples, result$number_passing,
                  result$number_not_passing, paste(round(result$percent_passing, 2), "%"),
                  paste(round(result$percent_not_passing, 2), "%"))
    )
    return(summary_data)
  })
  #Create scatter plot for filtered counts
  create_diagnosis_scatter_plot <- function(data) {
    req(input$CountsfileInput)
    result <- normalized_counts()
    data <- result$counts_data
    number_samples <- ncol(data)-1 #total number of samples in matrix
    
    gene_stats <- data %>%
      mutate(median_count = apply(select(., -Gene), 1, function(x) median(x, na.rm = TRUE)),# Gene median count
             variance = apply(select(., -Gene), 1, var), #Gene variance
             num_zeros = apply(select(., -Gene), 1, function(x) sum(x == 0, na.rm = TRUE))) #Number of Zero counts per gene
    
    #Median Count versus Variance
    plot1 <- ggplot(gene_stats, aes(x = log(median_count), y = log(variance), 
                                    color = ifelse(variance >= input$varianceSlider & (number_samples-num_zeros) >= input$nonZeroSlider, "Pass", "Fail"))) +
      geom_point() +
      labs(x = "Log Median Count", y = "Log Variance", title = "Gene Median Count versus Variance") +
      scale_color_manual(name = "Gene Status", values = c("Pass" = "darkgreen", "Fail" = "lightgrey")) +
      theme_minimal()
    
    #Median Count versus Number of Zeros
    plot2 <- ggplot(gene_stats, aes(x = median_count, y = num_zeros, 
                                    color = ifelse(variance >= input$varianceSlider & (number_samples-num_zeros) >= input$nonZeroSlider, "Pass", "Fail"))) +
      geom_point() +
      labs(x = "Median Count", y = "Number of Zeros", title = "Gene Median Count versus Number of Zeros") +
      scale_color_manual(name = "Gene Status", values = c("Pass" = "darkgreen", "Fail" = "lightgrey")) +
      theme_minimal()

    combined_plot <- cowplot::plot_grid(plot1, plot2, ncol = 1, align = "v") #plots stacked vertically
    return(combined_plot)
  }
  #Render diagnosis scatter plot
  output$diagnosis_scatter_plot <- renderPlot({
    result <- normalized_counts()
    data <- result$counts_data
    
    create_diagnosis_scatter_plot(data)
    })
  
  #Render heatmap
  output$heatmap <- renderPlot({
    result <- normalized_counts()
    data <- result$filtered_data #ensuring data is user with slider
    data <- data %>% select(-c("row_variance", "variance_percentile", "zero_count")) #removing non-counts rows
    data <- data %>% select_if(is.numeric)
    data_matrix <- as.matrix(data) #heatmap must have matrix data
    
    colors <- colorRampPalette(c("white", "blue", "red"))(50)
    heatmap_plot <- heatmap.2(data_matrix, col = colors, trace = "none")
    
    return(heatmap_plot)
  })
  
  #PCA for Counts Matrix
  create_pca <- function(data) {
    numeric_data <- data[, sapply(data, is.numeric)] #making data numeric
    expr_mat <- as.matrix(numeric_data) #pca plot requires numeric data
    rownames(expr_mat) <- data$"Gene"
    expr_mat_centered <- scale(expr_mat)
    pca <- prcomp(expr_mat_centered, center = FALSE, scale = TRUE)
    return(pca)
    }
  #Rendering PCA
  output$pca <- renderPlot({
    req(input$fileInput)
    metadata <- summary_table()
    groups <- metadata$Diagnosis #getting diagnosis to color points by
    result <- normalized_counts() 
    result <- result$counts_data #getting counts data (not filtered)
    pca <- create_pca(result) #creating pca object (will have PC1-69)
    
    x_pc <- input$x_axis #UI x-axis input
    x_data <- pca$x[, x_pc] #PC data for x
    y_pc <- input$y_axis #UI y-axis input
    y_data <- pca$x[, y_pc] #PC data for y
    pc_data <- as.data.frame(pca$x[, c(x_pc, y_pc)]) #merging into one data frame for plotting
    pc_data$Diagnosis <- groups
    print(pc_data)
    percent_variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100 #percent variance for all PCs
    percent_variance_x <- round(percent_variance[x_pc], 2) #variance for PC on x
    percent_variance_y <- round(percent_variance[y_pc], 2) #variance for PC on y
    #Creating plot
    ggplot(pc_data, aes(x = pc_data[, 1], y = pc_data[, 2])) +
      geom_point() +
      scale_color_manual(name = "Diagnosis", values = c("Neurologically normal" = "blue", "Huntington's Disease" = "red")) +
      labs(x = paste0("PC", x_pc, " (", percent_variance_x , "%)")
           , y = paste0("PC", y_pc, " (", percent_variance_y , "%)")
           , title = "PCA Plot")
  })
  #Getting DeSeq2 Results Matrix in a Filterable Table
  deseq_data <- reactive({
    req(input$DeseqfileInput)
    filepath <- input$DeseqfileInput$datapath # getting path from input button
    # Attempt to read as CSV
    data <- tryCatch({
      read.csv(file = filepath)
    }, error = function(e_csv) {
      # If reading as CSV fails, attempt to read as TSV
      tryCatch({
        read.delim(file = filepath, sep = "\t")
      }, error = function(e_tsv) {
        return("Please input a CSV or TSV file.") # if cannot be read as csv or tsv, return an error message 
      })
    })
    return(data)
  })
  #Output for deseq data
  output$deseq_table <- renderDT({
    data <- deseq_data()
    datatable(data, options = list(ordering = TRUE))
  })
  #Volcano Plot Function
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1 = "darkgreen", color2 = "lightgrey") {
      plot <- ggplot(dataf, aes(x = !!sym(x_name), y = !!sym(y_name), 
                                color = !!sym("padj") < 1*10^(slider))) +
        geom_point(size = 3, alpha = 0.7) +  
        labs(x = x_name,y = y_name) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        scale_color_manual(name = paste("p-adjusted", "< 1*10^",slider),
                           values = c("TRUE" = color1, "FALSE" = color2))
      
      return(plot)
    }
  #Rendering Volcano Plot for Deseq2 results
  output$deseq_volcano <- renderPlot({
    req(input$DeseqfileInput)
    data <- deseq_data()
    x_input <- input$x
    y_input <- input$y
    p <- volcano_plot(data, x_input, y_input, input$pAdjustMagnitude)
    print(p)
  })
  
  #Visualizing Gene Expression Tab:
  #Getting genes for search bar
  observe({
    x <- normalized_counts()
    data <- x$counts_data
    data_genes <- data$Gene #Getting column of gene symbols
    updateSelectizeInput(session, "selectedGene", choices = data_genes, server = TRUE)
  })
  
  #Getting data for geneplot:
  gene_data <- function(gene) {
    x <- normalized_counts()
    data <- as.data.frame(x$counts_data)
    summary_data <- as.data.frame(load_data()) %>% select("Sample Name", "Diagnosis")
    #category_data <- summary_data[summary_data$Diagnosis == sample_category, ] #selecting for category
    data <- t(data) #transposing data
    data <- as.data.frame(cbind(Sample = rownames(data), data))
    colnames(data) <- data[1,] #genes to column names for counts data
    merged <- left_join(summary_data, data, by = c("Sample Name" = "Gene")) #left join selects only category samples
    
    selected <- merged[, c("Sample Name", "Diagnosis", gene)]
    return(selected)
  }
  #Getting gene information to print above plots
  #gene_info <- function(gene){
    #data <- deseq_data()
    #selected <- data[data$X == gene_id, c("Symbol", "padj", "log2FoldChange")]
    #info_df <- data.frame(
      #"Gene Symbol"= gene_data$Symbol,
      #"Adjusted P-value" = gene_data$padj,
      #"Log2 Fold Change" = gene_data$log2FoldChange
    #)
  #}
  
  #output$genetable <- renderTable({
    #req(input$selectedGene)
    #gene <- input$selectedGene
    #df <- gene_info(gene)
    #return(df)
  #})
  
  observeEvent(input$submitButton, {
    selected_gene <- input$selectedGene
    category <- input$category
    type <- input$plot_type
    
    data <- gene_data(selected_gene)
    colnames(data)[3] = "Counts"
    data$Counts <- as.numeric(data$Counts) #setting to numeric values
    
    if (type == "Bar plot") {
      output$geneplot <- renderPlot({
        bar_plot <- ggplot(data, aes(x = reorder(data$"Sample Name", Counts), y = Counts, fill = Diagnosis)) +
          geom_bar(stat = "identity") +
          labs(title = paste("Gene Expression for ", selected_gene), x = "Samples", y = "Normalized Counts") +
          theme(axis.text.x = element_text(color = "transparent"),
                axis.ticks.x=element_blank(),
                legend.position = "right",
                legend.title = element_blank(),
                plot.margin = margin(.5, .5, .5, .5, "cm"))
        return(bar_plot)
      }, height = 400, width = 850)
    }
    if (type == "Boxplot") {
      output$geneplot <- renderPlot({
        box_plot <- ggplot(data, aes(x = Diagnosis, y = Counts)) +
          geom_boxplot(fill = "skyblue") +
          labs(title = paste("Gene Expression for ", selected_gene), x = "Samples", y = "Normalized Counts")
        return(box_plot)
      })
    }
    if (type == "Violin plot") {
      output$geneplot <- renderPlot({
        violin <- ggplot(data, aes(x = Diagnosis, y = Counts, fill = Diagnosis)) +
          geom_violin(alpha = 0.7) +
          labs(x = "Samples", y = "Normalized Counts", title = paste("Gene Expression for ", selected_gene)) +
          theme_minimal()
        return(violin)
      })
    }
    if (type == "Beeswarm plot") {
      output$geneplot <- renderPlot({
        beeswarm_plot <- ggplot(data, aes(x = Diagnosis, y = Counts)) +
          geom_beeswarm(fill = "lavender") +
          labs(title = paste("Gene Expression for ", selected_gene), x = "Samples", y = "Normalized Counts") +
          theme_minimal()
        return(beeswarm_plot)
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)