library(shiny)
library(stringr)

# User Interface (UI) definition
ui <- fluidPage(
  titlePanel("Species Identification"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("ref_file", "Reference file from DB"),
      fileInput("query_file", "Query file"),
      numericInput("num_peaks", "Number of top peaks selected per spectrum:", value = 20),
      numericInput("num_hits", "Number of hits presented:", value = 2),
      numericInput("tolerance", "Tolerance(kDa)", value = 0.5),
      actionButton("run_button", "Execute")
    ),
    
    mainPanel(
      plotOutput("plot"),
      plotOutput("plot1"),
      tableOutput("result")
    )
  )
)

# Server definition
server <- function(input, output) {
  source('~/Documents/Project_stage/mine/app_1/Identification_project.R')
  
  observeEvent(input$run_button, {
    # Check if files are selected
    if (is.null(input$ref_file) || is.null(input$query_file)) {
      return()  # Stop execution if files are not selected
    }
    
    # Get the paths of the selected files
    ref_file_path <- input$ref_file$datapath
    query_file_path <- input$query_file$datapath
    
    # Execute the species_identification function with user-provided parameters
    data <- species_identification(ref_file_path, query_file_path, input$num_peaks, input$num_hits, input$tolerance)
    output$result <- renderTable(data)
    
    # Show a notification message once the function has finished executing
    showNotification("Execution completed!", duration = 10)
    
    hits_table <- as.data.frame(data)
    
    # Extract specific parts from species names
    hits_table$Species_short <- str_extract(hits_table$query, "(?<=Culex_)[^_]+")
    hits_table$Species_found <- str_extract(hits_table[, 3], "(?<=:)[^[:space:]]+")
    
    # Count the occurrences of each species
    species_count_short <- table(hits_table$Species_short)
    species_count_found <- table(hits_table$Species_found)
    
    output$plot <- renderPlot({
      # Plot the number of specimens per species in the query file
      pie(species_count_short, labels = paste(names(species_count_short), " (", species_count_short, ")", sep = ""),
          col = rainbow(length(species_count_short)), cex = 0.8, main = "Number of specimens per species in query file", init.angle = 45)
    })
    
    output$plot1 <- renderPlot({
      # Plot the number of specimens per species found
      pie(species_count_found, labels = paste(names(species_count_found), " (", species_count_found, ")", sep = ""),
          col = rainbow(length(species_count_found)), cex = 0.8, main = "Number of specimens per species found", init.angle = 45)
    })
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
