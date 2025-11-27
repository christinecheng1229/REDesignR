# TODO Source: https://shiny.posit.co/r/gallery/application-layout/tabsets/
# NOTE: try to refrain from using single letter variable names (e.g., d, n) for input/output/buttons etc.

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Simulates, Optimizes, and Visualizes Restriction Enzyme Digestion"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # TODO detailed description of what the tool does and how the user can interact/use it
      tags$p("This Shiny application is dedicated to [...]. It can be used to calculate [...]."),
      tags$p("Select appropriate input and press 'Run'."),

      # Input: Select the random distribution type ----
      # TODO ensure input has detailed descriptions to guide users interacting with the tool
        # good to add default values too e.g. logL: -5080,
      # tags$p("Enter or select values required to perform analysis. Default
      #                   values are shown. Press 'Run' when done."),
      # textInput(inputId = "logL",
      #           label = "Enter loglikelihood value", "-5080"),

      radioButtons("dist", "Distribution type:",
                   c("Normal" = "norm",
                     "Uniform" = "unif",
                     "Log-normal" = "lnorm",
                     "Exponential" = "exp")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # Input: Slider for the number of observations to generate ----
      sliderInput("n",
                  "Number of observations:",
                  value = 500,
                  min = 1,
                  max = 1000)

    ),

        # if you only want action after the user finishes input (E.g., manually typing in input field): actionButton()


    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot")),
                  tabPanel("Summary", verbatimTextOutput("summary")),
                  tabPanel("Table", tableOutput("table"))
      )

    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression
  d <- reactive({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)

    dist(input$n)
  })

  # TODO also insure the output has a description so the user knows what the output is showing
  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  output$plot <- renderPlot({
    dist <- input$dist
    n <- input$n

    hist(d(),
         main = paste("r", dist, "(", n, ")", sep = ""),
         col = "#75AADB", border = "white")
  })

  # Generate a summary of the data ----
  output$summary <- renderPrint({
    summary(d())
  })

  # Generate an HTML table view of the data ----
  output$table <- renderTable({
    d()
  })

}

# Create Shiny app ----
shinyApp(ui, server)
