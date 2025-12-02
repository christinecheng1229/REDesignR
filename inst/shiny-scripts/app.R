# Purpose: R script for Shiny App.
# Author: Christine Cheng
# Date: December 1, 2025
# Version: 1.2
# Bugs and Issues: None known.

# ---------------- Define UI ----------------
ui <- fluidPage(
  titlePanel("REDesignR Interactive App"),

  tags$p("This Shiny App is a tool for simulating restriction enzyme
         co-digestion and generating restriction maps and gel simulations to
         aid in digestion experiment design optimization."),

  fluidRow(

    # ------------ Inputs panel ------------
    column(
      width = 4,

      # ------------ Description ------------
      h3("REDesignR Co-Digestion Simulation"),
      helpText(
        "Use this panel to specify the DNA sequence to digest and two",
        "restriction enzymes along with their recognition sequences.",
        "The results will update the visualization tabs on the right when",
        "you click 'Run simulation'."
      ),

      tags$hr(),

      # ------------ Inputs ------------
      h3("Co-digestion Inputs"),

      textAreaInput(
        inputId = "dna_sequence",
        label   = "DNA sequence",
        placeholder = "Paste a raw sequence (A/C/G/T characters only)...",
        rows = 4
      ),

      # ------------ Enzymes file upload ------------
      fileInput(
        inputId = "enzyme_csv",
        label   = "Upload 2-enzyme CSV (cols: Name, RecognitionSeq)",
        accept  = ".csv"
      ),

      # ------------ Download sample enzyme dataset ------------
      downloadButton(
        outputId = "download_sample_csv",
        label = "Download sample enzyme CSV"
      ),

      helpText("— OR — enter enzyme data manually:"),

      # ------------ Manual Enzyme Inputs ------------
      textInput(
        inputId = "enzyme_1_name",
        label   = "Enzyme 1 Name",
        placeholder = "AaaI"
      ),
      textInput(
        inputId = "enzyme_1_seq",
        label   = "Enzyme 1 Recognition Sequence",
        placeholder = "C/GGCCG"
      ),

      textInput(
        inputId = "enzyme_2_name",
        label   = "Enzyme 2 Name",
        placeholder = "AagI"
      ),
      textInput(
        inputId = "enzyme_2_seq",
        label   = "Enzyme 2 Recognition Sequence",
        placeholder = "AT/CGAT"
      ),

      tags$hr(),

      # Run simulation button
      actionButton(
        inputId = "run_codigest",
        label = "Run simulation",
        class = "btn-primary"
      )
    ),

    # ------------ Output Tabs ------------
    column(
      width = 8,

      tabsetPanel(

        # ------------ Restriction Map Tab ------------
        tabPanel(
          title = "Restriction Map",
          br(),
          h3("Restriction Map"),
          plotOutput("restriction_map_plot")
        ),

        # ------------ Gel Simulation Tab ------------
        tabPanel(
          title = "Gel Simulation",
          br(),
          h3("Simulated Agarose Gel"),
          plotOutput("gel_plot")
        )
      )
    )
  )
)


# ---------------- Define Server ----------------
server <- function(input, output, session) {

  # ----------------- Allow user to download sample dataset ------------------
  output$download_sample_csv <- downloadHandler(
    filename = function() {
      "REtype2_itype2_511.csv"
    },
    content = function(file) {
      sample_path <- system.file("extdata", "REtype2_itype2_511.csv",
                                 package = "REDesignR")

      if (sample_path == "") {
        stop("Sample CSV not found inside inst/extdata of the REDesignR package.")
      }

      file.copy(sample_path, file)
    }
  )

  # Construct enzyme dataframe from CSV or manual inputs
  enzyme_tbl <- reactive({
    if (!is.null(input$enzyme_csv)) {
      # User uploaded a CSV
      df <- read.csv(input$enzyme_csv$datapath, stringsAsFactors = FALSE)

      validate(
        need(
          nrow(df) == 2 &&
            all(c("Name", "RecognitionSeq") %in% names(df)),
          "CSV must contain exactly 2 rows and columns: Name, RecognitionSeq"
        )
      )

      return(as_tibble(df))
    }

    # Otherwise use manual inputs
    validate(
      need(
        input$enzyme_1_name != "" &&
          input$enzyme_1_seq  != "" &&
          input$enzyme_2_name != "" &&
          input$enzyme_2_seq  != "",
        "Please enter both enzymes (Name + RecognitionSeq) or upload a CSV."
      )
    )

    tibble::tibble(
      Name           = c(input$enzyme_1_name, input$enzyme_2_name),
      RecognitionSeq = c(input$enzyme_1_seq,  input$enzyme_2_seq)
    )
  })

  # Store co-digestion results
  codigest_out <- reactiveVal(NULL)

  observeEvent(input$run_codigest, {
    req(input$dna_sequence)

    # convert raw string to DNAString
    dna <- Biostrings::DNAString(gsub("\\s+", "", input$dna_sequence))

    enzymes <- enzyme_tbl()

    # Call your provided function EXACTLY as specified
    res <- simulateCoDigest(
      dnaSeq  = dna,
      enzymes = enzymes
    )

    # simulateCoDigest returns NULL if no digestion occurs
    if (is.null(res)) {
      codigest_out(NULL)
      return(NULL)
    }

    codigest_out(res)
  })


  # ------------ Restriction Map Visualization ------------
  output$restriction_map_plot <- renderPlot({
    req(codigest_out())

    df <- codigest_out()$digestDf
    multi <- codigest_out()$isCoDigest

    plotRestrictionMap(
      codigestDf = df,
      multiDigest = multi,
      showLengths = TRUE
    )
  })


  # ------------ Gel Simulation Visualization ------------
  output$gel_plot <- renderPlot({
    req(codigest_out())

    df <- codigest_out()$digestDf
    multi <- codigest_out()$isCoDigest

    simulateGel(
      codigestDf     = df,
      multiDigest    = multi,
      labelFragments = FALSE
    )
  })
}

# ----- Create Shiny App -----
shinyApp(ui = ui, server = server)

# [END]
