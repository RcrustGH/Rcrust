###############################
## Rcrust (ui.r)
###############################
# function-def: textInputRow(inputId, label, value = "")
textInputRow <<- function(inputId, label, value = "") {
  div(
    style = "display:inline-block",
    tags$label(label, `for` = inputId),
    tags$input(id = inputId, type = "text", value = value, class = "input-small")
  )
}
# load placeholder for working_file if there isnt one
if (!exists("working_file")) {
  working_file <- ""
}
shinyUI(
  fixedPage(
    # Add custom CSS & Javascript for Progress Indicator
    tagList(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
        tags$script(type = "text/javascript", src = "busy.js")
      )
    ),
    div(
      class = "busy",
      p("Busy.."),
      img(src = "35.gif")
    ),
    # Logo for Rcrust
    titlePanel(img(src = "Rcrust_logo.png", align = "left", width = 150, height = 60), windowTitle = "Rcrust"),
    fixedRow(textInputRow("projects_directory", "Projects Directory", value = paste0(substring(getwd(), 1, nchar(getwd()) - 4), "Projects")),
      textInputRow("working_file", "Working File", working_file),
      actionButton("save", "Save"),
      actionButton("load", "Load"),
      actionButton("run", "Run"),
      actionButton("clear", "Clear"),
      actionButton("console", "Console"),
      align = "right"
    ),
    fixedRow(verbatimTextOutput("print_message")),
    # Define Tabs
    tabsetPanel(
      id = "inTabset",
      # Input Parameters UI Tab
      tabPanel(
        "Input Parameters",
        # Size Panel
        wellPanel(
          "Size", actionButton("size_panel_minimiser", label = "..."),
          conditionalPanel(
            "(input.size_panel_minimiser)%2 == 0",
            fixedRow(
              column(3, textInput("x_n", "X")),
              column(3, textInput("y_n", "Y"))
            )
          )
        ),
        # PT Panel
        wellPanel(
          "Pressure and Temperature",
          actionButton("pt_panel_minimiser", label = "..."),
          actionButton("pt_import_panel_minimiser", label = "<-"),
          conditionalPanel(
            "(input.pt_panel_minimiser)%2 == 0",
            conditionalPanel(
              "(input.pt_import_panel_minimiser)%2 == 1",
              fileInput("file_pt", "Import P-T definitions"),
              actionButton("import_pt", "Import")
            ),
            textInput("n_pt_def", "Number of PT definitions"),
            # Dynamic PT definition input boxes
            uiOutput("pt")
          )
        ),
        # Bulk Composition Panel
        wellPanel(
          "Bulk Composition",
          actionButton("bulk_comp_panel_minimiser", label = "..."),
          actionButton("bulk_import_panel_minimiser", label = "<-"),
          conditionalPanel(
            "(input.bulk_comp_panel_minimiser)%2 == 0",
            conditionalPanel(
              "(input.bulk_import_panel_minimiser)%2 == 1",
              fileInput("file_bulk", "Import bulk composition definitions"),
              actionButton("import_bulk", "Import")
            ),
            textInput("n_comp_trans", "Number of Component Transformations"),
            # Dynamic transformation input boxes
            uiOutput("trans"),
            conditionalPanel(
              "input.bulk_def_file == false",
              uiOutput("maj"),
              checkboxInput("set_oxygen_fugacity", "Set oxygen fugacity", value = FALSE),
              checkboxInput("calculate_traces", "Partition traces above solidus", value = FALSE),
              conditionalPanel(
                "input.calculate_traces == true",
                selectInput("apply_trace_correction", "Apply saturation corrections",
                  c(
                    "None",
                    "Zircon Saturation (Watson & Harrison 1983)",
                    "Monazite saturation (Montel 1993)",
                    "Zircon and Monazite saturation",
                    "Apatite saturation",
                    "Apatite & Monazite Saturation"
                  ),
                  selected = "None",
                  selectize = TRUE
                ),
                conditionalPanel(
                  "input.apply_trace_correction == \"Apatite & Monazite Saturation\"",
                  fixedRow(
                    column(
                      4,
                      selectInput("monazite_saturation", "Monazite saturation options",
                        c("Stepanov et al. 2012"),
                        selected = "Stepanov et al. 2012",
                        selectize = TRUE
                      )
                    ),
                    column(4,
                      offset = 2,
                      selectInput("apatite_saturation_ApMnz", "Apatite saturation options",
                        c(
                          "Harrison & Watson 1984",
                          "H&W with Bea et al. 1992",
                          "H&W with Pichavant et al. 1992",
                          "Wolf & London 1994"
                        ),
                        selected = "Harrison & Watson 1984",
                        selectize = TRUE
                      )
                    )
                  ),
                  fixedRow(
                    column(
                      4,
                      textInput(inputId = "Xmz", label = "Xmz (mole ratio LREE)", value = "0.83"),
                    ),
                    column(4,
                      offset = 2,
                      textInput(inputId = "D_ApMelt_LREE", label = "D Ap/Melt LREE")
                    )
                  )
                ),
                conditionalPanel(
                  "input.apply_trace_correction == \"Apatite saturation\"",
                  selectInput("apatite_saturation_Ap", "Apatite saturation options",
                    c(
                      "Harrison & Watson 1984",
                      "H&W with Bea et al. 1992",
                      "H&W with Pichavant et al. 1992",
                      "Wolf & London 1994"
                    ),
                    selected = "Harrison & Watson 1984",
                    selectize = TRUE
                  )
                ),
                textInput("kd_file", "Kd file"),
                uiOutput("traces")
              ),
              textInput("n_bulk_def", "Number of bulk definitions"),
              # Dynamic bulk definition input boxes
              uiOutput("bulk")
            ),
            conditionalPanel(
              "input.bulk_def_file == true",
              textInput("bulk_file", "Bulk file")
            ),
            checkboxInput("bulk_def_file", "Import definitions from file", value = FALSE)
          )
        )
      ),
      # End of Input Parameters UI Tab
      # Phase Manipulations UI Tab
      tabPanel(
        "Phase Manipulations",
        # Phase Addition Panel
        wellPanel(
          "Phase Addition", actionButton("phase_addition_panel_minimiser", label = "..."), actionButton("phase_addition_import_panel_minimiser", label = "<-"),
          conditionalPanel(
            "(input.phase_addition_panel_minimiser)%2 == 0",
            conditionalPanel(
              "(input.phase_addition_import_panel_minimiser)%2 == 1",
              fileInput("file_ph_add", "Import phase addition definitions"),
              actionButton("import_ph_add", "Import")
            ),
            checkboxInput("ph_add", "Perform Phase Addition?", value = FALSE),
            conditionalPanel(
              "input.ph_add == true",
              textInput("n_ph_add_def", "Number of addition definitions"),
              # Dynamic phase addition definition input boxes
              uiOutput("ph_add")
            )
          )
        ),
        # Phase Extraction Panel
        wellPanel(
          "Phase Extraction", actionButton("phase_extraction_panel_minimiser", label = "..."), actionButton("phase_extraction_import_panel_minimiser", label = "<-"),
          conditionalPanel(
            "(input.phase_extraction_panel_minimiser)%2 == 0",
            conditionalPanel(
              "(input.phase_extraction_import_panel_minimiser)%2 == 1",
              fileInput("file_ph_extr", "Import phase extraction definitions"),
              actionButton("import_ph_extr", "Import")
            ),
            checkboxInput("ph_extr", "Perform Phase Extraction?", value = FALSE),
            conditionalPanel(
              "input.ph_extr == true",
              checkboxInput("reequilibrate_steps", "Re-equilibrate reactive subsystem after phase extraction?", value = TRUE),
              textInput("n_ph_extr_def", "Number of extraction definitions"),
              # Dynamic phase extraction input boxes
              uiOutput("ph_extr")
            )
          )
        )
      ),
      # End of Phase Manipulations UI Tab
      # Modelling Options UI Tab
      tabPanel(
        "Modelling Options",
        # Modelling data Panel
        wellPanel(
          "Modelling Data", actionButton("modelling_data_panel_minimiser", label = "..."),
          conditionalPanel(
            "(input.modelling_data_panel_minimiser)%2 == 0",
            textInput("meemum_path", "Meemum version", "meemum.exe"),
            textInput("perplex_option_file", "Perple_X Option File", "perplex_option.dat"),
            textInput("thermodynamic_data_file", "Thermodynamic Data File", "hp11ver.dat"),
            textInput("solution_models_file", "Solution Models File", "solution_model_673.dat"),
            # Dynamic solution models input boxes
            uiOutput("solution_models")
          )
        ),
        # Additional optional parameters Panel
        wellPanel(
          "Additional optional parameters", actionButton("additional_optional_parameters_panel_minimiser", label = "..."),
          conditionalPanel(
            "(input.additional_optional_parameters_panel_minimiser)%2 == 0",
            textInput("saturated_components", "Saturated components"),
            textInput("saturated_phase_components", "Saturated phase components"),
            textInput("independent_potential_fugacity_activity", "Independent potential/fugacity/activity"),
            checkboxInput("calculate_activities", "Calculate Activities?", value = FALSE),
            textInput("G_pure_phases", "Gibbs enthalpy of pure phases for activity calculation"),
            textInput("exclude_phases", "Exclude phases"),
            textInput("n_phases_to_rename", "Number of phases to rename"),
            uiOutput("rename_phases_inputs")
          )
        ),
        # Sean-tag
        # Component Packet Panel
        wellPanel(
          "Component Packet", actionButton("component_packet_panel_minimiser", label = "..."),
          conditionalPanel(
            "(input.component_packet_panel_minimiser)%2 == 0",
            checkboxInput("component_packet", "Component Packet", value = FALSE),
            conditionalPanel(
              "input.component_packet == true",
              tags$div(
                title = "Enter components separated by commas (case-sensitive).\nIf the component does not exist in the thermodynamic data file then it will be created in Major elements input.",
                textInput("cp_components", "Packet Components")
              ),
              # Dynamic component packet input boxes
              uiOutput("cp_packet_ui")
            )
          )
        ),
        # Extra Settings Panel
        wellPanel(
          "Extra Settings", actionButton("extra_settings_panel_minimiser", label = "..."),
          conditionalPanel(
            "(input.extra_settings_panel_minimiser)%2 == 0",
            checkboxInput("print_meem", "Print meemum outputs?", value = FALSE),
            checkboxInput("export_meemum_output", "Export meemum outputs?", value = FALSE),
            selectInput("end_of_calc", "When calculation is complete:", c("Return to Interface" = "Return to Interface", "Logout" = "Logout", "Shutdown" = "Shutdown"), selected = "Return to Interface", selectize = TRUE)
          )
        )
      ),
      # End of Modelling Options UI Tab
      # Outputs UI Tab
      tabPanel(
        "Outputs",
        sidebarLayout(
          sidebarPanel(
            textInput("phase_aliases", "Phase Aliases"),
            selectInput("output_type", "Select Output", c("Data File", "Grid", "Phase Abundance Along Path", "PAM"), selected = "Data File", selectize = TRUE),
            uiOutput("output_form_selection"),
            # Dynamic output selection boxes
            uiOutput("output_selection")
          ),
          mainPanel(
            h4(textOutput("output_header", container = span)),
            uiOutput("output_view")
          )
        )
      ),
      selected = "Input Parameters"
    )
    # End of tabs
  )
  # End of page
)
# End of Shiny
