################################################################################
#  mod_cpm_contour.R  --  CPM Contour Plotting (Shiny module)
################################################################################
#
#  Migrated from app_v10l.R:
#    UI:     lines  966 - 1090
#    Server: lines 2532 - 3261
#
#  Workflow:
#    1. Upload dF/dT data CSV (from GraphPad Prism). First column =
#       temperature, remaining columns = repeated sample names (one per
#       replicate). The tool averages replicates per condition.
#    2. Upload Tm data CSV. First column = condition labels, remaining
#       columns = replicate Tm measurements. The tool computes mean & SEM.
#    3. User pairs each Tm-file row to one of the dF/dT-file sample
#       names (since the two files don't have to use the same labels).
#    4. "Process Data" computes the mean / SEM dF/dT matrix and records
#       the pairing for the Tm plot.
#    5. Two visualisations are produced:
#       - dF/dT heatmap: temperature × sample, fill = normalised dF/dT.
#       - Tm vs. log[concentration] scatter, with the same heatmap shown
#         beneath the points as a fingerprint.
#
#  This tool uses reactiveValues (not reactiveVal) because the original
#  had several pieces of state that update in lockstep. I preserved that.
#
################################################################################

# -- UI -----------------------------------------------------------------------
cpm_contour_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::div(style = "display: none;",
      shiny::actionButton(ns("clear"), "\U0001f504  Clear All Data", class = "btn-clear")
    ),

    shiny::div(class = "sticky-tool",
      shiny::fluidRow(
        # ----- Left column: workflow (scrolls) --------------------------
        shiny::column(3,
          shiny::div(class = "workflow-col",
            lab_card(
              step_title(1, "Upload dF/dT Data"),
              info_box(paste("Upload CSV file exported from GraphPad Prism",
                             "containing dF/dT data with replicate samples.")),
              shiny::fileInput(ns("dfdt_file"), NULL, accept = c(".csv"),
                               buttonLabel = "Browse\u2026",
                               placeholder = "No dF/dT file selected"),
              shiny::uiOutput(ns("dfdt_status"))
            ),
            lab_card(
              step_title(2, "Upload Tm Data"),
              info_box("Upload CSV file containing Tm scatter data for the same samples."),
              shiny::fileInput(ns("tm_file"), NULL, accept = c(".csv"),
                               buttonLabel = "Browse\u2026",
                               placeholder = "No Tm file selected"),
              shiny::uiOutput(ns("tm_status"))
            ),
            lab_card(
              step_title(3, "Match Tm to Samples"),
              info_box("Match each Tm concentration to its corresponding dF/dT sample."),
              shiny::uiOutput(ns("tm_pairing_ui"))
            ),
            lab_card(
              step_title(4, "Plot Settings"),
              shiny::tags$label("dF/dT colour palette",
                                style = "color:var(--muted);font-size:0.78rem;"),
              shiny::selectInput(ns("palette"), NULL,
                choices = c(
                  "Grayscale (white \U2192 black)" = "grayscale",
                  "Inferno (dark \U2192 yellow)"   = "inferno",
                  "Magma (dark \U2192 white)"      = "magma",
                  "Viridis (dark \U2192 yellow)"   = "viridis",
                  "Plasma (dark \U2192 pink)"      = "plasma",
                  "Blue \U2192 White \U2192 Red"   = "RdBu",
                  "Yellow \U2192 Orange \U2192 Red"= "YlOrRd"),
                selected = "grayscale", width = "100%"),
              shiny::br(),
              shiny::tags$label("T\u2098 trace colour",
                                style = "color:var(--muted);font-size:0.78rem;"),
              info_box("Colour for the Tm line, error bars, and data points. Choose something that contrasts with your contour palette so they don't blend in."),
              shiny::selectInput(ns("tm_colour"), NULL,
                choices = c(
                  "Cyan"        = "#06B6D4",
                  "Bright cyan" = "#00E5FF",
                  "Magenta"     = "#E91E63",
                  "Yellow"      = "#FACC15",
                  "Lime"        = "#84CC16",
                  "Orange"      = "#F97316",
                  "White"       = "#FFFFFF",
                  "Red"         = "#EF4444"),
                selected = "#06B6D4", width = "100%"),
              shiny::br(),
              shiny::tags$label("Intensity threshold",
                                style = "color:var(--muted);font-size:0.78rem;"),
              info_box(paste("Values below this threshold are shown as background colour.",
                             "Remaining values are rescaled to use the full palette.")),
              shiny::sliderInput(ns("threshold"), NULL,
                                 min = 0, max = 0.95, value = 0, step = 0.05,
                                 width = "100%"),
              shiny::br(),
              shiny::tags$label("X-axis title",
                                style = "color:var(--muted);font-size:0.78rem;"),
              shiny::textInput(ns("x_title"), NULL,
                value = "Log [Concentration (\u03bcM)]",
                width = "100%"),
              shiny::br(),
              shiny::tags$button("\u2705  Temperature Range", class = "adv-toggle",
                onclick = sprintf("$('#%s').slideToggle(200)", ns("yrange_panel"))),
              shiny::div(id = ns("yrange_panel"), style = "display:none;",
                shiny::div(class = "settings-group",
                  shiny::div(class = "settings-group-title",
                             "Tm Plot Y-Axis Range (\u00b0C)"),
                  shiny::div(style = "font-size:0.75rem;color:var(--muted);margin-bottom:0.5rem;",
                    "Full data range shown by default. Adjust to zoom in on a specific region."),
                  shiny::fluidRow(
                    shiny::column(6, shiny::numericInput(ns("ymin"), "Min",
                                                         value = 25.5, step = 1, width = "100%")),
                    shiny::column(6, shiny::numericInput(ns("ymax"), "Max",
                                                         value = 89.5, step = 1, width = "100%")))
                )
              )
            ),
            lab_card(
              step_title(5, "Standard Panel (Optional)"),
              info_box(paste("Upload a dF/dT file for a control sample (e.g.",
                             "free nanobody) to display as a small reference",
                             "panel to the right of the main Tm plot. Useful",
                             "for showing that high-Tm peaks in titration data",
                             "match a known control's melt signature.")),
              shiny::fileInput(ns("standard_file"), NULL, accept = c(".csv"),
                               buttonLabel = "Browse\u2026",
                               placeholder = "No standard file (optional)"),
              shiny::uiOutput(ns("standard_status")),
              shiny::br(),
              shiny::selectInput(ns("standard_colour"), "Panel colour",
                choices = c(
                  "Auto (match main palette)" = "auto",
                  "Red"           = "#E41A1C",
                  "Blue"          = "#377EB8",
                  "Green"         = "#4DAF4A",
                  "Purple"        = "#984EA3",
                  "Orange"        = "#FF7F00",
                  "Teal"          = "#1B9E77",
                  "Magenta"       = "#E7298A",
                  "Yellow-Brown"  = "#E6AB02"),
                selected = "auto", width = "100%"),
              shiny::sliderInput(ns("standard_opacity"),
                "Panel opacity",
                min = 0.1, max = 1.0, value = 0.85, step = 0.05,
                width = "100%")
            ),
            # The Analyse box (last in the left column) combines what used
            # to be separate "Process Data" and "Export" cards. Putting
            # them together matches the pattern used by every other tool -
            # the run button sits above the export buttons since you can
            # only export after analysis succeeds.
            lab_card(
              step_title(6, "Analyse"),
              shiny::actionButton(ns("process"), "\u25b6  Run Analysis",
                                  class = "btn-run"),
              shiny::br(), shiny::br(),
              shiny::uiOutput(ns("sample_count")),
              shiny::br(),
              shiny::downloadButton(ns("export_stats"),   "\u2b07 Statistics CSV", class = "btn-download"),
              shiny::br(), shiny::br(),
              shiny::downloadButton(ns("export_heatmap"), "\u2b07 Heatmap PNG",    class = "btn-download"),
              shiny::br(), shiny::br(),
              shiny::downloadButton(ns("export_tmplot"),  "\u2b07 Tm Plot PNG",    class = "btn-download")
            )
          )  # close workflow-col
        ),

        # ----- Middle column: previews (sticky) ------------------------
        shiny::column(7,
          shiny::div(class = "preview-col",
            lab_card(
              shiny::div(class = "lab-card-title", "\U0001f4c8  dF/dT Heatmap"),
              shiny::plotOutput(ns("heatmap"), height = "450px")
            ),
            lab_card(
              shiny::div(class = "lab-card-title", "\U0001f321\ufe0f  Tm vs. Concentration"),
              shiny::plotOutput(ns("tm_plot"), height = "340px")
            )
          )  # close preview-col
        ),

        # ----- Right column: stats outputs (scrolls) ------------------
        shiny::column(2,
          shiny::div(class = "workflow-col",
            lab_card(
              shiny::div(class = "lab-card-title", "\U0001f4ca  Sample Statistics"),
              shiny::uiOutput(ns("stats_ui"))
            ),
            lab_card(
              shiny::div(class = "lab-card-title", "\U0001f321\ufe0f  Tm Summary"),
              shiny::uiOutput(ns("tm_summary"))
            )
          )  # close workflow-col
        )
      )
    )
  )
}


# -- Server -------------------------------------------------------------------
cpm_contour_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # All state is held in one reactiveValues bundle because several pieces
    # update together (e.g. clear must reset everything; process produces
    # both dfdt_processed and tm_processed in one go).
    state <- shiny::reactiveValues(
      dfdt_raw       = NULL,   # raw dF/dT (temps, sample_names, values)
      tm_raw         = NULL,   # raw Tm    (concentrations, values, means, sems)
      dfdt_processed = NULL,   # averaged-per-sample dF/dT
      tm_processed   = NULL,   # paired Tm summary
      sample_names   = NULL,
      temperatures   = NULL,
      tm_pairing     = NULL    # named character: i-th entry = sample paired to Tm row i
    )

    # =====================================================================
    # State sources of truth
    # =====================================================================
    # Mirrors AKTA/BCA/CPM-QC pattern: a current_*_file reactive per upload
    # input so the rest of the code reads from one place, populated either
    # by user upload OR by the example loader on session start.
    current_dfdt_file <- shiny::reactiveVal(NULL)
    current_tm_file   <- shiny::reactiveVal(NULL)
    current_standard_file <- shiny::reactiveVal(NULL)

    # =====================================================================
    # Parsing helpers - extracted so the upload observers and example
    # loader share identical parse logic.
    # =====================================================================
    .parse_dfdt <- function(path) {
      raw_data <- utils::read.csv(path,
                                  stringsAsFactors = FALSE,
                                  check.names      = FALSE,
                                  fileEncoding     = "UTF-8-BOM")
      temperatures <- as.numeric(raw_data[, 1])
      temperatures <- temperatures[!is.na(temperatures)]
      sample_names <- colnames(raw_data)[-1]
      dfdt_values <- as.data.frame(lapply(raw_data[, -1, drop = FALSE],
                                          as.numeric))
      dfdt_values <- dfdt_values[seq_len(length(temperatures)), ,
                                 drop = FALSE]
      colnames(dfdt_values) <- sample_names
      list(temperatures = temperatures,
           sample_names = sample_names,
           values       = dfdt_values)
    }

    .parse_tm <- function(path) {
      raw_data <- utils::read.csv(path,
                                  stringsAsFactors = FALSE,
                                  check.names      = FALSE,
                                  fileEncoding     = "UTF-8-BOM")
      concentrations <- as.character(raw_data[, 1])
      n_rows         <- nrow(raw_data)
      tm_values <- list()
      tm_means  <- numeric(n_rows)
      tm_sems   <- numeric(n_rows)
      for (i in seq_len(n_rows)) {
        values <- as.numeric(raw_data[i, -1])
        values <- values[!is.na(values)]
        tm_values[[i]] <- values
        tm_means[i]    <- if (length(values) > 0) mean(values) else NA
        tm_sems[i]     <- if (length(values) > 1)
                            stats::sd(values) / sqrt(length(values))
                          else 0
      }
      list(concentrations = concentrations,
           values         = tm_values,
           means          = tm_means,
           sems           = tm_sems,
           n_rows         = n_rows)
    }

    # =====================================================================
    # Internal: clear all derived state. Called on new uploads and
    # by the navbar Clear button.
    # =====================================================================
    .clear_state <- function() {
      state$dfdt_raw         <- NULL
      state$tm_raw           <- NULL
      state$dfdt_processed   <- NULL
      state$tm_processed     <- NULL
      state$sample_names     <- NULL
      state$temperatures     <- NULL
      state$tm_pairing       <- NULL
      state$standard_raw     <- NULL
      state$standard_processed <- NULL
    }

    # =====================================================================
    # File uploads route through current_*_file -> observer parses
    # =====================================================================

    shiny::observeEvent(input$dfdt_file, {
      # New dFdT file invalidates everything downstream - including any
      # processed matrices and the Tm pairing (sample names will differ).
      .clear_state()
      current_dfdt_file(input$dfdt_file)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$tm_file, {
      # New Tm file only invalidates the Tm side: keep the dfdt_raw if
      # one is already loaded. But clear the pairing because row count
      # may differ.
      state$tm_raw       <- NULL
      state$tm_processed <- NULL
      state$tm_pairing   <- NULL
      current_tm_file(input$tm_file)
    }, ignoreInit = TRUE)

    # ---- Standard overlay upload --------------------------------------
    # Optional control file (e.g. nanobody-only dF/dT). Uses the same
    # .parse_dfdt() helper as the main data since the format is identical;
    # the only difference is downstream presentation as a faded overlay
    # on the Tm plot.
    shiny::observeEvent(input$standard_file, {
      state$standard_raw       <- NULL
      state$standard_processed <- NULL
      current_standard_file(input$standard_file)
    }, ignoreInit = TRUE)

    # Parse observers - single path for both upload and example load.
    shiny::observeEvent(current_dfdt_file(), {
      cf <- current_dfdt_file()
      shiny::req(cf)
      tryCatch({
        parsed <- .parse_dfdt(cf$datapath)
        state$dfdt_raw <- parsed
        n_unique <- length(unique(parsed$sample_names))
        n_reps   <- as.integer(table(parsed$sample_names)[1])
        shiny::showNotification(
          sprintf("\u2713 dF/dT loaded: %d samples \u00d7 %d replicates, %d temps",
                  n_unique, n_reps, length(parsed$temperatures)),
          type = "message", duration = 3)
      }, error = function(e) {
        shiny::showNotification(paste("Error reading dF/dT file:", e$message),
                                type = "error", duration = 5)
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    shiny::observeEvent(current_tm_file(), {
      cf <- current_tm_file()
      shiny::req(cf)
      tryCatch({
        parsed <- .parse_tm(cf$datapath)
        state$tm_raw     <- parsed
        state$tm_pairing <- NULL
        shiny::showNotification(
          sprintf("\u2713 Tm loaded: %d conditions", parsed$n_rows),
          type = "message", duration = 3)
      }, error = function(e) {
        shiny::showNotification(paste("Error reading Tm file:", e$message),
                                type = "error", duration = 5)
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    shiny::observeEvent(current_standard_file(), {
      cf <- current_standard_file()
      shiny::req(cf)
      tryCatch({
        parsed <- .parse_dfdt(cf$datapath)
        state$standard_raw <- parsed
        n_unique <- length(unique(parsed$sample_names))
        shiny::showNotification(
          sprintf("\u2713 Standard loaded: %d sample(s)", n_unique),
          type = "message", duration = 3)
      }, error = function(e) {
        shiny::showNotification(paste("Error reading standard file:", e$message),
                                type = "error", duration = 5)
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # Compute the standard's normalized mean matrix the same way as the
    # main data: average replicates, normalise each column to [0,1].
    # Triggered whenever standard_raw changes (which is just on upload).
    shiny::observeEvent(state$standard_raw, {
      r <- state$standard_raw
      shiny::req(r)
      tryCatch({
        unique_samples <- unique(r$sample_names)
        n_temps        <- length(r$temperatures)
        n_samples      <- length(unique_samples)
        mean_matrix <- matrix(NA, nrow = n_temps, ncol = n_samples,
                              dimnames = list(NULL, unique_samples))
        for (i in seq_along(unique_samples)) {
          cols <- which(r$sample_names == unique_samples[i])
          for (j in seq_len(n_temps)) {
            v <- as.numeric(r$values[j, cols])
            v <- v[!is.na(v)]
            if (length(v) > 0) mean_matrix[j, i] <- mean(v)
          }
        }
        norm_matrix <- apply(mean_matrix, 2, function(col) {
          mn <- min(col, na.rm = TRUE); mx <- max(col, na.rm = TRUE)
          rng <- mx - mn
          if (is.na(rng) || rng == 0) rep(0, length(col)) else (col - mn) / rng
        })
        # R's apply() drops the dimension when there's only 1 column, so
        # a single-sample standard (which is the common case - one
        # NB-only control) would arrive here as a length-n_temps vector
        # instead of a 1-column matrix. Subsequent code subsets it as
        # norm[mask, col_idx] which fails on a vector. Reshape back.
        if (!is.matrix(norm_matrix)) {
          norm_matrix <- matrix(norm_matrix, ncol = n_samples,
                                dimnames = list(NULL, unique_samples))
        }
        colnames(norm_matrix) <- unique_samples
        state$standard_processed <- list(temperatures = r$temperatures,
                                         sample_names = unique_samples,
                                         norm         = norm_matrix)
      }, error = function(e) {
        shiny::showNotification(paste("Error processing standard:", e$message),
                                type = "error", duration = 5)
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # =====================================================================
    # Pairing UI - one selectInput per Tm row
    # =====================================================================
    output$tm_pairing_ui <- shiny::renderUI({
      shiny::req(state$tm_raw, state$dfdt_raw)
      tm_raw         <- state$tm_raw
      unique_samples <- unique(state$dfdt_raw$sample_names)

      pairing_inputs <- lapply(seq_len(tm_raw$n_rows), function(i) {
        conc_label <- tm_raw$concentrations[i]
        tm_mean    <- tm_raw$means[i]
        sel        <- if (!is.null(state$tm_pairing) &&
                          i <= length(state$tm_pairing) &&
                          !is.na(state$tm_pairing[i])) {
                        state$tm_pairing[i]
                      } else {
                        "-- Select sample --"
                      }

        shiny::div(
          style = paste("margin-bottom:10px;background:#1E2D45;",
                        "padding:8px;border-radius:4px;"),
          shiny::div(style = "color:#7A8FAD;font-size:0.85rem;margin-bottom:4px;",
            sprintf("Row %d: %s (Tm: %.1f\u00b0C)", i, conc_label, tm_mean)),
          shiny::selectInput(ns(paste0("tm_pair_", i)), NULL,
            choices  = c("-- Select sample --", unique_samples),
            selected = sel,
            width    = "100%")
        )
      })
      shiny::tagList(pairing_inputs)
    })

    # =====================================================================
    # Analysis core - extracted into a helper so we can call it from:
    #   - user clicks "Run Analysis" (input$process)
    #   - auto-run when dfdt/tm/pairing are all ready (e.g. example load)
    # =====================================================================
    .run_analysis <- function() {
      if (is.null(state$dfdt_raw)) return(invisible())
      tryCatch({
        temperatures <- state$dfdt_raw$temperatures
        sample_names <- state$dfdt_raw$sample_names
        dfdt_values  <- state$dfdt_raw$values

        unique_samples <- unique(sample_names)
        state$sample_names <- unique_samples
        state$temperatures <- temperatures

        n_temps   <- length(temperatures)
        n_samples <- length(unique_samples)

        mean_matrix <- matrix(NA, nrow = n_temps, ncol = n_samples,
                              dimnames = list(NULL, unique_samples))
        sem_matrix  <- matrix(NA, nrow = n_temps, ncol = n_samples,
                              dimnames = list(NULL, unique_samples))
        n_matrix    <- matrix(NA, nrow = n_temps, ncol = n_samples,
                              dimnames = list(NULL, unique_samples))

        for (i in seq_along(unique_samples)) {
          sample      <- unique_samples[i]
          sample_cols <- which(sample_names == sample)
          sample_data <- dfdt_values[, sample_cols, drop = FALSE]
          for (j in seq_len(n_temps)) {
            values <- as.numeric(sample_data[j, ])
            values <- values[!is.na(values)]
            if (length(values) > 0) {
              mean_matrix[j, i] <- mean(values)
              n_matrix[j, i]    <- length(values)
              sem_matrix[j, i]  <- if (length(values) > 1)
                                     stats::sd(values) / sqrt(length(values))
                                   else 0
            }
          }
        }

        # Per-sample normalisation: scale so each column's max=1, min=0.
        # Matches CPM Peak Picker convention.
        norm_matrix <- apply(mean_matrix, 2, function(col) {
          mn  <- min(col, na.rm = TRUE)
          mx  <- max(col, na.rm = TRUE)
          rng <- mx - mn
          if (is.na(rng) || rng == 0) rep(0, length(col))
          else (col - mn) / rng
        })
        colnames(norm_matrix) <- unique_samples

        state$dfdt_processed <- list(temperatures = temperatures,
                                     sample_names = unique_samples,
                                     mean         = mean_matrix,
                                     norm         = norm_matrix,
                                     sem          = sem_matrix,
                                     n            = n_matrix)

        # Collect Tm pairing from the selectInputs the user has touched.
        # For the example-load case, .auto_pair_example_samples() sets
        # state$tm_pairing programmatically so the selectInputs render
        # with the right defaults - by the time .run_analysis() reads
        # the inputs they reflect those defaults.
        if (!is.null(state$tm_raw)) {
          tm_raw  <- state$tm_raw
          pairing <- character(tm_raw$n_rows)
          for (i in seq_len(tm_raw$n_rows)) {
            v <- input[[paste0("tm_pair_", i)]]
            pairing[i] <- if (!is.null(v) && v != "-- Select sample --") v
                          else NA_character_
          }
          state$tm_pairing <- pairing

          valid <- !is.na(pairing)
          if (any(valid)) {
            state$tm_processed <- list(
              concentrations = tm_raw$concentrations[valid],
              sample_names   = pairing[valid],
              means          = tm_raw$means[valid],
              sems           = tm_raw$sems[valid])
          }
        }

        shiny::showNotification("\u2713 Analysis complete",
                                type = "message", duration = 3)
      }, error = function(e) {
        shiny::showNotification(paste("Error processing data:", e$message),
                                type = "error", duration = 5)
      })
    }

    # Trigger 1: user clicked "Run Analysis"
    shiny::observeEvent(input$process, .run_analysis())

    # Trigger 2: auto-run when dfdt_raw, tm_raw, and pairing inputs are all
    # ready. Same pending-flag pattern as CPM Peak / CPM QC - the pairing
    # UI renders only AFTER both raw datasets are loaded, so we can't
    # auto-run directly from observeEvent on either dataset.
    .contour_pending_run <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(state$dfdt_raw, {
      if (!is.null(state$dfdt_raw) && !is.null(state$tm_raw))
        .contour_pending_run(TRUE)
    }, ignoreInit = TRUE)
    shiny::observeEvent(state$tm_raw, {
      if (!is.null(state$dfdt_raw) && !is.null(state$tm_raw))
        .contour_pending_run(TRUE)
    }, ignoreInit = TRUE)
    shiny::observe({
      if (!isTRUE(.contour_pending_run())) return()
      shiny::req(state$dfdt_raw, state$tm_raw)
      # Wait for all the tm_pair_* inputs to be bound. The pairing UI
      # renders one selectInput per Tm row; we know how many rows there
      # are from state$tm_raw$n_rows.
      n_rows <- state$tm_raw$n_rows
      for (i in seq_len(n_rows)) {
        if (is.null(input[[paste0("tm_pair_", i)]])) return()
      }
      .contour_pending_run(FALSE)
      .run_analysis()
    })

    # =====================================================================
    # Clear
    # =====================================================================
    shiny::observeEvent(input$clear, {
      current_dfdt_file(NULL); current_tm_file(NULL)
      current_standard_file(NULL)
      .clear_state()
      # Reload examples so preview is never empty after Clear.
      tryCatch(.load_example_files(), error = function(e) NULL)
      shiny::showNotification("\u2713 All data cleared",
                              type = "message", duration = 2)
    })

    # =====================================================================
    # Trigger 3: session start with bundled examples
    # =====================================================================
    # Loads BOTH example files (dF/dT + Tm). The auto-pair observer
    # below handles the pairing; this just kicks off the file load.
    .load_example_files <- function() {
      cf_dfdt <- .cpm_contour_dfdt_example_file()
      cf_tm   <- .cpm_contour_tm_example_file()
      if (is.null(cf_dfdt) || is.null(cf_tm)) return(invisible())
      current_dfdt_file(cf_dfdt)
      current_tm_file(cf_tm)
    }

    # Auto-pair Tm rows to dF/dT samples by row order, whenever both
    # raw datasets are loaded AND no pairing has been set yet. This
    # fires for ALL loads (example AND user uploads), because lining up
    # row N of the Tm file with the Nth unique sample of the dF/dT
    # file is the standard layout - and if a user's data doesn't match
    # this convention they can still override individual rows via the
    # pairing dropdowns.
    #
    # The `is.null(state$tm_pairing)` guard means we don't trample over
    # manual edits: once the user has paired (or modified) anything,
    # this observer no-ops until pairing is reset (which happens on
    # new file uploads).
    shiny::observe({
      shiny::req(state$dfdt_raw, state$tm_raw)
      if (!is.null(state$tm_pairing)) return()  # don't overwrite manual pairing
      unique_samples <- unique(state$dfdt_raw$sample_names)
      n_rows         <- state$tm_raw$n_rows
      pairing <- character(n_rows)
      for (i in seq_len(n_rows)) {
        pairing[i] <- if (i <= length(unique_samples)) unique_samples[i]
                      else NA_character_
      }
      state$tm_pairing <- pairing
    })

    .example_loader_obs <- shiny::observe({
      .example_loader_obs$destroy()
      tryCatch(.load_example_files(), error = function(e)
        message("[CPM Contour] example load failed: ", conditionMessage(e)))
    })

    # =====================================================================
    # Status indicators
    # =====================================================================
    output$dfdt_status <- shiny::renderUI({
      if (!is.null(state$dfdt_raw)) {
        n_samples <- length(state$dfdt_raw$sample_names)
        n_temps   <- length(state$dfdt_raw$temperatures)
        status_pill("ready",
          sprintf("\u2713 Loaded: %d samples, %d temps", n_samples, n_temps))
      }
    })

    output$tm_status <- shiny::renderUI({
      if (!is.null(state$tm_raw)) {
        n_conc <- length(state$tm_raw$concentrations)
        status_pill("ready",
          sprintf("\u2713 Loaded: %d concentrations", n_conc))
      }
    })

    # ---- Standard panel status ----------------------------------------
    # The standard is now rendered as its own side panel rather than
    # overlaid on the main plot, so we just need a load-status pill -
    # no per-sample concentration input needed (all unique samples in
    # the standard file are shown as separate tile columns in file order).
    output$standard_status <- shiny::renderUI({
      if (!is.null(state$standard_raw)) {
        n <- length(unique(state$standard_raw$sample_names))
        status_pill("ready", sprintf("\u2713 Loaded: %d standard sample(s)", n))
      }
    })

    output$sample_count <- shiny::renderUI({
      if (!is.null(state$dfdt_processed)) {
        n_unique <- length(state$sample_names)
        status_pill("ready",
          sprintf("\u2713 Processed: %d unique samples", n_unique))
      }
    })

    # =====================================================================
    # Summary panels
    # =====================================================================
    output$tm_summary <- shiny::renderUI({
      shiny::req(state$tm_processed)
      tm <- state$tm_processed
      rows <- lapply(seq_along(tm$sample_names), function(i) {
        shiny::div(style = "background:#1E2D45;padding:10px;border-radius:6px;margin-bottom:10px;",
          shiny::div(style = "color:#00C2FF;font-weight:bold;margin-bottom:5px;",
                     tm$sample_names[i]),
          shiny::div(style = "color:#7A8FAD;font-size:0.85rem;",
                     sprintf("Condition: %s", tm$concentrations[i])),
          shiny::div(style = "color:#7A8FAD;font-size:0.85rem;",
                     sprintf("Tm: %.2f \u00b1 %.2f \u00b0C",
                             tm$means[i], tm$sems[i])))
      })
      shiny::tagList(rows)
    })

    output$stats_ui <- shiny::renderUI({
      shiny::req(state$dfdt_processed)
      p <- state$dfdt_processed
      rows <- lapply(seq_along(p$sample_names), function(i) {
        sample     <- p$sample_names[i]
        n_reps     <- p$n[1, i]
        peak_idx   <- which.max(p$mean[, i])
        peak_temp  <- p$temperatures[peak_idx]
        peak_value <- p$mean[peak_idx, i]
        shiny::div(style = "background:#1E2D45;padding:10px;border-radius:6px;margin-bottom:10px;",
          shiny::div(style = "color:#00C2FF;font-weight:bold;margin-bottom:5px;",
                     sample),
          shiny::div(style = "color:#7A8FAD;font-size:0.85rem;",
                     sprintf("Replicates: %d", n_reps)),
          shiny::div(style = "color:#7A8FAD;font-size:0.85rem;",
                     sprintf("Peak: %.2f\u00b0C (%.4f)",
                             peak_temp, peak_value)))
      })
      shiny::tagList(rows)
    })

    # =====================================================================
    # Plots
    # =====================================================================
    output$heatmap <- shiny::renderPlot({
      shiny::req(state$dfdt_processed)
      p   <- state$dfdt_processed
      pal <- .contour_palette(input$palette, 100)
      thr <- input$threshold

      df_long <- data.frame(
        Temperature = rep(p$temperatures, times = length(p$sample_names)),
        Sample      = rep(p$sample_names, each  = length(p$temperatures)),
        Norm_dFdT   = .contour_thresh(as.vector(p$norm), thr))
      df_long$Sample <- factor(df_long$Sample, levels = p$sample_names)

      ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Temperature,
                                            fill = Norm_dFdT)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(colours = pal, limits = c(0, 1),
                                      name = "Norm.\ndF/dT") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = NULL, y = "Temperature (\u00b0C)") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          plot.background   = ggplot2::element_rect(fill = "#0B1623", colour = NA),
          panel.background  = ggplot2::element_rect(fill = "#0B1623", colour = NA),
          plot.title        = ggplot2::element_text(colour = CG_PALETTE$txt, face = "bold", size = 13),
          axis.text.x       = ggplot2::element_text(colour = CG_PALETTE$muted, angle = 40, hjust = 1, size = 9),
          axis.text.y       = ggplot2::element_text(colour = CG_PALETTE$muted, size = 9),
          axis.title.y      = ggplot2::element_text(colour = CG_PALETTE$muted, size = 10),
          legend.background = ggplot2::element_rect(fill = "#0B1623", colour = NA),
          legend.text       = ggplot2::element_text(colour = CG_PALETTE$muted, size = 8),
          legend.title      = ggplot2::element_text(colour = CG_PALETTE$muted, size = 9),
          panel.grid        = ggplot2::element_blank(),
          panel.border      = ggplot2::element_rect(colour = CG_PALETTE$muted,
                                                    fill = NA, linewidth = 0.6))
    }, bg = "#0B1623")

    # The standard is now rendered as a SEPARATE side panel, so we just
    # pass through the processed object. No per-sample concentration
    # logic needed - the side panel shows every unique sample in the
    # file as its own tile column in file order.
    .standard_overlay_data <- shiny::reactive({
      state$standard_processed
    })

    output$tm_plot <- shiny::renderPlot({
      shiny::req(state$tm_processed, state$dfdt_processed)
      p_main <- .contour_tm_plot(
        tm           = state$tm_processed,
        processed    = state$dfdt_processed,
        tm_raw       = state$tm_raw,
        tm_pairing   = state$tm_pairing,
        palette_name = input$palette,
        threshold    = input$threshold,
        y_lo         = input$ymin,
        y_hi         = input$ymax,
        dark         = TRUE,
        tm_colour    = input$tm_colour %||% "#06B6D4",
        x_title      = input$x_title %||% "Log [Concentration (\u03bcM)]"
      )
      sp <- .standard_overlay_data()
      if (is.null(sp)) {
        return(p_main)
      }
      p_std <- .contour_standard_plot(
        sp           = sp,
        threshold    = input$threshold,
        y_lo         = input$ymin,
        y_hi         = input$ymax,
        dark         = TRUE,
        panel_colour = input$standard_colour,
        opacity      = input$standard_opacity,
        palette_name = input$palette,
        main_tile_w  = .contour_main_geom(state$tm_processed$concentrations)$tile_w,
        main_x_range = .contour_main_geom(state$tm_processed$concentrations)$x_range,
        layout_width_ratio = 4)
      # Compose the two plots. patchwork is preferred because it
      # panel-aligns the two plots (essential here - the user needs
      # to compare temperatures across both plots, so the panels MUST
      # have identical vertical extents). If patchwork isn't loaded
      # for some reason (install missing, R session needs restart),
      # fall back to gridExtra::grid.arrange so the tool still works,
      # just with slightly imperfect alignment.
      if (requireNamespace("patchwork", quietly = TRUE) &&
          "patchwork" %in% loadedNamespaces()) {
        p_main + p_std + patchwork::plot_layout(widths = c(4, 1))
      } else {
        warning("patchwork not loaded - falling back to gridExtra. ",
                "Install with install.packages('patchwork') and ",
                "restart R for properly aligned panels.")
        gridExtra::grid.arrange(p_main, p_std, ncol = 2, widths = c(4, 1))
      }
    }, bg = "#0B1623")

    # =====================================================================
    # Exports
    # =====================================================================
    output$export_stats <- shiny::downloadHandler(
      filename = function() ts_filename("cpm_contour_statistics", "csv"),
      content  = function(file) {
        shiny::req(state$dfdt_processed)
        p <- state$dfdt_processed

        export_df <- data.frame(Temperature = p$temperatures)
        for (i in seq_along(p$sample_names)) {
          s <- p$sample_names[i]
          export_df[[paste0(s, "_Mean")]] <- p$mean[, i]
          export_df[[paste0(s, "_SEM")]]  <- p$sem[,  i]
        }
        utils::write.csv(export_df, file, row.names = FALSE)

        if (!is.null(state$tm_processed)) {
          tm    <- state$tm_processed
          tm_df <- data.frame(Condition = tm$concentrations,
                              Sample    = tm$sample_names,
                              Tm_Mean   = round(tm$means, 3),
                              Tm_SEM    = round(tm$sems,  3))
          write("",            file, append = TRUE)
          write("# Tm Summary", file, append = TRUE)
          utils::write.table(tm_df, file, sep = ",",
                             row.names = FALSE, append = TRUE)
        }
        shiny::showNotification("\u2713 Statistics exported!",
                                type = "message", duration = 3)
      })

    output$export_heatmap <- shiny::downloadHandler(
      filename = function() ts_filename("cpm_heatmap", "png"),
      content  = function(file) {
        shiny::req(state$dfdt_processed)
        p   <- state$dfdt_processed
        pal <- .contour_palette(shiny::isolate(input$palette), 256)
        thr <- shiny::isolate(input$threshold)

        df_long <- data.frame(
          Temperature = rep(p$temperatures, times = length(p$sample_names)),
          Sample      = rep(p$sample_names, each  = length(p$temperatures)),
          Norm_dFdT   = .contour_thresh(as.vector(p$norm), thr))
        df_long$Sample <- factor(df_long$Sample, levels = p$sample_names)

        plt <- ggplot2::ggplot(df_long,
                               ggplot2::aes(x = Sample, y = Temperature,
                                            fill = Norm_dFdT)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradientn(colours = pal, limits = c(0, 1),
                                        name = "Norm.\ndF/dT") +
          ggplot2::scale_x_discrete(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::labs(x = NULL, y = "Temperature (\u00b0C)") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(
            plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
            panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
            plot.title        = ggplot2::element_text(face = "bold", size = 14),
            axis.text.x       = ggplot2::element_text(angle = 40, hjust = 1, size = 10),
            axis.text.y       = ggplot2::element_text(size = 10),
            axis.title.y      = ggplot2::element_text(size = 11),
            panel.grid        = ggplot2::element_blank(),
            panel.border      = ggplot2::element_rect(colour = "black",
                                                      fill = NA, linewidth = 0.6))
        ggplot2::ggsave(file, plt, width = 10, height = 7, dpi = 300, bg = "white")
        shiny::showNotification("\u2713 Heatmap exported!",
                                type = "message", duration = 3)
      })

    output$export_tmplot <- shiny::downloadHandler(
      filename = function() ts_filename("cpm_tmplot", "png"),
      content  = function(file) {
        shiny::req(state$tm_processed, state$dfdt_processed)
        p_main <- .contour_tm_plot(
          tm           = state$tm_processed,
          processed    = state$dfdt_processed,
          tm_raw       = state$tm_raw,
          tm_pairing   = state$tm_pairing,
          palette_name = shiny::isolate(input$palette),
          threshold    = shiny::isolate(input$threshold),
          y_lo         = shiny::isolate(input$ymin),
          y_hi         = shiny::isolate(input$ymax),
          dark         = FALSE,
          tm_colour    = shiny::isolate(input$tm_colour) %||% "#06B6D4",
          x_title      = shiny::isolate(input$x_title) %||% "Log [Concentration (\u03bcM)]"
        )
        sp <- shiny::isolate(.standard_overlay_data())
        if (is.null(sp)) {
          ggplot2::ggsave(file, p_main, width = 10, height = 7, dpi = 300, bg = "white")
        } else {
          p_std <- .contour_standard_plot(
            sp           = sp,
            threshold    = shiny::isolate(input$threshold),
            y_lo         = shiny::isolate(input$ymin),
            y_hi         = shiny::isolate(input$ymax),
            dark         = FALSE,
            panel_colour = shiny::isolate(input$standard_colour),
            opacity      = shiny::isolate(input$standard_opacity),
            palette_name = shiny::isolate(input$palette),
            main_tile_w  = .contour_main_geom(state$tm_processed$concentrations)$tile_w,
            main_x_range = .contour_main_geom(state$tm_processed$concentrations)$x_range,
            layout_width_ratio = 4)
          # patchwork preferred; gridExtra fallback - see notes in
          # output$tm_plot for the rationale.
          if (requireNamespace("patchwork", quietly = TRUE) &&
              "patchwork" %in% loadedNamespaces()) {
            combined <- p_main + p_std + patchwork::plot_layout(widths = c(4, 1))
          } else {
            combined <- gridExtra::arrangeGrob(p_main, p_std, ncol = 2,
                                               widths = c(4, 1))
          }
          ggplot2::ggsave(file, combined, width = 12.5, height = 7,
                          dpi = 300, bg = "white")
        }
        shiny::showNotification("\u2713 Tm plot exported!",
                                type = "message", duration = 3)
      })

    # ---- Public reactive for cross-tool export --------------------------
    shiny::reactive({
      if (is.null(state$dfdt_processed)) return(NULL)
      list(dfdt_processed = state$dfdt_processed,
           tm_processed   = state$tm_processed)
    })
  })
}


# -- Internal helpers ---------------------------------------------------------

# Resolve palette name to a vector of n colours. Same mapping the original
# used inline in two places; centralising avoids drift.
.contour_palette <- function(name, n = 256) {
  switch(name %||% "grayscale",
    grayscale = grDevices::colorRampPalette(c("#FFFFFF", "#000000"))(n),
    inferno   = scales::viridis_pal(option = "inferno")(n),
    magma     = scales::viridis_pal(option = "magma")(n),
    viridis   = scales::viridis_pal(option = "viridis")(n),
    plasma    = scales::viridis_pal(option = "plasma")(n),
    RdBu      = rev(grDevices::hcl.colors(n, "RdBu")),
    YlOrRd    = grDevices::hcl.colors(n, "YlOrRd"),
    grDevices::colorRampPalette(c("#FFFFFF", "#000000"))(n))
}

# Threshold normalised dF/dT: values below `thr` become 0; values >=thr are
# rescaled so [thr,1] -> [0,1]. This re-uses the full palette for the part
# of the data the user actually cares about.
.contour_thresh <- function(v, thr = 0) {
  if (is.null(thr) || thr <= 0) return(v)
  ifelse(v < thr, 0, (v - thr) / (1 - thr))
}

# Build the combined "Tm vs log[concentration] with heatmap fingerprint"
# plot. Used by both on-screen render and PNG export, with `dark` flipping
# the colour scheme.
# Compute the main Tm plot's tile width (in data units) and the full
# rendered x-axis range (including the 6% expand on each side). Both
# .contour_tm_plot and .contour_standard_plot need these numbers - the
# former to draw its own tiles, the latter to size its tiles so they
# render at the same pixel width.
.contour_main_geom <- function(concentrations) {
  conc_num <- suppressWarnings(as.numeric(concentrations))
  conc_log <- ifelse(is.na(conc_num) | conc_num == 0, 0.001, conc_num)
  log_x    <- log10(conc_log)

  tile_w <- if (length(unique(log_x)) > 1)
              min(diff(sort(unique(log_x)))) * 0.85
            else 0.4

  data_range <- if (length(unique(log_x)) > 1) max(log_x) - min(log_x) else 1
  # expand = expansion(mult = 0.06) adds 6% of the data range on each side
  x_range_full <- data_range * 1.12
  if (x_range_full <= 0) x_range_full <- 1

  list(tile_w = tile_w, x_range = x_range_full, log_x = log_x)
}

.contour_tm_plot <- function(tm, processed, tm_raw, tm_pairing,
                             palette_name, threshold,
                             y_lo, y_hi, dark = TRUE,
                             tm_colour = "#06B6D4",
                             x_title = "Log [Concentration (\u03bcM)]") {

  geom     <- .contour_main_geom(tm$concentrations)
  log_x    <- geom$log_x
  tile_w   <- geom$tile_w

  y_lo <- if (!is.null(y_lo) && !is.na(y_lo)) y_lo
          else min(processed$temperatures)
  y_hi <- if (!is.null(y_hi) && !is.na(y_hi)) y_hi
          else max(processed$temperatures)
  temp_mask <- processed$temperatures >= y_lo & processed$temperatures <= y_hi
  temps_sub <- processed$temperatures[temp_mask]
  temp_step <- if (length(temps_sub) > 1) mean(diff(temps_sub)) else 1.0

  tile_rows <- lapply(seq_along(tm$concentrations), function(i) {
    sname <- tm$sample_names[i]
    sidx  <- which(processed$sample_names == sname)
    if (length(sidx) == 0) return(NULL)
    nv <- .contour_thresh(processed$norm[temp_mask, sidx[1]], threshold)
    data.frame(log_x = log_x[i], Temperature = temps_sub, norm_dFdT = nv)
  })
  tile_df <- do.call(rbind, Filter(Negate(is.null), tile_rows))

  pal <- .contour_palette(palette_name, 256)

  # Optional: faint individual replicate points (dark mode only -
  # the export version omits these for clarity, matching the original)
  rep_df <- NULL
  if (dark && !is.null(tm_raw) && !is.null(tm_pairing)) {
    valid    <- !is.na(tm_pairing) & tm_pairing != ""
    raw_conc <- tm_raw$concentrations[valid]
    raw_vals <- tm_raw$values[valid]
    if (length(raw_conc) > 0) {
      rows <- lapply(seq_along(raw_conc), function(i) {
        cn <- suppressWarnings(as.numeric(raw_conc[i]))
        cx <- if (is.na(cn) || cn == 0) 0.001 else cn
        data.frame(log_x = log10(cx), Tm = raw_vals[[i]])
      })
      rep_df <- do.call(rbind, rows)
    }
  }

  df <- data.frame(log_x = log_x, Mean = tm$means, SEM = tm$sems)
  x_breaks <- seq(floor(min(log_x)), ceiling(max(log_x)))

  # Foreground colours - tm_colour is the user's choice (cyan default)
  # for the line/errorbars/points so they don't blend into the contour.
  # Same colour fills the data point circles too, giving a solid look.
  fg          <- tm_colour
  axis_line   <- if (dark) "#3A4D63" else "grey40"
  bg_colour   <- if (dark) "#0B1623" else "white"
  text_colour <- if (dark) CG_PALETTE$muted else "black"

  # Replicate dots also follow tm_colour, but more faded so they don't
  # compete with the main line
  rep_alpha <- 0.55

  p <- ggplot2::ggplot()
  if (!is.null(tile_df) && nrow(tile_df) > 0) {
    p <- p + ggplot2::geom_tile(
      data = tile_df,
      ggplot2::aes(x = log_x, y = Temperature, fill = norm_dFdT),
      width = tile_w, height = temp_step,
      alpha = if (dark) 0.65 else 1) +
      ggplot2::scale_fill_gradientn(colours = pal, limits = c(0, 1),
                                    name = "Norm.\ndF/dT")
  }
  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    p <- p + ggplot2::geom_point(data = rep_df,
      ggplot2::aes(x = log_x, y = Tm),
      colour = fg, size = 1.8, alpha = rep_alpha, shape = 16,
      inherit.aes = FALSE)
  }
  p <- p +
    ggplot2::geom_line(data = df, ggplot2::aes(x = log_x, y = Mean),
      colour = fg, linewidth = 0.9, inherit.aes = FALSE) +
    ggplot2::geom_errorbar(data = df,
      ggplot2::aes(x = log_x, ymin = Mean - SEM, ymax = Mean + SEM),
      width = 0.06, colour = fg, linewidth = 0.9, inherit.aes = FALSE) +
    ggplot2::geom_point(data = df, ggplot2::aes(x = log_x, y = Mean),
      colour = fg, fill = fg, size = 3.8, shape = 21, stroke = 1.8,
      inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = x_breaks,
                                labels = as.character(x_breaks),
                                expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi),
                                breaks = pretty(c(y_lo, y_hi), n = 8),
                                expand = c(0, 0)) +
    ggplot2::labs(x = x_title,
                  y = "Temperature (\u00b0C)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_colour, colour = NA),
      panel.background  = ggplot2::element_rect(fill = bg_colour, colour = NA),
      plot.title        = ggplot2::element_text(colour = text_colour,
                                                face = "bold", size = 13),
      axis.title        = ggplot2::element_text(colour = text_colour, size = 11),
      axis.text         = ggplot2::element_text(colour = text_colour, size = 10),
      axis.line         = ggplot2::element_line(colour = axis_line, linewidth = 0.6),
      panel.grid        = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = bg_colour, colour = NA),
      legend.text       = ggplot2::element_text(colour = text_colour, size = 8),
      legend.title      = ggplot2::element_text(colour = text_colour, size = 9))

  p
}


# --------------------------------------------------------------------------
# Side panel: standard-sample contour
# --------------------------------------------------------------------------
# Builds a narrow companion plot showing the standard's normalised dF/dT
# as one tile column per unique sample in the standard file, using a
# white->panel_colour ramp (or dark->panel_colour in dark mode). Y-axis
# matches the main plot's range; y-axis labels are suppressed since the
# main plot to the left already has them. Each sample's name appears as
# the x-axis tick label.
.contour_standard_plot <- function(sp, threshold, y_lo, y_hi, dark = TRUE,
                                   panel_colour = "auto", opacity = 0.85,
                                   palette_name = "grayscale",
                                   main_tile_w = NULL, main_x_range = NULL,
                                   layout_width_ratio = 4) {
  if (is.null(sp) || is.null(sp$norm) || ncol(sp$norm) < 1) return(NULL)

  y_lo <- if (!is.null(y_lo) && !is.na(y_lo)) y_lo else min(sp$temperatures)
  y_hi <- if (!is.null(y_hi) && !is.na(y_hi)) y_hi else max(sp$temperatures)
  mask <- sp$temperatures >= y_lo & sp$temperatures <= y_hi
  temps_sub <- sp$temperatures[mask]
  temp_step <- if (length(temps_sub) > 1) mean(diff(temps_sub)) else 1.0

  unique_samples <- sp$sample_names
  # Build long-format data.frame: one row per (sample, temperature)
  tile_rows <- lapply(seq_along(unique_samples), function(i) {
    nv <- .contour_thresh(sp$norm[mask, i], threshold)
    nv[is.na(nv)] <- 0
    data.frame(sample = factor(unique_samples[i], levels = unique_samples),
               x_idx  = i,
               Temperature = temps_sub,
               norm_dFdT   = nv)
  })
  tile_df <- do.call(rbind, Filter(Negate(is.null), tile_rows))

  # Build the colour ramp.
  #   - "auto" (default): use the same multi-stop palette as the main
  #     plot (.contour_palette). This makes the standard panel a direct
  #     visual extension of the main data - identical colour mapping,
  #     just in a separate column.
  #   - hex value: use the user's chosen single-hue ramp from the
  #     background colour up to that colour. Useful when the user wants
  #     the standard panel visually distinct from the main data.
  ramp <- if (identical(panel_colour, "auto") || is.null(panel_colour)) {
    .contour_palette(palette_name, 256)
  } else {
    anchor_colour <- if (dark) "#0B1623" else "white"
    grDevices::colorRampPalette(c(anchor_colour, panel_colour))(256)
  }

  bg_colour   <- if (dark) "#0B1623" else "white"
  axis_line   <- if (dark) "#3A4D63" else "grey40"
  text_colour <- if (dark) CG_PALETTE$muted else "black"

  # The standard plot's x-axis spans the sample indices plus expansion
  # padding on each side (add = 0.6).
  n_std <- length(unique_samples)
  std_x_range <- n_std + 2 * 0.6  # left and right expand

  # Compute the tile width in DATA UNITS so that its rendered PIXEL width
  # matches the main plot's tile pixel width. Reasoning:
  #
  #   main tile pixel width = (main_tile_w / main_x_range) * main_panel_pixels
  #   std panel pixels      = main_panel_pixels / layout_width_ratio
  #     (from patchwork's plot_layout(widths = c(4, 1)), the std panel
  #      gets 1/4 of the main panel's pixel width)
  #
  # For std tile pixels to equal main tile pixels:
  #   std_tile_w / std_x_range = layout_width_ratio * main_tile_w / main_x_range
  #
  # Falls back to a fixed 0.85 if main geometry isn't supplied (e.g. the
  # function is called standalone outside the renderPlot path).
  std_tile_w <- if (!is.null(main_tile_w) && !is.null(main_x_range) &&
                   main_x_range > 0) {
    val <- layout_width_ratio * main_tile_w / main_x_range * std_x_range
    # Clamp - don't let it exceed the natural per-sample spacing (1 unit)
    # which would cause adjacent tiles to overlap in multi-sample standards.
    min(val, 0.92)
  } else {
    0.85
  }

  ggplot2::ggplot(tile_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = x_idx, y = Temperature, fill = norm_dFdT),
      width = std_tile_w, height = temp_step, alpha = opacity) +
    ggplot2::scale_fill_gradientn(colours = ramp, limits = c(0, 1),
                                  guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(unique_samples),
      labels = unique_samples,
      expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi),
                                breaks = pretty(c(y_lo, y_hi), n = 8),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, title = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_colour, colour = NA),
      panel.background  = ggplot2::element_rect(fill = bg_colour, colour = NA),
      # Panel border draws a box around the standard's plotting area,
      # visually separating it from the main plot. Replaces axis.line.x
      # since the border serves as the bottom edge.
      panel.border      = ggplot2::element_rect(colour = axis_line,
                                                fill = NA, linewidth = 0.6),
      plot.title        = ggplot2::element_text(colour = text_colour,
                                                face = "bold", size = 11,
                                                hjust = 0.5),
      axis.title        = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(colour = text_colour, size = 9,
                                                angle = 45, hjust = 1),
      axis.text.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.line.x       = ggplot2::element_blank(),
      axis.line.y       = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(5, 8, 5, 4))
}

# ---- Example data loaders --------------------------------------------------
# Two small CSV files bundled in inst/examples/:
#   - cpm_contour_dfdt_example.csv (~10 KB): full dF/dT matrix with replicates
#   - cpm_contour_tm_example.csv (~150 B):  Tm scatter for the same samples
# Each loader copies its file to tempdir and returns a data.frame matching
# what shiny::fileInput would have produced for a single upload. Caches
# the temp path so revisits don't redo the copy.
.cpm_contour_dfdt_example_cache <- new.env(parent = emptyenv())
.cpm_contour_tm_example_cache   <- new.env(parent = emptyenv())

.cpm_contour_example_loader <- function(cache_env, basename_file, display_name) {
  if (!is.null(cache_env$path) && file.exists(cache_env$path)) {
    return(data.frame(name = display_name, datapath = cache_env$path,
                      stringsAsFactors = FALSE))
  }
  app_dir_local <- if (exists("app_dir", envir = globalenv())) {
    get("app_dir", envir = globalenv())
  } else getwd()
  candidates <- unique(c(
    file.path(app_dir_local, "inst", "examples", basename_file),
    file.path(getwd(),       "inst", "examples", basename_file),
    file.path("inst", "examples", basename_file)
  ))
  src <- candidates[file.exists(candidates)][1]
  if (is.na(src)) return(NULL)
  tryCatch({
    out <- file.path(tempdir(), display_name)
    file.copy(src, out, overwrite = TRUE)
    cache_env$path <- out
    data.frame(name = display_name, datapath = out, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

.cpm_contour_dfdt_example_file <- function() {
  .cpm_contour_example_loader(
    .cpm_contour_dfdt_example_cache,
    "cpm_contour_dfdt_example.csv",
    "NB65_Titre___HsUCP1_Apo_Low_Salt_All_Data_dFdT.csv"
  )
}

.cpm_contour_tm_example_file <- function() {
  .cpm_contour_example_loader(
    .cpm_contour_tm_example_cache,
    "cpm_contour_tm_example.csv",
    "NB65_Titre___HsUCP1_Apo_Low_Salt_Tm__Scatter_.csv"
  )
}
