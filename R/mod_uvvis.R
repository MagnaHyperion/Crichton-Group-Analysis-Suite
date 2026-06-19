################################################################################
#  mod_uvvis.R  --  UV-Vis Spectra + A280 Concentration (Shiny module)
################################################################################
#
#  Two related workflows in one tool:
#
#    1. SPECTRA MODE (multiple samples)
#       Overlay UV-Vis absorbance spectra from a DS-11 / NanoDrop wide-format
#       CSV export. Optional annotation of key-wavelength readings (280, 494,
#       260, 220 nm) using the paired "nm N" / "Abs N" metadata columns.
#
#    2. CONCENTRATION MODE (single sample)
#       Publication-style single-sample plot with a side panel showing the
#       calculated protein concentration:
#           Conc (mg/mL) = A280 x 100 / (E1% x pathlength_mm)
#       Users can override the file's E1% (e.g. to use a corrected value
#       from ProtParam) and see the concentration update live.
#
#  File format assumed: DS-11 CSV export with metadata columns first
#  (Sample Number, Sample Name, Concentration, E1%, A280, Pathlength (mm),
#  260/280, ...) followed by wavelength columns headed by their nanometre
#  value (e.g. "220", "221", ... "350" or "750").
#
################################################################################

# -- UI -----------------------------------------------------------------------
uvvis_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    # Hidden Clear button - fired by the global navbar Clear.
    shiny::div(style = "display: none;",
      shiny::actionButton(ns("clear"), "\U0001f504  Clear",
                          class = "btn-clear")
    ),

    shiny::div(class = "sticky-tool",
      shiny::fluidRow(

        # ============ Left column: workflow ==========================
        shiny::column(4,
          shiny::div(class = "workflow-col",

            lab_card(
              step_title(1, "Upload Data File"),
              info_box(paste("DS-11 / NanoDrop CSV export. First row of",
                             "columns is metadata (Sample Name, E1%, A280,",
                             "Pathlength, 260/280...) followed by",
                             "wavelength columns (220, 221, ...).")),
              shiny::fileInput(ns("file"), NULL,
                accept = c(".csv"),
                buttonLabel = "Browse\u2026",
                placeholder = "DS-11 CSV export"),
              shiny::uiOutput(ns("file_status"))
            ),

            lab_card(
              step_title(2, "Sample Selection"),
              info_box(paste("Select one sample for the concentration",
                             "panel view (recommended for figures).",
                             "Select multiple to overlay spectra.")),
              shiny::uiOutput(ns("sample_selector"))
            ),

            lab_card(
              step_title(3, "Plot Settings"),
              shiny::tags$label("Plot title",
                                style = "color:var(--muted);font-size:0.78rem;"),
              shiny::textInput(ns("plot_title"), NULL,
                value = "", placeholder = "e.g. Post-purification spectra",
                width = "100%"),
              shiny::checkboxInput(ns("show_title"),
                "Show plot title",
                value = FALSE),
              shiny::tags$div(style = "color:var(--muted);font-size:0.75rem;margin-top:-0.5rem;margin-bottom:0.5rem;",
                "Off by default - the sample name already appears in the",
                " concentration panel. Useful for multi-sample overlays",
                " where you want a group label."),
              shiny::tags$label("Wavelength window (nm)",
                                style = "color:var(--muted);font-size:0.78rem;"),
              shiny::fluidRow(
                shiny::column(6, shiny::numericInput(ns("wl_min"), "Min",
                                                     value = 220, min = 190,
                                                     max = 900, step = 5,
                                                     width = "100%")),
                shiny::column(6, shiny::numericInput(ns("wl_max"), "Max",
                                                     value = 350, min = 190,
                                                     max = 900, step = 5,
                                                     width = "100%"))
              ),
              shiny::br(),
              shiny::checkboxInput(ns("show_markers"),
                "Annotate key-wavelength readings (280, 494, 260, 220 nm)",
                value = FALSE),
              shiny::tags$div(style = "color:var(--muted);font-size:0.75rem;margin-top:-0.5rem;",
                "Uses the file's nm 1 / Abs 1 columns to place markers on the spectrum."),
              shiny::br(),
              shiny::tags$label("Line colour",
                                style = "color:var(--muted);font-size:0.78rem;"),
              shiny::selectInput(ns("line_colour"), NULL,
                choices = c(
                  "Blue"    = "#2E5CB8",
                  "Red"     = "#C0392B",
                  "Green"   = "#27AE60",
                  "Purple"  = "#8E44AD",
                  "Orange"  = "#E67E22",
                  "Teal"    = "#16A085",
                  "Charcoal"= "#2C3E50",
                  "Rust"    = "#D35400"),
                selected = "#2E5CB8", width = "100%"),
              shiny::tags$div(style = "color:var(--muted);font-size:0.75rem;margin-top:-0.5rem;",
                "Single-sample mode uses this colour. Multi-sample overlays cycle through a palette starting from it.")
            ),

            lab_card(
              step_title(4, "Concentration Panel"),
              info_box(paste("Available in single-sample mode only.",
                             "Displays A280, E1%, calculated concentration,",
                             "and 260/280 ratio in a side panel.",
                             "Adjust E1% below to recalculate concentration",
                             "using the formula:",
                             "Conc = A280 x 100 / (E1% x pathlength_mm).")),
              shiny::checkboxInput(ns("show_conc_panel"),
                "Show concentration side panel",
                value = TRUE),
              shiny::uiOutput(ns("conc_overrides_ui")),
              shiny::uiOutput(ns("conc_live_display"))
            ),

            lab_card(
              step_title(5, "Analyse"),
              shiny::actionButton(ns("run"), "\u25b6  Refresh Plot",
                                  class = "btn-run"),
              shiny::br(), shiny::br(),
              shiny::uiOutput(ns("download_buttons"))
            )
          )
        ),

        # ============ Right column: sticky preview =====================
        shiny::column(8,
          shiny::div(class = "preview-col",
            lab_card(
              shiny::div(class = "lab-card-title",
                         "\U0001f4c8  Absorbance Spectrum"),
              shiny::plotOutput(ns("plot"), height = "520px")
            )
          )
        )
      )
    )
  )
}


# -- Server -------------------------------------------------------------------
uvvis_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- State ---------------------------------------------------------
    current_file <- shiny::reactiveVal(NULL)
    parsed_data  <- shiny::reactiveVal(NULL)

    # ---- Internal: clear -----------------------------------------------
    .clear_state <- function() {
      parsed_data(NULL)
    }

    shiny::observeEvent(input$clear, {
      current_file(NULL)
      .clear_state()
      tryCatch(.load_example(), error = function(e) NULL)
      shiny::showNotification("UV-Vis data cleared",
                              type = "message", duration = 2)
    })

    # ---- File upload → current_file → parse ----------------------------
    shiny::observeEvent(input$file, {
      .clear_state()
      current_file(input$file)
    }, ignoreInit = TRUE)

    shiny::observeEvent(current_file(), {
      cf <- current_file()
      shiny::req(cf)
      tryCatch({
        parsed <- .parse_uvvis_csv(cf$datapath)
        parsed_data(parsed)
        shiny::showNotification(
          sprintf("\u2713 Loaded: %d sample(s), wavelengths %g-%g nm",
                  nrow(parsed$meta),
                  min(parsed$wavelengths), max(parsed$wavelengths)),
          type = "message", duration = 3)
      }, error = function(e) {
        shiny::showNotification(paste("Error reading file:", e$message),
                                type = "error", duration = 6)
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Session start: bundle example --------------------------------
    .load_example <- function() {
      cf <- .uvvis_example_file()
      if (!is.null(cf)) current_file(cf)
    }
    .example_loader_obs <- shiny::observe({
      .example_loader_obs$destroy()
      tryCatch(.load_example(), error = function(e)
        message("[UV-Vis] example load failed: ", conditionMessage(e)))
    })

    # ---- File status pill ---------------------------------------------
    output$file_status <- shiny::renderUI({
      p <- parsed_data()
      if (is.null(p)) return(NULL)
      status_pill("ready", sprintf("\u2713 Loaded: %d sample(s), %d wavelengths",
                                    nrow(p$meta), length(p$wavelengths)))
    })

    # ---- Sample selection area (checkboxes + renames in one output) ----
    #
    # Renders ONCE per file. Everything (checkboxes, toggle button,
    # rename text inputs) sits in a single uiOutput that depends on
    # parsed_data() alone - never on applied_renames. That way clicking
    # Refresh Plot doesn't destroy the rename text inputs.
    #
    # Rename application is button-gated via `applied_renames`. On
    # Refresh, we sweep the current input$rename_* values into
    # applied_renames AND call updateCheckboxGroupInput to refresh the
    # checkbox labels in place (without tearing them down).
    applied_renames <- shiny::reactiveVal(list())

    shiny::observeEvent(parsed_data(), {
      applied_renames(list())
    }, ignoreNULL = TRUE)

    output$sample_selector <- shiny::renderUI({
      p <- parsed_data()
      if (is.null(p)) {
        return(shiny::tags$div(style = "color:var(--muted);font-size:0.85rem;",
                               "Upload a file to see available samples."))
      }
      originals <- p$meta$sample_label

      # Build the rename text inputs. Wrap in an explicit div so Shiny
      # sees a proper container for the input bindings rather than a
      # bare list.
      rename_boxes <- lapply(seq_along(originals), function(i) {
        shiny::textInput(ns(paste0("rename_", i)),
          label       = originals[i],
          value       = "",
          placeholder = originals[i],
          width       = "100%")
      })
      # Prepend the intro text as an actual child so the div holds
      # info-then-inputs as one flat set of children
      rename_children <- c(
        list(shiny::div(
          style = "font-size:0.75rem;color:var(--muted);margin-bottom:0.5rem;",
          "Leave blank to keep the original name.")),
        rename_boxes
      )

      shiny::tagList(
        shiny::checkboxGroupInput(ns("selected_samples"), NULL,
          choices  = stats::setNames(originals, originals),  # labels = originals initially
          selected = originals[1]),
        shiny::tags$button("\u270f  Rename samples", class = "adv-toggle",
          onclick = sprintf("$('#%s').slideToggle(200)", ns("rename_panel"))),
        shiny::div(id = ns("rename_panel"), style = "display:none;",
          do.call(shiny::div, c(list(class = "settings-group"),
                                rename_children)))
      )
    })

    # Refresh Plot: capture current text-input values into
    # applied_renames AND live-update the checkbox labels so the user
    # sees their new names appear on the checkboxes without tearing
    # down the DOM (which would destroy the rename text inputs).
    shiny::observeEvent(input$run, {
      p <- parsed_data()
      if (is.null(p)) return()
      new_map <- list()
      for (i in seq_along(p$meta$sample_label)) {
        v <- input[[paste0("rename_", i)]]
        if (!is.null(v) && nzchar(trimws(v))) {
          new_map[[p$meta$sample_label[i]]] <- trimws(v)
        }
      }
      applied_renames(new_map)

      # Refresh the checkbox labels in place
      originals <- p$meta$sample_label
      displays  <- vapply(originals, function(o) {
        if (!is.null(new_map[[o]])) new_map[[o]] else o
      }, character(1), USE.NAMES = FALSE)
      shiny::updateCheckboxGroupInput(session, "selected_samples",
        choices  = stats::setNames(originals, displays),
        selected = shiny::isolate(input$selected_samples))
    })

    # Reactive mapping from original -> display name, sourced from
    # applied_renames (which only updates on Refresh Plot).
    .display_names <- shiny::reactive({
      p <- parsed_data()
      if (is.null(p)) return(character(0))
      ar <- applied_renames()
      vapply(p$meta$sample_label, function(orig) {
        v <- ar[[orig]]
        if (is.null(v) || !nzchar(v)) orig else v
      }, character(1), USE.NAMES = FALSE)
    })

    .display_for <- function(original) {
      p <- parsed_data()
      if (is.null(p)) return(original)
      idx <- match(original, p$meta$sample_label)
      if (is.na(idx)) return(original)
      dn <- .display_names()
      if (length(dn) < idx) return(original)
      dn[idx]
    }

    # ---- Concentration overrides UI (single-sample mode only) ---------
    # E1% and pathlength override boxes are only meaningful when exactly
    # one sample is selected. If the user has picked multiple samples we
    # collapse this section into a hint explaining why.
    output$conc_overrides_ui <- shiny::renderUI({
      p   <- parsed_data()
      sel <- input$selected_samples
      if (is.null(p) || length(sel) != 1) {
        return(shiny::tags$div(style = "color:var(--muted);font-size:0.78rem;margin-top:0.5rem;",
          "E1% override is available only when exactly one sample is selected."))
      }
      idx <- match(sel, p$meta$sample_label)
      shiny::tagList(
        shiny::tags$hr(style = "border-color:#2A3B52;margin:0.6rem 0;"),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput(ns("e1_override"),
            "E1% override",
            value = p$meta$e1[idx], min = 0.1, step = 0.1, width = "100%")),
          shiny::column(6, shiny::numericInput(ns("path_override"),
            "Pathlength (mm)",
            value = p$meta$pathlength[idx], min = 0.1, step = 0.1,
            width = "100%"))
        )
      )
    })

    # ---- Live concentration display -----------------------------------
    output$conc_live_display <- shiny::renderUI({
      cv <- .conc_values()
      if (is.null(cv)) return(NULL)
      shiny::tags$div(class = "settings-group", style = "margin-top:0.5rem;",
        shiny::tags$div(class = "settings-group-title", "Live values"),
        shiny::tags$div(
          style = "display:grid;grid-template-columns:1fr auto;gap:0.3rem 0.8rem;font-size:0.85rem;",
          shiny::tags$span("A280:", style = "color:var(--muted);"),
          shiny::tags$span(sprintf("%.3f", cv$a280)),
          shiny::tags$span("E1%:", style = "color:var(--muted);"),
          shiny::tags$span(sprintf("%.2f", cv$e1)),
          shiny::tags$span("Pathlength:", style = "color:var(--muted);"),
          shiny::tags$span(sprintf("%.1f mm", cv$pathlength)),
          shiny::tags$span("Conc:", style = "color:var(--muted);font-weight:bold;"),
          shiny::tags$span(sprintf("%.3f mg/mL", cv$conc),
                           style = "font-weight:bold;color:var(--accent);"),
          shiny::tags$span("260/280:", style = "color:var(--muted);"),
          shiny::tags$span(sprintf("%.3f", cv$ratio_260_280))
        )
      )
    })

    # ---- Concentration values reactive --------------------------------
    # Combines file-derived defaults with user overrides. Returns NULL
    # unless in single-sample mode with valid data.
    .conc_values <- shiny::reactive({
      p   <- parsed_data()
      sel <- input$selected_samples
      if (is.null(p) || length(sel) != 1) return(NULL)
      idx <- match(sel, p$meta$sample_label)
      m   <- p$meta[idx, , drop = FALSE]

      # Prefer user overrides when they're valid numbers; fall back to
      # the file's stored value otherwise
      e1        <- .num_or_default(input$e1_override,    m$e1)
      pathlen   <- .num_or_default(input$path_override,  m$pathlength)
      a280      <- m$a280
      ratio     <- m$ratio_260_280

      conc <- if (is.finite(a280) && is.finite(e1) && is.finite(pathlen) &&
                  e1 > 0 && pathlen > 0) {
        a280 * 100 / (e1 * pathlen)
      } else {
        NA_real_
      }
      list(sample_name = .display_for(m$sample_label),
           sample_number = m$sample_number,
           a280 = a280, e1 = e1, pathlength = pathlen,
           conc = conc, ratio_260_280 = ratio)
    })

    # ---- Reactive plot builder ----------------------------------------
    # Depends on parsed data, selected samples, and all plot settings.
    # Extracted into its own reactive so the render + export share exactly
    # the same plot object.
    .current_plot <- shiny::reactive({
      p   <- parsed_data()
      sel <- input$selected_samples
      shiny::req(p, sel)

      long <- .build_long_df(p, sel,
                             wl_min = input$wl_min, wl_max = input$wl_max)
      if (nrow(long) == 0) return(NULL)

      # Translate the Sample factor's values (originals) to their current
      # display names so the legend, plot title, and downstream renders
      # reflect any renames the user has entered. Keep the SAME order
      # the user selected in - selected_samples is already ordered by
      # click order.
      display_map <- vapply(sel, .display_for, character(1))
      long$Sample <- factor(display_map[as.character(long$Sample)],
                            levels = unique(display_map))

      show_panel <- isTRUE(input$show_conc_panel) && length(sel) == 1
      cv         <- if (show_panel) .conc_values() else NULL

      .render_uvvis_figure(
        long           = long,
        meta_selected  = p$meta[p$meta$sample_label %in% sel, , drop = FALSE],
        key_points     = if (isTRUE(input$show_markers)) p$key_points else NULL,
        wl_min         = input$wl_min,
        wl_max         = input$wl_max,
        base_colour    = input$line_colour %||% "#2E5CB8",
        show_panel     = show_panel,
        conc_values    = cv,
        single_sample  = length(sel) == 1,
        plot_title     = input$plot_title,
        show_title     = isTRUE(input$show_title)
      )
    })

    # ---- Plot render --------------------------------------------------
    output$plot <- shiny::renderPlot({
      shiny::req(parsed_data())
      pl <- .current_plot()
      if (is.null(pl)) return(NULL)
      pl
    }, bg = "white")

    # ---- Download buttons ---------------------------------------------
    output$download_buttons <- shiny::renderUI({
      shiny::req(parsed_data(), input$selected_samples)
      shiny::tagList(
        shiny::downloadButton(ns("download_png"), "\u2b07  PNG",
                              class = "btn-primary"),
        shiny::downloadButton(ns("download_csv"), "\u2b07  CSV")
      )
    })

    output$download_png <- shiny::downloadHandler(
      filename = function() ts_filename("uvvis_spectrum", "png"),
      content  = function(file) {
        pl <- shiny::isolate(.current_plot())
        shiny::req(pl)
        # Wider export when the concentration panel is displayed
        w <- if (isTRUE(shiny::isolate(input$show_conc_panel)) &&
                 length(shiny::isolate(input$selected_samples)) == 1) 11 else 9
        ggplot2::ggsave(file, pl, width = w, height = 6, dpi = 300,
                        bg = "white")
        shiny::showNotification("\u2713 PNG exported",
                                type = "message", duration = 3)
      })

    output$download_csv <- shiny::downloadHandler(
      filename = function() ts_filename("uvvis_data", "csv"),
      content  = function(file) {
        p   <- shiny::isolate(parsed_data())
        sel <- shiny::isolate(input$selected_samples)
        shiny::req(p, sel)
        long <- .build_long_df(p, sel, wl_min = shiny::isolate(input$wl_min),
                                       wl_max = shiny::isolate(input$wl_max))
        # Apply any user renames so the exported CSV matches what the
        # plot legend showed. Isolate the display_names lookup so it
        # doesn't create a reactive dependency in the handler.
        display_map <- shiny::isolate({
          vapply(sel, .display_for, character(1))
        })
        long$Sample <- factor(display_map[as.character(long$Sample)],
                              levels = unique(display_map))
        utils::write.csv(long, file, row.names = FALSE)
        shiny::showNotification("\u2713 CSV exported",
                                type = "message", duration = 3)
      })

    # Return an empty reactive - this module doesn't participate in the
    # cross-tool export bundle (which was removed anyway, but the app
    # server still calls uvvis_server(id) unconditionally).
    shiny::reactive(NULL)
  })
}


# =============================================================================
# Parsing & data helpers  (module-private, prefixed with a dot)
# =============================================================================

# Parse a DS-11 / NanoDrop CSV. Auto-detects the wavelength columns by
# treating any column whose header is a bare integer as a wavelength.
# Returns a list with:
#   $meta        - data.frame, one row per sample: sample_label, sample_number,
#                  a280, e1, pathlength, ratio_260_280, concentration_reported
#   $wavelengths - numeric vector
#   $abs_matrix  - matrix (samples x wavelengths) of absorbance values
#   $key_points  - data.frame of (sample_label, Wavelength, Absorbance) rows
#                  from the "nm N" / "Abs N" paired metadata columns
.parse_uvvis_csv <- function(path) {
  raw <- utils::read.csv(path, check.names = FALSE,
                         stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
  col_names <- names(raw)

  # Wavelength columns
  is_wl   <- grepl("^[0-9]+$", col_names)
  wl_cols <- col_names[is_wl]
  if (length(wl_cols) == 0)
    stop("No wavelength columns found (expected numeric headers like 220, 221...).")
  wavelengths <- as.numeric(wl_cols)

  # Absorbance matrix
  abs_matrix <- as.matrix(raw[, wl_cols, drop = FALSE])
  storage.mode(abs_matrix) <- "double"

  # Sample labels
  sample_label <- if ("Sample Name" %in% col_names) {
    ifelse(nzchar(trimws(raw[["Sample Name"]])),
           trimws(raw[["Sample Name"]]),
           paste("Sample", seq_len(nrow(raw))))
  } else {
    paste("Sample", seq_len(nrow(raw)))
  }

  # Duplicate labels are unusual but possible on the DS-11 - append
  # (2), (3) etc so the sample selector doesn't render ambiguous entries
  if (anyDuplicated(sample_label) > 0) {
    sample_label <- make.unique(sample_label, sep = " ")
  }

  # Sample metadata
  get_num <- function(col) {
    if (col %in% col_names) suppressWarnings(as.numeric(raw[[col]]))
    else                    rep(NA_real_, nrow(raw))
  }
  meta <- data.frame(
    sample_label           = sample_label,
    sample_number          = if ("Sample Number" %in% col_names)
                               raw[["Sample Number"]] else seq_len(nrow(raw)),
    a280                   = get_num("A280"),
    e1                     = get_num("E1%"),
    pathlength             = {
      pl <- get_num("Pathlength (mm)")
      ifelse(is.na(pl) | pl <= 0, 10, pl)  # default 10 mm if missing
    },
    ratio_260_280          = get_num("260/280"),
    concentration_reported = get_num("Concentration"),
    stringsAsFactors       = FALSE
  )

  # Key points from paired "nm N" / "Abs N" columns
  nm_cols  <- grep("^nm [0-9]+$",  col_names, value = TRUE)
  abs_hits <- grep("^Abs [0-9]+$", col_names, value = TRUE)
  key_points <- NULL
  if (length(nm_cols) > 0 && length(abs_hits) > 0) {
    nm_idx  <- as.integer(sub("nm ",  "", nm_cols))
    abs_idx <- as.integer(sub("Abs ", "", abs_hits))
    shared  <- intersect(nm_idx, abs_idx)
    rows <- lapply(seq_len(nrow(raw)), function(r) {
      do.call(rbind, lapply(shared, function(k) {
        wl <- suppressWarnings(as.numeric(raw[[paste0("nm ",  k)]][r]))
        ab <- suppressWarnings(as.numeric(raw[[paste0("Abs ", k)]][r]))
        if (is.na(wl) || is.na(ab)) return(NULL)
        data.frame(sample_label = sample_label[r],
                   Wavelength   = wl, Absorbance = ab,
                   stringsAsFactors = FALSE)
      }))
    })
    key_points <- do.call(rbind, Filter(Negate(is.null), rows))
  }

  list(meta = meta, wavelengths = wavelengths, abs_matrix = abs_matrix,
       key_points = key_points)
}

# Build a long-format data.frame for plotting: one row per (sample, wl)
# with sample as an ordered factor matching input selection order.
.build_long_df <- function(parsed, selected_samples, wl_min, wl_max) {
  meta_idx <- match(selected_samples, parsed$meta$sample_label)
  meta_idx <- meta_idx[!is.na(meta_idx)]
  if (length(meta_idx) == 0) return(data.frame())

  wl_mask <- rep(TRUE, length(parsed$wavelengths))
  if (!is.null(wl_min) && !is.na(wl_min))
    wl_mask <- wl_mask & parsed$wavelengths >= wl_min
  if (!is.null(wl_max) && !is.na(wl_max))
    wl_mask <- wl_mask & parsed$wavelengths <= wl_max

  wl_kept <- parsed$wavelengths[wl_mask]
  rows <- lapply(meta_idx, function(i) {
    data.frame(Sample     = parsed$meta$sample_label[i],
               Wavelength = wl_kept,
               Absorbance = as.numeric(parsed$abs_matrix[i, wl_mask]),
               stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, rows)
  long$Sample <- factor(long$Sample, levels = selected_samples)
  long
}

# Helpers -----------------------------------------------------------------
.num_or_default <- function(x, default) {
  if (is.null(x) || is.na(x) || !is.finite(as.numeric(x))) default
  else as.numeric(x)
}


# =============================================================================
# Figure builder
# =============================================================================
#
# Produces either:
#   - A single ggplot (spectrum only, multi- or single-sample), or
#   - A patchwork composition of a spectrum + concentration side panel
#     (single-sample mode with panel toggled on)
#
.render_uvvis_figure <- function(long, meta_selected, key_points,
                                 wl_min, wl_max, base_colour,
                                 show_panel, conc_values,
                                 single_sample,
                                 plot_title = "", show_title = FALSE) {

  # Palette: rotate through 8 known colours if multi-sample, otherwise
  # use just the user's picked colour
  palette_all <- c("#2E5CB8", "#C0392B", "#27AE60", "#8E44AD",
                   "#E67E22", "#16A085", "#2C3E50", "#D35400")
  n_samples <- length(unique(long$Sample))
  palette_use <- if (single_sample) {
    base_colour
  } else {
    # Put user's chosen colour first, then cycle the rest
    others <- setdiff(palette_all, base_colour)
    c(base_colour, others)[seq_len(n_samples)]
  }

  # -- Build spectrum plot ------------------------------------------------
  p <- ggplot2::ggplot(long,
    ggplot2::aes(x = Wavelength, y = Absorbance, colour = Sample)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = palette_use)

  # Optional key-wavelength markers (in-window only)
  if (!is.null(key_points) && nrow(key_points) > 0) {
    kp <- key_points[key_points$sample_label %in% levels(long$Sample), ,
                     drop = FALSE]
    if (!is.null(wl_min)) kp <- kp[kp$Wavelength >= wl_min, , drop = FALSE]
    if (!is.null(wl_max)) kp <- kp[kp$Wavelength <= wl_max, , drop = FALSE]
    if (nrow(kp) > 0) {
      kp$Sample <- factor(kp$sample_label, levels = levels(long$Sample))
      p <- p +
        ggplot2::geom_point(data = kp,
          ggplot2::aes(x = Wavelength, y = Absorbance, colour = Sample),
          size = 2.4, show.legend = FALSE, inherit.aes = FALSE)
      lbl <- sprintf("%g nm\nA=%.3f", kp$Wavelength, kp$Absorbance)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(data = kp,
          ggplot2::aes(x = Wavelength, y = Absorbance, label = lbl,
                       colour = Sample),
          size = 2.9, fontface = "bold", lineheight = 0.85,
          seed = 1, box.padding = 0.5, point.padding = 0.3,
          min.segment.length = 0, segment.size = 0.3,
          segment.colour = "grey55", show.legend = FALSE,
          inherit.aes = FALSE)
      } else {
        p <- p + ggplot2::geom_text(data = kp,
          ggplot2::aes(x = Wavelength, y = Absorbance, label = lbl,
                       colour = Sample),
          vjust = -0.5, size = 2.9, fontface = "bold", lineheight = 0.85,
          show.legend = FALSE, inherit.aes = FALSE)
      }
    }
  }

  # Y-axis label wording: use "Absorbance" when we can't derive a
  # confident pathlength; if all selected samples share a pathlength and
  # it's the classic 10 mm, use the reference figure's "10 mm Absorbance"
  # form
  y_axis_label <- if (!is.null(meta_selected) && nrow(meta_selected) > 0) {
    pl <- unique(meta_selected$pathlength)
    if (length(pl) == 1 && isTRUE(all.equal(pl, 10)))      "10 mm Absorbance"
    else if (length(pl) == 1)  sprintf("%g mm Absorbance", pl)
    else                        "Absorbance (a.u.)"
  } else "Absorbance (a.u.)"

  # X-axis breaks: sensible density based on span
  span <- (wl_max %||% max(long$Wavelength)) - (wl_min %||% min(long$Wavelength))
  x_step <- if (span <= 150) 20 else if (span <= 300) 50 else 100
  x_breaks <- seq(0, 1000, by = x_step)
  x_breaks <- x_breaks[x_breaks >= (wl_min %||% min(long$Wavelength)) &
                       x_breaks <= (wl_max %||% max(long$Wavelength))]

  p <- p +
    ggplot2::scale_x_continuous(breaks = x_breaks,
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(x = "Wavelength",  y = y_axis_label, colour = NULL) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title            = ggplot2::element_text(face = "bold", size = 15,
                                                    hjust = 0.5),
      axis.title            = ggplot2::element_text(face = "bold"),
      axis.line             = ggplot2::element_line(linewidth = 0.5),
      axis.ticks            = ggplot2::element_line(linewidth = 0.5),
      panel.grid            = ggplot2::element_blank())

  # Legend placement: overlay for multi-sample plots, suppress for
  # single-sample (each single-sample workflow has the sample name in
  # the concentration panel already, so the legend would be noise)
  if (single_sample) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(
      legend.position      = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background    = ggplot2::element_rect(fill = "white",
                                                    colour = "grey80"),
      legend.key           = ggplot2::element_blank())
  }

  # User-controlled title: applies whenever `show_title` is on and a
  # non-blank string was provided. Alignment stays centred (set in the
  # base theme). Without a title, the plot area starts higher and the
  # spectrum uses more vertical space.
  if (isTRUE(show_title) && !is.null(plot_title) &&
      nzchar(trimws(plot_title))) {
    p <- p + ggplot2::labs(title = plot_title)
  }

  # -- Compose with concentration panel if requested ---------------------
  if (isTRUE(show_panel) && !is.null(conc_values)) {
    p_conc <- .conc_side_panel(conc_values)
    if (requireNamespace("patchwork", quietly = TRUE) &&
        "patchwork" %in% loadedNamespaces()) {
      return(p + p_conc + patchwork::plot_layout(widths = c(3, 1)))
    } else {
      # gridExtra fallback (unaligned - but functional)
      return(gridExtra::arrangeGrob(p, p_conc, ncol = 2, widths = c(3, 1)))
    }
  }
  p
}

# =============================================================================
# Bundled example file resolver
# =============================================================================
# Same defensive path-resolution pattern as .bca_example_file: try the
# global app_dir first, then getwd(), then a relative path.
.uvvis_example_cache <- new.env(parent = emptyenv())

.uvvis_example_file <- function() {
  if (!is.null(.uvvis_example_cache$path) &&
      file.exists(.uvvis_example_cache$path)) {
    return(data.frame(name     = basename(.uvvis_example_cache$path),
                      datapath = .uvvis_example_cache$path,
                      stringsAsFactors = FALSE))
  }

  app_dir_local <- if (exists("app_dir", envir = globalenv())) {
    get("app_dir", envir = globalenv())
  } else getwd()

  candidates <- unique(c(
    file.path(app_dir_local, "inst", "examples", "uvvis_a280_example.csv"),
    file.path(getwd(),       "inst", "examples", "uvvis_a280_example.csv"),
    file.path("inst", "examples", "uvvis_a280_example.csv")
  ))
  src_path <- candidates[file.exists(candidates)][1]
  if (is.na(src_path)) return(NULL)

  tryCatch({
    out_path <- file.path(tempdir(), "uvvis_a280_example.csv")
    file.copy(src_path, out_path, overwrite = TRUE)
    .uvvis_example_cache$path <- out_path
    data.frame(name = "uvvis_a280_example.csv",
               datapath = out_path,
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
.conc_side_panel <- function(cv) {
  # Build the 4-row concentration side panel (Sample header + A280 / E1% /
  # Conc / 260-280 rows), styled to match the reference figure.
  header_bg     <- "#1F5D8C"
  label_bg      <- "#B7DBEE"
  cell_bg       <- "#FFFFFF"
  border_col    <- "#1F5D8C"

  # 4 data rows + 1 header row = 5 total. Rows numbered top to bottom
  # (row 5 at top for header, row 1 at bottom).
  # Use fixed cell heights and 2 columns (label / value).
  rows <- data.frame(
    y      = c(5, 4, 3, 2, 1),
    label  = c("", "A280", "E1%", "Conc.\n(mg/mL)", "260/280"),
    value  = c(cv$sample_name,
               sprintf("%.2f", cv$a280),
               sprintf("%.2f", cv$e1),
               sprintf("%.2f", cv$conc),
               sprintf("%.2f", cv$ratio_260_280)),
    header = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  # Draw header as full-width single cell; other rows as 2 cells (label,
  # value). Use geom_rect + geom_text so we get exact control.
  # x-range: labels in [0, 0.5], values in [0.5, 1]; header spans [0, 1].
  cell_h  <- 1
  label_w <- 0.5

  rects <- rbind(
    data.frame(xmin = 0,       xmax = 1, ymin = 5 - cell_h/2, ymax = 5 + cell_h/2,
               fill = header_bg, stringsAsFactors = FALSE),
    # Label cells (rows 1-4)
    data.frame(xmin = 0,       xmax = label_w,
               ymin = c(4,3,2,1) - cell_h/2,
               ymax = c(4,3,2,1) + cell_h/2,
               fill = label_bg, stringsAsFactors = FALSE),
    # Value cells (rows 1-4)
    data.frame(xmin = label_w, xmax = 1,
               ymin = c(4,3,2,1) - cell_h/2,
               ymax = c(4,3,2,1) + cell_h/2,
               fill = cell_bg, stringsAsFactors = FALSE)
  )

  # Text: labels white on blue (header + label cells), values black
  # (value cells)
  text_df <- rbind(
    data.frame(x = 0.5, y = 5, label = rows$value[1],
               colour = "white", face = "bold", size = 4.5,
               stringsAsFactors = FALSE),
    data.frame(x = 0.25, y = 4:1, label = rows$label[2:5],
               colour = header_bg, face = "bold", size = 3.8,
               stringsAsFactors = FALSE),
    data.frame(x = 0.75, y = 4:1, label = rows$value[2:5],
               colour = "black", face = "plain", size = 4.2,
               stringsAsFactors = FALSE)
  )

  ggplot2::ggplot() +
    ggplot2::geom_rect(data = rects,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   fill = fill),
      colour = border_col, linewidth = 0.4) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(data = text_df,
      ggplot2::aes(x = x, y = y, label = label, colour = colour,
                   fontface = face, size = size),
      lineheight = 0.85) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 8, 20, 8))
}
