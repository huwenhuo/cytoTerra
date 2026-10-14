# Internal helpers for the interactive Shiny gating workspace.
#'
#' @noRd


# Return named selectInput choices:
# display text = marker description [raw FCS channel name]
# selected value = raw FCS channel name
cytoTerra_workspace_channel_choices = function(
    flow_frame,
    available_channels = colnames(flowCore::exprs(flow_frame))
) {

  parameter_data <- tryCatch(
    flowCore::pData(flowCore::parameters(flow_frame)),
    error = function(e) NULL
  )

  if (is.null(parameter_data) ||
      nrow(parameter_data) == 0L) {
    return(stats::setNames(
      available_channels,
      available_channels
    ))
  }

  # FlowCore ordinarily uses `name` and `desc`; support conventional
  # FCS-style names too if they are present instead.
  name_column <- intersect(
    c("name", "$PnN"),
    colnames(parameter_data)
  )

  desc_column <- intersect(
    c("desc", "$PnS"),
    colnames(parameter_data)
  )

  if (length(name_column) == 0L) {
    return(stats::setNames(
      available_channels,
      available_channels
    ))
  }

  parameter_names <- as.character(
    parameter_data[[name_column[[1L]]]]
  )

  parameter_descriptions <- if (length(desc_column) > 0L) {
    as.character(parameter_data[[desc_column[[1L]]]])
  } else {
    rep(NA_character_, length(parameter_names))
  }

  parameter_index <- match(
    available_channels,
    parameter_names
  )

  channel_descriptions <- parameter_descriptions[
    parameter_index
  ]

  display_labels <- vapply(
    seq_along(available_channels),
    function(i) {
      channel_name <- available_channels[[i]]
      channel_description <- channel_descriptions[[i]]

      if (!is.na(channel_description) &&
          nzchar(trimws(channel_description))) {
        paste0(
          channel_description,
          " [",
          channel_name,
          "]"
        )
      } else {
        channel_name
      }
    },
    character(1L)
  )

  stats::setNames(
    available_channels,
    display_labels
  )
}


cytoTerra_workspace_transform_values = function(
    values,
    transform,
    cofactor
) {
  switch(
    transform,
    asinh = asinh(values / cofactor),
    log1p = log1p(values),
    none = values
  )
}


cytoTerra_workspace_transform_label = function(
    transform,
    cofactor
) {
  switch(
    transform,
    asinh = paste0(
      "asinh(x / ",
      signif(cofactor, 6),
      ")"
    ),
    log1p = "log1p(x)",
    none = "raw FCS scale"
  )
}


# Extract the newest Plotly-drawn SVG path from a relayout payload.
cytoTerra_workspace_get_path_from_relayout = function(relayout) {

  if (is.null(relayout) || length(relayout) == 0L) {
    return(NULL)
  }

  # Direct payload from Plotly, commonly:
  # list("shapes[0].path" = "M ... L ... Z")
  relayout_names <- names(relayout)

  if (!is.null(relayout_names)) {
    path_names <- grep(
      "^shapes\\[[0-9]+\\]\\.path$",
      relayout_names,
      value = TRUE
    )

    if (length(path_names) > 0L) {
      path <- relayout[[path_names[[length(path_names)]]]]

      if (is.character(path) &&
          length(path) == 1L &&
          !is.na(path) &&
          nzchar(path)) {
        return(path)
      }
    }
  }

  # A complete replacement shape collection may be sent instead.
  if (is.list(relayout) &&
      !is.null(relayout$shapes)) {

    shapes <- relayout$shapes

    # A single shape can occasionally arrive as a named list rather
    # than a list of shape lists.
    if (is.list(shapes) &&
        !is.null(shapes$path) &&
        is.character(shapes$path) &&
        length(shapes$path) == 1L &&
        !is.na(shapes$path) &&
        nzchar(shapes$path)) {
      return(shapes$path)
    }

    if (is.list(shapes)) {
      for (i in rev(seq_along(shapes))) {
        shape <- shapes[[i]]

        if (is.list(shape) &&
            !is.null(shape$path) &&
            is.character(shape$path) &&
            length(shape$path) == 1L &&
            !is.na(shape$path) &&
            nzchar(shape$path)) {
          return(shape$path)
        }
      }
    }
  }

  NULL
}

# Parse a Plotly SVG polygon path, for example:
# M 1,2 L 3,4 L 5,6 Z
cytoTerra_workspace_polygon_from_path = function(path) {

  if (is.null(path) ||
      !is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)) {
    return(NULL)
  }

  number_matches <- regmatches(
    path,
    gregexpr(
      "[-+]?(?:\\d*\\.?\\d+)(?:[eE][-+]?\\d+)?",
      path,
      perl = TRUE
    )
  )[[1L]]

  if (length(number_matches) < 6L ||
      length(number_matches) %% 2L != 0L) {
    return(NULL)
  }

  values <- suppressWarnings(
    as.numeric(number_matches)
  )

  if (any(!is.finite(values))) {
    return(NULL)
  }

  coordinates <- data.frame(
    x = values[seq(1L, length(values), by = 2L)],
    y = values[seq(2L, length(values), by = 2L)]
  )

  # Plotly sometimes repeats the first vertex as the final vertex.
  if (nrow(coordinates) >= 4L &&
      identical(
        coordinates$x[[1L]],
        coordinates$x[[nrow(coordinates)]]
      ) &&
      identical(
        coordinates$y[[1L]],
        coordinates$y[[nrow(coordinates)]]
      )) {
    coordinates <- coordinates[
      -nrow(coordinates),
      ,
      drop = FALSE
    ]
  }

  if (nrow(coordinates) < 3L) {
    return(NULL)
  }

  coordinates
}


cytoTerra_workspace_gate_tree_text = function(tree) {

  if (is.null(tree) || nrow(tree) == 0L) {
    return("No gates are defined.")
  }

  paste(
    paste0(
      strrep("  ", tree$depth),
      ifelse(tree$depth == 0L, "", "└─ "),
      tree$gate_name,
      " [",
      tree$type,
      "]"
    ),
    collapse = "\n"
  )
}


# Convert raw polygon-gate coordinates into a Plotly SVG shape on the
# currently selected display scale.
cytoTerra_workspace_polygon_to_plotly_shape = function(
    gate,
    x_channel,
    y_channel,
    transform,
    cofactor,
    line_color = "#D62728"
) {

  if (is.null(gate) ||
      !identical(gate$type, "polygon") ||
      is.null(gate$coordinates) ||
      !all(c("x", "y") %in% names(gate$coordinates)) ||
      !identical(gate$channels[["x"]], x_channel) ||
      !identical(gate$channels[["y"]], y_channel)) {
    return(NULL)
  }

  coordinates <- gate$coordinates

  x_display <- cytoTerra_workspace_transform_values(
    coordinates$x,
    transform = transform,
    cofactor = cofactor
  )

  y_display <- cytoTerra_workspace_transform_values(
    coordinates$y,
    transform = transform,
    cofactor = cofactor
  )

  if (any(!is.finite(x_display)) ||
      any(!is.finite(y_display))) {
    return(NULL)
  }

  path <- paste0(
    "M ",
    formatC(x_display[[1L]], digits = 12, format = "fg"),
    ",",
    formatC(y_display[[1L]], digits = 12, format = "fg"),
    paste0(
      " L ",
      formatC(
        x_display[-1L],
        digits = 12,
        format = "fg"
      ),
      ",",
      formatC(
        y_display[-1L],
        digits = 12,
        format = "fg"
      ),
      collapse = ""
    ),
    " Z"
  )

  list(
    type = "path",
    path = path,
    xref = "x",
    yref = "y",
    fillcolor = "rgba(0,0,0,0)",
    line = list(
      color = line_color,
      width = 2
    )
  )
}


cytoTerra_workspace_ui = function(
    sample_id,
    channel_choices,
    transform,
    cofactor
) {

  shiny::fluidPage(

    shiny::tags$head(
      shiny::tags$style(
        shiny::HTML("
          .cytoTerra-tree {
            font-family: monospace;
            font-size: 0.9em;
            white-space: pre-wrap;
            max-height: 340px;
            overflow-y: auto;
            border: 1px solid #d9d9d9;
            border-radius: 4px;
            background: #fafafa;
            padding: 10px;
          }

          .cytoTerra-status {
            margin-top: 10px;
            padding: 10px;
            border-radius: 4px;
            background: #f5f5f5;
          }
        ")
      )
    ),

    shiny::titlePanel(
      paste0(
        "cytoTerra manual gating workspace — ",
        sample_id
      )
    ),

    shiny::fluidRow(

      shiny::column(
        width = 3,

        shiny::h4("Gate hierarchy"),

        shiny::div(
          class = "cytoTerra-tree",
          shiny::verbatimTextOutput("gate_tree")
        ),

        shiny::h4("Current parent population"),

        shiny::uiOutput("parent_gate_ui"),

        shiny::actionButton(
          "set_parent_gate",
          "Set current parent",
          width = "100%"
        ),

        shiny::br(),
        shiny::br(),

        shiny::actionButton(
          "go_root",
          "Go to root",
          width = "100%"
        ),

        shiny::br(),
        shiny::br(),

        shiny::actionButton(
          "go_parent",
          "Up to parent",
          width = "100%"
        ),

        shiny::br(),
        shiny::br(),

        shiny::uiOutput("children_ui"),

        shiny::actionButton(
          "go_child",
          "Go to selected child",
          width = "100%"
        ),

        shiny::tags$hr(),

        shiny::h4("Delete current gate"),

        shiny::checkboxInput(
          "delete_descendants",
          "Also delete all descendants",
          value = FALSE
        ),

        shiny::actionButton(
          "delete_current_gate",
          "Delete current gate",
          class = "btn-danger",
          width = "100%"
        )
      ),

      shiny::column(
        width = 9,

        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::textInput(
              "gate_name",
              "New child-gate name"
            )
          ),

          shiny::column(
            width = 4,
            shiny::radioButtons(
              "gate_mode",
              "Gate type",
              choices = c(
                "2D polygon" = "polygon",
                "1D range" = "range"
              ),
              selected = "polygon",
              inline = TRUE
            )
          ),

          shiny::column(
            width = 4,
            shiny::selectInput(
              "plot_transform",
              "Display transformation",
              choices = c(
                "asinh" = "asinh",
                "log1p" = "log1p",
                "Raw FCS scale" = "none"
              ),
              selected = transform
            )
          )
        ),

        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::numericInput(
              "plot_cofactor",
              "asinh cofactor",
              value = cofactor,
              min = 0.000001,
              step = 10
            )
          ),

          shiny::column(
            width = 4,
            shiny::actionButton(
              "redraw",
              "Redraw plot",
              width = "100%"
            )
          ),

          shiny::column(
            width = 4,
            shiny::textOutput("current_gate_text")
          )
        ),

        shiny::uiOutput("parameter_controls"),

        shiny::div(
          style = "overflow-x: auto;",
          plotly::plotlyOutput(
            "gate_plot",
            width = "650px",
            height = "650px"
          )
        ),

        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::actionButton(
              "save_gate",
              "Save gate to project",
              class = "btn-primary",
              width = "100%"
            )
          ),

          shiny::column(
            width = 6,
            shiny::actionButton(
              "clear_drawing",
              "Clear current drawing",
              width = "100%"
            )
          )
        ),

        shiny::div(
          class = "cytoTerra-status",
          shiny::textOutput("status")
        )
      )
    )
  )
}


cytoTerra_workspace_server = function(
    project,
    sample_id,
    expression_matrix,
    available_channels,
    channel_choices,
    start_gate,
    n_plot,
    point_size
) {

  force(project)
  force(sample_id)
  force(expression_matrix)
  force(available_channels)
  force(channel_choices)
  force(start_gate)
  force(n_plot)
  force(point_size)

  function(input, output, session) {

    current_gate <- shiny::reactiveVal(start_gate)
    latest_polygon_path <- shiny::reactiveVal(NULL)
    project_revision <- shiny::reactiveVal(0L)

    status_message <- shiny::reactiveVal(
      "Choose a parent population, draw a gate, name it, and save it to the project."
    )

    current_parent_mask <- shiny::reactive({

      project_revision()

      gate_name <- current_gate()

      if (identical(gate_name, "root")) {
        return(rep(TRUE, nrow(expression_matrix)))
      }

      project$evaluate_gate(
        sample_id = sample_id,
        gate_name = gate_name
      )
    })

    output$current_gate_text <- shiny::renderText({
      project_revision()

      gate_name <- current_gate()

      if (identical(gate_name, "root")) {
        return("Parent: root (all events)")
      }

      gate <- project$get_gates(sample_id, gate_name)

      paste0(
        "Parent: ",
        gate_name,
        " [",
        gate$type,
        "]"
      )
    })

    output$gate_tree <- shiny::renderText({
      project_revision()

      tree <- project$get_gate_tree(
        sample_id = sample_id,
        include_root = TRUE
      )

      cytoTerra_workspace_gate_tree_text(tree)
    })

    output$parent_gate_ui <- shiny::renderUI({

      project_revision()

      tree <- project$get_gate_tree(
        sample_id = sample_id,
        include_root = TRUE
      )

      gate_names <- tree$gate_name

      gate_labels <- ifelse(
        gate_names == "root",
        "root (all events)",
        paste0(
          strrep("  ", tree$depth),
          gate_names,
          " [",
          tree$type,
          "]"
        )
      )

      shiny::selectInput(
        inputId = "parent_gate",
        label = "Parent population",
        choices = stats::setNames(gate_names, gate_labels),
        selected = current_gate()
      )
    })

    output$children_ui <- shiny::renderUI({
      project_revision()

      children <- project$get_child_gates(
        sample_id = sample_id,
        parent_gate = current_gate(),
        recursive = FALSE
      )

      shiny::selectInput(
        "child_gate",
        "Child gate",
        choices = children
      )
    })

    shiny::observeEvent(input$go_root, {
      current_gate("root")
      latest_polygon_path(NULL)
      status_message("Current parent changed to root.")
    })

    shiny::observeEvent(input$go_parent, {

      if (identical(current_gate(), "root")) {
        status_message("Already at root.")
        return()
      }

      parent_gate <- project$get_parent_gate(
        sample_id = sample_id,
        gate_name = current_gate()
      )

      current_gate(parent_gate)
      latest_polygon_path(NULL)

      status_message(
        paste0("Current parent changed to '", parent_gate, "'.")
      )
    })

    shiny::observeEvent(input$go_child, {

      if (is.null(input$child_gate) ||
          !nzchar(input$child_gate)) {
        status_message("Select a child gate first.")
        return()
      }

      current_gate(input$child_gate)
      latest_polygon_path(NULL)

      status_message(
        paste0(
          "Current parent changed to '",
          input$child_gate,
          "'."
        )
      )
    })

    shiny::observeEvent(input$set_parent_gate, {

      if (is.null(input$parent_gate) ||
          !nzchar(input$parent_gate)) {
        status_message("Select a parent population first.")
        return()
      }

      current_gate(input$parent_gate)
      latest_polygon_path(NULL)

      status_message(
        paste0(
          "Current parent changed to '",
          input$parent_gate,
          "'. New gates will be created under this population."
        )
      )
    })

    output$parameter_controls <- shiny::renderUI({

      if (identical(input$gate_mode, "polygon")) {
        return(
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::selectInput(
                "x_channel",
                "X parameter",
                choices = channel_choices,
                selected = available_channels[[1L]]
              )
            ),

            shiny::column(
              width = 6,
              shiny::selectInput(
                "y_channel",
                "Y parameter",
                choices = channel_choices,
                selected = if (length(available_channels) > 1L) {
                  available_channels[[2L]]
                } else {
                  available_channels[[1L]]
                }
              )
            )
          )
        )
      }

      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::selectInput(
            "range_channel",
            "Histogram parameter",
            choices = channel_choices,
            selected = available_channels[[1L]]
          )
        ),

        shiny::column(
          width = 4,
          shiny::numericInput(
            "range_lower",
            "Lower bound (display scale)",
            value = 0
          )
        ),

        shiny::column(
          width = 4,
          shiny::numericInput(
            "range_upper",
            "Upper bound (display scale)",
            value = 1
          )
        )
      )
    })

    plotted_data <- shiny::reactive({

      input$redraw
      project_revision()

      parent_index <- which(current_parent_mask())

      if (length(parent_index) == 0L) {
        return(NULL)
      }

      plot_index <- if (length(parent_index) <= n_plot) {
        parent_index
      } else {
        sample(
          parent_index,
          size = as.integer(n_plot),
          replace = FALSE
        )
      }

      list(
        index = plot_index,
        transform = input$plot_transform,
        cofactor = input$plot_cofactor
      )
    })

    output$gate_plot <- plotly::renderPlotly({

      plot_data <- plotted_data()
      shiny::req(plot_data)

      current_transform <- plot_data$transform
      current_cofactor <- plot_data$cofactor

      if (!is.numeric(current_cofactor) ||
          length(current_cofactor) != 1L ||
          is.na(current_cofactor) ||
          !is.finite(current_cofactor) ||
          current_cofactor <= 0) {
        current_cofactor <- 150
      }

      if (identical(input$gate_mode, "range")) {

        shiny::req(input$range_channel)

        values_raw <- expression_matrix[
          plot_data$index,
          input$range_channel
        ]

        valid <- is.finite(values_raw)

        if (identical(current_transform, "log1p")) {
          valid <- valid & values_raw >= 0
        }

        values_display <- cytoTerra_workspace_transform_values(
          values_raw[valid],
          transform = current_transform,
          cofactor = current_cofactor
        )

        return(
          plotly::plot_ly(
            x = values_display,
            type = "histogram",
            marker = list(color = "#1F77B4")
          ) |>
            plotly::layout(
              width = 650, 
              height = 650, 
              autosize = FALSE,
                
              title = paste0(
                sample_id,
                " — parent: ",
                current_gate()
              ),
              xaxis = list(
                title = paste0(
                  input$range_channel,
                  " (",
                  cytoTerra_workspace_transform_label(
                    current_transform,
                    current_cofactor
                  ),
                  ")"
                )
              ),
              yaxis = list(title = "Event count")
            ) |>
            plotly::config(displaylogo = FALSE)
        )
      }

      shiny::req(input$x_channel, input$y_channel)

      if (identical(input$x_channel, input$y_channel)) {
        return(
          plotly::plot_ly() |>
            plotly::layout(
              title = "Choose two different parameters."
            )
        )
      }

      x_raw <- expression_matrix[
        plot_data$index,
        input$x_channel
      ]

      y_raw <- expression_matrix[
        plot_data$index,
        input$y_channel
      ]

      valid <- is.finite(x_raw) & is.finite(y_raw)

      if (identical(current_transform, "log1p")) {
        valid <- valid & x_raw >= 0 & y_raw >= 0
      }

      x_display <- cytoTerra_workspace_transform_values(
        x_raw[valid],
        transform = current_transform,
        cofactor = current_cofactor
      )

      y_display <- cytoTerra_workspace_transform_values(
        y_raw[valid],
        transform = current_transform,
        cofactor = current_cofactor
      )

      # Display existing direct-child polygon gates, if they use the
      # currently selected X and Y channels.
      child_gates <- project$get_child_gates(
        sample_id = sample_id,
        parent_gate = current_gate(),
        recursive = FALSE
      )

      shapes <- lapply(
        child_gates,
        function(gate_name) {
          gate <- project$get_gates(sample_id, gate_name)

          cytoTerra_workspace_polygon_to_plotly_shape(
            gate = gate,
            x_channel = input$x_channel,
            y_channel = input$y_channel,
            transform = current_transform,
            cofactor = current_cofactor
          )
        }
      )

      shapes <- Filter(Negate(is.null), shapes)

      p <- plotly::plot_ly(
        data = data.frame(
          x = x_display,
          y = y_display
        ),
        x = ~x,
        y = ~y,
        type = "scattergl",
        mode = "markers",
        source = "gate_plot_source",
        marker = list(
          color = "#1F77B4",
          size = point_size,
          opacity = 0.35
        ),
        hoverinfo = "skip"
      ) |>
        plotly::layout(
          width = 650, 
          height = 650, 
          autosize = FALSE,
            
          title = paste0(
            sample_id,
            " — parent: ",
            current_gate(),
            "<br><sup>Draw a closed polygon, name the new gate, then save it.</sup>"
          ),
          xaxis = list(
            title = paste0(
              input$x_channel,
              " (",
              cytoTerra_workspace_transform_label(
                current_transform,
                current_cofactor
              ),
              ")"
            ),
            zeroline = FALSE,
            scaleanchor = 'y', ## suggest off
            scaleratio = 1
          ),
          yaxis = list(
            title = paste0(
              input$y_channel,
              " (",
              cytoTerra_workspace_transform_label(
                current_transform,
                current_cofactor
              ),
              ")"
            ),
            zeroline = FALSE
          ),
          shapes = shapes,
          dragmode = "drawclosedpath"
        ) |>
        plotly::config(
          displaylogo = FALSE,
          editable = FALSE,
          modeBarButtonsToAdd = c(
            "drawclosedpath",
            "eraseshape"
          )
        )

      # Required before event_data("plotly_relayout", ...) is used.
      cytoTerra_workspace_capture_relayout(p)
    })

    shiny::observeEvent(
      input$gate_plot_relayout,
      {

        relayout <- input$gate_plot_relayout

        path <- cytoTerra_workspace_get_path_from_relayout(
          relayout
        )

        if (is.null(path)) {
          return()
        }

        latest_polygon_path(path)

        status_message(
          "Polygon captured. Click Save gate to project."
        )
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$clear_drawing, {
      latest_polygon_path(NULL)

      status_message(
        "Current unsaved polygon was cleared. Draw a new polygon if needed."
      )
    })

    shiny::observeEvent(input$save_gate, {

      gate_name <- trimws(input$gate_name)

      if (!nzchar(gate_name)) {
        status_message("Enter a new child-gate name first.")
        return()
      }

      if (identical(gate_name, "root")) {
        status_message("The gate name 'root' is reserved.")
        return()
      }

      if (!is.null(project$get_gates(sample_id, gate_name))) {
        status_message(
          paste0(
            "A gate named '",
            gate_name,
            "' already exists. Choose a different name."
          )
        )
        return()
      }

      current_transform <- input$plot_transform
      current_cofactor <- input$plot_cofactor

      if (!is.numeric(current_cofactor) ||
          length(current_cofactor) != 1L ||
          is.na(current_cofactor) ||
          !is.finite(current_cofactor) ||
          current_cofactor <= 0) {
        status_message("Enter a positive asinh cofactor.")
        return()
      }

      if (identical(input$gate_mode, "polygon")) {

        if (is.null(input$x_channel) ||
            is.null(input$y_channel) ||
            identical(input$x_channel, input$y_channel)) {
          status_message("Choose two different polygon parameters.")
          return()
        }

        coordinates <- cytoTerra_workspace_polygon_from_path(
          latest_polygon_path()
        )

        if (is.null(coordinates)) {
          status_message(
            paste(
              "No readable polygon was found.",
              "Select Draw closed shape in the Plotly toolbar and draw a closed polygon."
            )
          )
          return()
        }

        project$add_polygon_gate(
          sample_id = sample_id,
          gate_name = gate_name,
          parent = current_gate(),
          x_channel = input$x_channel,
          y_channel = input$y_channel,
          coordinates = coordinates,
          coordinate_scale = "transformed",
          transform = current_transform,
          cofactor = current_cofactor,
          x_channel_description = input$x_channel,
          y_channel_description = input$y_channel,
          note = "Created in the cytoTerra Shiny gating workspace."
        )

        latest_polygon_path(NULL)

        status_message(
          paste0(
            "Saved polygon gate '",
            gate_name,
            "' under parent '",
            current_gate(),
            "'."
          )
        )
      } else {

        lower <- input$range_lower
        upper <- input$range_upper

        if (!is.numeric(lower) ||
            !is.numeric(upper) ||
            length(lower) != 1L ||
            length(upper) != 1L ||
            is.na(lower) ||
            is.na(upper) ||
            !is.finite(lower) ||
            !is.finite(upper) ||
            lower >= upper) {
          status_message(
            "Enter finite range boundaries where lower is less than upper."
          )
          return()
        }

        project$add_range_gate(
          sample_id = sample_id,
          gate_name = gate_name,
          parent = current_gate(),
          channel = input$range_channel,
          lower = lower,
          upper = upper,
          coordinate_scale = "transformed",
          transform = current_transform,
          cofactor = current_cofactor,
          channel_input = input$range_channel,
          channel_description = input$range_channel,
          note = "Created in the cytoTerra Shiny gating workspace."
        )

        status_message(
          paste0(
            "Saved range gate '",
            gate_name,
            "' under parent '",
            current_gate(),
            "'."
          )
        )
      }

      shiny::updateTextInput(
        session,
        "gate_name",
        value = ""
      )

      project_revision(project_revision() + 1L)
    })

    shiny::observeEvent(input$delete_current_gate, {

      gate_name <- current_gate()

      if (identical(gate_name, "root")) {
        status_message("The root population cannot be deleted.")
        return()
      }

      child_gates <- project$get_child_gates(
        sample_id = sample_id,
        parent_gate = gate_name,
        recursive = FALSE
      )

      if (length(child_gates) > 0L &&
          !isTRUE(input$delete_descendants)) {
        status_message(
          paste0(
            "'",
            gate_name,
            "' has child gates. Check “Also delete all descendants” to remove the branch."
          )
        )
        return()
      }

      shiny::showModal(
        shiny::modalDialog(
          title = "Delete gate?",
          paste0(
            "Delete gate '",
            gate_name,
            "'",
            if (isTRUE(input$delete_descendants)) {
              " and all descendant gates?"
            } else {
              "?"
            }
          ),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(
              "confirm_delete_current_gate",
              "Delete",
              class = "btn-danger"
            )
          )
        )
      )
    })

    shiny::observeEvent(input$confirm_delete_current_gate, {

      gate_name <- current_gate()

      shiny::removeModal()

      project$delete_gate(
        sample_id = sample_id,
        gate_name = gate_name,
        include_descendants = isTRUE(input$delete_descendants)
      )

      current_gate("root")
      latest_polygon_path(NULL)

      project_revision(project_revision() + 1L)

      status_message(
        paste0(
          "Deleted gate branch beginning at '",
          gate_name,
          "'. Current parent is now root."
        )
      )
    })

    output$status <- shiny::renderText({
      status_message()
    })
  }
}

cytoTerra_workspace_capture_relayout = function(plot) {

  htmlwidgets::onRender(
    plot,
    "
    function(el, x) {
      el.on('plotly_relayout', function(eventData) {
        Shiny.setInputValue(
          'gate_plot_relayout',
          eventData,
          { priority: 'event' }
        );
      });
    }
    "
  )
}
