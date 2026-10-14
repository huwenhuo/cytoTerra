# Interactive Shiny workspace for manual gating.
#'
#' @noRd
cytoTerra_methods_gating_workspace <- list(

  launch_gating_workspace = function(
      sample_id,
      start_gate = "root",
      n_plot = 30000L,
      transform = c("asinh", "log1p", "none"),
      cofactor = 150,
      point_size = 3,
      shiny_port = 3838L,
      ood_origin = "https://ood-bouchet.ycrc.yale.edu"
  ) {

    transform <- match.arg(transform)

    cytoTerra_validate_sample_ids(self, sample_id)

    if (!is.character(start_gate) ||
        length(start_gate) != 1L ||
        is.na(start_gate) ||
        !nzchar(start_gate)) {
      stop(
        "`start_gate` must be one non-empty gate name or 'root'.",
        call. = FALSE
      )
    }

    if (!identical(start_gate, "root") &&
        is.null(self$get_gates(sample_id, start_gate))) {
      stop(
        "No gate named '", start_gate,
        "' was found for sample '", sample_id, "'.",
        call. = FALSE
      )
    }

    if (!is.numeric(n_plot) ||
        length(n_plot) != 1L ||
        is.na(n_plot) ||
        !is.finite(n_plot) ||
        n_plot < 100L) {
      stop(
        "`n_plot` must be one finite number of at least 100.",
        call. = FALSE
      )
    }

    if (!is.numeric(cofactor) ||
        length(cofactor) != 1L ||
        is.na(cofactor) ||
        !is.finite(cofactor) ||
        cofactor <= 0) {
      stop(
        "`cofactor` must be one positive finite number.",
        call. = FALSE
      )
    }

    if (!is.numeric(point_size) ||
        length(point_size) != 1L ||
        is.na(point_size) ||
        !is.finite(point_size) ||
        point_size <= 0) {
      stop(
        "`point_size` must be one positive finite number.",
        call. = FALSE
      )
    }

    if (!is.numeric(shiny_port) ||
        length(shiny_port) != 1L ||
        is.na(shiny_port) ||
        !is.finite(shiny_port) ||
        shiny_port != as.integer(shiny_port) ||
        shiny_port < 1024L ||
        shiny_port > 65535L) {
      stop(
        "`shiny_port` must be an integer between 1024 and 65535.",
        call. = FALSE
      )
    }

    if (!is.character(ood_origin) ||
        length(ood_origin) != 1L ||
        is.na(ood_origin) ||
        !nzchar(ood_origin) ||
        !grepl("^https?://", ood_origin)) {
      stop(
        "`ood_origin` must be one non-empty http(s) URL.",
        call. = FALSE
      )
    }

    if (!requireNamespace("shiny", quietly = TRUE)) {
      stop(
        "`launch_gating_workspace()` requires the shiny package.",
        call. = FALSE
      )
    }

    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop(
        "`launch_gating_workspace()` requires the plotly package.",
        call. = FALSE
      )
    }

    shiny_port <- as.integer(shiny_port)
    ood_origin <- sub("/+$", "", ood_origin)

    flow_frame <- self$get_flow_frame(sample_id)
    expression_matrix <- flowCore::exprs(flow_frame)
    available_channels <- colnames(expression_matrix)

    if (length(available_channels) == 0L) {
      stop(
        "Sample '", sample_id,
        "' has no expression-matrix channels.",
        call. = FALSE
      )
    }

    channel_choices <- cytoTerra_workspace_channel_choices(
      flow_frame = flow_frame,
      available_channels = available_channels
    )

    compute_node <- Sys.getenv("host")
    jupyter_port <- Sys.getenv("port")

    if (!nzchar(compute_node) || !nzchar(jupyter_port)) {
      stop(
        paste(
          "Could not obtain the Open OnDemand compute node and Jupyter port",
          "from environment variables `host` and `port`."
        ),
        call. = FALSE
      )
    }

    shiny_url <- sprintf(
      "%s/node/%s/%s/proxy/%s/",
      ood_origin,
      compute_node,
      jupyter_port,
      shiny_port
    )

    message(
      "\n",
      "============================================================\n",
      "cytoTerra Shiny gating workspace\n",
      "Open this URL in a new browser tab:\n",
      shiny_url, "\n",
      "============================================================\n"
    )

    if (requireNamespace("IRdisplay", quietly = TRUE)) {
      IRdisplay::display_html(
        sprintf(
          paste0(
            '<p><a href="%s" target="_blank" rel="noopener noreferrer">',
            "Open cytoTerra manual-gating workspace",
            "</a></p>"
          ),
          shiny_url
        )
      )
    }

    ui <- cytoTerra_workspace_ui(
      sample_id = sample_id,
      channel_choices = channel_choices,
      transform = transform,
      cofactor = cofactor
    )

    server <- cytoTerra_workspace_server(
      project = self,
      sample_id = sample_id,
      expression_matrix = expression_matrix,
      available_channels = available_channels,
      channel_choices = channel_choices,
      start_gate = start_gate,
      n_plot = n_plot,
      point_size = point_size
    )

    shiny::runApp(
      shiny::shinyApp(ui = ui, server = server),
      host = "127.0.0.1",
      port = shiny_port,
      launch.browser = FALSE
    )

    invisible(self)
  }
)