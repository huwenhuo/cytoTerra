# Gate manually check

cytoTerra_methods_plot <- list(

  plot_2d_cells = function(
    sample_ids = NULL,
    parent_gate,
    x_channel,
    y_channel,
    layout = c(1, 1),
    max_cells = 20000,
    point_size = 0.15,
    alpha = 0.35,
    transform = c("asinh", "none"),
    cofactor = 150
  ) {

    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(
      project = self,
      sample_ids = sample_ids
    )

    if (missing(parent_gate) ||
        !is.character(parent_gate) ||
        length(parent_gate) != 1L ||
        is.na(parent_gate) ||
        !nzchar(parent_gate)) {
      stop("Please provide one non-empty parent_gate name.", call. = FALSE)
    }

    if (missing(x_channel) ||
        !is.character(x_channel) ||
        length(x_channel) != 1L ||
        is.na(x_channel) ||
        !nzchar(x_channel)) {
      stop("Please provide one non-empty x_channel.", call. = FALSE)
    }

    if (missing(y_channel) ||
        !is.character(y_channel) ||
        length(y_channel) != 1L ||
        is.na(y_channel) ||
        !nzchar(y_channel)) {
      stop("Please provide one non-empty y_channel.", call. = FALSE)
    }

    if (!is.numeric(layout) ||
        length(layout) != 2L ||
        any(!is.finite(layout)) ||
        any(layout < 1)) {
      stop(
        "layout must be a numeric vector of length 2: c(nrow, ncol).",
        call. = FALSE
      )
    }

    if (!is.numeric(max_cells) ||
        length(max_cells) != 1L ||
        is.na(max_cells) ||
        max_cells < 1L) {
      stop("max_cells must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(point_size) ||
        length(point_size) != 1L ||
        is.na(point_size) ||
        point_size <= 0) {
      stop("point_size must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(alpha) ||
        length(alpha) != 1L ||
        is.na(alpha) ||
        alpha < 0 ||
        alpha > 1) {
      stop("alpha must be one number between 0 and 1.", call. = FALSE)
    }

    transform <- match.arg(transform)

    if (!is.numeric(cofactor) ||
        length(cofactor) != 1L ||
        is.na(cofactor) ||
        cofactor <= 0) {
      stop("cofactor must be one positive numeric value.", call. = FALSE)
    }

    axis_label <- function(channel) {
      if (transform == "asinh") {
        paste0(channel, " [asinh(x / ", cofactor, ")]")
      } else {
        channel
      }
    }

    plot_data <- lapply(sample_ids, function(sample_id) {

      # Resolve marker description (e.g., "CD3") independently for
      # each sample, since FCS parameter names can vary among samples.
      flow_frame <- self$get_flow_frame(sample_id)

      x_channel_info <- cytoTerra_resolve_channel(
        flow_frame = flow_frame,
        channel = x_channel
      )

      y_channel_info <- cytoTerra_resolve_channel(
        flow_frame = flow_frame,
        channel = y_channel
      )

      x_fcs_channel <- as.character(x_channel_info$name)
      y_fcs_channel <- as.character(y_channel_info$name)

      if (length(x_fcs_channel) != 1L ||
          is.na(x_fcs_channel) ||
          !nzchar(x_fcs_channel)) {
        stop(
          "Could not resolve x_channel '", x_channel,
          "' to exactly one FCS parameter for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      if (length(y_fcs_channel) != 1L ||
          is.na(y_fcs_channel) ||
          !nzchar(y_fcs_channel)) {
        stop(
          "Could not resolve y_channel '", y_channel,
          "' to exactly one FCS parameter for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      # Retrieve events within parent_gate and all of its ancestor gates.
      cells <- cytoTerra_get_parent_gate_cells(
        project = self,
        sample_id = sample_id,
        parent_gate = parent_gate
      )

      if (!x_fcs_channel %in% colnames(cells)) {
        stop(
          "Resolved x channel '", x_fcs_channel,
          "' was not found in parent-gated cells for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      if (!y_fcs_channel %in% colnames(cells)) {
        stop(
          "Resolved y channel '", y_fcs_channel,
          "' was not found in parent-gated cells for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      x_values <- cells[[x_fcs_channel]]
      y_values <- cells[[y_fcs_channel]]

      # Display-only transform. Raw event data and gate membership
      # remain unchanged.
      if (transform == "asinh") {
        x_values <- asinh(x_values / cofactor)
        y_values <- asinh(y_values / cofactor)
      }

      sample_plot_data <- data.frame(
        sample_id = sample_id,
        x_value = x_values,
        y_value = y_values,
        stringsAsFactors = FALSE
      )

      sample_plot_data <- sample_plot_data[
        is.finite(sample_plot_data$x_value) &
          is.finite(sample_plot_data$y_value),
        ,
        drop = FALSE
      ]

      # Plot-only random downsampling.
      if (nrow(sample_plot_data) > max_cells) {
        sample_plot_data <- sample_plot_data[
          sample(seq_len(nrow(sample_plot_data)), max_cells),
          ,
          drop = FALSE
        ]
      }

      sample_plot_data
    })

    plot_data <- do.call(rbind, plot_data)

    if (is.null(plot_data) || nrow(plot_data) == 0L) {
      stop(
        "No finite cells were available to plot after applying parent_gate '",
        parent_gate, "'.",
        call. = FALSE
      )
    }

    plot_data$sample_id <- factor(
      plot_data$sample_id,
      levels = sample_ids
    )

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = x_value, y = y_value)
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha,
        color = "#1F78B4"
      ) +
      ggplot2::facet_wrap(
        ~sample_id,
        nrow = layout[1],
        ncol = layout[2]
      ) +
      ggplot2::labs(
        title = paste0(
          parent_gate, " cells: ",
          x_channel, " versus ", y_channel
        ),
        subtitle = "Parent-gated cells; display transformation only",
        x = axis_label(x_channel),
        y = axis_label(y_channel)
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        aspect.ratio = 1,
        strip.background = ggplot2::element_rect(
          fill = "grey95",
          color = "grey70"
        )
      )
  },


  plot_1d_cells = function(
    sample_ids = NULL,
    parent_gate,
    channel,
    layout = c(1, 1),
    bins = 100,
    max_cells = 50000,
    fill = "#2C7FB8",
    alpha = 0.75,
    transform = c("asinh", "none"),
    cofactor = 150
  ) {

    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(
      project = self,
      sample_ids = sample_ids
    )

    if (missing(parent_gate) ||
        !is.character(parent_gate) ||
        length(parent_gate) != 1L ||
        is.na(parent_gate) ||
        !nzchar(parent_gate)) {
      stop("Please provide one non-empty parent_gate name.", call. = FALSE)
    }

    if (missing(channel) ||
        !is.character(channel) ||
        length(channel) != 1L ||
        is.na(channel) ||
        !nzchar(channel)) {
      stop("Please provide one non-empty channel.", call. = FALSE)
    }

    if (!is.numeric(layout) ||
        length(layout) != 2L ||
        any(!is.finite(layout)) ||
        any(layout < 1)) {
      stop(
        "layout must be a numeric vector of length 2: c(nrow, ncol).",
        call. = FALSE
      )
    }

    if (!is.numeric(bins) ||
        length(bins) != 1L ||
        is.na(bins) ||
        bins < 1L) {
      stop("bins must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(max_cells) ||
        length(max_cells) != 1L ||
        is.na(max_cells) ||
        max_cells < 1L) {
      stop("max_cells must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(alpha) ||
        length(alpha) != 1L ||
        is.na(alpha) ||
        alpha < 0 ||
        alpha > 1) {
      stop("alpha must be one number between 0 and 1.", call. = FALSE)
    }

    transform <- match.arg(transform)

    if (!is.numeric(cofactor) ||
        length(cofactor) != 1L ||
        is.na(cofactor) ||
        cofactor <= 0) {
      stop("cofactor must be one positive numeric value.", call. = FALSE)
    }

    x_axis_label <- if (transform == "asinh") {
      paste0(channel, " [asinh(x / ", cofactor, ")]")
    } else {
      channel
    }

    plot_data <- lapply(sample_ids, function(sample_id) {

      flow_frame <- self$get_flow_frame(sample_id)

      channel_info <- cytoTerra_resolve_channel(
        flow_frame = flow_frame,
        channel = channel
      )

      fcs_channel <- as.character(channel_info$name)

      if (length(fcs_channel) != 1L ||
          is.na(fcs_channel) ||
          !nzchar(fcs_channel)) {
        stop(
          "Could not resolve channel '", channel,
          "' to exactly one FCS parameter for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      cells <- cytoTerra_get_parent_gate_cells(
        project = self,
        sample_id = sample_id,
        parent_gate = parent_gate
      )

      if (!fcs_channel %in% colnames(cells)) {
        stop(
          "Resolved channel '", fcs_channel,
          "' was not found in parent-gated cells for sample '",
          sample_id, "'.",
          call. = FALSE
        )
      }

      values <- cells[[fcs_channel]]

      # Display-only transform.
      if (transform == "asinh") {
        values <- asinh(values / cofactor)
      }

      sample_plot_data <- data.frame(
        sample_id = sample_id,
        value = values,
        stringsAsFactors = FALSE
      )

      sample_plot_data <- sample_plot_data[
        is.finite(sample_plot_data$value),
        ,
        drop = FALSE
      ]

      # Correct plotting-only downsampling.
      if (nrow(sample_plot_data) > max_cells) {
        sample_plot_data <- sample_plot_data[
          sample(seq_len(nrow(sample_plot_data)), max_cells),
          ,
          drop = FALSE
        ]
      }

      sample_plot_data
    })

    plot_data <- do.call(rbind, plot_data)

    if (is.null(plot_data) || nrow(plot_data) == 0L) {
      stop(
        "No finite cells were available to plot after applying parent_gate '",
        parent_gate, "'.",
        call. = FALSE
      )
    }

    plot_data$sample_id <- factor(
      plot_data$sample_id,
      levels = sample_ids
    )

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = value)
    ) +
      ggplot2::geom_histogram(
        bins = bins,
        fill = fill,
        color = "white",
        alpha = alpha
      ) +
      ggplot2::facet_wrap(
        ~sample_id,
        nrow = layout[1],
        ncol = layout[2],
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(parent_gate, ": ", channel),
        subtitle = "Parent-gated cells; display transformation only",
        x = x_axis_label,
        y = "Cell count"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(
          fill = "grey95",
          color = "grey70"
        ),
        panel.spacing = grid::unit(0.8, "lines")
      )
  }
)
