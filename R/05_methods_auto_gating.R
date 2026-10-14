# Automatic gating methods.

cytoTerra_methods_auto_gating <- list(

  # Automatically create an FSC-A/FSC-H singlet gate.
  auto_gate_singlets = function(
      sample_ids = NULL,
      gate_name = "Singlets",
      parent = "root",
      channels = c(x = "FSC-A", y = "FSC-H"),
      n_max = 4000L,
      mad_multiplier = 4,
      save_as_default = FALSE,
      verbose = TRUE
  ) {

    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(self, sample_ids)

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
	  stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
	}
	
    if (isTRUE(verbose)) {
      message(
        "Automatic singlet gating: ",
        length(sample_ids),
        " sample",
        if (length(sample_ids) == 1L) "" else "s",
        "; gate name = '",
        gate_name,
        "'."
      )
    }

    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop(
        "`gate_name` must be one non-empty character value.",
        call. = FALSE
      )
    }

    if (!is.character(parent) ||
        length(parent) != 1L ||
        is.na(parent) ||
        !nzchar(parent)) {
      stop(
        "`parent` must be one non-empty character value.",
        call. = FALSE
      )
    }

    if (!is.character(channels) ||
        length(channels) != 2L ||
        !identical(names(channels), c("x", "y"))) {
      stop(
        "`channels` must be named exactly c(x = 'FSC-A', y = 'FSC-H').",
        call. = FALSE
      )
    }

    if (!is.numeric(n_max) || length(n_max) != 1L || n_max < 100L) {
      stop(
        "`n_max` must be one numeric value of at least 100.",
        call. = FALSE
      )
    }

    if (!is.numeric(mad_multiplier) ||
        length(mad_multiplier) != 1L ||
        mad_multiplier <= 0) {
      stop(
        "`mad_multiplier` must be one positive numeric value.",
        call. = FALSE
      )
    }

    # Parent gating cannot yet be applied until apply_gating() exists.
    if (!identical(parent, "root")) {
      stop(
        "Currently, `auto_gate_singlets()` supports only `parent = 'root'`. ",
        "Parent-event membership will be supported after gate application ",
        "is implemented.",
        call. = FALSE
      )
    }

     for (sample_index in seq_along(sample_ids)) {
       sample_id <- sample_ids[[sample_index]]
     
       if (isTRUE(verbose)) {
         message(
           "[",
           sample_index,
           "/",
           length(sample_ids),
           "] Fitting singlet gate for sample '",
           sample_id,
           "'..."
         )
       }
     
       started_at <- Sys.time()
     
       flow_frame <- self$get_flow_frame(sample_id)

      available_channels <- colnames(flowCore::exprs(flow_frame))
      missing_channels <- setdiff(unname(channels), available_channels)

      if (length(missing_channels) > 0L) {
        stop(
          "Sample '", sample_id, "' does not contain: ",
          paste(missing_channels, collapse = ", "),
          ". Available channel names should be inspected with ",
          "`colnames(flowCore::exprs(project$get_flow_frame('",
          sample_id, "')))`. ",
          call. = FALSE
        )
      }

      event_index <- cytoTerra_draw_event_sample(
        flow_frame = flow_frame,
        n_max = n_max
      )

      # Extract only two channels and only the sampled events.
      event_data <- flowCore::exprs(flow_frame)[
        event_index,
        unname(channels),
        drop = FALSE
      ]

      singlet_fit <- cytoTerra_fit_singlet_polygon(
        x = event_data[, channels[["x"]]],
        y = event_data[, channels[["y"]]],
        mad_multiplier = mad_multiplier
      )

      gate_definition <- list(
        type = "polygon",
        parent = parent,
        channels = channels,
        coordinates = singlet_fit$coordinates,
        method = "automatic_singlet",
        algorithm = "iterative_fsc_diagonal_mad_band",
        parameters = list(
          n_events_sampled = length(event_index),
          mad_multiplier = mad_multiplier,
          x_quantiles = c(0.01, 0.99)
        ),
        diagnostics = singlet_fit$fit
      )

      self$update_gate(
        sample_ids = sample_id,
        gate_name = gate_name,
        gate_definition = gate_definition,
        save_as_default = save_as_default
      )

       if (isTRUE(verbose)) {
        elapsed_seconds <- round(
          as.numeric(difftime(Sys.time(), started_at, units = "secs")),
          digits = 1
        )

        message(
          "    Saved '",
          gate_name,
          "' for sample '",
          sample_id,
          "' using ",
          singlet_fit$fit$events_used,
          " retained events (",
          elapsed_seconds,
          " sec)."
        )
      }

    }

    if (isTRUE(verbose)) {
      message("Automatic singlet gating completed.")
    }

    invisible(self)
  }, 

  plot_singlet_gates = function(
      sample_ids = NULL,
      gate_name = "Singlets",
      n_plot = 10000L,
      layout = c(2L, 2L),
      file = NULL,
      width = 10,
      height = 8,
      point_alpha = 0.20,
      point_cex = 0.35
  ) {
    # Default: plot all samples in the project.
    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(self, sample_ids)

    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop("`gate_name` must be one non-empty character value.", call. = FALSE)
    }

    if (!is.numeric(n_plot) || length(n_plot) != 1L || n_plot < 100L) {
      stop("`n_plot` must be one number of at least 100.", call. = FALSE)
    }

    if (!is.numeric(layout) ||
        length(layout) != 2L ||
        any(layout < 1L)) {
      stop("`layout` must be c(rows, columns), with positive values.", call. = FALSE)
    }

    if (!is.null(file) &&
        (!is.character(file) || length(file) != 1L || !nzchar(file))) {
      stop("`file` must be NULL or one non-empty file path.", call. = FALSE)
    }

    if (!is.numeric(point_alpha) ||
        length(point_alpha) != 1L ||
        point_alpha < 0 ||
        point_alpha > 1) {
      stop("`point_alpha` must be between 0 and 1.", call. = FALSE)
    }

    # A PDF is recommended for batches. Each page holds `layout[1] * layout[2]`
    # samples, and base R automatically starts a new page as needed.
    pdf_opened <- FALSE

    if (!is.null(file)) {
      output_directory <- dirname(file)

      if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE)
      }

      grDevices::pdf(
        file = file,
        width = width,
        height = height,
        onefile = TRUE
      )

      pdf_opened <- TRUE

      on.exit(
        grDevices::dev.off(),
        add = TRUE
      )

      message(
        "Writing singlet-gate QC plots for ",
        length(sample_ids),
        " sample",
        if (length(sample_ids) == 1L) "" else "s",
        " to: ",
        normalizePath(file, winslash = "/", mustWork = FALSE)
      )
    }

    old_par <- graphics::par(no.readonly = TRUE)

    on.exit(
      graphics::par(old_par),
      add = TRUE
    )

    graphics::par(
      mfrow = as.integer(layout),
      mar = c(4.2, 4.2, 3.2, 1.2),
      oma = c(0, 0, 0, 0)
    )

    plot_summary <- vector("list", length(sample_ids))
    names(plot_summary) <- sample_ids

    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[[i]]

      # This assumes get_gates(sample_id, gate_name) returns the gate
      # definition saved by auto_gate_singlets().
      gate <- self$get_gates(
        sample_id = sample_id,
        gate_name = gate_name
      )

      if (is.null(gate)) {
        stop(
          "No gate named '", gate_name,
          "' was found for sample '", sample_id,
          "'. Run `auto_gate_singlets()` first.",
          call. = FALSE
        )
      }

      if (!identical(gate$type, "polygon")) {
        stop(
          "Gate '", gate_name,
          "' for sample '", sample_id,
          "' is not a polygon gate.",
          call. = FALSE
        )
      }

      channels <- gate$channels

      if (is.null(channels) ||
          !identical(names(channels), c("x", "y"))) {
        stop(
          "Gate '", gate_name,
          "' for sample '", sample_id,
          "' has invalid channel metadata.",
          call. = FALSE
        )
      }

      coordinates <- gate$coordinates

      if (!is.data.frame(coordinates) ||
          !all(c("x", "y") %in% names(coordinates))) {
        stop(
          "Gate '", gate_name,
          "' for sample '", sample_id,
          "' has invalid polygon coordinates.",
          call. = FALSE
        )
      }

      flow_frame <- self$get_flow_frame(sample_id)

      available_channels <- colnames(flowCore::exprs(flow_frame))
      missing_channels <- setdiff(unname(channels), available_channels)

      if (length(missing_channels) > 0L) {
        stop(
          "Sample '", sample_id,
          "' is missing gate channel(s): ",
          paste(missing_channels, collapse = ", "),
          call. = FALSE
        )
      }

      event_index <- cytoTerra_draw_event_sample(
        flow_frame = flow_frame,
        n_max = n_plot
      )

      event_data <- flowCore::exprs(flow_frame)[
        event_index,
        unname(channels),
        drop = FALSE
      ]

      x <- event_data[, channels[["x"]]]
      y <- event_data[, channels[["y"]]]

      valid <- is.finite(x) & is.finite(y)

      x <- x[valid]
      y <- y[valid]

      if (length(x) == 0L) {
        graphics::plot.new()
        graphics::title(
          main = paste0(sample_id, ": no finite FSC events")
        )
        next
      }

      point_color <- grDevices::adjustcolor(
        col = "black",
        alpha.f = point_alpha
      )

      graphics::plot(
        x = x,
        y = y,
        pch = 16,
        cex = point_cex,
        col = point_color,
        xlab = channels[["x"]],
        ylab = channels[["y"]],
        main = paste0(
          sample_id,
          " — ",
          gate_name,
          "\n",
          length(x),
          " sampled events"
        )
      )

      graphics::polygon(
        x = coordinates$x,
        y = coordinates$y,
        border = "#D62728",
        lwd = 2,
        col = grDevices::adjustcolor(
          col = "#D62728",
          alpha.f = 0.12
        )
      )

      graphics::legend(
        "topleft",
        legend = paste0("Automatic gate: ", gate_name),
        border = NA,
        bty = "n",
        lwd = 2,
        col = "#D62728",
        cex = 0.8
      )

      plot_summary[[sample_id]] <- data.frame(
        sample_id = sample_id,
        sample_name = self$get_sample_name(sample_id),
        gate_name = gate_name,
        channels = paste(unname(channels), collapse = " / "),
        events_plotted = length(x),
        events_used_to_fit = gate$diagnostics$events_used,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(file)) {
      message("Finished writing singlet-gate QC plots.")
    }

    invisible(do.call(rbind, plot_summary))
  },

  set_active_gate = function(sample_id, gate_name) {
     cytoTerra_validate_sample_ids(self, sample_id)

     gate <- self$get_gates(
       sample_id = sample_id,
       gate_name = gate_name
     )

     if (is.null(gate)) {
       stop(
         "No gate named '", gate_name,
         "' exists for sample '", sample_id, "'.",
         call. = FALSE
       )
     }

     self$active$sample_id <- sample_id
     self$active$gate_name <- gate_name

     message(
       "Active gate set to '", gate_name,
       "' for sample '", sample_id, "'."
     )

     invisible(self)
  },

# Automatically create a 1D positive or negative threshold gate.
# The threshold is fitted within the selected parent-gate population.
#
# If automatic mixture fitting or automatic threshold selection fails,
# `manual_threshold` is used when supplied.

auto_gate_1d = function(
    sample_ids = NULL,
    gate_name = NULL,
    parent = "root",
    channel,
    selection = c("positive", "negative"),
    probability = 0.90,
    n_max = 4000L,
    transform = c("log1p", "asinh", "none"),
    cofactor = 150,
    min_ashman_d = 2,
    manual_threshold = NULL,
    save_as_default = FALSE,
    on_threshold_failure = c("valley", "stop"),
    verbose = TRUE
) {

  selection <- match.arg(selection)
  transform <- match.arg(transform)
  on_threshold_failure <- match.arg(on_threshold_failure)

  # Default: apply to all samples.
  if (is.null(sample_ids)) {
    sample_ids <- self$get_sample_ids()
  }

  cytoTerra_validate_sample_ids(self, sample_ids)

  # Validate channel.
  if (missing(channel) ||
      !is.character(channel) ||
      length(channel) != 1L ||
      is.na(channel) ||
      !nzchar(channel)) {
    stop(
      "`channel` must be one non-empty channel name.",
      call. = FALSE
    )
  }

  # Default gate name.
  if (is.null(gate_name)) {
    gate_name <- paste0(
      channel,
      "_",
      if (identical(selection, "positive")) "Positive" else "Negative"
    )
  }

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  # Validate parent.
  if (!is.character(parent) ||
      length(parent) != 1L ||
      is.na(parent) ||
      !nzchar(parent)) {
    stop(
      "`parent` must be one non-empty gate name or 'root'.",
      call. = FALSE
    )
  }

  if (identical(parent, gate_name)) {
    stop(
      "`parent` cannot be the same as `gate_name`.",
      call. = FALSE
    )
  }

  # Validate probability.
  if (!is.numeric(probability) ||
      length(probability) != 1L ||
      is.na(probability) ||
      probability <= 0.5 ||
      probability >= 1) {
    stop(
      "`probability` must be strictly between 0.5 and 1.",
      call. = FALSE
    )
  }

  # Validate maximum number of fitting events.
  if (!is.numeric(n_max) ||
      length(n_max) != 1L ||
      is.na(n_max) ||
      n_max < 100L) {
    stop(
      "`n_max` must be one number of at least 100.",
      call. = FALSE
    )
  }

  n_max <- as.integer(n_max)

  # Validate cofactor for asinh fitting.
  if (identical(transform, "asinh")) {
    if (!is.numeric(cofactor) ||
        length(cofactor) != 1L ||
        is.na(cofactor) ||
        cofactor <= 0) {
      stop(
        "`cofactor` must be one positive number when ",
        "`transform = 'asinh'`.",
        call. = FALSE
      )
    }
  }

  # Validate minimum Ashman's D.
  if (!is.numeric(min_ashman_d) ||
      length(min_ashman_d) != 1L ||
      is.na(min_ashman_d) ||
      min_ashman_d <= 0) {
    stop(
      "`min_ashman_d` must be one positive number.",
      call. = FALSE
    )
  }

  # manual_threshold is defined on the transformed fitting scale:
  # - asinh: asinh(raw / cofactor)
  # - log1p: log1p(raw)
  # - none: raw FCS scale
  if (!is.null(manual_threshold)) {
    if (!is.numeric(manual_threshold) ||
        length(manual_threshold) != 1L ||
        is.na(manual_threshold) ||
        !is.finite(manual_threshold)) {
      stop(
        "`manual_threshold` must be NULL or one finite numeric value.",
        call. = FALSE
      )
    }
  }

  if (!is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(
      "Automatic 1D ", selection, " gating for ",
      length(sample_ids), " sample",
      if (length(sample_ids) == 1L) "" else "s",
      "; parent = '", parent,
      "'; channel = '", channel, "'."
    )

    if (!is.null(manual_threshold)) {
      message(
        "Manual fallback threshold = ", manual_threshold,
        " on the ", transform, " fitting scale."
      )
    }
  }

  for (sample_index in seq_along(sample_ids)) {

    sample_id <- sample_ids[[sample_index]]

    if (isTRUE(verbose)) {
      message(
        "[", sample_index, "/", length(sample_ids),
        "] Processing sample '", sample_id, "'..."
      )
    }

    flow_frame <- self$get_flow_frame(sample_id)
    expression_matrix <- flowCore::exprs(flow_frame)

    channel_info <- cytoTerra_resolve_channel(
      flow_frame = flow_frame,
      channel = channel
    )

    if (isTRUE(channel_info$is_scatter_or_time)) {
      warning(
        "Sample '", sample_id,
        "': '", channel_info$name,
        "' appears to be FSC, SSC, or Time rather than fluorescence.",
        call. = FALSE
      )
    }

    # Get the recursively evaluated parent population.
    if (identical(parent, "root")) {

      parent_mask <- rep(TRUE, nrow(expression_matrix))

    } else {

      parent_gate <- self$get_gates(
        sample_id = sample_id,
        gate_name = parent
      )

      if (is.null(parent_gate)) {
        stop(
          "Sample '", sample_id,
          "' has no parent gate named '", parent, "'.",
          call. = FALSE
        )
      }

      parent_mask <- self$evaluate_gate(
        sample_id = sample_id,
        gate_name = parent
      )
    }

    parent_event_index <- which(parent_mask)

    if (length(parent_event_index) == 0L) {
      stop(
        "Parent gate '", parent,
        "' contains zero events in sample '", sample_id, "'.",
        call. = FALSE
      )
    }

    # Sample only events that belong to the parent population.
    event_index <- if (length(parent_event_index) <= n_max) {
      parent_event_index
    } else {
      sample(parent_event_index, size = n_max, replace = FALSE)
    }

    raw_values <- expression_matrix[
      event_index,
      channel_info$name
    ]

    # Transform values before mixture fitting.
    if (identical(transform, "log1p")) {

      valid <- is.finite(raw_values) & raw_values >= 0
      model_values <- log1p(raw_values[valid])

    } else if (identical(transform, "asinh")) {

      valid <- is.finite(raw_values)
      model_values <- asinh(raw_values[valid] / cofactor)

    } else {

      valid <- is.finite(raw_values)
      model_values <- raw_values[valid]
    }

    if (length(model_values) < 100L) {
      stop(
        "Fewer than 100 valid events remain for sample '", sample_id,
        "' in parent gate '", parent, "'.",
        call. = FALSE
      )
    }

    # --------------------------------------------------------------
    # Try automatic fitting; use manual fallback when it fails.
    # --------------------------------------------------------------
    mixture <- tryCatch(
      cytoTerra_fit_bimodal_mixture(model_values),
      error = function(error) error
    )

    use_manual_fallback <- FALSE
    fallback_reason <- NULL
    threshold_result <- NULL

    if (inherits(mixture, "error")) {

      use_manual_fallback <- TRUE
      fallback_reason <- paste(
        "mixture-model fitting failed:",
        conditionMessage(mixture)
      )

    } else if (!is.finite(mixture$ashman_d) ||
               mixture$ashman_d < min_ashman_d) {

      use_manual_fallback <- TRUE
      fallback_reason <- paste0(
        "insufficient bimodal separation; Ashman's D = ",
        round(mixture$ashman_d, 2),
        ", required >= ", min_ashman_d
      )
    }

    # Try automatic threshold selection only if mixture quality was acceptable.
    if (!use_manual_fallback) {

      threshold_result <- tryCatch(
        cytoTerra_find_mixture_threshold(
          mixture = mixture,
          selection = selection,
          probability = probability,
          on_failure = on_threshold_failure
        ),
        error = function(error) error
      )

      if (inherits(threshold_result, "error")) {
        use_manual_fallback <- TRUE
        fallback_reason <- paste(
          "automatic threshold selection failed:",
          conditionMessage(threshold_result)
        )
      }
    }

    # Choose automatic or manual threshold.
    if (use_manual_fallback) {

      if (is.null(manual_threshold)) {
        stop(
          "Automatic gating failed for sample '", sample_id,
          "' in channel '", channel,
          "' within parent '", parent, "' (",
          fallback_reason, "). ",
          "Provide `manual_threshold` to use a fallback gate.",
          call. = FALSE
        )
      }

      threshold_transformed <- manual_threshold
      threshold_method <- "manual_fallback"
      probability_achieved <- NA_real_

    } else {

      threshold_transformed <- threshold_result$threshold
      threshold_method <- threshold_result$threshold_method
      probability_achieved <- threshold_result$probability_achieved
    }

    # Convert threshold from the fitting scale to raw FCS scale.
    threshold_raw <- if (identical(transform, "log1p")) {

      expm1(threshold_transformed)

    } else if (identical(transform, "asinh")) {

      cofactor * sinh(threshold_transformed)

    } else {

      threshold_transformed
    }

    # Store gate threshold in raw FCS coordinates.
    gate_definition <- list(
      type = "threshold",
      parent = parent,

      channels = channel_info$name,
      channel_input = channel_info$input,
      channel_description = channel_info$description,

      threshold = threshold_raw,
      direction = if (identical(selection, "positive")) {
        "above"
      } else {
        "below"
      },

      method = if (use_manual_fallback) {
        "manual_fallback_after_automatic_failure"
      } else {
        "automatic_bimodal_mixture"
      },

      selection = selection,
      probability = probability,
      transform = transform,
      cofactor = if (identical(transform, "asinh")) cofactor else NULL,

      parameters = list(
        n_events_parent = length(parent_event_index),
        n_events_sampled = length(event_index),
        n_events_valid = length(model_values),
        min_ashman_d = min_ashman_d,
        manual_threshold_available = !is.null(manual_threshold)
      ),

      diagnostics = list(
        threshold = threshold_raw,
        threshold_transformed = threshold_transformed,
        threshold_method = threshold_method,
        probability_achieved = probability_achieved,
        used_manual_fallback = use_manual_fallback,
        fallback_reason = fallback_reason,

        component_weights = if (inherits(mixture, "error")) {
          NA
        } else {
          mixture$weights
        },

        component_means = if (inherits(mixture, "error")) {
          NA
        } else {
          mixture$means
        },

        component_sds = if (inherits(mixture, "error")) {
          NA
        } else {
          mixture$sds
        },

        ashman_d = if (inherits(mixture, "error")) {
          NA_real_
        } else {
          mixture$ashman_d
        },

        iterations = if (inherits(mixture, "error")) {
          NA_integer_
        } else {
          mixture$iterations
        }
      )
    )

    self$update_gate(
      sample_ids = sample_id,
      gate_name = gate_name,
      gate_definition = gate_definition,
      save_as_default = save_as_default
    )

    if (isTRUE(verbose)) {

      if (use_manual_fallback) {
        message(
          "    Saved '", gate_name,
          "' using MANUAL fallback; threshold = ",
          signif(threshold_raw, 5),
          " raw (", threshold_transformed, " on ", transform,
          " scale); reason: ", fallback_reason, "."
        )
      } else {
        message(
          "    Saved '", gate_name,
          "' automatically; threshold = ",
          signif(threshold_raw, 5),
          " raw; Ashman's D = ",
          round(mixture$ashman_d, 2), "."
        )
      }
    }
  }

  if (isTRUE(verbose)) {
    message("Automatic 1D gating completed.")
  }

  invisible(self)
},

plot_1d_gates = function(
    sample_ids = NULL,
    gate_name,
    n_plot = 20000L,
    breaks = 100L,
    layout = c(1L, 1L),
    point_threshold_color = "#D62728",
    transform = c("asinh", "none"),
    cofactor = 150,
    verbose = TRUE
) {
    if (is.null(sample_ids)) {
        sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(self, sample_ids)

    if (!is.character(gate_name) || length(gate_name) != 1L ||
        is.na(gate_name) || !nzchar(gate_name)) {
        stop("`gate_name` must be one non-empty character value.",
             call. = FALSE)
    }

    if (!is.numeric(n_plot) || length(n_plot) != 1L || n_plot < 100L) {
        stop("`n_plot` must be one number of at least 100.",
             call. = FALSE)
    }

    if (!is.numeric(breaks) || length(breaks) != 1L || breaks < 10L) {
        stop("`breaks` must be one number of at least 10.",
             call. = FALSE)
    }

    if (!is.numeric(layout) || length(layout) != 2L || any(layout < 1L)) {
        stop("`layout` must be c(rows, columns), with positive values.",
             call. = FALSE)
    }

    transform <- match.arg(transform)

    if (!is.numeric(cofactor) || length(cofactor) != 1L ||
        is.na(cofactor) || cofactor <= 0) {
        stop("`cofactor` must be one positive numeric value.",
             call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    graphics::par(
        mfrow = as.integer(layout),
        mar = c(4.5, 4.5, 3.5, 1.5)
    )

    results <- vector("list", length(sample_ids))
    names(results) <- sample_ids

    for (sample_index in seq_along(sample_ids)) {
        sample_id <- sample_ids[[sample_index]]

        if (isTRUE(verbose)) {
            message(
                "[", sample_index, "/", length(sample_ids),
                "] Plotting gate '", gate_name,
                "' for sample '", sample_id, "'..."
            )
        }

        gate <- self$get_gates(
            sample_id = sample_id,
            gate_name = gate_name
        )

        if (is.null(gate)) {
            stop(
                "No gate named '", gate_name,
                "' was found for sample '", sample_id, "'.",
                call. = FALSE
            )
        }

        if (!identical(gate$type, "threshold")) {
            stop(
                "Gate '", gate_name, "' for sample '", sample_id,
                "' is type '", gate$type,
                "', but `plot_1d_gates()` requires type 'threshold'.",
                call. = FALSE
            )
        }

        if (!is.character(gate$channels) || length(gate$channels) != 1L) {
            stop(
                "Threshold gate '", gate_name, "' for sample '",
                sample_id,
                "' must contain exactly one character `channels` value.",
                call. = FALSE
            )
        }

        if (!is.numeric(gate$threshold) || length(gate$threshold) != 1L ||
            !is.finite(gate$threshold)) {
            stop(
                "Threshold gate '", gate_name, "' for sample '",
                sample_id, "' has an invalid `threshold`.",
                call. = FALSE
            )
        }

        flow_frame <- self$get_flow_frame(sample_id)
        exprs <- flowCore::exprs(flow_frame)
        channel <- gate$channels[[1]]

        if (!channel %in% colnames(exprs)) {
            stop(
                "Sample '", sample_id,
                "' does not contain channel '", channel, "'.",
                call. = FALSE
            )
        }

        # Restrict plotting to the parent population.
        parent_mask <- if (identical(gate$parent, "root")) {
            rep(TRUE, nrow(exprs))
        } else {
            self$evaluate_gate(
                sample_id = sample_id,
                gate_name = gate$parent
            )
        }

        parent_indices <- which(parent_mask)

        if (length(parent_indices) == 0L) {
            graphics::plot.new()
            graphics::title(
                main = paste0(sample_id, ": no cells in parent gate")
            )
            next
        }

        # Randomly subsample parent-gated events for display.
        event_index <- if (length(parent_indices) <= n_plot) {
            parent_indices
        } else {
            sample(parent_indices, size = as.integer(n_plot))
        }

        raw_values <- exprs[event_index, channel]
        raw_values <- raw_values[is.finite(raw_values)]

        if (length(raw_values) == 0L) {
            graphics::plot.new()
            graphics::title(
                main = paste0(
                    sample_id, ": no finite values in parent gate for ", channel
                )
            )
            next
        }

        # Display-only transformation, matching plot_1d_cells().
        if (identical(transform, "asinh")) {
            values <- asinh(raw_values / cofactor)
            threshold_for_plot <- asinh(gate$threshold / cofactor)
            x_axis_label <- paste0(
                channel, " expression [asinh(x / ", cofactor, ")]"
            )
        } else {
            values <- raw_values
            threshold_for_plot <- gate$threshold
            x_axis_label <- paste0(channel, " expression [raw scale]")
        }

        direction <- if (!is.null(gate$direction)) {
            gate$direction
        } else {
            NA_character_
        }

        selection <- if (!is.null(gate$selection)) {
            gate$selection
        } else if (identical(direction, "above")) {
            "positive"
        } else if (identical(direction, "below")) {
            "negative"
        } else {
            "unknown"
        }

        gate_mask <- self$evaluate_gate(
            sample_id = sample_id,
            gate_name = gate_name
        )

        events_in_parent <- sum(parent_mask)
        events_selected <- sum(gate_mask & parent_mask)

        selected_fraction <- if (events_in_parent > 0L) {
            events_selected / events_in_parent
        } else {
            NA_real_
        }

        marker <- if (!is.null(gate$channel_description) &&
            !is.na(gate$channel_description) &&
            nzchar(gate$channel_description)) {
            gate$channel_description
        } else {
            channel
        }

        graphics::hist(
            values,
            breaks = breaks,
            col = "#BDBDBD",
            border = "white",
            main = paste0(
                sample_id, " — ", marker, " (", selection, " gate)"
            ),
            xlab = x_axis_label,
            ylab = "Event count"
        )

        graphics::abline(
            v = threshold_for_plot,
            col = point_threshold_color,
            lwd = 3,
            lty = 2
        )

        threshold_label <- if (identical(transform, "asinh")) {
            c(
                paste0("Raw threshold: ", signif(gate$threshold, 5)),
                paste0(
                    "asinh threshold: ",
                    signif(threshold_for_plot, 5)
                )
            )
        } else {
            paste0("Threshold: ", signif(gate$threshold, 5))
        }

        graphics::legend(
            "topright",
            legend = c(
                threshold_label,
                paste0(selection, ": ", round(100 * selected_fraction, 1), "%")
            ),
            col = c(
                point_threshold_color,
                rep(NA, length(threshold_label))
            ),
            lwd = c(3, rep(NA, length(threshold_label))),
            lty = c(2, rep(NA, length(threshold_label))),
            bty = "n",
            cex = 0.85
        )

        results[[sample_id]] <- data.frame(
            sample_id = sample_id,
            sample_name = self$get_sample_name(sample_id),
            gate_name = gate_name,
            channel = channel,
            marker = marker,
            selection = selection,
            threshold = gate$threshold,
            display_transform = transform,
            cofactor = if (identical(transform, "asinh")) cofactor else NA_real_,
            events_plotted = length(values),
            events_in_parent = events_in_parent,
            selected_fraction = selected_fraction,
            stringsAsFactors = FALSE
        )
    }

    invisible(do.call(rbind, results))
},

  # Automatically create a 2D polygon gate from a Gaussian-mixture cluster.
  #'
  #' The selected cluster is converted to a regular cytoTerra polygon gate
  #' for each sample. Gates are stored in raw FCS parameter units and can be
  #' retrieved, modified, copied, reset, and evaluated like any other
  #' cytoTerra gate.
  #'
  #' @param sample_ids Sample IDs to process. Defaults to all samples.
  #' @param gate_name Name assigned to the resulting polygon gate.
  #' @param parent Parent gate name, or "root".
  #' @param x_channel Channel plotted/modelled on the x axis.
  #' @param y_channel Channel plotted/modelled on the y axis.
  #' @param selection_rule Positional rule used to choose a cluster.
  #' @param components Candidate number(s) of Gaussian-mixture components.
  #' @param transform Transformation used for clustering only.
  #' @param cofactor Asinh cofactor when `transform = "asinh"`.
  #' @param confidence_level Central cluster fraction used for polygon creation.
  #' @param n_fit Maximum parent events used to fit each mixture model.
  #' @param seed Optional random seed for reproducible fitting subsamples.
  #' @param save_as_default Whether to save generated gates as defaults.
  #' @param verbose Whether to report progress.
  #'
  #' @return Invisibly returns the CytoTerraProject object.
  auto_gate_2d_cluster = function(
      sample_ids = NULL,
      gate_name = NULL,
      parent = "root",
      x_channel,
      y_channel,
      selection_rule = c(
        "high_x_low_y",
        "high_x_high_y",
        "low_x_low_y",
        "low_x_high_y",
        "high_x",
        "low_x",
        "high_y",
        "low_y"
      ),
      components = 3L,
      transform = c("asinh", "log1p", "none"),
      cofactor = 150,
      confidence_level = 0.95,
      n_fit = 20000L,
      seed = NULL,
      save_as_default = FALSE,
      verbose = TRUE
  ) {

    selection_rule <- match.arg(selection_rule)
    transform <- match.arg(transform)

    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(self, sample_ids)

    if (missing(x_channel) ||
        !is.character(x_channel) ||
        length(x_channel) != 1L ||
        is.na(x_channel) ||
        !nzchar(x_channel)) {
      stop(
        "`x_channel` must be one non-empty channel name.",
        call. = FALSE
      )
    }

    if (missing(y_channel) ||
        !is.character(y_channel) ||
        length(y_channel) != 1L ||
        is.na(y_channel) ||
        !nzchar(y_channel)) {
      stop(
        "`y_channel` must be one non-empty channel name.",
        call. = FALSE
      )
    }

    if (identical(x_channel, y_channel)) {
      stop(
        "`x_channel` and `y_channel` must be different channels.",
        call. = FALSE
      )
    }

    if (is.null(gate_name)) {
      gate_name <- paste0(
        x_channel,
        "_",
        y_channel,
        "_",
        selection_rule
      )
    }

    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop(
        "`gate_name` must be one non-empty character value.",
        call. = FALSE
      )
    }

    if (!is.character(parent) ||
        length(parent) != 1L ||
        is.na(parent) ||
        !nzchar(parent)) {
      stop(
        "`parent` must be one non-empty gate name or 'root'.",
        call. = FALSE
      )
    }

    if (!is.numeric(confidence_level) ||
        length(confidence_level) != 1L ||
        is.na(confidence_level) ||
        confidence_level <= 0 ||
        confidence_level >= 1) {
      stop(
        "`confidence_level` must be strictly between 0 and 1.",
        call. = FALSE
      )
    }

    if (!is.numeric(n_fit) ||
        length(n_fit) != 1L ||
        is.na(n_fit) ||
        n_fit < 100L) {
      stop(
        "`n_fit` must be one number of at least 100.",
        call. = FALSE
      )
    }

    if (!is.numeric(cofactor) ||
        length(cofactor) != 1L ||
        is.na(cofactor) ||
        cofactor <= 0) {
      stop(
        "`cofactor` must be one positive number.",
        call. = FALSE
      )
    }

    if (!is.numeric(components) ||
        length(components) < 1L ||
        anyNA(components) ||
        any(components < 2L) ||
        any(components != as.integer(components))) {
      stop(
        "`components` must contain one or more integer values of at least 2.",
        call. = FALSE
      )
    }

    if (!is.null(seed) &&
        (!is.numeric(seed) ||
         length(seed) != 1L ||
         is.na(seed))) {
      stop(
        "`seed` must be NULL or one numeric value.",
        call. = FALSE
      )
    }

    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose)) {
      stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message(
        "Automatic 2D cluster gating for ",
        length(sample_ids), " sample",
        if (length(sample_ids) == 1L) "" else "s",
        "; x = '", x_channel,
        "'; y = '", y_channel,
        "'; selection rule = '", selection_rule, "'."
      )
    }

    for (sample_index in seq_along(sample_ids)) {
      sample_id <- sample_ids[[sample_index]]

      if (isTRUE(verbose)) {
        message(
          "[", sample_index, "/", length(sample_ids),
          "] Fitting 2D mixture for sample '", sample_id, "'..."
        )
      }

      flow_frame <- self$get_flow_frame(sample_id)

      x_info <- cytoTerra_resolve_channel(
        flow_frame = flow_frame,
        channel = x_channel
      )

      y_info <- cytoTerra_resolve_channel(
        flow_frame = flow_frame,
        channel = y_channel
      )

      if (identical(x_info$name, y_info$name)) {
        stop(
          "Sample '", sample_id,
          "': `x_channel` and `y_channel` resolved to the same parameter '",
          x_info$name, "'.",
          call. = FALSE
        )
      }

      if (isTRUE(x_info$is_scatter_or_time)) {
        warning(
          "Sample '", sample_id,
          "': x channel '", x_info$name,
          "' appears to be FSC, SSC, or Time. Continuing because it was ",
          "explicitly requested.",
          call. = FALSE
        )
      }

      if (isTRUE(y_info$is_scatter_or_time)) {
        warning(
          "Sample '", sample_id,
          "': y channel '", y_info$name,
          "' appears to be FSC, SSC, or Time. Continuing because it was ",
          "explicitly requested.",
          call. = FALSE
        )
      }

      # Parent-aware fitting: clustering is performed only on events that
      # already pass the declared parent gate.
      parent_mask <- if (identical(parent, "root")) {
        rep(TRUE, nrow(flowCore::exprs(flow_frame)))
      } else {
        self$evaluate_gate(
          sample_id = sample_id,
          gate_name = parent
        )
      }

      parent_event_count <- sum(parent_mask)

      if (parent_event_count < 100L) {
        stop(
          "Sample '", sample_id,
          "' has only ", parent_event_count,
          " events in parent '", parent,
          "'. At least 100 parent events are required for 2D automatic gating.",
          call. = FALSE
        )
      }

      expression_matrix <- flowCore::exprs(flow_frame)

      raw_values <- expression_matrix[
        parent_mask,
        c(x_info$name, y_info$name),
        drop = FALSE
      ]

      # Ensure the utility and diagnostics preserve the actual parameter names.
      colnames(raw_values) <- c(x_info$name, y_info$name)

      # Use a distinct deterministic seed per sample when a seed is supplied.
      # This retains reproducibility while avoiding identical subsamples across
      # different samples.
      sample_seed <- if (is.null(seed)) {
        NULL
      } else {
        as.integer(seed) + sample_index - 1L
      }

      mixture_result <- cytoTerra_fit_2d_mixture(
        raw_values = raw_values,
        transform = transform,
        cofactor = cofactor,
        components = components,
        n_fit = n_fit,
        seed = sample_seed
      )

      cluster_result <- cytoTerra_select_cluster(
        transformed_values = mixture_result$transformed_values,
        classification = mixture_result$classification,
        selection_rule = selection_rule
      )

      polygon_matrix <- cytoTerra_cluster_to_polygon(
        raw_values = mixture_result$raw_values,
        transformed_values = mixture_result$transformed_values,
        classification = mixture_result$classification,
        component = cluster_result$component,
        level = confidence_level
      )

      # cytoTerra polygon gates use generic x/y coordinate column names.
      # The actual FCS channels are recorded separately in `channels`.
      coordinates <- data.frame(
        x = polygon_matrix[, 1L],
        y = polygon_matrix[, 2L],
        stringsAsFactors = FALSE
      )

      selected_cluster_events <- sum(
        mixture_result$classification == cluster_result$component
      )

      gate_definition <- list(
        type = "polygon",
        parent = parent,

        channels = c(
          x = x_info$name,
          y = y_info$name
        ),

        channel_input = c(
          x = x_info$input,
          y = y_info$input
        ),

        channel_description = c(
          x = x_info$description,
          y = y_info$description
        ),

        coordinates = coordinates,

        method = "automatic_2d_gaussian_mixture",
        selection_rule = selection_rule,
        transform = transform,
        cofactor = cofactor,

        parameters = list(
          components_requested = as.integer(components),
          n_events_parent = parent_event_count,
          n_events_valid = nrow(mixture_result$raw_values),
          n_events_fit = length(mixture_result$fit_index),
          confidence_level = confidence_level
        ),

        diagnostics = list(
          selected_component = cluster_result$component,
          component_ids = cluster_result$component_ids,
          cluster_means_transformed = cluster_result$cluster_means,
          cluster_sizes = cluster_result$cluster_sizes,
          cluster_scores = cluster_result$scores,
          selected_cluster_events = selected_cluster_events,
          polygon_vertices = nrow(coordinates),
          transform = transform,
          cofactor = cofactor,
          model_name = mixture_result$fit$modelName,
          components_selected = mixture_result$fit$G,
          bic = mixture_result$fit$bic
        )
      )

      self$update_gate(
        sample_ids = sample_id,
        gate_name = gate_name,
        gate_definition = gate_definition,
        save_as_default = save_as_default
      )

      if (isTRUE(verbose)) {
        message(
          "    Saved '", gate_name,
          "'; x = '", x_info$name,
          "'; y = '", y_info$name,
          "'; component = ", cluster_result$component,
          "; selected-cluster events = ", selected_cluster_events,
          "; polygon vertices = ", nrow(coordinates), "."
        )
      }
    }

    if (isTRUE(verbose)) {
      message("Automatic 2D cluster gating completed.")
    }

    invisible(self)
  },

    # Plot saved cytoTerra polygon gates on parent-population events.
  #'
  #' Events are plotted in raw FCS parameter units because cytoTerra polygon
  #' coordinates are stored and evaluated in raw parameter units.
  #'
  #' @param sample_ids Sample IDs to plot. Defaults to all samples.
  #' @param gate_name Name of an existing polygon gate.
  #' @param n_plot Maximum number of parent events plotted per sample.
  #' @param layout Plot layout as c(rows, columns).
  #' @param file Optional PDF output path.
  #' @param width PDF width in inches when `file` is supplied.
  #' @param height PDF height in inches when `file` is supplied.
  #' @param parent_color Color for parent events outside the gate.
  #' @param gated_color Color for events passing the polygon gate.
  #' @param gate_color Color for the polygon boundary.
  #' @param parent_alpha Transparency for parent events.
  #' @param gated_alpha Transparency for gated events.
  #' @param point_cex Point size.
  #' @param verbose Whether to report progress.
  #'
  #' @return Invisibly returns a sample-level plotting summary data frame.
  plot_2d_gates = function(
      sample_ids = NULL,
      gate_name,
      n_plot = 20000L,
      layout = c(2L, 2L),
      file = NULL,
      width = 10,
      height = 8,
      parent_color = "#7F7F7F",
      gated_color = "#1F77B4",
      gate_color = "#D62728",
      parent_alpha = 0.18,
      gated_alpha = 0.45,
      point_cex = 0.35,
      transform = c("gate", "asinh", "log1p", "none"),
      cofactor = NULL,
      verbose = TRUE
  ) {
    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    transform <- match.arg(transform)

    cytoTerra_validate_sample_ids(self, sample_ids)

    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop("`gate_name` must be one non-empty character value.", call. = FALSE)
    }

    if (!is.numeric(n_plot) ||
        length(n_plot) != 1L ||
        is.na(n_plot) ||
        n_plot < 100L) {
      stop("`n_plot` must be one number of at least 100.", call. = FALSE)
    }

    if (!is.numeric(layout) ||
        length(layout) != 2L ||
        any(is.na(layout)) ||
        any(layout < 1L)) {
      stop(
        "`layout` must be c(rows, columns), with positive values.",
        call. = FALSE
      )
    }

    if (!is.null(file) &&
        (!is.character(file) ||
         length(file) != 1L ||
         is.na(file) ||
         !nzchar(file))) {
      stop("`file` must be NULL or one non-empty file path.", call. = FALSE)
    }

    if (!is.numeric(width) ||
        length(width) != 1L ||
        is.na(width) ||
        width <= 0) {
      stop("`width` must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(height) ||
        length(height) != 1L ||
        is.na(height) ||
        height <= 0) {
      stop("`height` must be one positive number.", call. = FALSE)
    }

    if (!is.numeric(parent_alpha) ||
        length(parent_alpha) != 1L ||
        is.na(parent_alpha) ||
        parent_alpha < 0 ||
        parent_alpha > 1) {
      stop("`parent_alpha` must be between 0 and 1.", call. = FALSE)
    }

    if (!is.numeric(gated_alpha) ||
        length(gated_alpha) != 1L ||
        is.na(gated_alpha) ||
        gated_alpha < 0 ||
        gated_alpha > 1) {
      stop("`gated_alpha` must be between 0 and 1.", call. = FALSE)
    }

    if (!is.numeric(point_cex) ||
        length(point_cex) != 1L ||
        is.na(point_cex) ||
        point_cex <= 0) {
      stop("`point_cex` must be one positive number.", call. = FALSE)
    }

    if (!is.null(cofactor) &&
        (!is.numeric(cofactor) ||
         length(cofactor) != 1L ||
         is.na(cofactor) ||
         !is.finite(cofactor) ||
         cofactor <= 0)) {
      stop(
        "`cofactor` must be NULL or one positive finite number.",
        call. = FALSE
      )
    }

    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose)) {
      stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(file)) {
      output_directory <- dirname(file)

      if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE)
      }

      grDevices::pdf(
        file = file,
        width = width,
        height = height,
        onefile = TRUE
      )

      on.exit(grDevices::dev.off(), add = TRUE)

      if (isTRUE(verbose)) {
        message(
          "Writing 2D gate QC plots for ",
          length(sample_ids),
          " sample",
          if (length(sample_ids) == 1L) "" else "s",
          " to: ",
          normalizePath(file, winslash = "/", mustWork = FALSE)
        )
      }
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    graphics::par(
      mfrow = as.integer(layout),
      mar = c(4.5, 4.5, 4.0, 1.5)
    )

    parent_point_color <- grDevices::adjustcolor(
      parent_color,
      alpha.f = parent_alpha
    )

    gated_point_color <- grDevices::adjustcolor(
      gated_color,
      alpha.f = gated_alpha
    )

    results <- vector("list", length(sample_ids))
    names(results) <- sample_ids

    for (sample_index in seq_along(sample_ids)) {
      sample_id <- sample_ids[[sample_index]]

      if (isTRUE(verbose)) {
        message(
          "[", sample_index, "/", length(sample_ids),
          "] Plotting 2D gate '", gate_name,
          "' for sample '", sample_id, "'..."
        )
      }

      gate <- self$get_gates(
        sample_id = sample_id,
        gate_name = gate_name
      )

      if (is.null(gate)) {
        stop(
          "No gate named '", gate_name,
          "' was found for sample '", sample_id, "'.",
          call. = FALSE
        )
      }

      if (!identical(gate$type, "polygon")) {
        stop(
          "Gate '", gate_name,
          "' for sample '", sample_id,
          "' is type '", gate$type,
          "', but `plot_2d_gates()` requires type 'polygon'.",
          call. = FALSE
        )
      }

      cytoTerra_validate_polygon_gate(gate)

      x_channel <- unname(gate$channels[["x"]])
      y_channel <- unname(gate$channels[["y"]])

      flow_frame <- self$get_flow_frame(sample_id)
      expression_matrix <- flowCore::exprs(flow_frame)

      missing_channels <- setdiff(
        c(x_channel, y_channel),
        colnames(expression_matrix)
      )

      if (length(missing_channels) > 0L) {
        stop(
          "Sample '", sample_id,
          "' does not contain polygon-gate channel(s): ",
          paste(missing_channels, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      # Full event masks are used for correct gate statistics.
      gate_mask <- self$evaluate_gate(
        sample_id = sample_id,
        gate_name = gate_name
      )

      parent_mask <- if (identical(gate$parent, "root")) {
        rep(TRUE, nrow(expression_matrix))
      } else {
        self$evaluate_gate(
          sample_id = sample_id,
          gate_name = gate$parent
        )
      }

      parent_index <- which(parent_mask)
      events_in_parent <- length(parent_index)
      events_in_gate <- sum(gate_mask)

      if (events_in_parent == 0L) {
        graphics::plot.new()
        graphics::title(
          main = paste0(
            sample_id,
            " — ",
            gate_name,
            "\nNo events in parent '",
            gate$parent,
            "'."
          )
        )

        results[[sample_id]] <- data.frame(
          sample_id = sample_id,
          sample_name = self$get_sample_name(sample_id),
          gate_name = gate_name,
          parent = gate$parent,
          x_channel = x_channel,
          y_channel = y_channel,
          events_in_parent = 0L,
          events_in_gate = 0L,
          fraction_of_parent = NA_real_,
          events_plotted = 0L,
          gated_events_plotted = 0L,
          stringsAsFactors = FALSE
        )

        next
      }

      # Sample only parent events for display. Event membership remains based
      # on the full gate evaluation, rather than being recalculated from the
      # displayed subset.
      plot_index <- if (events_in_parent <= n_plot) {
        parent_index
      } else {
        sample(
          parent_index,
          size = as.integer(n_plot),
          replace = FALSE
        )
      }

      # Extract raw values. Gate membership is always evaluated in raw units;
      # transformation below affects visualization only.
      x_raw <- expression_matrix[plot_index, x_channel]
      y_raw <- expression_matrix[plot_index, y_channel]
      plotted_gate_mask <- gate_mask[plot_index]

      # "gate" uses the transformation originally used when this automatic
      # gate was fitted. This is normally asinh for 2D fluorescence gating.
      plot_transform <- if (identical(transform, "gate")) {
        if (!is.null(gate$transform)) {
          gate$transform
        } else {
          "none"
        }
      } else {
        transform
      }

      plot_cofactor <- if (!is.null(cofactor)) {
        cofactor
      } else if (!is.null(gate$cofactor)) {
        gate$cofactor
      } else {
        150
      }

      if (!plot_transform %in% c("asinh", "log1p", "none")) {
        stop(
          "Gate '", gate_name,
          "' for sample '", sample_id,
          "' has unsupported transform '", plot_transform, "'.",
          call. = FALSE
        )
      }

      # log1p is only valid for non-negative fluorescence values under this
      # package's auto-gating convention.
      valid <- is.finite(x_raw) & is.finite(y_raw)

      if (identical(plot_transform, "log1p")) {
        valid <- valid & x_raw >= 0 & y_raw >= 0
      }

      x_raw <- x_raw[valid]
      y_raw <- y_raw[valid]
      plotted_gate_mask <- plotted_gate_mask[valid]

      transform_values <- function(values) {
        switch(
          plot_transform,
          asinh = asinh(values / plot_cofactor),
          log1p = log1p(values),
          none = values
        )
      }

      # Plot coordinates are transformed, but saved gate coordinates remain
      # raw values in the underlying gate definition.
      x <- transform_values(x_raw)
      y <- transform_values(y_raw)

      polygon_x <- transform_values(gate$coordinates$x)
      polygon_y <- transform_values(gate$coordinates$y)

      if (length(x) == 0L) {
        graphics::plot.new()
        graphics::title(
          main = paste0(
            sample_id,
            " — ",
            gate_name,
            "\nNo finite parent events for plotting."
          )
        )

        results[[sample_id]] <- data.frame(
          sample_id = sample_id,
          sample_name = self$get_sample_name(sample_id),
          gate_name = gate_name,
          parent = gate$parent,
          x_channel = x_channel,
          y_channel = y_channel,
          events_in_parent = events_in_parent,
          events_in_gate = events_in_gate,
          fraction_of_parent = events_in_gate / events_in_parent,
          events_plotted = 0L,
          gated_events_plotted = 0L,
          stringsAsFactors = FALSE
        )

        next
      }

      x_marker <- if (!is.null(gate$channel_description[["x"]]) &&
                      !is.na(gate$channel_description[["x"]]) &&
                      nzchar(gate$channel_description[["x"]])) {
        gate$channel_description[["x"]]
      } else {
        x_channel
      }

      y_marker <- if (!is.null(gate$channel_description[["y"]]) &&
                      !is.na(gate$channel_description[["y"]]) &&
                      nzchar(gate$channel_description[["y"]])) {
        gate$channel_description[["y"]]
      } else {
        y_channel
      }

      transform_label <- switch(
        plot_transform,
        asinh = paste0("asinh(x / ", signif(plot_cofactor, 5), ")"),
        log1p = "log1p(x)",
        none = "raw scale"
      )

      x_label <- paste0(
        x_marker,
        "\n(",
        x_channel,
        "; ",
        transform_label,
        ")"
      )

      y_label <- paste0(
        y_marker,
        "\n(",
        y_channel,
        "; ",
        transform_label,
        ")"
      ) 

      fraction_of_parent <- events_in_gate / events_in_parent

      outside_gate_mask <- !plotted_gate_mask

      graphics::plot(
        x = x[outside_gate_mask],
        y = y[outside_gate_mask],
        pch = 16,
        cex = point_cex,
        col = parent_point_color,
        xlab = x_label,
        ylab = y_label,
        main = paste0(
          sample_id,
          " — ",
          gate_name,
          "\n",
          events_in_gate, " / ", events_in_parent,
          " parent events (",
          formatC(
            100 * fraction_of_parent,
            format = "f",
            digits = 1
          ),
          "%)"
        ),
        xlim = range(c(x, polygon_x), finite = TRUE),
        ylim = range(c(y, polygon_y), finite = TRUE)
      )

      # Draw gate-positive parent events in a separate color.
      if (any(plotted_gate_mask)) {
        graphics::points(
          x = x[plotted_gate_mask],
          y = y[plotted_gate_mask],
          pch = 16,
          cex = point_cex,
          col = gated_point_color
        )
      }

      graphics::polygon(
        x = polygon_x,
        y = polygon_y,
        border = gate_color,
        lwd = 2,
        col = grDevices::adjustcolor(
          gate_color,
          alpha.f = 0.08
        )
      )

      graphics::legend(
        "topright",
        legend = c(
          paste0(
            "Parent outside gate: ",
            events_in_parent - events_in_gate
          ),
          paste0(
            "In gate: ",
            events_in_gate,
            " (",
            formatC(
              100 * fraction_of_parent,
              format = "f",
              digits = 1
            ),
            "% of parent)"
          ),
          paste0(
            "Display scale: ",
            transform_label
          ),
          "Saved polygon"
        ),
        pch = c(16, 16, NA, NA),
        col = c(parent_color, gated_color, NA, gate_color),
        lwd = c(NA, NA, NA, 2),
        pt.cex = c(0.8, 0.8, NA, NA),
        bty = "n",
        cex = 0.75
      )

      results[[sample_id]] <- data.frame(
        sample_id = sample_id,
        sample_name = self$get_sample_name(sample_id),
        gate_name = gate_name,
        parent = gate$parent,
        x_channel = x_channel,
        y_channel = y_channel,
        events_in_parent = events_in_parent,
        events_in_gate = events_in_gate,
        fraction_of_parent = fraction_of_parent,
        events_plotted = length(x),
        gated_events_plotted = sum(plotted_gate_mask),
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(file) && isTRUE(verbose)) {
      message("Finished writing 2D gate QC plots.")
    }

    invisible(do.call(rbind, results))
  }

)

