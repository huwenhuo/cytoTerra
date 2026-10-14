# Internal helpers for gate validation.
#'
#' @noRd
cytoTerra_validate_sample_ids <- function(project, sample_ids) {
  if (!is.character(sample_ids) ||
      length(sample_ids) == 0L ||
      anyNA(sample_ids) ||
      any(!nzchar(sample_ids))) {
    print(sample_ids)
    print("\n")
    stop(
      "`sample_ids` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  unknown_samples <- setdiff(sample_ids, names(project$samples))

  if (length(unknown_samples) > 0L) {
    stop(
      "Unknown sample ID(s): ",
      paste(unknown_samples, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' @noRd
cytoTerra_validate_gate_definition <- function(
    gate_definition,
    gate_name
) {

  if (!is.list(gate_definition)) {
    stop(
      "`gate_definition` must be a list.",
      call. = FALSE
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

  # ---------------------------------------------------------------
  # Validate gate type.
  # ---------------------------------------------------------------
  if (is.null(gate_definition$type) ||
      !is.character(gate_definition$type) ||
      length(gate_definition$type) != 1L ||
      is.na(gate_definition$type) ||
      !nzchar(gate_definition$type)) {
    stop(
      "Gate '", gate_name,
      "' must include one non-empty character `type` value.",
      call. = FALSE
    )
  }

  allowed_types <- c("polygon", "threshold", "range")

  if (!gate_definition$type %in% allowed_types) {
    stop(
      "Gate '", gate_name,
      "' has unsupported type '", gate_definition$type,
      "'. Supported types are: ",
      paste(allowed_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Every non-root gate must have one parent.
  # ---------------------------------------------------------------
  if (is.null(gate_definition$parent) ||
      !is.character(gate_definition$parent) ||
      length(gate_definition$parent) != 1L ||
      is.na(gate_definition$parent) ||
      !nzchar(gate_definition$parent)) {
    stop(
      "Gate '", gate_name,
      "' must include one non-empty character `parent` value.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # All gate types require channel metadata.
  # Exact requirements are checked below by gate type.
  # ---------------------------------------------------------------
  if (is.null(gate_definition$channels) ||
      !is.character(gate_definition$channels) ||
      length(gate_definition$channels) == 0L ||
      anyNA(gate_definition$channels) ||
      any(!nzchar(gate_definition$channels))) {
    stop(
      "Gate '", gate_name,
      "' must include non-empty character `channels`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Polygon gate: named x/y channels and at least three vertices.
  # Coordinates must already be in raw FCS units when saved.
  # ---------------------------------------------------------------
  if (identical(gate_definition$type, "polygon")) {

    if (length(gate_definition$channels) != 2L ||
        !identical(names(gate_definition$channels), c("x", "y"))) {
      stop(
        "Polygon gate '", gate_name,
        "' requires exactly two named channels: ",
        "`c(x = ..., y = ...)`.",
        call. = FALSE
      )
    }

    if (identical(
      unname(gate_definition$channels[["x"]]),
      unname(gate_definition$channels[["y"]])
    )) {
      stop(
        "Polygon gate '", gate_name,
        "' requires different x and y channels.",
        call. = FALSE
      )
    }

    coordinates <- gate_definition$coordinates

    if (!is.data.frame(coordinates) ||
        !all(c("x", "y") %in% names(coordinates)) ||
        nrow(coordinates) < 3L) {
      stop(
        "Polygon gate '", gate_name,
        "' requires `coordinates`: a data frame with numeric ",
        "`x` and `y` columns and at least three vertices.",
        call. = FALSE
      )
    }

    if (!is.numeric(coordinates$x) ||
        !is.numeric(coordinates$y) ||
        any(!is.finite(coordinates$x)) ||
        any(!is.finite(coordinates$y))) {
      stop(
        "Polygon gate '", gate_name,
        "' has non-finite or non-numeric polygon coordinates.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------
  # Threshold gate: one channel, one finite threshold, and a direction.
  # ---------------------------------------------------------------
  if (identical(gate_definition$type, "threshold")) {

    if (length(gate_definition$channels) != 1L) {
      stop(
        "Threshold gate '", gate_name,
        "' requires exactly one channel.",
        call. = FALSE
      )
    }

    if (is.null(gate_definition$threshold) ||
        !is.numeric(gate_definition$threshold) ||
        length(gate_definition$threshold) != 1L ||
        is.na(gate_definition$threshold) ||
        !is.finite(gate_definition$threshold)) {
      stop(
        "Threshold gate '", gate_name,
        "' requires one finite numeric `threshold` value.",
        call. = FALSE
      )
    }

    if (is.null(gate_definition$direction) ||
        !is.character(gate_definition$direction) ||
        length(gate_definition$direction) != 1L ||
        is.na(gate_definition$direction) ||
        !gate_definition$direction %in% c("above", "below")) {
      stop(
        "Threshold gate '", gate_name,
        "' requires `direction = \"above\"` or `direction = \"below\"`.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------
  # Range gate: one channel and raw-scale lower/upper boundaries.
  #
  # lower may be -Inf, and upper may be Inf. This allows:
  #   - negative-only / low range:  -Inf to threshold
  #   - positive-only / high range: threshold to Inf
  # ---------------------------------------------------------------
  if (identical(gate_definition$type, "range")) {

    if (length(gate_definition$channels) != 1L) {
      stop(
        "Range gate '", gate_name,
        "' requires exactly one channel.",
        call. = FALSE
      )
    }

    if (is.null(gate_definition$lower) ||
        !is.numeric(gate_definition$lower) ||
        length(gate_definition$lower) != 1L ||
        is.na(gate_definition$lower) ||
        (!is.finite(gate_definition$lower) &&
         !identical(gate_definition$lower, -Inf))) {
      stop(
        "Range gate '", gate_name,
        "' requires a finite numeric `lower` bound or `-Inf`.",
        call. = FALSE
      )
    }

    if (is.null(gate_definition$upper) ||
        !is.numeric(gate_definition$upper) ||
        length(gate_definition$upper) != 1L ||
        is.na(gate_definition$upper) ||
        (!is.finite(gate_definition$upper) &&
         !identical(gate_definition$upper, Inf))) {
      stop(
        "Range gate '", gate_name,
        "' requires a finite numeric `upper` bound or `Inf`.",
        call. = FALSE
      )
    }

    if (gate_definition$lower >= gate_definition$upper) {
      stop(
        "Range gate '", gate_name,
        "' requires `lower < upper`.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' @noRd
cytoTerra_validate_gate_channels <- function(
    project,
    sample_id,
    gate_definition
) {
  expression_matrix <- flowCore::exprs(
    project$samples[[sample_id]]$flow_frame
  )

  missing_channels <- setdiff(
    unname(gate_definition$channels),
    colnames(expression_matrix)
  )

  if (length(missing_channels) > 0L) {
    stop(
      "Gate channels not found in sample '", sample_id, "': ",
      paste(missing_channels, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Test whether points lie inside a 2D polygon.
# Points exactly on a polygon edge are considered inside.
#'
#' @param x Numeric x coordinates.
#' @param y Numeric y coordinates.
#' @param polygon Numeric matrix/data.frame with two columns: x then y.
#'
#' @return A logical vector with one value per input point.
#' @noRd
cytoTerra_points_in_polygon <- function(x, y, polygon) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    stop(
      "`x` and `y` must be numeric vectors of equal length.",
      call. = FALSE
    )
  }

  polygon <- as.matrix(polygon)

  if (!is.numeric(polygon) ||
      ncol(polygon) != 2L ||
      nrow(polygon) < 3L ||
      any(!is.finite(polygon))) {
    stop(
      "`polygon` must be a finite numeric matrix with at least three ",
      "rows and exactly two columns.",
      call. = FALSE
    )
  }

  valid_points <- is.finite(x) & is.finite(y)
  inside <- rep(FALSE, length(x))

  if (!any(valid_points)) {
    return(inside)
  }

  x_valid <- x[valid_points]
  y_valid <- y[valid_points]

  # The polygon does not need to be explicitly closed. This loop connects
  # the final vertex to the first vertex automatically.
  vertex_count <- nrow(polygon)
  previous_vertex <- vertex_count

  inside_valid <- rep(FALSE, length(x_valid))
  on_boundary <- rep(FALSE, length(x_valid))

  coordinate_scale <- max(
    1,
    abs(polygon),
    abs(x_valid),
    abs(y_valid)
  )

  tolerance <- sqrt(.Machine$double.eps) * coordinate_scale

  for (current_vertex in seq_len(vertex_count)) {
    x1 <- polygon[previous_vertex, 1L]
    y1 <- polygon[previous_vertex, 2L]
    x2 <- polygon[current_vertex, 1L]
    y2 <- polygon[current_vertex, 2L]

    # First identify events lying directly on an edge. Those should be
    # retained by a gate, rather than being unstable due to floating point
    # rounding near the boundary.
    cross_product <- (x_valid - x1) * (y2 - y1) -
      (y_valid - y1) * (x2 - x1)

    within_x <- x_valid >= min(x1, x2) - tolerance &
      x_valid <= max(x1, x2) + tolerance

    within_y <- y_valid >= min(y1, y2) - tolerance &
      y_valid <= max(y1, y2) + tolerance

    on_boundary <- on_boundary |
      (abs(cross_product) <= tolerance & within_x & within_y)

    # Ray-casting algorithm for interior points.
    crosses_y <- (y1 > y_valid) != (y2 > y_valid)

    x_at_crossing <- x1 +
      (y_valid - y1) * (x2 - x1) / (y2 - y1)

    crossing <- crosses_y & x_valid < x_at_crossing

    inside_valid <- xor(inside_valid, crossing)

    previous_vertex <- current_vertex
  }

  inside[valid_points] <- inside_valid | on_boundary
  inside
}


# Validate a native cytoTerra polygon gate.
#'
#' @param gate Gate-definition list.
#' @noRd
cytoTerra_validate_polygon_gate <- function(gate) {
  if (!is.character(gate$channels) ||
      length(gate$channels) != 2L ||
      !all(c("x", "y") %in% names(gate$channels))) {
    stop(
      "A polygon gate requires named channels: c(x = ..., y = ...).",
      call. = FALSE
    )
  }

  if (anyNA(gate$channels) || any(!nzchar(gate$channels))) {
    stop(
      "Polygon-gate channel names must be non-empty.",
      call. = FALSE
    )
  }

  coordinates <- gate$coordinates

  if (!is.data.frame(coordinates) ||
      !all(c("x", "y") %in% names(coordinates)) ||
      nrow(coordinates) < 3L) {
    stop(
      "A polygon gate requires `coordinates`: a data frame with ",
      "`x` and `y` columns and at least three vertices.",
      call. = FALSE
    )
  }

  if (!is.numeric(coordinates$x) ||
      !is.numeric(coordinates$y) ||
      any(!is.finite(coordinates$x)) ||
      any(!is.finite(coordinates$y))) {
    stop(
      "Polygon coordinates must be finite numeric values.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Evaluate one cytoTerra gate within its parent population.
#'
#' The returned mask always has one value per original event in the sample.
#' A child gate is automatically restricted to events passing its parent.
#'
#' @param project A CytoTerraProject object.
#' @param sample_id Sample identifier.
#' @param gate_name Gate name to evaluate.
#' @param .ancestry Internal recursion state used to detect circular parents.
#'
#' @return Logical event-membership vector.
#' @noRd
cytoTerra_evaluate_gate <- function(
    project,
    sample_id,
    gate_name,
    .ancestry = character()
) {

  # ---------------------------------------------------------------
  # Root contains every event in the FCS file.
  # ---------------------------------------------------------------
  if (identical(gate_name, "root")) {
    flow_frame <- project$get_flow_frame(sample_id)

    return(
      rep(
        TRUE,
        nrow(flowCore::exprs(flow_frame))
      )
    )
  }

  # ---------------------------------------------------------------
  # Guard against invalid circular parent relationships.
  # ---------------------------------------------------------------
  if (gate_name %in% .ancestry) {
    stop(
      "Circular gate-parent relationship detected for sample '",
      sample_id,
      "': ",
      paste(c(.ancestry, gate_name), collapse = " -> "),
      ".",
      call. = FALSE
    )
  }

  gate <- project$get_gates(
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

  # ---------------------------------------------------------------
  # Validate generic gate metadata.
  # ---------------------------------------------------------------
  if (is.null(gate$type) ||
      !is.character(gate$type) ||
      length(gate$type) != 1L ||
      is.na(gate$type) ||
      !nzchar(gate$type)) {
    stop(
      "Gate '", gate_name,
      "' for sample '", sample_id,
      "' must contain one non-empty `type` value.",
      call. = FALSE
    )
  }

  allowed_types <- c("polygon", "threshold", "range")

  if (!gate$type %in% allowed_types) {
    stop(
      "Gate '", gate_name,
      "' has unsupported type '", gate$type,
      "'. Allowed types are: ",
      paste(allowed_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (is.null(gate$parent) ||
      !is.character(gate$parent) ||
      length(gate$parent) != 1L ||
      is.na(gate$parent) ||
      !nzchar(gate$parent)) {
    stop(
      "Gate '", gate_name,
      "' for sample '", sample_id,
      "' must contain one non-empty `parent` value.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Retrieve raw FCS event values and recursively evaluate parent.
  # Gates themselves are evaluated on raw FCS scale.
  # ---------------------------------------------------------------
  flow_frame <- project$get_flow_frame(sample_id)
  expression_matrix <- flowCore::exprs(flow_frame)
  event_count <- nrow(expression_matrix)

  parent_mask <- if (identical(gate$parent, "root")) {
    rep(TRUE, event_count)
  } else {
    cytoTerra_evaluate_gate(
      project = project,
      sample_id = sample_id,
      gate_name = gate$parent,
      .ancestry = c(.ancestry, gate_name)
    )
  }

  if (!is.logical(parent_mask) ||
      length(parent_mask) != event_count) {
    stop(
      "Parent gate '", gate$parent,
      "' returned an invalid membership mask for sample '",
      sample_id,
      "'.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Evaluate gate-specific membership.
  # ---------------------------------------------------------------
  gate_mask <- switch(

    gate$type,

    # -------------------------------------------------------------
    # One-dimensional threshold gate:
    #
    # direction = "above": channel >= threshold
    # direction = "below": channel <= threshold
    # -------------------------------------------------------------
    threshold = {

      if (is.null(gate$channels) ||
          !is.character(gate$channels) ||
          length(gate$channels) != 1L ||
          is.na(gate$channels[[1L]]) ||
          !nzchar(gate$channels[[1L]])) {
        stop(
          "Threshold gate '", gate_name,
          "' must contain exactly one non-empty FCS channel in `channels`.",
          call. = FALSE
        )
      }

      if (!is.numeric(gate$threshold) ||
          length(gate$threshold) != 1L ||
          is.na(gate$threshold) ||
          !is.finite(gate$threshold)) {
        stop(
          "Threshold gate '", gate_name,
          "' has an invalid finite `threshold`.",
          call. = FALSE
        )
      }

      if (!identical(gate$direction, "above") &&
          !identical(gate$direction, "below")) {
        stop(
          "Threshold gate '", gate_name,
          "' must have `direction = \"above\"` or `direction = \"below\"`.",
          call. = FALSE
        )
      }

      channel <- gate$channels[[1L]]

      if (!channel %in% colnames(expression_matrix)) {
        stop(
          "Sample '", sample_id,
          "' does not contain channel '", channel,
          "' required by threshold gate '", gate_name, "'.",
          call. = FALSE
        )
      }

      values <- expression_matrix[, channel]

      if (identical(gate$direction, "above")) {
        is.finite(values) & values >= gate$threshold
      } else {
        is.finite(values) & values <= gate$threshold
      }
    },

    # -------------------------------------------------------------
    # One-dimensional range gate:
    #
    # lower <= channel <= upper
    #
    # `lower` may be -Inf and `upper` may be Inf.
    # -------------------------------------------------------------
    range = {

      if (is.null(gate$channels) ||
          !is.character(gate$channels) ||
          length(gate$channels) != 1L ||
          is.na(gate$channels[[1L]]) ||
          !nzchar(gate$channels[[1L]])) {
        stop(
          "Range gate '", gate_name,
          "' must contain exactly one non-empty FCS channel in `channels`.",
          call. = FALSE
        )
      }

      if (!is.numeric(gate$lower) ||
          length(gate$lower) != 1L ||
          is.na(gate$lower) ||
          (!is.finite(gate$lower) &&
           !identical(gate$lower, -Inf))) {
        stop(
          "Range gate '", gate_name,
          "' must contain a finite `lower` bound or `-Inf`.",
          call. = FALSE
        )
      }

      if (!is.numeric(gate$upper) ||
          length(gate$upper) != 1L ||
          is.na(gate$upper) ||
          (!is.finite(gate$upper) &&
           !identical(gate$upper, Inf))) {
        stop(
          "Range gate '", gate_name,
          "' must contain a finite `upper` bound or `Inf`.",
          call. = FALSE
        )
      }

      if (gate$lower >= gate$upper) {
        stop(
          "Range gate '", gate_name,
          "' requires `lower < upper`.",
          call. = FALSE
        )
      }

      channel <- gate$channels[[1L]]

      if (!channel %in% colnames(expression_matrix)) {
        stop(
          "Sample '", sample_id,
          "' does not contain channel '", channel,
          "' required by range gate '", gate_name, "'.",
          call. = FALSE
        )
      }

      values <- expression_matrix[, channel]

      # Excludes NA, NaN, Inf, and -Inf event values.
      is.finite(values) &
        values >= gate$lower &
        values <= gate$upper
    },

    # -------------------------------------------------------------
    # Two-dimensional polygon gate.
    # -------------------------------------------------------------
    polygon = {

      cytoTerra_validate_polygon_gate(gate)

      x_channel <- unname(gate$channels[["x"]])
      y_channel <- unname(gate$channels[["y"]])

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

      cytoTerra_points_in_polygon(
        x = expression_matrix[, x_channel],
        y = expression_matrix[, y_channel],
        polygon = as.matrix(
          gate$coordinates[, c("x", "y"), drop = FALSE]
        )
      )
    }
  )

  if (!is.logical(gate_mask) ||
      length(gate_mask) != event_count) {
    stop(
      "Gate '", gate_name,
      "' returned an invalid membership mask for sample '",
      sample_id,
      "'.",
      call. = FALSE
    )
  }

  # A child gate can only include events already included by its parent.
  parent_mask & gate_mask
}

cytoTerra_get_parent_gate_cells <- function(
  project,
  sample_id,
  parent_gate
) {

  # Confirm that this gate exists for this particular sample.
  gate <- project$get_gates(
    sample_id = sample_id,
    gate_name = parent_gate
  )

  if (is.null(gate)) {
    stop(
      "Gate '", parent_gate,
      "' was not found for sample '", sample_id, "'.",
      call. = FALSE
    )
  }

  # Get raw FCS event data in the same event order used for gating.
  flow_frame <- project$get_flow_frame(sample_id)

  events <- as.data.frame(
    flowCore::exprs(flow_frame)
  )

  # Existing CytoTerra gate evaluator should resolve the full path:
  # raw events -> all ancestors -> requested parent_gate.
  gate_mask <- project$evaluate_gate(
    sample_id = sample_id,
    gate_name = parent_gate
  )

  # Accommodate either a logical membership vector or integer indices.
  if (is.logical(gate_mask)) {

    if (length(gate_mask) != nrow(events)) {
      stop(
        "The gate evaluation result for '", parent_gate,
        "' in sample '", sample_id,
        "' does not match the number of FCS events.",
        call. = FALSE
      )
    }

    return(events[which(gate_mask), , drop = FALSE])
  }

  if (is.numeric(gate_mask) || is.integer(gate_mask)) {
    return(events[gate_mask, , drop = FALSE])
  }

  stop(
    "Gate evaluation for '", parent_gate,
    "' in sample '", sample_id,
    "' did not return logical event membership or event indices.",
    call. = FALSE
  )
}

