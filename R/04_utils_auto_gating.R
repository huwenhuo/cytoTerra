# Internal helpers for automatic gating.
#'
#' @noRd
cytoTerra_draw_event_sample <- function(flow_frame, n_max = 200000L) {
  n_events <- nrow(flow_frame)

  if (is.null(n_events) || n_events < 1L) {
    stop("The flowFrame contains no events.", call. = FALSE)
  }

  n_take <- min(as.integer(n_max), n_events)

  if (n_take == n_events) {
    return(seq_len(n_events))
  }

  sort(sample.int(n_events, size = n_take, replace = FALSE))
}


#' Fit a diagonal FSC-A/FSC-H singlet band.
#'
#' @noRd
cytoTerra_fit_singlet_polygon <- function(
    x,
    y,
    x_quantiles = c(0.01, 0.99),
    mad_multiplier = 4,
    iterations = 3L
) {
  valid <- is.finite(x) & is.finite(y) & x > 0 & y > 0

  x <- x[valid]
  y <- y[valid]

  if (length(x) < 100L) {
    stop(
      "Too few finite, positive events are available to fit a singlet gate.",
      call. = FALSE
    )
  }

  x_limits <- as.numeric(
    stats::quantile(
      x,
      probs = x_quantiles,
      names = FALSE,
      na.rm = TRUE
    )
  )

  keep <- x >= x_limits[[1]] & x <= x_limits[[2]]

  if (sum(keep) < 100L) {
    stop(
      "Too few events remain after FSC range filtering.",
      call. = FALSE
    )
  }

  # Iteratively fit the FSC-A/FSC-H diagonal and remove outliers.
  for (i in seq_len(iterations)) {
    fit_data <- data.frame(
      x = x[keep],
      y = y[keep]
    )

    fit <- stats::lm(
      y ~ x,
      data = fit_data
    )

    predicted_y <- stats::predict(
      fit,
      newdata = data.frame(x = x)
    )

    residuals <- y - predicted_y

    residual_mad <- stats::mad(
      residuals[keep],
      center = 0,
      constant = 1,
      na.rm = TRUE
    )

    if (!is.finite(residual_mad) || residual_mad <= 0) {
      stop(
        "Could not estimate a stable FSC-A/FSC-H singlet relationship.",
        call. = FALSE
      )
    }

    keep <- x >= x_limits[[1]] &
      x <= x_limits[[2]] &
      abs(residuals) <= mad_multiplier * residual_mad

    if (sum(keep) < 100L) {
      stop(
        "Automatic singlet fitting retained too few events.",
        call. = FALSE
      )
    }
  }

  # Refit using the final retained events.
  fit <- stats::lm(
    y ~ x,
    data = data.frame(
      x = x[keep],
      y = y[keep]
    )
  )

  predicted_y <- stats::predict(
    fit,
    newdata = data.frame(x = x)
  )

  residuals <- y - predicted_y

  residual_mad <- stats::mad(
    residuals[keep],
    center = 0,
    constant = 1,
    na.rm = TRUE
  )

  band_width <- mad_multiplier * residual_mad

  x_min <- x_limits[[1]]
  x_max <- x_limits[[2]]

  predicted_limits <- stats::predict(
    fit,
    newdata = data.frame(x = c(x_min, x_max))
  )

  list(
    coordinates = data.frame(
      x = c(x_min, x_max, x_max, x_min),
      y = c(
        predicted_limits[[1]] - band_width,
        predicted_limits[[2]] - band_width,
        predicted_limits[[2]] + band_width,
        predicted_limits[[1]] + band_width
      )
    ),
    fit = list(
      intercept = unname(stats::coef(fit)[["(Intercept)"]]),
      slope = unname(stats::coef(fit)[["x"]]),
      residual_mad = residual_mad,
      band_width = band_width,
      events_used = sum(keep)
    )
  )
}

# Fit a two-component univariate Gaussian mixture using EM.
#'
#' @noRd
cytoTerra_fit_bimodal_mixture <- function(
    values,
    max_iterations = 200L,
    tolerance = 1e-6,
    min_sd = 1e-8
) {
  values <- values[is.finite(values)]

  if (length(values) < 100L) {
    stop(
      "At least 100 finite events are required for automatic 1D gating.",
      call. = FALSE
    )
  }

  # Initial component centers from the lower and upper quartiles.
  means <- as.numeric(
    stats::quantile(values, probs = c(0.25, 0.75), names = FALSE)
  )

  overall_sd <- stats::sd(values)

  if (!is.finite(overall_sd) || overall_sd <= min_sd) {
    stop(
      "The selected channel has insufficient variation for bimodal gating.",
      call. = FALSE
    )
  }

  sds <- c(overall_sd / 2, overall_sd / 2)
  weights <- c(0.5, 0.5)

  previous_log_likelihood <- -Inf

  for (iteration in seq_len(max_iterations)) {
    density_1 <- weights[[1]] *
      stats::dnorm(values, mean = means[[1]], sd = sds[[1]])

    density_2 <- weights[[2]] *
      stats::dnorm(values, mean = means[[2]], sd = sds[[2]])

    total_density <- density_1 + density_2

    # Avoid division by zero for extreme numeric cases.
    total_density[total_density <= 0] <- .Machine$double.xmin

    posterior_1 <- density_1 / total_density
    posterior_2 <- density_2 / total_density

    responsibilities <- cbind(posterior_1, posterior_2)

    effective_n <- colSums(responsibilities)

    weights <- effective_n / length(values)

    means <- colSums(responsibilities * values) / effective_n

    sds <- sqrt(
      colSums(
        responsibilities *
          (values - rep(means, each = length(values)))^2
      ) / effective_n
    )

    sds <- pmax(sds, min_sd)

    log_likelihood <- sum(log(total_density))

    if (abs(log_likelihood - previous_log_likelihood) < tolerance) {
      break
    }

    previous_log_likelihood <- log_likelihood
  }

  # Always order the model as lower component, then upper component.
  component_order <- order(means)

  weights <- weights[component_order]
  means <- means[component_order]
  sds <- sds[component_order]

  # Ashman's D: a practical measure of separation for two populations.
  ashman_d <- abs(means[[2]] - means[[1]]) /
    sqrt((sds[[1]]^2 + sds[[2]]^2) / 2)

  list(
    weights = weights,
    means = means,
    sds = sds,
    ashman_d = ashman_d,
    iterations = iteration
  )
}


# Find a 1D mixture threshold for positive or negative selection.
#'
#' Attempts to find the requested posterior-probability boundary. If none can
#' be found between the component means, it can optionally fall back to the
#' minimum fitted mixture density ("the valley") between those means.
#'
#' @noRd
cytoTerra_find_mixture_threshold <- function(
    mixture,
    selection = c("positive", "negative"),
    probability = 0.90,
    on_failure = c("stop", "valley")
) {
  selection <- match.arg(selection)
  on_failure <- match.arg(on_failure)

  if (!is.numeric(probability) ||
      length(probability) != 1L ||
      is.na(probability) ||
      probability <= 0.5 ||
      probability >= 1) {
    stop(
      "`probability` must be one number strictly between 0.5 and 1.",
      call. = FALSE
    )
  }

  selected_component <- if (identical(selection, "positive")) 2L else 1L

  posterior_selected <- function(x) {
    density_1 <- mixture$weights[[1]] *
      stats::dnorm(x, mean = mixture$means[[1]], sd = mixture$sds[[1]])

    density_2 <- mixture$weights[[2]] *
      stats::dnorm(x, mean = mixture$means[[2]], sd = mixture$sds[[2]])

    total_density <- density_1 + density_2

    if (total_density <= 0) {
      return(NA_real_)
    }

    c(density_1, density_2)[selected_component] / total_density
  }

  mixture_density <- function(x) {
    mixture$weights[[1]] *
      stats::dnorm(x, mean = mixture$means[[1]], sd = mixture$sds[[1]]) +
      mixture$weights[[2]] *
      stats::dnorm(x, mean = mixture$means[[2]], sd = mixture$sds[[2]])
  }

  # The conventional gate boundary / density valley is searched between
  # the lower- and higher-expression component means.
  interval <- sort(mixture$means)

  root_function <- function(x) {
    posterior_selected(x) - probability
  }

  left_value <- root_function(interval[[1]])
  right_value <- root_function(interval[[2]])

  # Use the requested posterior-probability threshold when it is bracketed.
  if (is.finite(left_value) &&
      is.finite(right_value) &&
      left_value * right_value <= 0) {
    return(list(
      threshold = stats::uniroot(
        f = root_function,
        interval = interval
      )$root,
      threshold_method = "posterior_probability",
      probability_achieved = probability
    ))
  }

  if (identical(on_failure, "stop")) {
    stop(
      "Could not find a probability threshold between the two fitted ",
      "populations. Try a lower `probability`, use ",
      "`on_threshold_failure = 'valley'`, or inspect the channel plot.",
      call. = FALSE
    )
  }

  # Fallback: find the local low-density point between the populations.
  # This is often the most interpretable cutoff for a borderline mixture.
  valley <- stats::optimize(
    f = mixture_density,
    interval = interval
  )$minimum

  warning(
    "Could not find the requested posterior probability of ",
    probability,
    " between the fitted component means. Using the fitted density valley ",
    "instead. Review this gate visually.",
    call. = FALSE
  )

  list(
    threshold = valley,
    threshold_method = "density_valley_fallback",
    probability_achieved = posterior_selected(valley)
  )
}

# Resolve an FCS parameter from either its channel name ($PnN) or its
# marker/antibody description ($PnS).
#'
#' @noRd
cytoTerra_resolve_channel <- function(flow_frame, channel) {
  if (!is.character(channel) ||
      length(channel) != 1L ||
      is.na(channel) ||
      !nzchar(channel)) {
    stop(
      "`channel` must be one non-empty character value.",
      call. = FALSE
    )
  }

  expression_channels <- colnames(flowCore::exprs(flow_frame))

  parameter_data <- Biobase::pData(
    flowCore::parameters(flow_frame)
  )

  # The `name` column ordinarily holds the FCS $PnN channel name.
  # The `desc` column ordinarily holds the FCS $PnS description/marker.
  parameter_names <- parameter_data$name

  if (is.null(parameter_names)) {
    parameter_names <- expression_channels
  }

  parameter_descriptions <- parameter_data$desc

  if (is.null(parameter_descriptions)) {
    parameter_descriptions <- rep(NA_character_, length(parameter_names))
  }

  parameter_names <- as.character(parameter_names)
  parameter_descriptions <- as.character(parameter_descriptions)

  # First prefer an exact detector/channel-name match.
  name_matches <- which(parameter_names == channel)

  # Then permit an exact antibody/description match.
  desc_matches <- which(
    !is.na(parameter_descriptions) &
      nzchar(parameter_descriptions) &
      parameter_descriptions == channel
  )

  matches <- unique(c(name_matches, desc_matches))

  # If case differs, allow a case-insensitive match.
  if (length(matches) == 0L) {
    name_matches <- which(tolower(parameter_names) == tolower(channel))

    desc_matches <- which(
      !is.na(parameter_descriptions) &
        nzchar(parameter_descriptions) &
        tolower(parameter_descriptions) == tolower(channel)
    )

    matches <- unique(c(name_matches, desc_matches))
  }

  if (length(matches) == 0L) {
    available <- data.frame(
      channel = parameter_names,
      description = ifelse(
        is.na(parameter_descriptions) | !nzchar(parameter_descriptions),
        "<none>",
        parameter_descriptions
      ),
      stringsAsFactors = FALSE
    )

    available_text <- paste0(
      "  - ",
      available$channel,
      "  [",
      available$description,
      "]"
    )

    stop(
      "Could not find '", channel, "' as either an FCS channel name ",
      "or parameter description.\n\nAvailable parameters:\n",
      paste(available_text, collapse = "\n"),
      call. = FALSE
    )
  }

  if (length(matches) > 1L) {
    matching_text <- paste0(
      "  - ",
      parameter_names[matches],
      "  [",
      ifelse(
        is.na(parameter_descriptions[matches]) |
          !nzchar(parameter_descriptions[matches]),
        "<none>",
        parameter_descriptions[matches]
      ),
      "]"
    )

    stop(
      "'",
      channel,
      "' matched more than one FCS parameter. Use the exact channel name ",
      "instead:\n",
      paste(matching_text, collapse = "\n"),
      call. = FALSE
    )
  }

  match_index <- matches[[1]]
  resolved_name <- parameter_names[[match_index]]

  # Confirm it is truly available in the expression matrix.
  if (!resolved_name %in% expression_channels) {
    stop(
      "Resolved parameter '", resolved_name,
      "' is not present in the flowFrame expression matrix.",
      call. = FALSE
    )
  }

  # This is only informational. Do not use it as a strict rule because
  # fluorescence channel naming differs substantially across instruments.
  is_scatter_or_time <- grepl(
    pattern = "^(FSC|SSC|TIME)",
    x = resolved_name,
    ignore.case = TRUE
  )

  list(
    input = channel,
    name = resolved_name,
    description = parameter_descriptions[[match_index]],
    is_scatter_or_time = is_scatter_or_time
  )
}


# Fit a Gaussian mixture model to two-channel cytometry data.
#'
#' The model is fit on transformed values, but raw values are retained so
#' the resulting gate can later be stored and evaluated in raw FCS units.
#'
#' @param raw_values Numeric matrix with exactly two columns. Rows are events.
#' @param transform One of "asinh", "log1p", or "none".
#' @param cofactor Positive asinh cofactor when transform = "asinh".
#' @param components Integer vector of candidate mixture-component counts.
#' @param n_fit Maximum number of events used to fit the mixture model.
#' @param seed Optional integer random seed used when subsampling.
#'
#' @return A list containing raw/transformed valid data, the fitted Mclust
#'   object, and one cluster classification per valid input event.
#' @noRd
cytoTerra_fit_2d_mixture <- function(
    raw_values,
    transform = c("asinh", "log1p", "none"),
    cofactor = 150,
    components = 3L,
    n_fit = 20000L,
    seed = NULL
) {
  transform <- match.arg(transform)

  raw_values <- as.matrix(raw_values)

  if (!is.numeric(raw_values) || ncol(raw_values) != 2L) {
    stop(
      "`raw_values` must be a numeric matrix with exactly two columns.",
      call. = FALSE
    )
  }

  if (is.null(colnames(raw_values)) ||
      anyNA(colnames(raw_values)) ||
      any(!nzchar(colnames(raw_values)))) {
    stop(
      "`raw_values` must have two non-empty channel column names.",
      call. = FALSE
    )
  }

  if (!is.numeric(cofactor) ||
      length(cofactor) != 1L ||
      !is.finite(cofactor) ||
      cofactor <= 0) {
    stop("`cofactor` must be one finite number greater than zero.", call. = FALSE)
  }

  if (!is.numeric(components) ||
      length(components) < 1L ||
      any(!is.finite(components)) ||
      any(components < 2L) ||
      any(components != as.integer(components))) {
    stop(
      "`components` must contain one or more integer values of at least 2.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_fit) ||
      length(n_fit) != 1L ||
      !is.finite(n_fit) ||
      n_fit < 100L) {
    stop("`n_fit` must be one number of at least 100.", call. = FALSE)
  }

  if (!is.null(seed) &&
      (!is.numeric(seed) ||
       length(seed) != 1L ||
       !is.finite(seed))) {
    stop("`seed` must be NULL or one finite number.", call. = FALSE)
  }

  valid <- is.finite(raw_values[, 1L]) &
    is.finite(raw_values[, 2L])

  # log1p is not defined for values less than -1. We require non-negative
  # values here, which is generally appropriate for untransformed channels.
  if (identical(transform, "log1p")) {
    valid <- valid &
      raw_values[, 1L] >= 0 &
      raw_values[, 2L] >= 0
  }

  raw_valid <- raw_values[valid, , drop = FALSE]

  if (nrow(raw_valid) < 100L) {
    stop(
      "At least 100 valid events are required for 2D mixture gating; found ",
      nrow(raw_valid),
      ".",
      call. = FALSE
    )
  }

  transformed_valid <- switch(
    transform,
    asinh = asinh(raw_valid / cofactor),
    log1p = log1p(raw_valid),
    none = raw_valid
  )

  candidate_components <- sort(unique(as.integer(components)))
  candidate_components <- candidate_components[
    candidate_components < nrow(transformed_valid)
  ]

  if (length(candidate_components) == 0L) {
    stop(
      "No requested component count is smaller than the number of valid events.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  fit_index <- if (nrow(transformed_valid) <= n_fit) {
    seq_len(nrow(transformed_valid))
  } else {
    sample.int(
      n = nrow(transformed_valid),
      size = as.integer(n_fit),
      replace = FALSE
    )
  }

  fit_values <- transformed_valid[fit_index, , drop = FALSE]

  fit <- tryCatch(
    mclust::Mclust(
      data = fit_values,
      G = candidate_components,
      verbose = FALSE
    ),
    error = function(error) {
      stop(
        "2D Gaussian-mixture fitting failed: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  if (is.null(fit) || is.null(fit$classification)) {
    stop(
      "2D Gaussian-mixture fitting did not return event classifications.",
      call. = FALSE
    )
  }

  # The model may have been fitted on a subset, but gates need membership for
  # every valid parent event. Predict component membership for all valid rows.
  classification <- if (length(fit_index) == nrow(transformed_valid)) {
    fit$classification
  } else {
    prediction <- tryCatch(
      stats::predict(fit, newdata = transformed_valid),
      error = function(error) {
        stop(
          "Could not classify all events using the fitted 2D mixture: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )

    prediction$classification
  }

  if (length(classification) != nrow(transformed_valid) ||
      anyNA(classification)) {
    stop(
      "2D mixture classification is incomplete.",
      call. = FALSE
    )
  }

  list(
    raw_values = raw_valid,
    transformed_values = transformed_valid,
    valid = valid,
    fit_index = fit_index,
    fit = fit,
    classification = as.integer(classification),
    transform = transform,
    cofactor = cofactor
  )
}


# Select one mixture component using a transparent positional rule.
#'
#' @param transformed_values Numeric two-column matrix used for fitting.
#' @param classification Integer component assignment for each row.
#' @param selection_rule Cluster-selection rule.
#'
#' @return A list containing the selected component and cluster diagnostics.
#' @noRd
cytoTerra_select_cluster <- function(
    transformed_values,
    classification,
    selection_rule = c(
      "high_x_low_y",
      "high_x_high_y",
      "low_x_low_y",
      "low_x_high_y",
      "high_x",
      "low_x",
      "high_y",
      "low_y"
    )
) {
  selection_rule <- match.arg(selection_rule)

  transformed_values <- as.matrix(transformed_values)

  if (!is.numeric(transformed_values) ||
      ncol(transformed_values) != 2L ||
      any(!is.finite(transformed_values))) {
    stop(
      "`transformed_values` must be a finite numeric matrix with two columns.",
      call. = FALSE
    )
  }

  if (!is.numeric(classification) ||
      length(classification) != nrow(transformed_values) ||
      anyNA(classification) ||
      any(classification < 1L)) {
    stop(
      "`classification` must contain one positive component ID per event.",
      call. = FALSE
    )
  }

  component_ids <- sort(unique(as.integer(classification)))

  cluster_means <- do.call(
    rbind,
    lapply(component_ids, function(component_id) {
      colMeans(
        transformed_values[
          classification == component_id,
          ,
          drop = FALSE
        ]
      )
    })
  )

  rownames(cluster_means) <- as.character(component_ids)
  colnames(cluster_means) <- colnames(transformed_values)

  cluster_sizes <- vapply(
    component_ids,
    function(component_id) sum(classification == component_id),
    integer(1L)
  )

  # Rank-based scores avoid a numerical-range problem between channels and
  # remain stable if cluster means on one axis are nearly identical.
  x_rank <- rank(cluster_means[, 1L], ties.method = "average")
  y_rank <- rank(cluster_means[, 2L], ties.method = "average")

  score <- switch(
    selection_rule,
    high_x_low_y = x_rank - y_rank,
    high_x_high_y = x_rank + y_rank,
    low_x_low_y = -x_rank - y_rank,
    low_x_high_y = -x_rank + y_rank,
    high_x = x_rank,
    low_x = -x_rank,
    high_y = y_rank,
    low_y = -y_rank
  )

  # Resolve ties deterministically. For the intended CD45-high/live rule,
  # prefer higher x, then lower y, then the larger cluster.
  ordering <- switch(
    selection_rule,

    high_x_low_y = order(
      -score,
      -cluster_means[, 1L],
      cluster_means[, 2L],
      -cluster_sizes
    ),

    high_x_high_y = order(
      -score,
      -cluster_means[, 1L],
      -cluster_means[, 2L],
      -cluster_sizes
    ),

    low_x_low_y = order(
      -score,
      cluster_means[, 1L],
      cluster_means[, 2L],
      -cluster_sizes
    ),

    low_x_high_y = order(
      -score,
      cluster_means[, 1L],
      -cluster_means[, 2L],
      -cluster_sizes
    ),

    high_x = order(-score, -cluster_sizes),
    low_x = order(-score, -cluster_sizes),
    high_y = order(-score, -cluster_sizes),
    low_y = order(-score, -cluster_sizes)
  )

  selected_row <- ordering[[1L]]
  selected_component <- component_ids[[selected_row]]

  list(
    component = selected_component,
    selection_rule = selection_rule,
    component_ids = component_ids,
    cluster_means = cluster_means,
    cluster_sizes = cluster_sizes,
    scores = score
  )
}


# Convert a selected 2D cluster to a raw-scale convex-hull polygon.
#'
#' The cluster is first trimmed using Mahalanobis distance on the transformed
#' scale. The retained events are then converted to a convex hull in raw FCS
#' coordinates, so the resulting gate is directly evaluable on raw data.
#'
#' @param raw_values Finite raw two-channel event matrix.
#' @param transformed_values Corresponding transformed two-channel matrix.
#' @param classification Cluster assignment for each row.
#' @param component Selected component ID.
#' @param level Probability retained from the central fitted cluster.
#'
#' @return A numeric two-column matrix of raw polygon coordinates.
#' @noRd
cytoTerra_cluster_to_polygon <- function(
    raw_values,
    transformed_values,
    classification,
    component,
    level = 0.95
) {
  raw_values <- as.matrix(raw_values)
  transformed_values <- as.matrix(transformed_values)

  if (!is.numeric(raw_values) ||
      !is.numeric(transformed_values) ||
      ncol(raw_values) != 2L ||
      ncol(transformed_values) != 2L ||
      nrow(raw_values) != nrow(transformed_values)) {
    stop(
      "`raw_values` and `transformed_values` must be equally sized ",
      "numeric matrices with two columns.",
      call. = FALSE
    )
  }

  if (!is.numeric(classification) ||
      length(classification) != nrow(raw_values)) {
    stop(
      "`classification` must have one value per event.",
      call. = FALSE
    )
  }

  if (!is.numeric(component) ||
      length(component) != 1L ||
      is.na(component)) {
    stop("`component` must be one component ID.", call. = FALSE)
  }

  if (!is.numeric(level) ||
      length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "`level` must be one finite number strictly between zero and one.",
      call. = FALSE
    )
  }

  selected <- classification == component

  if (sum(selected) < 10L) {
    stop(
      "The selected 2D cluster contains fewer than 10 events.",
      call. = FALSE
    )
  }

  cluster_transformed <- transformed_values[selected, , drop = FALSE]
  cluster_raw <- raw_values[selected, , drop = FALSE]

  center <- colMeans(cluster_transformed)
  covariance <- stats::cov(cluster_transformed)

  # A small regularization keeps Mahalanobis distance usable for compact or
  # nearly collinear cytometry clusters.
  axis_scale <- apply(cluster_transformed, 2L, stats::sd)
  regularization <- max(axis_scale, 1e-6)^2 * 1e-8

  covariance <- covariance +
    diag(regularization, nrow = 2L, ncol = 2L)

  distance <- stats::mahalanobis(
    x = cluster_transformed,
    center = center,
    cov = covariance
  )

  keep <- distance <= stats::qchisq(level, df = 2L)

  # Do not fail solely because a very compact cluster had too few points
  # within the theoretical ellipse. In that case, use all cluster members.
  if (sum(keep) < 3L) {
    keep <- rep(TRUE, nrow(cluster_raw))
  }

  polygon_candidates <- cluster_raw[keep, , drop = FALSE]

  # Remove duplicate vertices before making a hull.
  polygon_candidates <- unique(polygon_candidates)

  if (nrow(polygon_candidates) < 3L) {
    stop(
      "The selected cluster does not contain enough distinct points to form ",
      "a polygon gate.",
      call. = FALSE
    )
  }

  hull_index <- grDevices::chull(
    polygon_candidates[, 1L],
    polygon_candidates[, 2L]
  )

  if (length(hull_index) < 3L) {
    stop(
      "The selected cluster is effectively one-dimensional and cannot form ",
      "a 2D polygon gate.",
      call. = FALSE
    )
  }

  coordinates <- polygon_candidates[hull_index, , drop = FALSE]
  colnames(coordinates) <- colnames(raw_values)

  coordinates
}
