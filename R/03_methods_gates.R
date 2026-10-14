# Gate storage and editing methods.
#'
#' @noRd
cytoTerra_methods_gates <- list(

update_gate = function(
    sample_ids = NULL,
    gate_name,
    gate_definition,
    save_as_default = FALSE
) {

  if (is.null(sample_ids)) {
    sample_ids <- self$get_sample_ids()
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

  cytoTerra_validate_sample_ids(self, sample_ids)

  cytoTerra_validate_gate_definition(gate_definition, gate_name)

  for (sample_id in sample_ids) {
    cytoTerra_validate_gate_channels(
      project = self,
      sample_id = sample_id,
      gate_definition = gate_definition
    )

    sample <- self$samples[[sample_id]]
    old_gate <- sample$gates[[gate_name]]

    # Preserve a prior version only if this gate already exists.
    if (!is.null(old_gate)) {
      sample$gate_history[[gate_name]] <- c(
        sample$gate_history[[gate_name]],
        list(old_gate)
      )
    }

    new_gate <- gate_definition
    new_gate$name <- gate_name
    new_gate$updated_at <- Sys.time()

    if (is.null(new_gate$created_at)) {
      new_gate$created_at <- Sys.time()
    }

    if (is.null(new_gate$method)) {
      new_gate$method <- "manual"
    }

    sample$gates[[gate_name]] <- new_gate

    if (isTRUE(save_as_default)) {
      sample$default_gates[[gate_name]] <- new_gate
    }

    sample$gate_membership <- list()
    self$samples[[sample_id]] <- sample
  }

  self$metadata$updated_at <- Sys.time()

  invisible(self)
},


  # Return all gates or one named gate for a sample.
  get_gates = function(sample_id, gate_name = NULL) {
    cytoTerra_validate_sample_ids(self, sample_id)

    gates <- self$samples[[sample_id]]$gates

    if (is.null(gate_name)) {
      return(gates)
    }

    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop(
        "`gate_name` must be NULL or one non-empty character value.",
        call. = FALSE
      )
    }

    gates[[gate_name]]
  },


  # Copy one gate definition from one sample to selected target samples.
  copy_gate = function(
      from_sample,
      gate_name,
      to_samples,
      overwrite = FALSE
  ) {
    cytoTerra_validate_sample_ids(self, from_sample)
    cytoTerra_validate_sample_ids(self, to_samples)

    source_gate <- self$samples[[from_sample]]$gates[[gate_name]]

    if (is.null(source_gate)) {
      stop(
        "Gate '", gate_name, "' does not exist in sample '",
        from_sample, "'.",
        call. = FALSE
      )
    }

    for (sample_id in to_samples) {
      existing_gate <- self$samples[[sample_id]]$gates[[gate_name]]

      if (!is.null(existing_gate) && !isTRUE(overwrite)) {
        warning(
          "Skipping sample '", sample_id,
          "': gate '", gate_name,
          "' already exists. Use `overwrite = TRUE` to replace it.",
          call. = FALSE
        )
        next
      }

      copied_gate <- source_gate
      copied_gate$method <- "copied"
      copied_gate$copied_from <- from_sample
      copied_gate$updated_at <- Sys.time()

      self$update_gate(
        sample_ids = sample_id,
        gate_name = gate_name,
        gate_definition = copied_gate,
        save_as_default = FALSE
      )
    }

    invisible(self)
  },


  # Save an existing sample gate as its reset/default gate.
  set_default_gate = function(sample_ids, gate_name) {
    cytoTerra_validate_sample_ids(self, sample_ids)

    for (sample_id in sample_ids) {
      sample <- self$samples[[sample_id]]
      gate <- sample$gates[[gate_name]]

      if (is.null(gate)) {
        warning(
          "Sample '", sample_id,
          "' has no gate named '", gate_name, "'.",
          call. = FALSE
        )
        next
      }

      sample$default_gates[[gate_name]] <- gate
      self$samples[[sample_id]] <- sample
    }

    self$metadata$updated_at <- Sys.time()

    invisible(self)
  },


  # Restore a gate from the saved default gate definition.
  reset_gate = function(sample_ids, gate_name) {
    cytoTerra_validate_sample_ids(self, sample_ids)

    for (sample_id in sample_ids) {
      default_gate <- self$samples[[sample_id]]$default_gates[[gate_name]]

      if (is.null(default_gate)) {
        warning(
          "Sample '", sample_id,
          "' has no default gate named '", gate_name, "'.",
          call. = FALSE
        )
        next
      }

      default_gate$method <- "reset_to_default"
      default_gate$updated_at <- Sys.time()

      self$update_gate(
        sample_ids = sample_id,
        gate_name = gate_name,
        gate_definition = default_gate,
        save_as_default = FALSE
      )
    }

    invisible(self)
  },


  # Set the gate currently selected in an interactive workflow.
  set_current_gate = function(gate_name) {
    if (!is.character(gate_name) ||
        length(gate_name) != 1L ||
        is.na(gate_name) ||
        !nzchar(gate_name)) {
      stop(
        "`gate_name` must be one non-empty character value.",
        call. = FALSE
      )
    }

    self$active$gate_name <- gate_name
    self$metadata$updated_at <- Sys.time()

    invisible(self)
  },


  # Get the currently selected gate name.
  get_current_gate = function() {
    self$active$gate_name
  },

  evaluate_gate = function(sample_id, gate_name) {
     cytoTerra_validate_sample_ids(self, sample_id)

     if (!is.character(gate_name) ||
         length(gate_name) != 1L ||
         is.na(gate_name) ||
         !nzchar(gate_name)) {
       stop("`gate_name` must be one non-empty character value.", call. = FALSE)
     }

     cytoTerra_evaluate_gate(
       project = self,
       sample_id = sample_id,
       gate_name = gate_name
     )
  },

get_gate_stats = function(
    sample_ids = NULL,
    gate_name
) {

  # Default: calculate statistics for every project sample.
  if (is.null(sample_ids)) {
    sample_ids <- self$get_sample_ids()
  }

  cytoTerra_validate_sample_ids(self, sample_ids)

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  gate_stats <- lapply(sample_ids, function(sample_id) {

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

    # Child mask is recursively gated: it includes membership in every
    # ancestor gate as well as membership in gate_name itself.
    child_mask <- self$evaluate_gate(
      sample_id = sample_id,
      gate_name = gate_name
    )

    # Evaluate the immediate parent to calculate fraction of parent.
    parent_mask <- if (identical(gate$parent, "root")) {

      flow_frame <- self$get_flow_frame(sample_id)

      rep(
        TRUE,
        nrow(flowCore::exprs(flow_frame))
      )

    } else {

      self$evaluate_gate(
        sample_id = sample_id,
        gate_name = gate$parent
      )
    }

    events_in_parent <- sum(parent_mask)
    events_in_gate <- sum(child_mask)
    events_total <- length(child_mask)

    data.frame(
      sample_id = sample_id,
      sample_name = self$get_sample_name(sample_id),
      gate_name = gate_name,
      parent = gate$parent,
      events_in_gate = events_in_gate,
      events_in_parent = events_in_parent,
      events_total = events_total,

      fraction_of_parent = if (events_in_parent > 0L) {
        events_in_gate / events_in_parent
      } else {
        NA_real_
      },

      fraction_of_all_events = if (events_total > 0L) {
        events_in_gate / events_total
      } else {
        NA_real_
      },

      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, gate_stats)
},

# -------------------------------------------------------------------
# Create a new polygon child gate.
# Coordinates are always stored internally on the raw FCS scale.
# -------------------------------------------------------------------
add_polygon_gate = function(
    sample_id,
    gate_name,
    parent = "root",
    x_channel,
    y_channel,
    coordinates,
    coordinate_scale = c("raw", "transformed"),
    transform = c("asinh", "log1p", "none"),
    cofactor = 150,
    x_channel_description = NULL,
    y_channel_description = NULL,
    note = NULL,
    save_as_default = FALSE
) {

  coordinate_scale <- match.arg(coordinate_scale)
  transform <- match.arg(transform)

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (identical(gate_name, "root")) {
    stop(
      "`gate_name` cannot be 'root'; that name is reserved.",
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

  if (!is.character(x_channel) ||
      length(x_channel) != 1L ||
      is.na(x_channel) ||
      !nzchar(x_channel)) {
    stop(
      "`x_channel` must be one non-empty FCS channel name.",
      call. = FALSE
    )
  }

  if (!is.character(y_channel) ||
      length(y_channel) != 1L ||
      is.na(y_channel) ||
      !nzchar(y_channel)) {
    stop(
      "`y_channel` must be one non-empty FCS channel name.",
      call. = FALSE
    )
  }

  if (identical(x_channel, y_channel)) {
    stop(
      "`x_channel` and `y_channel` must be different.",
      call. = FALSE
    )
  }

  if (!is.data.frame(coordinates) ||
      !all(c("x", "y") %in% names(coordinates)) ||
      nrow(coordinates) < 3L) {
    stop(
      "`coordinates` must be a data frame with numeric x/y columns ",
      "and at least three vertices.",
      call. = FALSE
    )
  }

  coordinates <- coordinates[, c("x", "y"), drop = FALSE]

  if (!is.numeric(coordinates$x) ||
      !is.numeric(coordinates$y) ||
      any(!is.finite(coordinates$x)) ||
      any(!is.finite(coordinates$y))) {
    stop(
      "`coordinates$x` and `coordinates$y` must be finite numeric values.",
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

  if (!is.null(note) &&
      (!is.character(note) ||
       length(note) != 1L ||
       is.na(note) ||
       !nzchar(note))) {
    stop(
      "`note` must be NULL or one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.logical(save_as_default) ||
      length(save_as_default) != 1L ||
      is.na(save_as_default)) {
    stop(
      "`save_as_default` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  existing_gates <- self$get_gates(sample_id)

  if (!is.null(existing_gates[[gate_name]])) {
    stop(
      "A gate named '", gate_name,
      "' already exists for sample '", sample_id,
      "'. Use `replace_polygon_gate()` to revise it.",
      call. = FALSE
    )
  }

  if (!identical(parent, "root") &&
      is.null(existing_gates[[parent]])) {
    stop(
      "Parent gate '", parent,
      "' does not exist for sample '", sample_id,
      "'.",
      call. = FALSE
    )
  }

  # Convert coordinates drawn on a transformed display back to raw FCS units.
  raw_coordinates <- if (identical(coordinate_scale, "raw")) {
    coordinates
  } else {
    switch(
      transform,
      asinh = data.frame(
        x = cofactor * sinh(coordinates$x),
        y = cofactor * sinh(coordinates$y)
      ),
      log1p = data.frame(
        x = expm1(coordinates$x),
        y = expm1(coordinates$y)
      ),
      none = coordinates
    )
  }

  if (any(!is.finite(raw_coordinates$x)) ||
      any(!is.finite(raw_coordinates$y))) {
    stop(
      "Coordinate conversion produced non-finite raw FCS values.",
      call. = FALSE
    )
  }

  if (is.null(x_channel_description)) {
    x_channel_description <- x_channel
  }

  if (is.null(y_channel_description)) {
    y_channel_description <- y_channel
  }

  new_gate <- list(
    type = "polygon",
    parent = parent,

    channels = c(
      x = x_channel,
      y = y_channel
    ),

    channel_input = c(
      x = x_channel,
      y = y_channel
    ),

    channel_description = c(
      x = x_channel_description,
      y = y_channel_description
    ),

    coordinates = raw_coordinates,

    # Display metadata. Gate membership remains evaluated in raw units.
    transform = transform,
    cofactor = cofactor,

    method = "manual_polygon",

    manual_creation = list(
      created_at = Sys.time(),
      coordinate_scale_supplied = coordinate_scale,
      transform_supplied = if (identical(
        coordinate_scale,
        "transformed"
      )) {
        transform
      } else {
        NULL
      },
      cofactor_supplied = if (identical(
        coordinate_scale,
        "transformed"
      ) && identical(transform, "asinh")) {
        cofactor
      } else {
        NULL
      },
      note = note
    )
  )

  self$update_gate(
    sample_ids = sample_id,
    gate_name = gate_name,
    gate_definition = new_gate,
    save_as_default = save_as_default
  )

  message(
    "Created polygon gate '", gate_name,
    "' under parent '", parent,
    "' for sample '", sample_id, "'."
  )

  invisible(self)
},


# -------------------------------------------------------------------
# Create a new 1D range gate.
#
# A range gate includes events satisfying:
#     lower <= channel value <= upper
#
# Bounds are stored internally on the raw FCS scale.
# -------------------------------------------------------------------
add_range_gate = function(
    sample_id,
    gate_name,
    parent = "root",
    channel,
    lower = -Inf,
    upper = Inf,
    coordinate_scale = c("raw", "transformed"),
    transform = c("asinh", "log1p", "none"),
    cofactor = 150,
    channel_input = NULL,
    channel_description = NULL,
    note = NULL,
    save_as_default = FALSE
) {

  coordinate_scale <- match.arg(coordinate_scale)
  transform <- match.arg(transform)

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (identical(gate_name, "root")) {
    stop(
      "`gate_name` cannot be 'root'; that name is reserved.",
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

  if (!is.character(channel) ||
      length(channel) != 1L ||
      is.na(channel) ||
      !nzchar(channel)) {
    stop(
      "`channel` must be one non-empty FCS channel name.",
      call. = FALSE
    )
  }

  if (!is.numeric(lower) ||
      length(lower) != 1L ||
      is.na(lower)) {
    stop(
      "`lower` must be one numeric value; -Inf is allowed.",
      call. = FALSE
    )
  }

  if (!is.numeric(upper) ||
      length(upper) != 1L ||
      is.na(upper)) {
    stop(
      "`upper` must be one numeric value; Inf is allowed.",
      call. = FALSE
    )
  }

  if (!is.finite(lower) && !identical(lower, -Inf)) {
    stop(
      "`lower` must be finite or -Inf.",
      call. = FALSE
    )
  }

  if (!is.finite(upper) && !identical(upper, Inf)) {
    stop(
      "`upper` must be finite or Inf.",
      call. = FALSE
    )
  }

  if (lower >= upper) {
    stop(
      "`lower` must be strictly less than `upper`.",
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

  if (!is.null(note) &&
      (!is.character(note) ||
       length(note) != 1L ||
       is.na(note) ||
       !nzchar(note))) {
    stop(
      "`note` must be NULL or one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.logical(save_as_default) ||
      length(save_as_default) != 1L ||
      is.na(save_as_default)) {
    stop(
      "`save_as_default` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  existing_gates <- self$get_gates(sample_id)

  if (!is.null(existing_gates[[gate_name]])) {
    stop(
      "A gate named '", gate_name,
      "' already exists for sample '", sample_id, "'.",
      call. = FALSE
    )
  }

  if (!identical(parent, "root") &&
      is.null(existing_gates[[parent]])) {
    stop(
      "Parent gate '", parent,
      "' does not exist for sample '", sample_id,
      "'.",
      call. = FALSE
    )
  }

  convert_to_raw <- function(value) {
    if (identical(coordinate_scale, "raw") ||
        identical(transform, "none")) {
      return(value)
    }

    if (identical(value, -Inf) ||
        identical(value, Inf)) {
      return(value)
    }

    switch(
      transform,
      asinh = cofactor * sinh(value),
      log1p = expm1(value),
      none = value
    )
  }

  raw_lower <- convert_to_raw(lower)
  raw_upper <- convert_to_raw(upper)

  if ((!is.finite(raw_lower) && !identical(raw_lower, -Inf)) ||
      (!is.finite(raw_upper) && !identical(raw_upper, Inf)) ||
      raw_lower >= raw_upper) {
    stop(
      "Range-bound conversion produced invalid raw FCS limits.",
      call. = FALSE
    )
  }

  if (is.null(channel_input)) {
    channel_input <- channel
  }

  if (is.null(channel_description)) {
    channel_description <- channel
  }

  new_gate <- list(
    type = "range",
    parent = parent,

    # `channels` follows the convention used by threshold gates.
    channels = channel,

    channel_input = channel_input,
    channel_description = channel_description,

    lower = raw_lower,
    upper = raw_upper,

    transform = transform,
    cofactor = cofactor,

    method = "manual_range",

    manual_creation = list(
      created_at = Sys.time(),
      coordinate_scale_supplied = coordinate_scale,
      transform_supplied = if (identical(
        coordinate_scale,
        "transformed"
      )) {
        transform
      } else {
        NULL
      },
      cofactor_supplied = if (identical(
        coordinate_scale,
        "transformed"
      ) && identical(transform, "asinh")) {
        cofactor
      } else {
        NULL
      },
      supplied_lower = lower,
      supplied_upper = upper,
      note = note
    )
  )

  self$update_gate(
    sample_ids = sample_id,
    gate_name = gate_name,
    gate_definition = new_gate,
    save_as_default = save_as_default
  )

  message(
    "Created range gate '", gate_name,
    "' under parent '", parent,
    "' for sample '", sample_id, "'."
  )

  invisible(self)
},


# -------------------------------------------------------------------
# Return the immediate parent gate name.
# Returns NULL when gate_name = "root".
# -------------------------------------------------------------------
get_parent_gate = function(
    sample_id,
    gate_name
) {

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (identical(gate_name, "root")) {
    return(NULL)
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

  gate$parent
},


# -------------------------------------------------------------------
# Return immediate child gates of a parent gate.
# Set recursive = TRUE to return all descendants.
# -------------------------------------------------------------------
get_child_gates = function(
    sample_id,
    parent_gate = "root",
    recursive = FALSE
) {

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.character(parent_gate) ||
      length(parent_gate) != 1L ||
      is.na(parent_gate) ||
      !nzchar(parent_gate)) {
    stop(
      "`parent_gate` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.logical(recursive) ||
      length(recursive) != 1L ||
      is.na(recursive)) {
    stop(
      "`recursive` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  gates <- self$get_gates(sample_id)

  if (!identical(parent_gate, "root") &&
      is.null(gates[[parent_gate]])) {
    stop(
      "No parent gate named '", parent_gate,
      "' was found for sample '", sample_id, "'.",
      call. = FALSE
    )
  }

  if (length(gates) == 0L) {
    return(character(0))
  }

  direct_children <- names(gates)[vapply(
    gates,
    function(gate) {
      identical(gate$parent, parent_gate)
    },
    logical(1L)
  )]

  direct_children <- sort(direct_children)

  if (!isTRUE(recursive) || length(direct_children) == 0L) {
    return(direct_children)
  }

  descendants <- direct_children

  for (child_gate in direct_children) {
    descendants <- c(
      descendants,
      self$get_child_gates(
        sample_id = sample_id,
        parent_gate = child_gate,
        recursive = TRUE
      )
    )
  }

  unique(descendants)
},


# -------------------------------------------------------------------
# Return a gate hierarchy as a data frame.
# Useful for an eventual FlowJo-like sidebar/tree widget.
# -------------------------------------------------------------------
get_gate_tree = function(
    sample_id,
    include_root = TRUE
) {

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.logical(include_root) ||
      length(include_root) != 1L ||
      is.na(include_root)) {
    stop(
      "`include_root` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  gates <- self$get_gates(sample_id)

  output <- list()
  visited <- character(0)

  add_branch <- function(
      current_gate,
      parent_gate,
      depth,
      path
  ) {

    if (current_gate %in% visited) {
      stop(
        "A cycle was detected in the gate hierarchy at gate '",
        current_gate,
        "'.",
        call. = FALSE
      )
    }

    visited <<- c(visited, current_gate)

    current_definition <- gates[[current_gate]]

    output[[length(output) + 1L]] <<- data.frame(
      gate_name = current_gate,
      parent = parent_gate,
      type = current_definition$type,
      depth = depth,
      path = path,
      stringsAsFactors = FALSE
    )

    children <- self$get_child_gates(
      sample_id = sample_id,
      parent_gate = current_gate,
      recursive = FALSE
    )

    for (child in children) {
      add_branch(
        current_gate = child,
        parent_gate = current_gate,
        depth = depth + 1L,
        path = paste(path, child, sep = "/")
      )
    }
  }

  root_children <- self$get_child_gates(
    sample_id = sample_id,
    parent_gate = "root",
    recursive = FALSE
  )

  if (isTRUE(include_root)) {
    output[[length(output) + 1L]] <- data.frame(
      gate_name = "root",
      parent = NA_character_,
      type = "root",
      depth = 0L,
      path = "root",
      stringsAsFactors = FALSE
    )
  }

  for (child in root_children) {
    add_branch(
      current_gate = child,
      parent_gate = "root",
      depth = 1L,
      path = paste("root", child, sep = "/")
    )
  }

  tree <- if (length(output) == 0L) {
    data.frame(
      gate_name = character(0),
      parent = character(0),
      type = character(0),
      depth = integer(0),
      path = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, output)
  }

  # Identify definitions whose parent is missing or which are disconnected
  # from root. The tree is still returned, but the issue is made explicit.
  orphaned_gates <- setdiff(names(gates), visited)

  if (length(orphaned_gates) > 0L) {
    warning(
      "The following gate(s) are disconnected from root or have a missing ",
      "parent and are absent from the returned tree: ",
      paste(orphaned_gates, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  rownames(tree) <- NULL
  tree
},


export_gate_tree <- function(
    project,
    sample_id,
    file
) {

  tree <- project$get_gate_tree(
    sample_id = sample_id,
    include_root = FALSE
  )

  gate_definitions <- stats::setNames(
    lapply(
      tree$gate_name,
      function(gate_name) {
        project$get_gates(
          sample_id = sample_id,
          gate_name = gate_name
        )
      }
    ),
    tree$gate_name
  )

  gate_export <- list(
    format = "cytoTerra_gate_tree",
    format_version = 1L,
    exported_at_utc = format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    ),
    source_sample = sample_id,
    tree = tree,
    gates = gate_definitions
  )

  saveRDS(gate_export, file = file)

  invisible(file)
},

export_all_gate_trees <- function(
    project,
    file
) {

  sample_ids <- project$get_sample_ids()

  samples <- stats::setNames(
    lapply(
      sample_ids,
      function(sample_id) {

        tree <- project$get_gate_tree(
          sample_id = sample_id,
          include_root = FALSE
        )

        gates <- stats::setNames(
          lapply(
            tree$gate_name,
            function(gate_name) {
              project$get_gates(
                sample_id = sample_id,
                gate_name = gate_name
              )
            }
          ),
          tree$gate_name
        )

        list(
          tree = tree,
          gates = gates
        )
      }
    ),
    sample_ids
  )

  gate_export <- list(
    format = "cytoTerra_all_gate_trees",
    format_version = 1L,
    exported_at_utc = format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    ),
    samples = samples
  )

  saveRDS(gate_export, file = file)

  invisible(file)
},

import_gate_tree = function(
    file,
    to_samples = NULL,
    overwrite = FALSE
) {

  if (!is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      !nzchar(file)) {
    stop(
      "`file` must be one non-empty file path.",
      call. = FALSE
    )
  }

  if (!file.exists(file)) {
    stop(
      "The gate-template file does not exist: ",
      file,
      call. = FALSE
    )
  }

  if (!is.logical(overwrite) ||
      length(overwrite) != 1L ||
      is.na(overwrite)) {
    stop(
      "`overwrite` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  # Deliberately conservative for the first implementation.
  # Deleting/replacing a gate can also affect its descendants, so that
  # behavior should be added separately and explicitly.
  if (isTRUE(overwrite)) {
    stop(
      paste(
        "`overwrite = TRUE` is not yet supported for gate-tree import.",
        "Use `overwrite = FALSE`, rename/delete conflicting target gates,",
        "then import again."
      ),
      call. = FALSE
    )
  }

  gate_export <- readRDS(file)

  if (!is.list(gate_export) ||
      !identical(gate_export$format, "cytoTerra_gate_tree") ||
      is.null(gate_export$gates) ||
      is.null(gate_export$tree)) {
    stop(
      paste(
        "This does not appear to be a valid cytoTerra single-sample",
        "gate-tree export."
      ),
      call. = FALSE
    )
  }

  tree <- gate_export$tree
  gates <- gate_export$gates

  required_tree_columns <- c(
    "gate_name",
    "depth",
    "type"
  )

  missing_tree_columns <- setdiff(
    required_tree_columns,
    colnames(tree)
  )

  if (length(missing_tree_columns) > 0L) {
    stop(
      "The exported gate tree is missing column(s): ",
      paste(missing_tree_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (nrow(tree) == 0L) {
    stop(
      "The exported template contains no gates.",
      call. = FALSE
    )
  }

  gate_names <- as.character(tree$gate_name)

  if (anyDuplicated(gate_names)) {
    stop(
      "The exported gate tree contains duplicate gate names.",
      call. = FALSE
    )
  }

  if (!all(gate_names %in% names(gates))) {
    missing_gates <- setdiff(gate_names, names(gates))

    stop(
      "The exported template is missing definition(s) for: ",
      paste(missing_gates, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Preserve parent-before-child ordering even if the input data frame
  # was reordered externally.
  tree <- tree[order(tree$depth), , drop = FALSE]
  gate_names <- as.character(tree$gate_name)

  if (is.null(to_samples)) {

    to_samples <- self$get_sample_ids()

    # A template normally originated from one current-project sample.
    # Do not import it back onto itself by default.
    if (!is.null(gate_export$source_sample) &&
        gate_export$source_sample %in% to_samples) {
      to_samples <- setdiff(
        to_samples,
        gate_export$source_sample
      )
    }
  }

  cytoTerra_validate_sample_ids(
    self,
    to_samples
  )

  if (length(to_samples) == 0L) {
    stop(
      "No target samples remain after excluding the template source sample.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # Preflight validation:
  # validate every target first, so an invalid target does not result
  # in a partially imported set of gates.
  # ------------------------------------------------------------------

  for (target_sample in to_samples) {

    flow_frame <- self$get_flow_frame(target_sample)

    target_channels <- colnames(
      flowCore::exprs(flow_frame)
    )

    for (gate_name in gate_names) {

      gate <- gates[[gate_name]]

      if (!is.list(gate) ||
          is.null(gate$type) ||
          is.null(gate$parent)) {
        stop(
          "Gate '", gate_name,
          "' has an invalid or incomplete saved definition.",
          call. = FALSE
        )
      }

      if (!gate$type %in% c("polygon", "range")) {
        stop(
          "Gate '", gate_name,
          "' has unsupported type '", gate$type,
          "'. Only polygon and range gates can currently be imported.",
          call. = FALSE
        )
      }

      # Existing same-name gates are blocked before import begins.
      if (!is.null(self$get_gates(target_sample, gate_name))) {
        stop(
          "Target sample '", target_sample,
          "' already contains a gate named '", gate_name,
          "'. Import was cancelled; no gates were added.",
          call. = FALSE
        )
      }

      required_channels <- if (identical(gate$type, "polygon")) {

        if (is.null(gate$channels) ||
            !all(c("x", "y") %in% names(gate$channels))) {
          stop(
            "Polygon gate '", gate_name,
            "' does not contain valid X/Y channel information.",
            call. = FALSE
          )
        }

        c(
          gate$channels[["x"]],
          gate$channels[["y"]]
        )
      } else {

        range_channel <- NULL

        if (!is.null(gate$channels) &&
            "channel" %in% names(gate$channels)) {
          range_channel <- gate$channels[["channel"]]
        }

        if (is.null(range_channel) &&
            !is.null(gate$channel)) {
          range_channel <- gate$channel
        }

        if (is.null(range_channel) ||
            !is.character(range_channel) ||
            length(range_channel) != 1L ||
            !nzchar(range_channel)) {
          stop(
            "Range gate '", gate_name,
            "' does not contain valid channel information.",
            call. = FALSE
          )
        }

        range_channel
      }

      missing_channels <- setdiff(
        required_channels,
        target_channels
      )

      if (length(missing_channels) > 0L) {
        stop(
          "Target sample '", target_sample,
          "' does not contain required channel(s) for gate '",
          gate_name,
          "': ",
          paste(missing_channels, collapse = ", "),
          ". Import was cancelled; no gates were added.",
          call. = FALSE
        )
      }
    }
  }

  # ------------------------------------------------------------------
  # Import gates. Since tree order is parent before child, each parent
  # exists before its child is added.
  # ------------------------------------------------------------------

  for (target_sample in to_samples) {

    for (gate_name in gate_names) {

      gate <- gates[[gate_name]]

      coordinate_scale <- if (!is.null(gate$coordinate_scale)) {
        gate$coordinate_scale
      } else {
        "raw"
      }

      gate_transform <- if (!is.null(gate$transform)) {
        gate$transform
      } else {
        "none"
      }

      gate_cofactor <- if (!is.null(gate$cofactor)) {
        gate$cofactor
      } else {
        150
      }

      import_note <- paste0(
        "Imported from gate-tree template",
        if (!is.null(gate_export$source_sample)) {
          paste0(
            " created from sample '",
            gate_export$source_sample,
            "'"
          )
        } else {
          ""
        },
        "."
      )

      if (identical(gate$type, "polygon")) {

        if (is.null(gate$coordinates) ||
            !all(c("x", "y") %in% names(gate$coordinates)) ||
            nrow(gate$coordinates) < 3L) {
          stop(
            "Polygon gate '", gate_name,
            "' has invalid saved coordinates.",
            call. = FALSE
          )
        }

        self$add_polygon_gate(
          sample_id = target_sample,
          gate_name = gate_name,
          parent = gate$parent,
          x_channel = gate$channels[["x"]],
          y_channel = gate$channels[["y"]],
          coordinates = gate$coordinates,
          coordinate_scale = coordinate_scale,
          transform = gate_transform,
          cofactor = gate_cofactor,
          note = import_note
        )
      }

      if (identical(gate$type, "range")) {

        range_channel <- if (!is.null(gate$channels) &&
                             "channel" %in% names(gate$channels)) {
          gate$channels[["channel"]]
        } else {
          gate$channel
        }

        if (is.null(gate$lower) ||
            is.null(gate$upper) ||
            !is.finite(gate$lower) ||
            !is.finite(gate$upper) ||
            gate$lower >= gate$upper) {
          stop(
            "Range gate '", gate_name,
            "' has invalid saved boundaries.",
            call. = FALSE
          )
        }

        self$add_range_gate(
          sample_id = target_sample,
          gate_name = gate_name,
          parent = gate$parent,
          channel = range_channel,
          lower = gate$lower,
          upper = gate$upper,
          coordinate_scale = coordinate_scale,
          transform = gate_transform,
          cofactor = gate_cofactor,
          note = import_note
        )
      }
    }
  }

  message(
    "Imported ",
    length(gate_names),
    " gate(s) into ",
    length(to_samples),
    " sample(s)."
  )

  invisible(self)
},
    
print_gate_tree = function(project, sample_id) {
  tree <- project$get_gate_tree(sample_id)

  cat(
    paste0(
      strrep("  ", tree$depth),
      ifelse(tree$depth == 0L, "", "└─ "),
      tree$gate_name,
      " [",
      tree$type,
      "]"
    ),
    sep = "\n"
  )

  invisible(tree)
},

# -------------------------------------------------------------------
# Delete one gate, optionally with every descendant gate.
#
# Default behavior protects against accidentally orphaning child gates.
# -------------------------------------------------------------------
delete_gate = function(
    sample_id,
    gate_name,
    include_descendants = FALSE,
    remove_default = TRUE
) {

  cytoTerra_validate_sample_ids(self, sample_id)

  if (!is.character(gate_name) ||
      length(gate_name) != 1L ||
      is.na(gate_name) ||
      !nzchar(gate_name)) {
    stop(
      "`gate_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (identical(gate_name, "root")) {
    stop(
      "The reserved root population cannot be deleted.",
      call. = FALSE
    )
  }

  if (!is.logical(include_descendants) ||
      length(include_descendants) != 1L ||
      is.na(include_descendants)) {
    stop(
      "`include_descendants` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(remove_default) ||
      length(remove_default) != 1L ||
      is.na(remove_default)) {
    stop(
      "`remove_default` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  sample <- self$samples[[sample_id]]
  gate <- sample$gates[[gate_name]]

  if (is.null(gate)) {
    stop(
      "No gate named '", gate_name,
      "' was found for sample '", sample_id, "'.",
      call. = FALSE
    )
  }

  child_gates <- self$get_child_gates(
    sample_id = sample_id,
    parent_gate = gate_name,
    recursive = FALSE
  )

  if (length(child_gates) > 0L &&
      !isTRUE(include_descendants)) {
    stop(
      "Cannot delete gate '", gate_name,
      "' because it has child gate(s): ",
      paste(child_gates, collapse = ", "),
      ". Set `include_descendants = TRUE` to delete the full branch.",
      call. = FALSE
    )
  }

  gates_to_remove <- if (isTRUE(include_descendants)) {
    c(
      gate_name,
      self$get_child_gates(
        sample_id = sample_id,
        parent_gate = gate_name,
        recursive = TRUE
      )
    )
  } else {
    gate_name
  }

  # Preserve deleted definitions in gate_history.
  for (current_gate_name in gates_to_remove) {
    old_definition <- sample$gates[[current_gate_name]]

    if (!is.null(old_definition)) {
      sample$gate_history[[current_gate_name]] <- c(
        sample$gate_history[[current_gate_name]],
        list(old_definition)
      )
    }

    sample$gates[[current_gate_name]] <- NULL

    if (isTRUE(remove_default)) {
      sample$default_gates[[current_gate_name]] <- NULL
    }
  }

  # Cached memberships are no longer valid after deleting any gate.
  sample$gate_membership <- list()

  self$samples[[sample_id]] <- sample
  self$metadata$updated_at <- Sys.time()

  message(
    "Deleted ",
    length(gates_to_remove),
    " gate",
    if (length(gates_to_remove) == 1L) "" else "s",
    " from sample '", sample_id,
    "': ",
    paste(gates_to_remove, collapse = ", "),
    "."
  )

  invisible(self)
}
    
)


