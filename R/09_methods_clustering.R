## some util function for plotting after clustering
cytoTerra_get_clustering_data = function(
    project,
    result_name = NULL
) {

  if (is.null(result_name)) {
    result_name <- project$active$clustering_result_name
  }

  if (is.null(result_name) ||
      !is.character(result_name) ||
      length(result_name) != 1L ||
      !nzchar(result_name)) {
    stop(
      paste0(
        "No active clustering result. Supply `result_name` or set ",
        "`project$active$clustering_result_name`."
      ),
      call. = FALSE
    )
  }

  result <- project$clustering_results[[result_name]]

  if (is.null(result)) {
    stop(
      "Clustering result '", result_name, "' does not exist.",
      call. = FALSE
    )
  }

  # Expected cytoTerra structure:
  # project$clustering_results[[result_name]] is a cell-level data.frame.
  if (is.data.frame(result)) {
    return(result)
  }

  # Also permit a list containing a data.frame in $data, if you use
  # that structure for some saved result types.
  if (is.list(result) && is.data.frame(result$data)) {
    return(result$data)
  }

  stop(
    "Clustering result '", result_name,
    "' must be a data.frame or a list containing a data.frame in `$data`.",
    call. = FALSE
  )
}

cytoTerra_add_sample_metadata = function(
    project,
    data,
    group_by
) {

  # A column already present in per-cell clustering data.
  if (group_by %in% colnames(data)) {
    return(data)
  }

  # Otherwise, look in project-level sample metadata.
  if (is.null(project$sample_meta)) {
    stop(
      "`group_by = '", group_by,
      "' is not in clustering data, and `project$sample_meta` is NULL.",
      call. = FALSE
    )
  }

  if (!"sample_id" %in% colnames(data)) {
    stop(
      "Clustering data must contain `sample_id` to use `sample_meta`.",
      call. = FALSE
    )
  }

  sample_meta <- as.data.frame(project$sample_meta)

  if (!group_by %in% colnames(sample_meta)) {
    stop(
      "`group_by = '", group_by,
      "' is not a column in clustering data or `project$sample_meta`.",
      call. = FALSE
    )
  }

  sample_meta$sample_id <- rownames(sample_meta)

  data[[group_by]] <- sample_meta[
    match(data$sample_id, sample_meta$sample_id),
    group_by
  ]

  data
}


cytoTerra_validate_markers = function(
    data,
    markers
) {

  if (is.null(markers) ||
      !is.character(markers) ||
      length(markers) == 0L) {
    stop(
      "`markers` must contain one or more marker names.",
      call. = FALSE
    )
  }

  missing_markers <- setdiff(markers, colnames(data))

  if (length(missing_markers) > 0L) {
    stop(
      "Marker(s) absent from clustering data: ",
      paste(missing_markers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

cytoTerra_downsample_by_metacluster = function(
    data,
    max_cells_per_cluster = 5000L,
    cluster_column = "FastPG_metacluster",
    seed = 1L
) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (!cluster_column %in% colnames(data)) {
    stop(
      "`cluster_column` is not present in `data`: ",
      cluster_column,
      ".",
      call. = FALSE
    )
  }

  if (is.null(max_cells_per_cluster)) {
    return(data)
  }

  if (!is.numeric(max_cells_per_cluster) ||
      length(max_cells_per_cluster) != 1L ||
      is.na(max_cells_per_cluster) ||
      !is.finite(max_cells_per_cluster) ||
      max_cells_per_cluster < 1) {
    stop(
      "`max_cells_per_cluster` must be NULL or a positive number.",
      call. = FALSE
    )
  }

  max_cells_per_cluster <- as.integer(
    max_cells_per_cluster
  )

  data <- data[
    !is.na(data[[cluster_column]]),
    ,
    drop = FALSE
  ]

  set.seed(seed)

  cluster_indices <- split(
    seq_len(nrow(data)),
    as.character(data[[cluster_column]])
  )

  selected_indices <- unlist(
    lapply(cluster_indices, function(indices) {

      if (length(indices) <= max_cells_per_cluster) {
        return(indices)
      }

      sample(
        indices,
        size = max_cells_per_cluster,
        replace = FALSE
      )
    }),
    use.names = FALSE
  )

  # Randomize plotting order so one cluster is not systematically
  # drawn on top of the others.
  selected_indices <- sample(selected_indices)

  data[
    selected_indices,
    ,
    drop = FALSE
  ]
}

## list 
cytoTerra_methods_clustering <- list(

  run_harmony_umap_flowsom = function(
      result_name,
      data = NULL,
      gate_names = NULL,
      sample_ids = NULL,
      markers,
      batch_col = "sample_id",
      scale = c("asinh", "raw", "log1p", "gate"),
      cofactor = 150,
      max_events = NULL,
      seed = 42,
      exclude_metaclusters = NULL,
      metacluster_col = "FlowSOM_metacluster",
      umap_neighbours = 30,
      umap_min_dist = 0.3,
      flowsom_xdim = 10,
      flowsom_ydim = 10,
      verbose = TRUE
  ) {

    scale <- match.arg(scale)

    # ---------------------------------------------------------------
    # Validate arguments
    # ---------------------------------------------------------------
    if (!is.character(result_name) ||
        length(result_name) != 1L ||
        is.na(result_name) ||
        !nzchar(result_name)) {
      stop(
        "`result_name` must be one non-empty character value.",
        call. = FALSE
      )
    }

    if (!is.character(markers) ||
        length(markers) < 1L ||
        anyNA(markers) ||
        any(!nzchar(markers))) {
      stop(
        "`markers` must be a non-empty character vector of marker names.",
        call. = FALSE
      )
    }

    if (anyDuplicated(markers)) {
      stop(
        "`markers` cannot contain duplicated marker names.",
        call. = FALSE
      )
    }

    if (!is.character(batch_col) ||
        length(batch_col) != 1L ||
        is.na(batch_col) ||
        !nzchar(batch_col)) {
      stop(
        "`batch_col` must be one non-empty column name.",
        call. = FALSE
      )
    }

    if (!is.character(metacluster_col) ||
        length(metacluster_col) != 1L ||
        is.na(metacluster_col) ||
        !nzchar(metacluster_col)) {
      stop(
        "`metacluster_col` must be one non-empty column name.",
        call. = FALSE
      )
    }

    if (!is.null(gate_names)) {
      if (!is.character(gate_names) ||
          length(gate_names) == 0L ||
          anyNA(gate_names) ||
          any(!nzchar(gate_names))) {
        stop(
          "`gate_names` must be NULL or a non-empty character vector.",
          call. = FALSE
        )
      }

      gate_names <- unique(gate_names)
    }

    if (!is.null(sample_ids)) {
      cytoTerra_validate_sample_ids(self, sample_ids)
    }

    if (!is.null(exclude_metaclusters)) {
      if (!is.atomic(exclude_metaclusters) ||
          length(exclude_metaclusters) == 0L ||
          anyNA(exclude_metaclusters)) {
        stop(
          "`exclude_metaclusters` must be NULL or a non-missing vector ",
          "of FlowSOM metacluster IDs.",
          call. = FALSE
        )
      }

      # With independent analyses, exclusions can only be applied to
      # data supplied explicitly by the caller.
      if (is.null(data)) {
        stop(
          "`exclude_metaclusters` requires explicitly supplied `data` ",
          "containing prior FlowSOM metacluster labels. ",
          "Each saved clustering run is independent.",
          call. = FALSE
        )
      }
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

    if (!is.null(max_events)) {
      if (!is.numeric(max_events) ||
          length(max_events) != 1L ||
          is.na(max_events) ||
          !is.finite(max_events) ||
          max_events < 1L) {
        stop(
          "`max_events` must be NULL or one positive number.",
          call. = FALSE
        )
      }

      max_events <- as.integer(max_events)
    }

    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed)) {
      stop(
        "`seed` must be one finite numeric value.",
        call. = FALSE
      )
    }

    if (!is.numeric(umap_neighbours) ||
        length(umap_neighbours) != 1L ||
        is.na(umap_neighbours) ||
        umap_neighbours < 2L) {
      stop(
        "`umap_neighbours` must be at least 2.",
        call. = FALSE
      )
    }

    if (!is.numeric(umap_min_dist) ||
        length(umap_min_dist) != 1L ||
        is.na(umap_min_dist) ||
        !is.finite(umap_min_dist) ||
        umap_min_dist < 0) {
      stop(
        "`umap_min_dist` must be one non-negative finite number.",
        call. = FALSE
      )
    }

    if (!is.numeric(flowsom_xdim) ||
        length(flowsom_xdim) != 1L ||
        is.na(flowsom_xdim) ||
        !is.finite(flowsom_xdim) ||
        flowsom_xdim < 1L) {
      stop(
        "`flowsom_xdim` must be a positive integer.",
        call. = FALSE
      )
    }

    if (!is.numeric(flowsom_ydim) ||
        length(flowsom_ydim) != 1L ||
        is.na(flowsom_ydim) ||
        !is.finite(flowsom_ydim) ||
        flowsom_ydim < 1L) {
      stop(
        "`flowsom_ydim` must be a positive integer.",
        call. = FALSE
      )
    }

    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose)) {
      stop(
        "`verbose` must be TRUE or FALSE.",
        call. = FALSE
      )
    }

    # Do not silently overwrite a distinct saved analysis.
    if (is.null(self$clustering_results)) {
      self$clustering_results <- list()
    }

    if (!is.null(self$clustering_results[[result_name]])) {
      stop(
        "A clustering result named '", result_name,
        "' already exists. Choose a new `result_name` to preserve ",
        "the existing independent analysis.",
        call. = FALSE
      )
    }

    # ---------------------------------------------------------------
    # Obtain input data
    #
    # Priority:
    # 1. Explicit `data`, if supplied.
    # 2. Fresh extraction from gate_names.
    # ---------------------------------------------------------------
    if (!is.null(data)) {

      dat <- data.table::as.data.table(
        data.table::copy(data)
      )

      if (isTRUE(verbose)) {
        message("Using explicitly supplied `data` for clustering.")
      }

    } else {

      if (is.null(gate_names)) {
        stop(
          "Supply `gate_names` when `data = NULL`.",
          call. = FALSE
        )
      }

      if (isTRUE(verbose)) {
        message(
          "Extracting events from gate",
          if (length(gate_names) == 1L) "" else "s",
          ": ",
          paste(shQuote(gate_names), collapse = ", "),
          "."
        )
      }

      dat <- self$extract_gated_events(
        gate_names = gate_names,
        sample_ids = sample_ids,
        scale = scale,
        cofactor = cofactor,
        max_events = max_events,
        seed = seed,
        include_sample_id = TRUE,
        include_event_id = TRUE,
        include_matched_gates = TRUE
      )

      dat <- data.table::as.data.table(
        data.table::copy(dat)
      )
    }

    if (nrow(dat) == 0L) {
      stop(
        "No events are available for clustering.",
        call. = FALSE
      )
    }

    # ---------------------------------------------------------------
    # Validate needed columns
    # ---------------------------------------------------------------
    missing_markers <- setdiff(markers, names(dat))

    if (length(missing_markers) > 0L) {
      stop(
        "The following requested marker column(s) are absent from the ",
        "clustering data: ",
        paste(missing_markers, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!batch_col %in% names(dat)) {
      stop(
        "The batch column '", batch_col,
        "' is absent from the clustering data.",
        call. = FALSE
      )
    }

    # ---------------------------------------------------------------
    # Optionally remove events assigned to old metaclusters.
    #
    # This only operates on explicit input data, not saved results.
    # ---------------------------------------------------------------
    excluded_event_count <- 0L

    if (!is.null(exclude_metaclusters)) {

      if (!metacluster_col %in% names(dat)) {
        stop(
          "Cannot exclude metaclusters because column '",
          metacluster_col,
          "' is not present in supplied `data`.",
          call. = FALSE
        )
      }

      exclude_metaclusters <- as.character(
        unique(exclude_metaclusters)
      )

      n_events_before_exclusion <- nrow(dat)

      dat <- dat[
        !as.character(get(metacluster_col)) %in%
          exclude_metaclusters
      ]

      excluded_event_count <- n_events_before_exclusion - nrow(dat)

      if (nrow(dat) == 0L) {
        stop(
          "No events remain after excluding FlowSOM metacluster(s): ",
          paste(exclude_metaclusters, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (isTRUE(verbose)) {
        message(
          "Excluded ", excluded_event_count, " event",
          if (excluded_event_count == 1L) "" else "s",
          " from metacluster(s): ",
          paste(exclude_metaclusters, collapse = ", "),
          "."
        )
      }
    }

    # ---------------------------------------------------------------
    # Validate batch information
    # ---------------------------------------------------------------
    if (anyNA(dat[[batch_col]])) {
      stop(
        "The batch column '", batch_col,
        "' contains missing values.",
        call. = FALSE
      )
    }

    dat[, (batch_col) := as.factor(get(batch_col))]

    n_batches <- data.table::uniqueN(dat[[batch_col]])

    if (n_batches < 1L) {
      stop(
        "No valid batches remain for clustering.",
        call. = FALSE
      )
    }

    # ---------------------------------------------------------------
    # Create and validate marker matrix
    # ---------------------------------------------------------------
    marker_matrix <- as.matrix(dat[, ..markers])
    storage.mode(marker_matrix) <- "double"

    if (any(!is.finite(marker_matrix))) {
      stop(
        "The marker matrix contains NA, Inf, or -Inf values. ",
        "Clean or filter the data before clustering.",
        call. = FALSE
      )
    }

    if (nrow(marker_matrix) < 3L) {
      stop(
        "At least three events are required for PCA, UMAP, and FlowSOM.",
        call. = FALSE
      )
    }

    marker_sds <- apply(marker_matrix, 2L, stats::sd)

    zero_variance_markers <- markers[
      !is.finite(marker_sds) | marker_sds == 0
    ]

    if (length(zero_variance_markers) > 0L) {
      stop(
        "The following marker(s) have zero variance and cannot be used ",
        "for scaled PCA: ",
        paste(zero_variance_markers, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      message(
        "Running PCA",
        if (n_batches > 1L) ", Harmony" else "",
        ", UMAP, and FlowSOM on ",
        nrow(dat), " events across ",
        n_batches, " batch",
        if (n_batches == 1L) "" else "es",
        "."
      )
    }

    # ---------------------------------------------------------------
    # 1. PCA
    # ---------------------------------------------------------------
    pca_rank <- min(
      30L,
      length(markers),
      nrow(marker_matrix) - 1L
    )

    pca_result <- stats::prcomp(
      marker_matrix,
      center = TRUE,
      scale. = TRUE,
      rank. = pca_rank
    )

    # ---------------------------------------------------------------
    # 2. Harmony batch correction
    # ---------------------------------------------------------------
    if (n_batches > 1L) {

      harmony_metadata <- data.frame(
        batch = as.factor(dat[[batch_col]])
      )

      names(harmony_metadata) <- batch_col

      harmony_pcs <- harmony::RunHarmony(
        data_mat = pca_result$x,
        meta_data = harmony_metadata,
        vars_use = batch_col,
        verbose = verbose
      )

    } else {

      warning(
        "Only one batch is present; Harmony was skipped and PCA ",
        "coordinates will be used for UMAP and FlowSOM.",
        call. = FALSE
      )

      harmony_pcs <- pca_result$x
    }

    harmony_colnames <- paste0(
      "Harmony_PC_",
      seq_len(ncol(harmony_pcs))
    )

    dat[, (harmony_colnames) := data.table::as.data.table(harmony_pcs)]

    # ---------------------------------------------------------------
    # 3. UMAP
    # ---------------------------------------------------------------
    effective_umap_neighbours <- min(
      as.integer(umap_neighbours),
      nrow(dat) - 1L
    )

    if (effective_umap_neighbours < 2L) {
      stop(
        "At least three events are required to run UMAP.",
        call. = FALSE
      )
    }

    set.seed(seed)

    dat_umap <- Spectre::run.umap(
      dat = dat,
      use.cols = harmony_colnames,
      umap.x.name = "UMAP_X",
      umap.y.name = "UMAP_Y",
      umap.seed = seed,
      neighbours = effective_umap_neighbours,
      min_dist = umap_min_dist,
      metric = "euclidean",
      verbose = verbose
    )

    # ---------------------------------------------------------------
    # 4. FlowSOM
    # ---------------------------------------------------------------
    set.seed(seed)

    dat_flowsom <- Spectre::run.flowsom(
      dat = dat_umap,
      use.cols = harmony_colnames,
      xdim = as.integer(flowsom_xdim),
      ydim = as.integer(flowsom_ydim)
    )

    # ---------------------------------------------------------------
    # Save the independent named clustering result
    # ---------------------------------------------------------------
    clustering_manifest <- list(
      result_name = result_name,
      gate_names = gate_names,

      sample_ids = if (is.null(sample_ids)) {
        sort(unique(as.character(dat_flowsom[[batch_col]])))
      } else {
        sample_ids
      },

      markers = markers,
      batch_col = batch_col,
      scale = scale,
      cofactor = cofactor,
      max_events = max_events,
      seed = seed,

      n_events_input = nrow(dat_flowsom) + excluded_event_count,
      n_events_excluded = excluded_event_count,
      excluded_metaclusters = exclude_metaclusters,
      prior_metacluster_column = if (is.null(exclude_metaclusters)) {
        NULL
      } else {
        metacluster_col
      },

      n_events_clustered = nrow(dat_flowsom),
      n_batches = n_batches,
      harmony_applied = n_batches > 1L,
      harmony_columns = harmony_colnames,

      umap_neighbours = effective_umap_neighbours,
      umap_min_dist = umap_min_dist,

      flowsom_xdim = as.integer(flowsom_xdim),
      flowsom_ydim = as.integer(flowsom_ydim),

      created_at = Sys.time()
    )

    self$clustering_results[[result_name]] <- list(
      name = result_name,
      method = "Harmony_PCA_UMAP_FlowSOM",

      # This contains the original event columns plus PCA/Harmony,
      # UMAP, FlowSOM cluster, and FlowSOM metacluster outputs.
      data = dat_flowsom,

      # Settings and provenance for the run.
      manifest = clustering_manifest,

      created_at = Sys.time(),
      updated_at = Sys.time()
    )

    self$active$clustering_result_name <- result_name
    self$metadata$updated_at <- Sys.time()

    if (isTRUE(verbose)) {
      message(
        "Clustering completed and saved as: ",
        "self$clustering_results[['", result_name, "']]"
      )
    }

    invisible(dat_flowsom)
  },

run_harmony_umap_flowsom_fastpg = function(
    result_name,
    data = NULL,
    gate_names = NULL,
    sample_ids = NULL,
    markers,
    batch_col = "sample_id",
    scale = c("asinh", "raw", "log1p", "gate"),
    cofactor = 150,
    max_events = NULL,
    seed = 42,

    # Optional: remove cells carrying these labels in supplied `data`.
    exclude_metaclusters = NULL,
    metacluster_col = "FlowSOM_metacluster",

    umap_neighbours = 30,
    umap_min_dist = 0.3,

    flowsom_xdim = 10,
    flowsom_ydim = 10,

    fastpg_k = 30,

    # Number of consolidated FastPG populations.
    # Set to NULL to skip consensus metaclustering.
    fastpg_metaclusters = 10,
    fastpg_metacluster_col = "FastPG_metacluster",

    num_threads = 8,
    verbose = TRUE
) {

  scale <- match.arg(scale)

  # ---------------------------------------------------------------
  # Validate arguments
  # ---------------------------------------------------------------
  if (!is.character(result_name) ||
      length(result_name) != 1L ||
      is.na(result_name) ||
      !nzchar(result_name)) {
    stop(
      "`result_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.character(markers) ||
      length(markers) < 1L ||
      anyNA(markers) ||
      any(!nzchar(markers))) {
    stop(
      "`markers` must be a non-empty character vector of marker names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(markers)) {
    stop(
      "`markers` cannot contain duplicated marker names.",
      call. = FALSE
    )
  }

  if (!is.character(batch_col) ||
      length(batch_col) != 1L ||
      is.na(batch_col) ||
      !nzchar(batch_col)) {
    stop(
      "`batch_col` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.character(metacluster_col) ||
      length(metacluster_col) != 1L ||
      is.na(metacluster_col) ||
      !nzchar(metacluster_col)) {
    stop(
      "`metacluster_col` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.character(fastpg_metacluster_col) ||
      length(fastpg_metacluster_col) != 1L ||
      is.na(fastpg_metacluster_col) ||
      !nzchar(fastpg_metacluster_col)) {
    stop(
      "`fastpg_metacluster_col` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!is.null(gate_names)) {
    if (!is.character(gate_names) ||
        length(gate_names) == 0L ||
        anyNA(gate_names) ||
        any(!nzchar(gate_names))) {
      stop(
        "`gate_names` must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }

    gate_names <- unique(gate_names)
  }

  if (!is.null(sample_ids)) {
    cytoTerra_validate_sample_ids(self, sample_ids)
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

  if (!is.null(max_events)) {
    if (!is.numeric(max_events) ||
        length(max_events) != 1L ||
        is.na(max_events) ||
        !is.finite(max_events) ||
        max_events < 1L) {
      stop(
        "`max_events` must be NULL or one positive finite number.",
        call. = FALSE
      )
    }

    max_events <- as.integer(max_events)
  }

  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed)) {
    stop(
      "`seed` must be one finite numeric value.",
      call. = FALSE
    )
  }

  if (!is.numeric(umap_neighbours) ||
      length(umap_neighbours) != 1L ||
      is.na(umap_neighbours) ||
      !is.finite(umap_neighbours) ||
      umap_neighbours < 2L) {
    stop(
      "`umap_neighbours` must be a finite number of at least 2.",
      call. = FALSE
    )
  }

  if (!is.numeric(umap_min_dist) ||
      length(umap_min_dist) != 1L ||
      is.na(umap_min_dist) ||
      !is.finite(umap_min_dist) ||
      umap_min_dist < 0) {
    stop(
      "`umap_min_dist` must be one non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(flowsom_xdim) ||
      length(flowsom_xdim) != 1L ||
      is.na(flowsom_xdim) ||
      !is.finite(flowsom_xdim) ||
      flowsom_xdim < 1L) {
    stop(
      "`flowsom_xdim` must be a positive integer.",
      call. = FALSE
    )
  }

  if (!is.numeric(flowsom_ydim) ||
      length(flowsom_ydim) != 1L ||
      is.na(flowsom_ydim) ||
      !is.finite(flowsom_ydim) ||
      flowsom_ydim < 1L) {
    stop(
      "`flowsom_ydim` must be a positive integer.",
      call. = FALSE
    )
  }

  if (!is.numeric(fastpg_k) ||
      length(fastpg_k) != 1L ||
      is.na(fastpg_k) ||
      !is.finite(fastpg_k) ||
      fastpg_k < 2L) {
    stop(
      "`fastpg_k` must be a finite number of at least 2.",
      call. = FALSE
    )
  }

  if (!is.null(fastpg_metaclusters)) {
    if (!is.numeric(fastpg_metaclusters) ||
        length(fastpg_metaclusters) != 1L ||
        is.na(fastpg_metaclusters) ||
        !is.finite(fastpg_metaclusters) ||
        fastpg_metaclusters < 2L ||
        fastpg_metaclusters != as.integer(fastpg_metaclusters)) {
      stop(
        "`fastpg_metaclusters` must be NULL or one integer of at least 2.",
        call. = FALSE
      )
    }

    fastpg_metaclusters <- as.integer(fastpg_metaclusters)
  }

  if (!is.numeric(num_threads) ||
      length(num_threads) != 1L ||
      is.na(num_threads) ||
      !is.finite(num_threads) ||
      num_threads < 1L) {
    stop(
      "`num_threads` must be a positive integer.",
      call. = FALSE
    )
  }

  if (!is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)) {
    stop(
      "`verbose` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.null(exclude_metaclusters)) {
    if (!is.atomic(exclude_metaclusters) ||
        length(exclude_metaclusters) == 0L ||
        anyNA(exclude_metaclusters)) {
      stop(
        "`exclude_metaclusters` must be NULL or a non-missing vector.",
        call. = FALSE
      )
    }

    if (is.null(data)) {
      stop(
        "`exclude_metaclusters` requires explicitly supplied `data` ",
        "that contains previous cluster labels.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------
  # Initialize result storage and prevent accidental overwrite
  # ---------------------------------------------------------------
  if (is.null(self$clustering_results)) {
    self$clustering_results <- list()
  }

  if (!is.null(self$clustering_results[[result_name]])) {
    stop(
      "A clustering result named '", result_name,
      "' already exists. Use a new `result_name`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Obtain data
  # ---------------------------------------------------------------
  if (!is.null(data)) {

    dat <- data.table::as.data.table(
      data.table::copy(data)
    )

    if (isTRUE(verbose)) {
      message("Using explicitly supplied `data` for clustering.")
    }

  } else {

    if (is.null(gate_names)) {
      stop(
        "Supply `gate_names` when `data = NULL`.",
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      message(
        "Extracting events from gate",
        if (length(gate_names) == 1L) "" else "s",
        ": ",
        paste(shQuote(gate_names), collapse = ", "),
        "."
      )
    }

    dat <- self$extract_gated_events(
      gate_names = gate_names,
      sample_ids = sample_ids,
      scale = scale,
      cofactor = cofactor,
      max_events = max_events,
      seed = seed,
      include_sample_id = TRUE,
      include_event_id = TRUE,
      include_matched_gates = TRUE
    )

    dat <- data.table::as.data.table(
      data.table::copy(dat)
    )
  }

  if (nrow(dat) == 0L) {
    stop(
      "No events are available for clustering. If multiple sibling ",
      "gate names were supplied, they may have been treated as an ",
      "intersection. Cluster their shared parent gate or provide ",
      "an explicit union through `data`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Validate columns
  # ---------------------------------------------------------------
  missing_markers <- setdiff(markers, names(dat))

  if (length(missing_markers) > 0L) {
    stop(
      "Requested marker column(s) are absent from the data: ",
      paste(missing_markers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!batch_col %in% names(dat)) {
    stop(
      "Batch column '", batch_col,
      "' is absent from the clustering data.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # Optionally remove cells with existing metacluster labels
  # ---------------------------------------------------------------
  excluded_event_count <- 0L

  if (!is.null(exclude_metaclusters)) {

    if (!metacluster_col %in% names(dat)) {
      stop(
        "Cannot exclude metaclusters: column '",
        metacluster_col,
        "' is absent from supplied `data`.",
        call. = FALSE
      )
    }

    exclude_metaclusters <- as.character(
      unique(exclude_metaclusters)
    )

    n_before_exclusion <- nrow(dat)

    dat <- dat[
      !as.character(get(metacluster_col)) %in%
        exclude_metaclusters
    ]

    excluded_event_count <- n_before_exclusion - nrow(dat)

    if (nrow(dat) == 0L) {
      stop(
        "No events remain after excluding metacluster(s): ",
        paste(exclude_metaclusters, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------
  # Validate batch and marker matrices
  # ---------------------------------------------------------------
  if (anyNA(dat[[batch_col]])) {
    stop(
      "Batch column '", batch_col,
      "' contains missing values.",
      call. = FALSE
    )
  }

  dat[, (batch_col) := as.factor(get(batch_col))]

  n_batches <- data.table::uniqueN(dat[[batch_col]])

  marker_matrix <- as.matrix(dat[, ..markers])
  storage.mode(marker_matrix) <- "double"

  if (any(!is.finite(marker_matrix))) {
    stop(
      "Marker data contain NA, Inf, or -Inf values.",
      call. = FALSE
    )
  }

  n_events <- nrow(marker_matrix)

  if (n_events < 3L) {
    stop(
      "At least three events are required for clustering.",
      call. = FALSE
    )
  }

  marker_sds <- apply(marker_matrix, 2L, stats::sd)

  zero_variance_markers <- markers[
    !is.finite(marker_sds) | marker_sds == 0
  ]

  if (length(zero_variance_markers) > 0L) {
    stop(
      "Markers with zero variance cannot be used: ",
      paste(zero_variance_markers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------
  # 1. PCA
  # ---------------------------------------------------------------
  pca_rank <- min(
    30L,
    length(markers),
    n_events - 1L
  )

  if (isTRUE(verbose)) {
    message(
      "Running PCA",
      if (n_batches > 1L) ", Harmony" else "",
      ", FlowSOM, FastPG, FastPG metaclustering, and UMAP on ",
      n_events, " events across ", n_batches, " batch",
      if (n_batches == 1L) "" else "es",
      "."
    )
  }

  set.seed(seed)

  pca_result <- stats::prcomp(
    marker_matrix,
    center = TRUE,
    scale. = TRUE,
    rank. = pca_rank
  )

  # ---------------------------------------------------------------
  # 2. Harmony correction
  # ---------------------------------------------------------------
  if (n_batches > 1L) {

    harmony_metadata <- data.frame(
      batch = as.factor(dat[[batch_col]])
    )

    names(harmony_metadata) <- batch_col

    set.seed(seed)

    harmony_pcs <- harmony::RunHarmony(
      data_mat = pca_result$x,
      meta_data = harmony_metadata,
      vars_use = batch_col,
      verbose = verbose
    )

  } else {

    warning(
      "Only one batch is present; Harmony was skipped and PCA ",
      "coordinates will be used.",
      call. = FALSE
    )

    harmony_pcs <- pca_result$x
  }

  harmony_pcs <- as.matrix(harmony_pcs)
  storage.mode(harmony_pcs) <- "double"

  if (nrow(harmony_pcs) != n_events ||
      ncol(harmony_pcs) < 1L ||
      any(!is.finite(harmony_pcs))) {
    stop(
      "Harmony/PCA produced invalid coordinates.",
      call. = FALSE
    )
  }

  harmony_colnames <- paste0(
    "Harmony_PC_",
    seq_len(ncol(harmony_pcs))
  )

  colnames(harmony_pcs) <- harmony_colnames

  dat[, (harmony_colnames) := data.table::as.data.table(harmony_pcs)]

  # ---------------------------------------------------------------
  # Determine safe event-count-dependent settings
  # ---------------------------------------------------------------
  effective_umap_neighbours <- min(
    as.integer(umap_neighbours),
    n_events - 1L
  )

  effective_fastpg_k <- min(
    as.integer(fastpg_k),
    n_events - 1L
  )

  if (effective_umap_neighbours < 2L ||
      effective_fastpg_k < 2L) {
    stop(
      "At least three events are required for FastPG and UMAP.",
      call. = FALSE
    )
  }

  effective_flowsom_xdim <- min(
    as.integer(flowsom_xdim),
    n_events
  )

  effective_flowsom_ydim <- min(
    as.integer(flowsom_ydim),
    max(1L, floor(n_events / effective_flowsom_xdim))
  )

  effective_flowsom_nodes <- effective_flowsom_xdim *
    effective_flowsom_ydim

  # ---------------------------------------------------------------
  # 3. FlowSOM clustering of Harmony/PCA coordinates
  # ---------------------------------------------------------------
  if (isTRUE(verbose)) {
    message("Running FlowSOM.")
  }

  set.seed(seed)

  harmony_flow_frame <- flowCore::flowFrame(
    exprs = harmony_pcs
  )

  flowsom_result <- FlowSOM::FlowSOM(
    input = harmony_flow_frame,
    colsToUse = seq_len(ncol(harmony_pcs)),
    xdim = effective_flowsom_xdim,
    ydim = effective_flowsom_ydim
  )

  flowsom_mapping <- flowsom_result$map$mapping

  if (is.null(flowsom_mapping) ||
      nrow(flowsom_mapping) != n_events) {
    stop(
      "FlowSOM did not return one assignment per event.",
      call. = FALSE
    )
  }

  dat[, FlowSOM_cluster := factor(flowsom_mapping[, 1L])]

  if (ncol(flowsom_mapping) >= 2L) {
    dat[, FlowSOM_metacluster := factor(flowsom_mapping[, 2L])]
  }

  # ---------------------------------------------------------------
  # 4. FastPG clustering of Harmony/PCA coordinates
  # ---------------------------------------------------------------
  if (isTRUE(verbose)) {
    message(
      "Running FastPG with k = ", effective_fastpg_k,
      " using ", as.integer(num_threads), " thread",
      if (as.integer(num_threads) == 1L) "" else "s",
      "."
    )
  }

  set.seed(seed)

  fastpg_result <- FastPG::fastCluster(
    data = harmony_pcs,
    k = effective_fastpg_k,
    num_threads = as.integer(num_threads)
  )

  if (is.null(fastpg_result$communities) ||
      length(fastpg_result$communities) != n_events ||
      anyNA(fastpg_result$communities)) {
    stop(
      "FastPG did not return a valid cluster assignment for each event.",
      call. = FALSE
    )
  }

  dat[, FastPG_cluster := factor(fastpg_result$communities)]

  # ---------------------------------------------------------------
  # 5. Consensus metaclustering of FastPG clusters
  #
  # Calculate average original marker expression for every FastPG
  # cluster, then use FlowSOM consensus metaclustering to assign
  # each fine-grained FastPG cluster to a larger population.
  # ---------------------------------------------------------------
  fastpg_metacluster_summary <- NULL
  effective_fastpg_metaclusters <- NULL

  if (!is.null(fastpg_metaclusters)) {

    if (isTRUE(verbose)) {
      message(
        "Running FlowSOM consensus metaclustering of FastPG clusters."
      )
    }

    fastpg_metacluster_summary <- dat[
      ,
      lapply(.SD, mean, na.rm = TRUE),
      by = FastPG_cluster,
      .SDcols = markers
    ]

    fastpg_metacluster_summary[, cluster_id := as.character(FastPG_cluster)]

    n_fastpg_clusters <- nrow(fastpg_metacluster_summary)

    if (n_fastpg_clusters == 1L) {

      effective_fastpg_metaclusters <- 1L

      fastpg_metacluster_summary[
        ,
        (fastpg_metacluster_col) := "Meta_1"
      ]

      dat[
        ,
        (fastpg_metacluster_col) := factor("Meta_1")
      ]

      warning(
        "FastPG produced only one cluster; assigning all events to Meta_1.",
        call. = FALSE
      )

    } else {

      effective_fastpg_metaclusters <- min(
        fastpg_metaclusters,
        n_fastpg_clusters
      )

      fastpg_marker_matrix <- as.matrix(
        fastpg_metacluster_summary[, ..markers]
      )

      storage.mode(fastpg_marker_matrix) <- "double"

      if (any(!is.finite(fastpg_marker_matrix))) {
        stop(
          "FastPG cluster mean marker values contain NA, Inf, or -Inf.",
          call. = FALSE
        )
      }

      rownames(fastpg_marker_matrix) <- paste0(
        "FastPG_cluster_",
        fastpg_metacluster_summary$cluster_id
      )

      set.seed(seed)

      consensus_ids <- FlowSOM::metaClustering_consensus(
        data = fastpg_marker_matrix,
        k = effective_fastpg_metaclusters,
        seed = as.integer(seed)
      )

      if (length(consensus_ids) != n_fastpg_clusters ||
          anyNA(consensus_ids)) {
        stop(
          "Consensus metaclustering did not return one label per ",
          "FastPG cluster.",
          call. = FALSE
        )
      }

      fastpg_metacluster_summary[
        ,
        (fastpg_metacluster_col) := paste0(
          "Meta_",
          consensus_ids
        )
      ]

      cluster_to_metacluster <- setNames(
        fastpg_metacluster_summary[[fastpg_metacluster_col]],
        fastpg_metacluster_summary$cluster_id
      )

      dat[
        ,
        (fastpg_metacluster_col) := unname(
          cluster_to_metacluster[
            as.character(FastPG_cluster)
          ]
        )
      ]

      meta_levels <- paste0(
        "Meta_",
        sort(unique(as.integer(consensus_ids)))
      )

      dat[
        ,
        (fastpg_metacluster_col) := factor(
          get(fastpg_metacluster_col),
          levels = meta_levels
        )
      ]
    }

    if (anyNA(dat[[fastpg_metacluster_col]])) {
      stop(
        "Failed to map one or more FastPG metacluster labels back ",
        "to event-level data.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------
  # 6. UMAP on Harmony/PCA coordinates
  # ---------------------------------------------------------------
  if (isTRUE(verbose)) {
    message("Running UMAP.")
  }

  set.seed(seed)

  umap_result <- uwot::umap(
    X = harmony_pcs,
    n_neighbors = effective_umap_neighbours,
    min_dist = umap_min_dist,
    metric = "euclidean",
    n_threads = as.integer(num_threads),
    verbose = verbose,
    ret_model = FALSE
  )

  if (!is.matrix(umap_result) ||
      nrow(umap_result) != n_events ||
      ncol(umap_result) < 2L ||
      any(!is.finite(umap_result[, 1:2, drop = FALSE]))) {
    stop(
      "UMAP did not return valid two-dimensional coordinates.",
      call. = FALSE
    )
  }

  dat[, UMAP_X := umap_result[, 1L]]
  dat[, UMAP_Y := umap_result[, 2L]]

  # ---------------------------------------------------------------
  # Save result and provenance
  # ---------------------------------------------------------------
  clustering_manifest <- list(
    result_name = result_name,
    method = "Harmony_PCA_FlowSOM_FastPG_Consensus_UMAP",

    gate_names = gate_names,

    sample_ids = if (is.null(sample_ids)) {
      sort(unique(as.character(dat[[batch_col]])))
    } else {
      sample_ids
    },

    markers = markers,
    batch_col = batch_col,

    scale = scale,
    cofactor = cofactor,
    max_events = max_events,
    seed = seed,

    n_events_input = n_events + excluded_event_count,
    n_events_excluded = excluded_event_count,
    excluded_metaclusters = exclude_metaclusters,
    prior_metacluster_column = if (is.null(exclude_metaclusters)) {
      NULL
    } else {
      metacluster_col
    },

    n_events_clustered = n_events,
    n_batches = n_batches,

    pca_rank = pca_rank,
    harmony_applied = n_batches > 1L,
    harmony_columns = harmony_colnames,

    flowsom_xdim_requested = as.integer(flowsom_xdim),
    flowsom_ydim_requested = as.integer(flowsom_ydim),
    flowsom_xdim_used = effective_flowsom_xdim,
    flowsom_ydim_used = effective_flowsom_ydim,
    flowsom_n_nodes = effective_flowsom_nodes,

    n_flowsom_clusters = data.table::uniqueN(dat$FlowSOM_cluster),

    fastpg_k_requested = as.integer(fastpg_k),
    fastpg_k_used = effective_fastpg_k,
    n_fastpg_clusters = data.table::uniqueN(dat$FastPG_cluster),

    fastpg_metaclusters_requested = fastpg_metaclusters,
    fastpg_metaclusters_used = effective_fastpg_metaclusters,
    fastpg_metacluster_column = if (is.null(fastpg_metaclusters)) {
      NULL
    } else {
      fastpg_metacluster_col
    },

    n_fastpg_metaclusters = if (is.null(fastpg_metaclusters)) {
      NULL
    } else {
      data.table::uniqueN(dat[[fastpg_metacluster_col]])
    },

    umap_neighbours_requested = as.integer(umap_neighbours),
    umap_neighbours_used = effective_umap_neighbours,
    umap_min_dist = umap_min_dist,

    num_threads = as.integer(num_threads),
    created_at = Sys.time()
  )

  self$clustering_results[[result_name]] <- list(
    name = result_name,
    method = "Harmony_PCA_FlowSOM_FastPG_Consensus_UMAP",

    # Cell-level data: original markers, Harmony PCs, UMAP positions,
    # FlowSOM_cluster, FlowSOM_metacluster (if returned),
    # FastPG_cluster, and FastPG_metacluster.
    data = dat,

    # Full model/results objects.
    flowsom_output = flowsom_result,
    fastpg_output = fastpg_result,

    # One row per FastPG cluster, containing average marker expression
    # and its consensus-derived FastPG metacluster label.
    fastpg_metacluster_summary = fastpg_metacluster_summary,

    manifest = clustering_manifest,

    created_at = Sys.time(),
    updated_at = Sys.time()
  )

  self$active$clustering_result_name <- result_name
  self$metadata$updated_at <- Sys.time()

  if (isTRUE(verbose)) {
    message(
      "Clustering completed and saved as: ",
      "self$clustering_results[['", result_name, "']]"
    )
  }

  invisible(dat)
},


set_active_clustering_result = function(
    result_name
) {

  if (!is.character(result_name) ||
      length(result_name) != 1L ||
      !nzchar(result_name)) {
    stop(
      "`result_name` must be one non-empty character string.",
      call. = FALSE
    )
  }

  if (is.null(self$clustering_results[[result_name]])) {
    stop(
      "Clustering result '", result_name, "' does not exist.",
      call. = FALSE
    )
  }

  self$active$clustering_result_name <- result_name

  invisible(self)
},

plot_umap = function(
    result_name = NULL,
    group_by = "FastPG_metacluster",
    max_cells_per_cluster = 5000L,
    point_size = 0.25,
    alpha = 0.70,
    shuffle = TRUE,
    seed = 1L
) {

  data <- cytoTerra_get_clustering_data(
    project = self,
    result_name = result_name
  )

  required_columns <- c("UMAP_X", "UMAP_Y")

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Clustering data is missing: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  data <- cytoTerra_downsample_by_metacluster(
    data = data,
    max_cells_per_cluster = max_cells_per_cluster,
    cluster_column = "FastPG_metacluster",
    seed = seed
  )

  data <- cytoTerra_add_sample_metadata(
    project = self,
    data = data,
    group_by = group_by
  )

  if (shuffle) {
    set.seed(seed)
    data <- data[sample.int(nrow(data)), , drop = FALSE]
  }

  data$plot_group <- as.factor(data[[group_by]])

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$UMAP_X,
      y = .data$UMAP_Y,
      colour = .data$plot_group
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = alpha
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(
        "UMAP — ",
        if (is.null(result_name)) {
          self$active$clustering_result_name
        } else {
          result_name
        }
      ),
      x = "UMAP 1",
      y = "UMAP 2",
      colour = group_by
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
},

plot_feature_expression = function(
    markers,
    result_name = NULL,
    max_cells_per_cluster = 5000L,
    point_size = 0.25,
    alpha = 0.70,
    shuffle = TRUE,
    seed = 1L,
    min_cutoff = NULL,
    max_cutoff = NULL
) {

  data <- cytoTerra_get_clustering_data(
    project = self,
    result_name = result_name
  )

  cytoTerra_validate_markers(
    data = data,
    markers = markers
  )

  required_columns <- c("UMAP_X", "UMAP_Y")

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Clustering data is missing: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  data <- cytoTerra_downsample_by_metacluster(
    data = data,
    max_cells_per_cluster = max_cells_per_cluster,
    cluster_column = "FastPG_metacluster",
    seed = seed
  )

  if (shuffle) {
    set.seed(seed)
    data <- data[sample.int(nrow(data)), , drop = FALSE]
  }

  plots <- lapply(markers, function(marker) {

    plot_data <- data.frame(
      UMAP_X = data$UMAP_X,
      UMAP_Y = data$UMAP_Y,
      expression = data[[marker]]
    )

    if (!is.null(min_cutoff)) {
      plot_data$expression <- pmax(
        plot_data$expression,
        min_cutoff
      )
    }

    if (!is.null(max_cutoff)) {
      plot_data$expression <- pmin(
        plot_data$expression,
        max_cutoff
      )
    }

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$UMAP_X,
        y = .data$UMAP_Y,
        colour = .data$expression
      )
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha
      ) +
      ggplot2::scale_colour_viridis_c(
        option = "magma",
        name = marker
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = paste0(marker, " expression"),
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.line = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  })

  names(plots) <- markers

  # Return one ggplot for one marker and a named list for multiple markers.
  if (length(plots) == 1L) {
    return(plots[[1]])
  }

  plots
},

plot_dotplot = function(
    markers,
    result_name = NULL,
    group_by = "FastPG_metacluster",
    expression_threshold = 0,
    scale = TRUE,
    cluster_rows = FALSE
) {

  data <- cytoTerra_get_clustering_data(
    project = self,
    result_name = result_name
  )

  cytoTerra_validate_markers(
    data = data,
    markers = markers
  )

  data <- cytoTerra_add_sample_metadata(
    project = self,
    data = data,
    group_by = group_by
  )

  data$plot_group <- as.character(data[[group_by]])

  if (anyNA(data$plot_group)) {
    warning(
      "Some cells have missing values for `group_by`; they were omitted.",
      call. = FALSE
    )

    data <- data[
      !is.na(data$plot_group),
      ,
      drop = FALSE
    ]
  }

  long_data <- do.call(
    rbind,
    lapply(markers, function(marker) {
      data.frame(
        group = data$plot_group,
        marker = marker,
        expression = data[[marker]],
        stringsAsFactors = FALSE
      )
    })
  )

  summary_data <- stats::aggregate(
    expression ~ group + marker,
    data = long_data,
    FUN = mean,
    na.rm = TRUE
  )

  names(summary_data)[
    names(summary_data) == "expression"
  ] <- "average_expression"

  percent_data <- stats::aggregate(
    I(expression > expression_threshold) ~ group + marker,
    data = long_data,
    FUN = function(x) {
      100 * mean(x, na.rm = TRUE)
    }
  )

  names(percent_data)[
    names(percent_data) == "I(expression > expression_threshold)"
  ] <- "percent_expressing"

  summary_data <- merge(
    summary_data,
    percent_data,
    by = c("group", "marker"),
    sort = FALSE
  )

  if (scale) {

    summary_data$scaled_expression <- ave(
      summary_data$average_expression,
      summary_data$marker,
      FUN = function(x) {
        if (length(unique(x)) <= 1L) {
          rep(0, length(x))
        } else {
          as.numeric(scale(x))
        }
      }
    )

    colour_column <- "scaled_expression"
    colour_label <- "Scaled average expression"
  } else {
    colour_column <- "average_expression"
    colour_label <- "Average expression"
  }

  group_levels <- unique(summary_data$group)

  if (cluster_rows && length(group_levels) > 2L) {

    average_matrix <- xtabs(
      average_expression ~ group + marker,
      data = summary_data
    )

    group_levels <- rownames(
      average_matrix
    )[stats::hclust(
      stats::dist(average_matrix)
    )$order]
  }

  summary_data$group <- factor(
    summary_data$group,
    levels = rev(group_levels)
  )

  summary_data$marker <- factor(
    summary_data$marker,
    levels = markers
  )

  ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x = .data$marker,
      y = .data$group,
      size = .data$percent_expressing,
      colour = .data[[colour_column]]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size(
      range = c(0, 9),
      limits = c(0, 100),
      name = "% expressing"
    ) +
    ggplot2::scale_colour_viridis_c(
      option = "magma",
      name = colour_label
    ) +
    ggplot2::labs(
      title = paste0("Marker expression by ", group_by),
      x = NULL,
      y = group_by
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
},

plot_expression_heatmap = function(
    markers,
    result_name = NULL,
    group_by = "FastPG_metacluster",
    scale = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE
) {

  data <- cytoTerra_get_clustering_data(
    project = self,
    result_name = result_name
  )

  cytoTerra_validate_markers(
    data = data,
    markers = markers
  )

  data <- cytoTerra_add_sample_metadata(
    project = self,
    data = data,
    group_by = group_by
  )

  data$plot_group <- as.character(data[[group_by]])

  data <- data[
    !is.na(data$plot_group),
    ,
    drop = FALSE
  ]

  average_expression <- stats::aggregate(
    data[, markers, drop = FALSE],
    by = list(group = data$plot_group),
    FUN = mean,
    na.rm = TRUE
  )

  expression_matrix <- as.matrix(
    average_expression[, markers, drop = FALSE]
  )

  rownames(expression_matrix) <- average_expression$group

  if (scale) {
    expression_matrix <- scale(expression_matrix)
  }

  if (cluster_rows && nrow(expression_matrix) > 2L) {
    row_order <- stats::hclust(
      stats::dist(expression_matrix)
    )$order

    expression_matrix <- expression_matrix[
      row_order,
      ,
      drop = FALSE
    ]
  }

  if (cluster_columns && ncol(expression_matrix) > 2L) {
    column_order <- stats::hclust(
      stats::dist(t(expression_matrix))
    )$order

    expression_matrix <- expression_matrix[
      ,
      column_order,
      drop = FALSE
    ]
  }

  heatmap_data <- as.data.frame(
    as.table(expression_matrix),
    stringsAsFactors = FALSE
  )

  colnames(heatmap_data) <- c(
    "group",
    "marker",
    "expression"
  )

  heatmap_data$group <- factor(
    heatmap_data$group,
    levels = rev(rownames(expression_matrix))
  )

  heatmap_data$marker <- factor(
    heatmap_data$marker,
    levels = colnames(expression_matrix)
  )

  ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(
      x = .data$marker,
      y = .data$group,
      fill = .data$expression
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = if (scale) {
        "Scaled mean\nexpression"
      } else {
        "Mean\nexpression"
      }
    ) +
    ggplot2::labs(
      title = paste0("Marker expression by ", group_by),
      x = NULL,
      y = group_by
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid = ggplot2::element_blank()
    )
},

get_metacluster_counts = function(
    result_name = NULL
) {

  data <- cytoTerra_get_clustering_data(
    project = self,
    result_name = result_name
  )

  cytoTerra_metacluster_counts(data)
}

)

## end of list

