# Internal mutable cytometry project object.
#'
#' @noRd
CytoTerraProject <- R6::R6Class(
  classname = "CytoTerraProject",

  public = c(
    list(
      
      metadata = NULL,
      sample_meta = NULL, 
      samples = NULL,
      sample_ids = NULL,
      sample_names = NULL, 
      gating_trees = NULL,
      channel_labels = NULL,
      active = NULL,
      clustering_results = NULL,

      initialize = function(
          flow_set,
          fcs_files,
          sample_ids,
          sample_names, 
          sample_meta = NULL,
          channel_labels = NULL
      ) {
        if (!inherits(flow_set, "flowSet")) {
          stop("`flow_set` must be a flowCore::flowSet.", call. = FALSE)
        }

        if (length(sample_ids) != length(flow_set)) {
          stop(
            "`sample_ids` must have one ID for each flowFrame.",
            call. = FALSE
          )
        }

        if (length(sample_names) != length(flow_set)) {
          stop(
            "`sample_names` must have one name for each flowFrame.",
            call. = FALSE
          )
        }

        if (anyDuplicated(sample_ids)) {
          stop("`sample_ids` must be unique.", call. = FALSE)
        }

        if (length(fcs_files) != length(flow_set)) {
          stop(
            "`fcs_files` must have one path for each flowFrame.",
            call. = FALSE
          )
        }

        self$metadata <- list(
          class_version = "0.1.0",
          created_at = Sys.time(),
          updated_at = Sys.time()
        )

        if (!is.null(sample_meta)) {
          if (!is.data.frame(sample_meta)) {
            stop("`sample_meta` must be a data.frame.", call. = FALSE)
          }

          if (is.null(rownames(sample_meta)) ||
              !all(sample_ids %in% rownames(sample_meta))) {
            stop(
              "`sample_meta` row names must contain all `sample_ids`.",
              call. = FALSE
            )
          }

          self$sample_meta <- sample_meta[sample_ids, , drop = FALSE]
        } else {
          self$sample_meta <- NULL
        }

        self$samples <- list()

        for (i in seq_along(sample_ids)) {
          sample_id <- sample_ids[[i]]
          sample_name <- sample_names[[i]]
          fcs_path <- fcs_files[[i]]

          self$samples[[sample_id]] <- list(
            sample_id = sample_id,
            sample_name = sample_name,
            flow_frame = flow_set[[i]],
            fcs = list(
              sample_name = sample_name,
              file_name = basename(fcs_path),
              file_path = normalizePath(
                fcs_path,
                winslash = "/",
                mustWork = FALSE
              )
            ),
            gates = list(),
            default_gates = list(),
            gate_history = list(),
            gate_membership = list()
          )
        }

        self$sample_ids <- sample_ids
        self$sample_names <- stats::setNames(sample_names, sample_ids)
          
        self$gating_trees <- list()
        self$channel_labels <- channel_labels

        self$active <- list(
          sample_ids = sample_ids[[1]],
          gate_name = NULL,
          gating_tree_name = NULL,
          clustering_result_name = NULL
        )

        self$clustering_results <- list()

        invisible(self)
      },

      print = function(...) {
        cat("<CytoTerraProject>\n")
        cat("  Samples: ", length(self$samples), "\n", sep = "")
        cat(
          "  Sample IDs: ",
          paste(names(self$samples), collapse = ", "),
          "\n",
          sep = ""
        )
        cat(
          "  Gating trees: ",
          length(self$gating_trees),
          "\n",
          sep = ""
        )

        invisible(self)
      }
    ),

    # Methods defined in separate R files:
    cytoTerra_methods_samples, 
    cytoTerra_methods_gates,
    cytoTerra_methods_auto_gating,
    cytoTerra_methods_manual_gates,
    cytoTerra_methods_dataset,
    cytoTerra_methods_plot,
    cytoTerra_methods_clustering,
    cytoTerra_methods_gating_workspace
  )
)
