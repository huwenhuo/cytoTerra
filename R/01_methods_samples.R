# Methods related to project samples.
#'
#' @noRd
cytoTerra_methods_samples <- list(

  get_sample_ids = function() {
    names(self$samples)
  },


  get_sample_names = function() {
    vapply(
      self$samples,
      FUN = function(sample) sample$sample_name,
      FUN.VALUE = character(1)
    )
  },

  get_sample_name = function(sample_id) {
    if (!is.character(sample_id) ||
        length(sample_id) != 1L ||
        !sample_id %in% names(self$samples)) {
      stop(
        "`sample_id` must be the ID of a sample in this project.",
        call. = FALSE
      )
    }

    self$samples[[sample_id]]$sample_name
  },


  get_flow_frame = function(sample_id) {
    if (!is.character(sample_id) ||
        length(sample_id) != 1L ||
        !sample_id %in% names(self$samples)) {
      stop(
        "`sample_id` must be the ID of a sample in this project.",
        call. = FALSE
      )
    }

    self$samples[[sample_id]]$flow_frame
  },

  get_fcs_info = function(sample_id) {
    if (!is.character(sample_id) ||
        length(sample_id) != 1L ||
        !sample_id %in% names(self$samples)) {
      stop(
        "`sample_id` must be the ID of a sample in this project.",
        call. = FALSE
      )
    }

    self$samples[[sample_id]]$fcs
  },

  get_channel_info = function(sample_ids = NULL) {
    # Default: return channel metadata for every sample.
    if (is.null(sample_ids)) {
      sample_ids <- self$get_sample_ids()
    }

    cytoTerra_validate_sample_ids(self, sample_ids)

    channel_info <- lapply(sample_ids, function(sample_id) {
      flow_frame <- self$get_flow_frame(sample_id)

      expression_channels <- colnames(flowCore::exprs(flow_frame))

      parameter_data <- Biobase::pData(
        flowCore::parameters(flow_frame)
      )

      # FCS $PnN: instrument parameter / expression-matrix channel name.
      parameter_names <- parameter_data$name

      if (is.null(parameter_names)) {
        parameter_names <- expression_channels
      }

      parameter_names <- as.character(parameter_names)

      # FCS $PnS: marker or antibody description, when supplied.
      parameter_descriptions <- parameter_data$desc

      if (is.null(parameter_descriptions)) {
        parameter_descriptions <- rep(
          NA_character_,
          length(parameter_names)
        )
      }

      parameter_descriptions <- as.character(parameter_descriptions)

      # If FCS metadata has no useful description, preserve NA rather than
      # incorrectly treating the detector name as an antibody marker.
      parameter_descriptions[
        is.na(parameter_descriptions) |
          !nzchar(parameter_descriptions)
      ] <- NA_character_

      # These classifications are intentionally conservative. A detector
      # channel not named FSC, SSC, or Time is "fluorescence_or_other";
      # it may still require visual/metadata confirmation.
      parameter_class <- ifelse(
        grepl("^FSC", parameter_names, ignore.case = TRUE),
        "scatter_forward",
        ifelse(
          grepl("^SSC", parameter_names, ignore.case = TRUE),
          "scatter_side",
          ifelse(
            grepl("^(TIME|EVENT)", parameter_names, ignore.case = TRUE),
            "time_or_event",
            "fluorescence_or_other"
          )
        )
      )

      # Normally `parameter_names` and `expression_channels` are identical.
      # This explicitly reports whether the named FCS parameter is available
      # for expression-value extraction.
      in_expression_matrix <- parameter_names %in% expression_channels

      data.frame(
        sample_id = sample_id,
        sample_name = self$get_sample_name(sample_id),
        parameter_index = seq_along(parameter_names),
        channel = parameter_names,
        description = parameter_descriptions,
        parameter_class = parameter_class,
        is_fluorescence_candidate =
          parameter_class == "fluorescence_or_other",
        in_expression_matrix = in_expression_matrix,
        stringsAsFactors = FALSE
      )
    })

    channel_info <- do.call(rbind, channel_info)
    rownames(channel_info) <- NULL

    channel_info
  }
)
