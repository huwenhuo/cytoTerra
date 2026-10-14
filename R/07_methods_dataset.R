# data extract et al

cytoTerra_methods_dataset <- list(

  #' Extract sequentially gated events for downstream analysis.
  #'
  #' Applies the requested gate and every ancestor gate in its hierarchy for
  #' each sample. The returned object contains the original flow-expression
  #' channels, with columns named from FCS channel descriptions when available.
  #'
  #' This method does not write files and does not modify gates or flowFrames.
  #' It is intended for downstream workflows such as FlowSOM.
  #'
  #' @param gate_name Name of the final gate to extract.
  #' @param sample_ids Optional character vector of sample IDs. If `NULL`,
  #'   all project samples are considered.
  #' @param scale Expression scale returned in the result: `"raw"`, `"asinh"`,
  #'   `"log1p"`, or `"gate"`. `"gate"` uses each sample's final gate transform.
  #' @param cofactor Positive asinh cofactor. Required for `"asinh"` unless the
  #'   gate has a stored cofactor when `scale = "gate"`.
  #' @param max_events Optional maximum number of gated events retained per
  #'   sample. If `NULL`, all gated events are returned.
  #' @param seed Optional random seed used when `max_events` requires sampling.
  #' @param include_sample_id Whether to include a `sample_id` column.
  #' @param include_event_id Whether to include an `event_id` column containing
  #'   original row positions in each sample's flowFrame.
  #'
  #' @return A `data.table::data.table`. Provenance is available through
  #'   `attr(result, "cytoTerra_manifest")` and
  #'   `attr(result, "cytoTerra_channel_map")`.

extract_gated_events = function(
    gate_names,
    sample_ids = NULL,
    scale = c("raw", "asinh", "log1p", "gate"),
    cofactor = NULL,
    max_events = NULL,
    seed = NULL,
    include_sample_id = TRUE,
    include_sample_name = TRUE,
    include_event_id = TRUE,
    include_matched_gates = TRUE
) {
    scale <- match.arg(scale)

    if (!is.character(gate_names) || length(gate_names) < 1L ||
        anyNA(gate_names) || any(!nzchar(gate_names))) {
        stop(
            "`gate_names` must be a non-empty character vector of gate names.",
            call. = FALSE
        )
    }
    gate_names <- unique(gate_names)

    if (is.null(sample_ids)) {
        sample_ids <- self$sample_ids
    }

    if (!is.character(sample_ids) || length(sample_ids) == 0L ||
        anyNA(sample_ids) || any(!nzchar(sample_ids))) {
        stop(
            "`sample_ids` must be NULL or a non-empty character vector.",
            call. = FALSE
        )
    }

    sample_ids <- unique(sample_ids)
    cytoTerra_validate_sample_ids(self, sample_ids)

    if (!is.null(cofactor) &&
        (!is.numeric(cofactor) || length(cofactor) != 1L ||
         is.na(cofactor) || !is.finite(cofactor) || cofactor <= 0)) {
        stop(
            "`cofactor` must be NULL or one positive finite number.",
            call. = FALSE
        )
    }

    if (!is.null(max_events) &&
        (!is.numeric(max_events) || length(max_events) != 1L ||
         is.na(max_events) || !is.finite(max_events) ||
         max_events < 1L || max_events != as.integer(max_events))) {
        stop(
            "`max_events` must be NULL or one positive integer.",
            call. = FALSE
        )
    }

    if (!is.null(seed) &&
        (!is.numeric(seed) || length(seed) != 1L ||
         is.na(seed) || !is.finite(seed) ||
         seed != as.integer(seed))) {
        stop("`seed` must be NULL or one integer.", call. = FALSE)
    }

    if (!is.logical(include_sample_id) || length(include_sample_id) != 1L ||
        is.na(include_sample_id)) {
        stop("`include_sample_id` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.logical(include_sample_name) || length(include_sample_name) != 1L ||
        is.na(include_sample_name)) {
        stop("`include_sample_name` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.logical(include_event_id) || length(include_event_id) != 1L ||
        is.na(include_event_id)) {
        stop("`include_event_id` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.logical(include_matched_gates) ||
        length(include_matched_gates) != 1L ||
        is.na(include_matched_gates)) {
        stop("`include_matched_gates` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(seed)) {
        set.seed(as.integer(seed))
    }

    if (!is.null(max_events)) {
        max_events <- as.integer(max_events)
    }

    get_gate_path <- function(current_sample_id, target_gate_name) {
        target_gate <- self$get_gates(
            sample_id = current_sample_id,
            gate_name = target_gate_name
        )

        if (is.null(target_gate)) {
            return(NULL)
        }

        gate_path <- target_gate_name
        visited <- target_gate_name
        current_gate <- target_gate

        while (!is.null(current_gate$parent) &&
               !identical(current_gate$parent, "root")) {

            parent_gate_name <- current_gate$parent

            if (!is.character(parent_gate_name) ||
                length(parent_gate_name) != 1L ||
                is.na(parent_gate_name) ||
                !nzchar(parent_gate_name)) {
                stop(
                    "Gate '", target_gate_name, "' in sample '",
                    current_sample_id, "' has an invalid parent definition.",
                    call. = FALSE
                )
            }

            if (parent_gate_name %in% visited) {
                stop(
                    "A cycle was found in the gate hierarchy for sample '",
                    current_sample_id, "': ",
                    paste(c(visited, parent_gate_name), collapse = " > "),
                    ".",
                    call. = FALSE
                )
            }

            parent_gate <- self$get_gates(
                sample_id = current_sample_id,
                gate_name = parent_gate_name
            )

            if (is.null(parent_gate)) {
                stop(
                    "Gate '", target_gate_name, "' in sample '",
                    current_sample_id, "' refers to parent gate '",
                    parent_gate_name, "', but that parent gate was not found.",
                    call. = FALSE
                )
            }

            gate_path <- c(parent_gate_name, gate_path)
            visited <- c(visited, parent_gate_name)
            current_gate <- parent_gate
        }

        list(
            target_gate = target_gate,
            gate_path = gate_path
        )
    }

    get_channel_map <- function(flow_frame, expression_matrix) {
        parameter_data <- tryCatch(
            Biobase::pData(flowCore::parameters(flow_frame)),
            error = function(error) NULL
        )

        parameter_names <- colnames(expression_matrix)
        channel_descriptions <- rep(NA_character_, length(parameter_names))

        if (!is.null(parameter_data) && "name" %in% colnames(parameter_data)) {
            parameter_match <- match(
                parameter_names,
                as.character(parameter_data$name)
            )

            if ("desc" %in% colnames(parameter_data)) {
                channel_descriptions <- as.character(
                    parameter_data$desc[parameter_match]
                )
            }
        }

        usable_description <- !is.na(channel_descriptions) &
            nzchar(channel_descriptions)

        exported_names <- parameter_names
        exported_names[usable_description] <-
            channel_descriptions[usable_description]

        reserved_names <- c(
            "sample_id",
            "sample_name",
            "event_id",
            "matched_gates"
        )

        exported_names[exported_names %in% reserved_names] <- paste0(
            exported_names[exported_names %in% reserved_names],
            "_channel"
        )

        exported_names <- make.unique(exported_names)

        data.table::data.table(
            fcs_parameter_name = parameter_names,
            channel_description = channel_descriptions,
            exported_column_name = exported_names
        )
    }

    transform_expression <- function(
        expression_matrix,
        requested_scale,
        current_gate,
        current_sample_id
    ) {
        effective_scale <- requested_scale
        effective_cofactor <- cofactor

        if (identical(requested_scale, "gate")) {
            effective_scale <- current_gate$transform

            if (is.null(effective_scale) ||
                !effective_scale %in% c("raw", "none", "asinh", "log1p")) {
                effective_scale <- "none"
            }

            if (is.null(effective_cofactor) &&
                !is.null(current_gate$cofactor) &&
                is.numeric(current_gate$cofactor) &&
                length(current_gate$cofactor) == 1L &&
                is.finite(current_gate$cofactor) &&
                current_gate$cofactor > 0) {
                effective_cofactor <- current_gate$cofactor
            }
        }

        if (identical(effective_scale, "raw") ||
            identical(effective_scale, "none")) {
            return(list(
                values = expression_matrix,
                scale = "raw",
                cofactor = NA_real_
            ))
        }

        if (identical(effective_scale, "asinh")) {
            if (is.null(effective_cofactor)) {
                stop(
                    "`cofactor` must be supplied when `scale = \"asinh\"`.",
                    call. = FALSE
                )
            }

            return(list(
                values = asinh(expression_matrix / effective_cofactor),
                scale = "asinh",
                cofactor = effective_cofactor
            ))
        }

        if (identical(effective_scale, "log1p")) {
            if (any(expression_matrix < 0, na.rm = TRUE)) {
                stop(
                    "Cannot apply `log1p` to sample '", current_sample_id,
                    "' because one or more channels contain negative values.",
                    call. = FALSE
                )
            }

            return(list(
                values = log1p(expression_matrix),
                scale = "log1p",
                cofactor = NA_real_
            ))
        }

        stop(
            "Unsupported scale '", effective_scale,
            "' for sample '", current_sample_id, "'.",
            call. = FALSE
        )
    }

    extracted_tables <- list()
    manifest_rows <- list()
    channel_maps <- list()

    for (current_sample_id in sample_ids) {
        flow_frame <- self$get_flow_frame(current_sample_id)
        expression_matrix <- flowCore::exprs(flow_frame)
        total_events <- nrow(expression_matrix)

        # Retrieve sample_name directly from self$sample_names using sample_id
        current_sample_name <- self$sample_names[[current_sample_id]]

        if (is.null(current_sample_name) || length(current_sample_name) == 0L) {
            current_sample_name <- tryCatch(
                Biobase::pData(flow_frame)$sample_name,
                error = function(e) NULL
            )
        }

        # Get hierarchy information for every requested gate.
        gate_information <- lapply(
            gate_names,
            function(current_gate_name) {
                get_gate_path(
                    current_sample_id = current_sample_id,
                    target_gate_name = current_gate_name
                )
            }
        )
        names(gate_information) <- gate_names

        found_gate_names <- names(gate_information)[
            !vapply(gate_information, is.null, logical(1))
        ]

        missing_gate_names <- setdiff(gate_names, found_gate_names)

        if (length(found_gate_names) == 0L) {
            manifest_rows[[current_sample_id]] <- data.table::data.table(
                sample_id = current_sample_id,
                sample_name = current_sample_name,
                gate_names = paste(gate_names, collapse = "; "),
                total_events = total_events,
                union_gated_events = 0L,
                returned_events = 0L,
                scale = scale,
                cofactor = if (is.null(cofactor)) NA_real_ else cofactor,
                status = "skipped: none of the requested gates were found"
            )
            next
        }

        if (identical(scale, "gate") && length(found_gate_names) > 1L) {
            stop(
                "`scale = \"gate\"` cannot be used with multiple gates because ",
                "the gates may have different transformation settings. Use ",
                "`scale = \"raw\"`, `\"asinh\"`, or `\"log1p\"` instead.",
                call. = FALSE
            )
        }

        # Obtain one membership mask for each requested, available gate.
        gate_masks <- lapply(found_gate_names, function(current_gate_name) {
            current_information <- gate_information[[current_gate_name]]

            path_masks <- lapply(
                current_information$gate_path,
                function(path_gate_name) {
                    current_mask <- self$evaluate_gate(
                        sample_id = current_sample_id,
                        gate_name = path_gate_name
                    )

                    if (!is.logical(current_mask) ||
                        length(current_mask) != total_events) {
                        stop(
                            "Gate '", path_gate_name,
                            "' returned an invalid membership mask for sample '",
                            current_sample_id, "'.",
                            call. = FALSE
                        )
                    }

                    current_mask[is.na(current_mask)] <- FALSE
                    current_mask
                }
            )

            Reduce(`&`, path_masks)
        })
        names(gate_masks) <- found_gate_names

        # UNION of all gate masks: overlapping events occur only once.
        union_gate_mask <- Reduce(`|`, gate_masks)
        gated_event_ids <- which(union_gate_mask)
        union_gated_events <- length(gated_event_ids)

        if (union_gated_events == 0L) {
            manifest_rows[[current_sample_id]] <- data.table::data.table(
                sample_id = current_sample_id,
                sample_name = current_sample_name,
                gate_names = paste(gate_names, collapse = "; "),
                gates_found = paste(found_gate_names, collapse = "; "),
                gates_missing = if (length(missing_gate_names)) {
                    paste(missing_gate_names, collapse = "; ")
                } else {
                    NA_character_
                },
                total_events = total_events,
                union_gated_events = 0L,
                returned_events = 0L,
                scale = scale,
                cofactor = if (is.null(cofactor)) NA_real_ else cofactor,
                status = "skipped: no gated events"
            )
            next
        }

        # Record all gates that matched each unique event.
        matched_gates <- vapply(gated_event_ids, function(event_id) {
            paste(
                names(gate_masks)[vapply(
                    gate_masks,
                    function(mask) mask[[event_id]],
                    logical(1)
                )],
                collapse = "; "
            )
        }, character(1))

        # Sampling happens after union/deduplication.
        if (!is.null(max_events) &&
            length(gated_event_ids) > max_events) {

            selected_index <- sort(sample.int(
                length(gated_event_ids),
                size = max_events,
                replace = FALSE
            ))

            gated_event_ids <- gated_event_ids[selected_index]
            matched_gates <- matched_gates[selected_index]
        }

        selected_expression <- expression_matrix[
            gated_event_ids,
            ,
            drop = FALSE
        ]

        transform_gate <- gate_information[[found_gate_names[[1]]]]$target_gate

        transformed <- transform_expression(
            expression_matrix = selected_expression,
            requested_scale = scale,
            current_gate = transform_gate,
            current_sample_id = current_sample_id
        )

        current_channel_map <- get_channel_map(
            flow_frame = flow_frame,
            expression_matrix = expression_matrix
        )

        current_channel_map <- data.table::as.data.table(current_channel_map)
        current_channel_map[, sample_id := current_sample_id]

        data.table::setcolorder(
            current_channel_map,
            c(
                "sample_id",
                "fcs_parameter_name",
                "channel_description",
                "exported_column_name"
            )
        )

        channel_maps[[current_sample_id]] <- current_channel_map

        extracted_sample <- data.table::as.data.table(transformed$values)

        data.table::setnames(
            extracted_sample,
            current_channel_map$exported_column_name
        )

        if (isTRUE(include_matched_gates)) {
            extracted_sample[, matched_gates := matched_gates]
        }

        if (isTRUE(include_event_id)) {
            extracted_sample[, event_id := as.integer(gated_event_ids)]
        }

        if (isTRUE(include_sample_name)) {
            extracted_sample[, sample_name := current_sample_name]
        }

        if (isTRUE(include_sample_id)) {
            extracted_sample[, sample_id := current_sample_id]
        }

        preferred_columns <- c(
            if (isTRUE(include_sample_id)) "sample_id",
            if (isTRUE(include_sample_name)) "sample_name",
            if (isTRUE(include_event_id)) "event_id",
            if (isTRUE(include_matched_gates)) "matched_gates"
        )

        data.table::setcolorder(
            extracted_sample,
            c(preferred_columns, setdiff(names(extracted_sample), preferred_columns))
        )

        extracted_tables[[current_sample_id]] <- extracted_sample

        manifest_rows[[current_sample_id]] <- data.table::data.table(
            sample_id = current_sample_id,
            sample_name = current_sample_name,
            gate_names = paste(gate_names, collapse = "; "),
            gates_found = paste(found_gate_names, collapse = "; "),
            gates_missing = if (length(missing_gate_names)) {
                paste(missing_gate_names, collapse = "; ")
            } else {
                NA_character_
            },
            total_events = total_events,
            union_gated_events = union_gated_events,
            returned_events = nrow(extracted_sample),
            scale = transformed$scale,
            cofactor = transformed$cofactor,
            status = if (length(gated_event_ids) < union_gated_events) {
                "included; sampled to max_events"
            } else if (length(missing_gate_names) > 0L) {
                "included; some requested gates were not found"
            } else {
                "included"
            }
        )
    }

    result <- data.table::rbindlist(
        extracted_tables,
        use.names = TRUE,
        fill = TRUE
    )

    manifest <- data.table::rbindlist(
        manifest_rows,
        use.names = TRUE,
        fill = TRUE
    )

    channel_map <- data.table::rbindlist(
        channel_maps,
        use.names = TRUE,
        fill = TRUE
    )

    data.table::setattr(result, "cytoTerra_manifest", manifest)
    data.table::setattr(result, "cytoTerra_channel_map", channel_map)
    data.table::setattr(result, "cytoTerra_gate_names", gate_names)
    data.table::setattr(result, "cytoTerra_requested_scale", scale)
    data.table::setattr(result, "cytoTerra_max_events", max_events)

    result
}
)
