# Internal constructor for a CytoTerraProject.
#'
#' @noRd
new_cyto_project <- function(
    flow_set,
    fcs_files,
    sample_ids,
    sample_names,
    sample_meta = NULL, 
    channel_labels = NULL
) {
  CytoTerraProject$new(
    flow_set = flow_set,
    fcs_files = fcs_files,
    sample_ids = sample_ids,
    sample_names = sample_names,
    sample_meta = sample_meta,
    channel_labels = channel_labels
  )
}
