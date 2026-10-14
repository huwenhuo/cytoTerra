#' Create a cytometry analysis project from FCS files
#'
#' Reads one or more FCS files and creates a mutable `CytoTerraProject`
#' object. The project retains each original `flowFrame` without modification.
#'
#' @param fcs_files A non-empty character vector of FCS file paths.
#' @param sample_ids Optional unique sample IDs. If `NULL`, IDs are derived
#'   from FCS filenames without their extensions.
#' @param transformation Logical. Whether to apply FCS transformations while
#'   reading. Defaults to `FALSE`.
#' @param truncate_max_range Logical. Whether to truncate expression values
#'   to the recorded FCS parameter range. Defaults to `FALSE`.
#' @param channel_labels Optional data frame containing channel annotations.
#'
#' @return A `CytoTerraProject` R6 object.
#' @export
#' @import data.table
cytoTerra <- function(
    fcs_files,
    sample_ids = NULL,
    sample_names = NULL, 
    sample_meta = NULL,
    transformation = FALSE,
    truncate_max_range = FALSE,
    channel_labels = NULL
) {
  if (!is.character(fcs_files) || length(fcs_files) == 0L) {
    stop(
      "`fcs_files` must be a non-empty character vector of FCS paths.",
      call. = FALSE
    )
  }

  if (anyNA(fcs_files) || any(!nzchar(fcs_files))) {
    stop(
      "`fcs_files` cannot contain missing or empty paths.",
      call. = FALSE
    )
  }

  missing_files <- fcs_files[!file.exists(fcs_files)]

  if (length(missing_files) > 0L) {
    stop(
      "These FCS file(s) do not exist:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  # Full original FCS filename, excluding ".fcs".
  if (is.null(sample_names)) {
    sample_names <- tools::file_path_sans_ext(basename(fcs_files))
  }

  if (!is.character(sample_names) ||
      length(sample_names) != length(fcs_files) ||
      anyNA(sample_names) ||
      any(!nzchar(sample_names))) {
    stop(
      "`sample_names` must be a non-empty character vector with ",
      "one name per FCS file.",
      call. = FALSE
    )
  }

  # Simple ID: everything before the first space.
  # "A1 421_Myeloids" becomes "A1".
  if (is.null(sample_ids)) {
    sample_ids <- gsub(" .*", "", sample_names)
  }

  if (!is.character(sample_ids) ||
      length(sample_ids) != length(fcs_files) ||
      anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop(
      "`sample_ids` must be a unique, non-empty character vector with ",
      "one ID per FCS file.",
      call. = FALSE
    )
  }

  flow_set <- flowCore::read.flowSet(
    files = fcs_files,
    transformation = transformation,
    truncate_max_range = truncate_max_range
  )

  flowCore::sampleNames(flow_set) <- sample_ids
  Biobase::pData(flow_set)$sample_id <- sample_ids
  Biobase::pData(flow_set)$sample_name <- sample_names

  new_cyto_project(
    flow_set = flow_set,
    fcs_files = fcs_files,
    sample_ids = sample_ids,
    sample_names = sample_names,
    sample_meta = sample_meta,
    channel_labels = channel_labels
  )
}
