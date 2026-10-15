# cytoTerra

**cytoTerra** provides scalable, reproducible tools for flow-cytometry analysis
and gating of large FCS data sets in Yale HPC workflows.

It supports:

- FCS import, validation, and channel inspection
- Reproducible storage of gate definitions
- Automated gating and clustering workflows
- Manual/interactive gating through Shiny
- Use in JupyterLab notebooks with an R (`IRkernel`) kernel
- Secure Open OnDemand access to interactive Shiny workspaces

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")

remotes::install_github("huwenhuo/cytoTerra")
```

## Quick start

### Basic operation
```r
library(cytoTerra)

# `fcs_files` is a character vector of FCS file paths.
# `df` contains sample-level metadata in the same order as `fcs_files`.
# my typical file name is 'A1 421_Plate_001.fcs
# I use A1 as sample_id
# sample_name is usually a meaning full name such as 'drug A treat rep 1'

project <- cytoTerra(
  fcs_files = fcs_files[1:3],
  sample_ids = df$sample_id[1:3],
  sample_names = df$Tx[1:3]
)

# Automatically create an initial singlet gate for each sample.
project$auto_gate_singlets()

# Launch the interactive workspace for a reference sample.
# In the workspace, select `Singlets` as the parent population and
# manually create downstream gates (for example, live CD45+ cells,
# CD4+CD8− cells, and CD4−CD8+ cells).
project$launch_gating_workspace(
  sample_id = "H1"
)

# Copy the complete gate tree from H1 to every other sample.
# Existing gates are retained because overwrite = FALSE by default.
project$copy_gate_tree(
  from_sample = "H1",
  overwrite = FALSE
)

project$get_gate_tree(
  sample_id = "A1",
  include_root = TRUE
)

gating_file <- "~/workspace/flow/lymphoid_all_samples_gate_tree.rds"

project$export_all_gate_trees(
  file = gating_file
)

saveRDS(
  project,
  file = "lymphoid_cytoTerra_project.rds"
)

```

### Importing gates into a new project

Keep this as a separate section, because import is intended for a **new or different project**, not the project from which the gates were just exported.

If importing one reference-sample template:

### Reuse a gating template

Export a gate tree from a manually gated reference sample:

```r
project$export_gate_tree(
  sample_id = "H1",
  file = "H1_gate_template.rds"
)

new_project$import_gate_tree(
  file = "H1_gate_template.rds",
  to_samples = c("A1", "A2", "A3"),
  overwrite = FALSE
)

```


Use the all-sample export/import API only if you have implemented a matching method such as:

```r
new_project$import_all_gate_trees(
  file = "lymphoid_all_samples_gate_tree.rds",
  overwrite = FALSE
)

```

### Clustering gated populations

Use selected gated populations for integrated dimensionality reduction
and clustering.

```r
# channel names extract:
project$get_channel_info(sample_id = 'A1')

```

```r
# marker names can be these typed on the machine such as CD3
# alternatively, the channel name also works such as BUV395-A
# cytoTerra will automatically match the marker name or channel name

cluster_markers <- c(
  "CD69", "MHCII", "TCF7-TCF-1", "CD62L", "CD44",
  "IFNg", "CD25", "CD103", "GranzymeB", "CD8a",
  "GATA3", "PD-1", "CD5", "TOX", "Nkp46", "CD127",
  "Ki-67", "CCR7", "Tim-3", "FOXP3", "LAG-3"
)


project$run_harmony_umap_flowsom_fastpg(
  result_name = "run1",
  gate_names = c(
    "CD4pos_CD8neg",
    "CD4neg_CD8pos"
  ),
  markers = cluster_markers
)

clustering_data <- project$clustering_results[["run1"]]

str(clustering_data)

```

### Visualization 

Built-in visualization functions for UMAPs, marker feature plots,
cluster-composition dot plots, and related summaries are available.

```r
project$plot_umap()
project$plot_feature_expression("PD-1")
project$plot_dotplot(cluster_markers)

```

In addition, clustering results can be extracted directly for custom plot:

```r
clustering_data <- project$clustering_results[["run1"]]
```

## R dependencies

cytoTerra installs its required R dependencies automatically. For interactive
manual gating, the R environment used by your Jupyter R kernel must have
`shiny` available:

```r
install.packages("shiny")
```

## JupyterLab and Open OnDemand setup

Interactive gating starts a local Shiny application from the R notebook kernel.
To make that local application accessible in the browser through JupyterLab or
Open OnDemand, install
[`jupyter-server-proxy`](https://github.com/jupyterhub/jupyter-server-proxy)
in the **Python environment that runs Jupyter Server**.

For example:

```bash
conda activate py310
python -m pip install jupyter-server-proxy
```

Verify that the extension is enabled:

```bash
jupyter server extension list
```

The output should include something similar to:

```text
jupyter_server_proxy enabled
jupyter_server_proxy ... OK
```

If it is not enabled, run:

```bash
jupyter server extension enable --py jupyter_server_proxy
```

After installing or enabling the extension, restart the JupyterLab/Open
OnDemand session.

> `jupyter-server-proxy` is a Python/Jupyter dependency, not an R package.
> Therefore it belongs in these setup instructions rather than cytoTerra's
> `DESCRIPTION` file.

## Launching an interactive Shiny workspace

cytoTerra's interactive-gating function should print a URL similar to:

```text
https://ood-bouchet.ycrc.yale.edu/node/<compute-node>/<jupyter-port>/proxy/<shiny-port>/
```

The compute node and Jupyter port are assigned dynamically for each Open
OnDemand session. The function determines them automatically from the current
environment.

For security, the Shiny application runs only on the local compute node:

```r
host = "127.0.0.1"
```

Do not expose the Shiny process directly to the network. Access it using the
Open OnDemand `/proxy/` URL printed by cytoTerra.

## Planned Preprocessing Features

### Compensation

The package will support conventional fluorescence compensation using a lazy-processing approach in the future. 
Compensation matrices will be stored with the project and applied only when data are accessed for visualization, gating, or downstream analysis.

- Each sample may retain its own compensation matrix, such as a matrix imported from FCS metadata.
- The project may also define a shared compensation matrix to apply consistently across selected samples.
- The active analysis setting will be stored in `project$active$compensate`, indicating whether plotting and gating use:
  - raw, uncompensated data; or
  - compensated data.
- Raw event data will remain unchanged; compensation will be applied dynamically when needed.

Planned compensation tools include:

1. Estimating a compensation matrix from unstained and single-stained controls.
2. Reviewing and editing compensation matrices interactively through the Shiny interface.
3. Comparing raw and compensated channel-pair plots for quality control.
4. Saving compensation settings and matrix provenance with the project.

### Spectral Unmixing

Spectral unmixing is generally performed in instrument-vendor software, such as Cytek SpectroFlo or Sony ID7000 software, before FCS export. The package will support analysis of already unmixed data.

A future capability may support **post-acquisition re-unmixing** of recorded spectral data. This may be useful when reference controls, unmixing settings, or analysis decisions need to be revisited after acquisition.

Planned re-unmixing tools include:

1. Estimating an unmixing/reference-spectra table from reference controls.
2. Performing spectral unmixing using the estimated or user-supplied reference table.
3. Reviewing unmixing quality interactively.
4. Exporting re-unmixed events as FCS files for use in other analysis software.

> Native spectral unmixing is a longer-term feature and will require validation with appropriate single-stained reference controls, autofluorescence controls, and instrument-specific detector configurations.

## License

MIT License. See [`LICENSE`](LICENSE).
