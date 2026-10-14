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

remotes::install_github("YOUR-GITHUB-USERNAME/cytoTerra")
```

Replace `YOUR-GITHUB-USERNAME` with your GitHub username or organization.

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

## License

MIT License. See [`LICENSE`](LICENSE).
