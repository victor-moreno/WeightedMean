# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

This project is a jamovi module developed using the `jmvcore` R package.

- **Install/Update Module**: From an R session (specifically jamovi's R environment), use:
  ```r
  jmvtools::install()
  ```
- **R Dependencies**: Ensure `jmvcore`, `R6`, and `jsonlite` are installed.

## Architecture

The module implements a weighted mean analysis where weights can be sourced from an external CSV or edited interactively.

### Project Structure
- `R/`: Contains the analysis implementation.
    - `weightedMean.b.R`: Contains the R6 class definitions for options, results, and the analysis base.
    - `weightedmean.h.R`: Contains the main analysis logic (`weightedMeanClass`) and the `weightedMean` entry point.
- `jamovi/`: Contains YAML configuration files that define the jamovi UI (options and results) and integration metadata.
- `DESCRIPTION`: Standard R package description file.
- `NAMESPACE`: R namespace definitions.
- `inst/`: Installation files.
- `data-raw/`: Test data (CSV/XLSX).

### Logic Flow
1. **UI Definition**: Jamovi reads `jamovi/*.yaml` to render the options panel.
2. **Execution**: When the analysis runs, `weightedMean()` is called, which instantiates `weightedMeanClass`.
3. **Weight Handling**: 
    - It attempts to read a CSV file specified in the options.
    - It applies overrides provided as a JSON string in the options (`tableOverrides`).
4. **Interactive Loop**: 
    - The analysis renders an HTML table in the results panel.
    - A JavaScript function captures user edits in this table.
    - These edits are sent back to jamovi via `postMessage` or `window.jamovi.setOption`, updating the `tableOverrides` option.
    - This triggers a re-run of the analysis with the new weights.
5. **Calculation**: Computes unweighted and weighted means for selected variables and an overall weighted mean.
