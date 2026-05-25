# legacy/ — archived pre-API-rename code

This directory holds code that is no longer part of the active
`causalBKMR` package. It is preserved for provenance only and is
excluded from the R package build via `.Rbuildignore`.

## R-old/

The previous-generation `run_gbkmr_panel()`, `prepare_gbkmr_data()`,
and related helpers. The signatures used the legacy `mediator_*`
naming (e.g., `mediator_basenames`, `fit_mediators`), which were
renamed to `confounder_*` in the active code under `R/`.

**Do not source these files.** They will shadow the active API and
cause "unused argument" errors at call sites that pass the new
argument names. If you genuinely need to reproduce a historical run,
do so in an isolated R session and revert afterwards.
