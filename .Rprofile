source("renv/activate.R")
if (!interactive()) {
  options(error = quote(dump.frames(
    dumpto = sprintf("%s_crash_dump", format(Sys.time(), "%Y-%m-%d-%H%M")),
    to.file = TRUE,
    include.GlobalEnv = TRUE
  )))
}
