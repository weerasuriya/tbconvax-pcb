# Startup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(JuliaCall)
  library(EasyABC)
  library(data.table)
  library(bayesplot)
  library(fs)
  library(yaml)
  library(glue)
  library(git2r)
  library(qs)
  library(ragg)
  library(ggplot2)
  library(gridExtra)
})

# Set working directory
setwd(here())

# Empty list for metadata
metadata <- list()

# Check if workspace is clean
if (!interactive()) {
  git_status <- git2r::status()
  stopifnot(
    length(git_status$staged) == 0,
    length(git_status$unstaged) == 0
  )
  metadata$job_commit <- git2r::last_commit()$sha
}

parameter_names <- read_yaml("data/param_names.yml")
CMA <- read_yaml("data/CMA_indicator.yml")$CMA
metadata$CMA <- CMA

# JULIA -------------------------------------------------------------------

julia <- julia_setup()

# Call function library
julia$command("using Revise")
julia$command("using tbconvax")

# Initialise fixed inputs
init_object <- julia$call("init_mat", "data")

# Julia functions
target_extract <- function(x) julia$call("target_extract", x)
pst <- function(x) julia$call("pst", x)
btc <- function(x) julia$call("btc", x)
pobgen <- function(x) julia$call("pobgen", x)
main <- function(x) julia$call("main", variable_input = pst(x), fixed_input = init_object)

n_targ <- 9

mw <- function(x) {
  tryCatch(
    {
      julia$call("wam", x)
    },
    error = function(e) {
      return(rep(1e5, n_targ))
    }
  )
}

# ABC priors --------------------------------------------------------------

pr <- rep(
  list(c("unif", 0, 1)),
  25
)

# ABC seed ----------------------------------------------------------------

seed <- T
if (seed) {
  meta_init <- read_yaml(path("output", CMA, "params", "seeds", "latest.yml"))
  ip <- meta_init$ts
  metadata$seed_id <- meta_init$id
} else {
  ip <- NULL
  metadata$seed_id <- NA
}

# ABC run -----------------------------------------------------------------

metadata$start_time <- format(Sys.time(), "%Y-%m-%d-%H%M")
metadata$JN <- glue("{metadata$start_time}_{CMA}")
cat(sprintf("start: %s", metadata$start_time, "\n"))

abc_op <- ABC_mcmc(
  method = "Marjoram_original",
  model = mw,
  prior = pr,
  summary_stat_target = rep(0, n_targ),
  n_rec = 100000,
  dist_max = n_targ,
  progress_bar = TRUE,
  tab_normalization = rep(1, n_targ),
  verbose = TRUE,
  n_between_sampling = 100,
  init_param = ip
)

metadata$end_time <- format(Sys.time(), "%Y-%m-%d-%H%M")
cat(sprintf("start: %s", metadata$end_time, "\n"))

# Post-process results ----------------------------------------------------

post <- t(apply(X = abc_op$param, MARGIN = 1, FUN = pst))
colnames(post) <- parameter_names
colnames(abc_op$param) <- parameter_names

tgt <- abc_op$stats
colnames(tgt) <- paste0("T", 1:n_targ)

## Interactive analysis block ---------------------------------------------

ANALYSE <- TRUE

if (ANALYSE) {
  try({
    trace_plot <- mcmc_trace(post)
    dens_plot <- mcmc_dens(post)
    ggsave(
      filename = path("output", CMA, "params", "figs", metadata$JN, "trace", ext = "png"),
      plot = trace_plot,
      device = ragg::agg_png, units = "in", width = 50, height = 25, limitsize = F
    )
    ggsave(
      filename = path("output", CMA, "params", "figs", metadata$JN, "dens", ext = "png"),
      plot = dens_plot,
      device = ragg::agg_png, units = "in", width = 50, height = 50, limitsize = F
    )
  })
}

# QS dump -----------------------------------------------------------------

try({
  dir_create("assets")
  qs::qsave(x = abc_op, file = sprintf("assets/abc_op_%s.qs", metadata$JN))
})

# Write out ---------------------------------------------------------------

post_dt <- as.data.table(post)
psid <- apply(post_dt, MARGIN = 1, FUN = digest::digest)
post_dt[, PSID := psid]

opp <- dir_create(here("output", CMA, "params", "sets"))

path_jn_raw <- path(opp, paste(metadata$JN, "raw", sep = "_"), ext = "csv")
path_jn_raw_gz <- path(opp, paste(metadata$JN, "raw", sep = "_"), ext = "csv.gz")
path_jn <- path(opp, metadata$JN, ext = "csv")
u_post_dt <- unique(post_dt)
set.seed(123)
subsample <- u_post_dt[sample(1:nrow(u_post_dt), size = 1000), ]

WRITE <- TRUE

if (WRITE) {
  fwrite(x = post_dt, file = path_jn_raw)
  fwrite(x = post_dt, file = path_jn_raw_gz)
  fwrite(x = subsample, file = path_jn)
  metadata$paramset_hash <- digest::digest(path_jn_raw, file = TRUE)

  param_latest_path <- here("output", CMA, "metadata", paste0(metadata$JN, ".yml"))
  param_unique_path <- here("output", CMA, "metadata", "latest.yml")

  write_yaml(x = metadata, file = param_latest_path)
  write_yaml(x = metadata, file = param_unique_path)

  git2r::add(path = c(
    path_jn,
    param_latest_path,
    param_unique_path
  ), repo = "../tbconvax-output")

  git2r::commit(message = glue("added: new parameter set from ABC - {metadata$JN}"), repo = "../tbconvax-output")
}
