library(socialmixr)
data("polymod")
cm <- contact_matrix(survey = polymod, age.limits = c(0, 15, 65), symmetric = TRUE, split = TRUE)

cm_raw <- contact_matrix(survey = polymod, age.limits = c(0, 15, 65), symmetric = TRUE)

cm_raw$matrix
cm_recalc <- sweep(
  x = cm$mean.contacts * cm$normalisation * cm$matrix * cm$contacts,
  MARGIN = 2, FUN = "*", STATS = cm$demography$proportion
)

stopifnot(all.equal(cm_recalc, cm_raw$matrix, check.attributes = FALSE))

cm_gamma <- cm$mean.contacts * cm$normalisation * cm$matrix * cm$contacts

data.table::fwrite(x = cm_gamma, file = here::here("data/proc_data/cm_gamma.csv"), col.names = F)
