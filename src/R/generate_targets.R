## WHO Raw data processing script
library(data.table)
library(here)
setwd(here())

TB_burden <- fread("data/raw_data/TB_burden_countries_2021-04-01.csv")
Agewise_TB_burden <- fread("data/raw_data/TB_burden_age_sex_2021-04-01.csv")

IN_burden <- TB_burden[country == "India", .(
  Year = year,
  IncidenceRate_mid = e_inc_100k,
  IncidenceRate_lo = e_inc_100k_lo,
  IncidenceRate_hi = e_inc_100k_hi,
  MortRate_mid = e_mort_100k,
  MortRate_lo = e_mort_100k_lo,
  MortRate_hi = e_mort_100k_hi
)]

IN_incidence <- dcast(melt(IN_burden[Year %in% range(Year)], id.vars = "Year", measure.vars = patterns("^IncidenceRate"), value.name = c("Indicator")), Year ~ variable)[, Indicator := "IncidenceRate"]
setnames(IN_incidence, new = c("Year", "mid", "lo", "hi", "Indicator"))

IN_Mort <- dcast(melt(IN_burden[Year %in% range(Year)], id.vars = "Year", measure.vars = patterns("^MortRate"), value.name = c("Indicator")), Year ~ variable)[, Indicator := "MortRate"]
setnames(IN_Mort, new = c("Year", "mid", "lo", "hi", "Indicator"))

IN_prev <- data.table(
  Year = 2015,
  mid = 253,
  lo = 195,
  hi = 312,
  Indicator = "PrevalenceRate"
)

IN_targets <- rbind(IN_prev, IN_Mort, IN_incidence)[, AgeGrp := "[0,99]"]

## Age-wise incidence

IN_agewise <- Agewise_TB_burden[country == "India" & sex == "a" & risk_factor == "all"][, .(Year = year, age_group, mid = best, lo, hi)]
Agewise_TB_burden[country == "India" & age_group == "65plus", .(age_group, best = sum(best), lo = sum(lo), hi = sum(hi))]
IN_inc_65p <- Agewise_TB_burden[country == "India" & age_group == "65plus", .(Year = year, age_group, mid = sum(best), lo = sum(lo), hi = sum(hi))][1, ]
IN_agewise <- rbind(IN_agewise, IN_inc_65p)
IN_agewise[, AgeGrp := fcase(
  age_group == "0-14", "[0,15)",
  age_group == "15plus", "[15,99]",
  age_group == "all", "[0,99]",
  age_group == "65plus", "[65,99]"
)][, age_group := NULL]

MIN_agewise <- melt(IN_agewise, id.vars = c("Year", "AgeGrp"))

# Population
IN_raw_pop <- rbind(
  fread("data/raw_data/WPP2019_PopulationBySingleAgeSex_1950-2019.csv"),
  fread("data/raw_data/WPP2019_PopulationBySingleAgeSex_2020-2100.csv")
)[Location == "India", .(
  Year = Time,
  AgeGrp,
  Pop = PopTotal
)]
tpop <- fread("data/raw_data/data_total_population_1900-2099.csv")
setnames(tpop, old = "99plus", new = "99")
mtpop <- melt(tpop, id.vars = "Year")[, variable := as.numeric(variable) - 1]
bag <- c(14, 99)
binned <- mtpop[, .(Year, AgeGrp = cut(variable, c(0, 15, 99), right = FALSE, include.lowest = TRUE), PopTotal = value)][, .(PopTotal = sum(PopTotal)), by = .(Year, AgeGrp)]
binned_65p <- mtpop[variable %between% c(65, 99), .(Year, AgeGrp = "[65,99]", PopTotal = value)][, .(PopTotal = sum(PopTotal)), by = .(Year, AgeGrp)]

total <- mtpop[, .(PopTotal = sum(value), AgeGrp = "[0,99]"), by = Year]
IN_pop <- rbind(binned, binned_65p, total)

IN_inc_agewise <- merge(IN_agewise, IN_pop[Year == 2019])[, .(Year,
  AgeGrp,
  mid = mid / (PopTotal * 1e3) * 1e5,
  lo = lo / (PopTotal * 1e3) * 1e5,
  hi = hi / (PopTotal * 1e3) * 1e5
)][, Indicator := "IncidenceRate"]

IN_final <- rbind(IN_targets, IN_inc_agewise)
setcolorder(IN_final, c("Year", "Indicator", "AgeGrp"))
IN_final <- IN_final[c(1:6, 8, 9)]

fwrite(IN_final, "data/proc_data/targets.csv")
