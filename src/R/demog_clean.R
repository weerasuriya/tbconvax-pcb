# 3 age class difference model
library(data.table)
library(here)
setwd(here())

# Total population matrix
tpop <- fread("data/raw_data/data_total_population_1900-2099.csv")

# Temporary max_age corrrection
tpop[, `99plus` := `98` * 0.9]

# Births
births <- tpop[, .(Year, B = `0`)]

# Convert to long
mtpop <- melt(tpop, id.vars = "Year")[, variable := as.numeric(variable) - 1]

# Bin into age groups / create border age groups - 14, 64, 99
bag <- c(14, 64, 99)
binned <- mtpop[, AgeGrp := cut(variable, c(0, bag + 1), right = FALSE, include.lowest = TRUE)]

# Total population by group
summed <- binned[, .(Total = sum(value)), by = .(Year, AgeGrp)]

# Select border age groups
border <- mtpop[variable %in% (bag + 1)]

# Line up next-year exit groups
summed[, ex_grp := fcase(
  AgeGrp == "[0,15)", 1,
  AgeGrp == "[15,65)", 2,
  AgeGrp == "[65,100]", 3
)]

border[, ex_grp := fcase(
  AgeGrp == "[15,65)", 1,
  AgeGrp == "[65,100]", 2
)]

# Adjustment for 2101 age_out group
border <- rbind(
  border,
  list(Year = 2101, variable = 15, value = 0, AgeGrp = "[15,65)", ex_grp = 1),
  list(Year = 2101, variable = 65, value = 0, AgeGrp = "[65,100]", ex_grp = 2)
)

# Exit group 3
exg3 <- data.table(
  Year = 1950:2101,
  ex_grp = 3,
  AgeGrp = "[65,100]",
  variable = 99,
  value = 0
)

border <- rbind(border, exg3)

# Exit dataframe, with age-out quantity
exit <- merge(summed[, .(Year, AgeGrp, ex_grp, Total)],
  border[, .(Year = Year - 1, variable, value, ex_grp)],
  by = c("Year", "ex_grp")
)[order(ex_grp, Year)]

# Mortality calculation
# Calculate survival matrices
pm <- data.matrix(tpop[, -1])
pma <- pm[2:151, 2:100] / pm[1:150, 1:99]
pma <- cbind(pma, 0)
pma <- rbind(pma, 0)

mort <- data.table(
  Year = 1950:2100,
  pma
)

setnames(mort, new = c("Year", 0:99))

lmort <- melt(mort, id.vars = "Year", variable.name = "variable", value.name = "surv")
lmort[, variable := as.numeric(as.character(variable))]
sym <- merge(mtpop, lmort)

# Calculate deaths
sym[, deaths := value * (1 - surv)]

mort_summed <- sym[, .(pop = sum(value), deaths = sum(deaths)), by = .(Year, AgeGrp)]
mort_summed[deaths < 0, deaths := 0]
mort_summed[, mort_rate := deaths / pop]

# Annual final
final <- merge(mort_summed, exit)[, .(Year, AgeGrp, pop, mort_rate, ao_rate = value / pop)]

# AO adjustment to allow aging only in one step
final_2s <- merge(mort_summed, exit)[, .(Year, AgeGrp, pop, mort_rate, ao_rate = value / (pop * exp(log(1 - mort_rate) * 0.5)))]
final_2s[Year == 2100, ao_rate := 0]

#### Testing ####
wfinal <- dcast(final, Year ~ AgeGrp, value.var = list("pop", "mort_rate", "ao_rate"))
wfinal_2s <- dcast(final_2s, Year ~ AgeGrp, value.var = list("pop", "mort_rate", "ao_rate"))

# Duplicate rows for 2-step
dupfinal <- wfinal_2s[rep(1:nrow(wfinal_2s), each = 2), ]
dupfinal[, Index := 1:302]

# Convert survival -> mortality rate
dupfinal[, `mort_rate_[0,15)` := 1 - exp(log(1 - `mort_rate_[0,15)`) * 0.5)]
dupfinal[, `mort_rate_[15,65)` := 1 - exp(log(1 - `mort_rate_[15,65)`) * 0.5)]
dupfinal[, `mort_rate_[65,100]` := 1 - exp(log(1 - `mort_rate_[65,100]`) * 0.5)]

# Zero out aging in first annual step
dupfinal[seq(1, 301, 2), `ao_rate_[0,15)` := 0]
dupfinal[seq(1, 301, 2), `ao_rate_[15,65)` := 0]
dupfinal[seq(1, 301, 2), `ao_rate_[65,100]` := 0]

# Duplicate births
dbirths <- births[rep(1:151, each = 2)]
# Zero out births in second annual step
dbirths[seq(2, 302, 2), B := 0]

#### Testing Loop ####
mat <- matrix(0, 3, 3)
tmd <- matrix(0, nrow = 302, ncol = 3)
tmd[1, 1] <- final_2s[1, pop]
tmd[1, 2] <- final_2s[2, pop]
tmd[1, 3] <- final_2s[3, pop]

for (i in 2:302) {
  mat[1, 1] <- (1 - dupfinal[i - 1, `mort_rate_[0,15)`] - dupfinal[i - 1, `ao_rate_[0,15)`])
  mat[2, 1] <- dupfinal[i - 1, `ao_rate_[0,15)`]
  mat[2, 2] <- (1 - dupfinal[i - 1, `mort_rate_[15,65)`] - dupfinal[i - 1, `ao_rate_[15,65)`])
  mat[3, 2] <- dupfinal[i - 1, `ao_rate_[15,65)`]
  mat[3, 3] <- (1 - dupfinal[i - 1, `mort_rate_[65,100]`])

  tmd[i, ] <- mat %*% tmd[i - 1, ]
  tmd[i, 1] <- tmd[i, 1] + dbirths[i, B]
}

# Simulated total population
sim <- data.table(
  Year = rep(1950:2100, each = 2),
  tmd
)

sim[, V4 := V1 + V2 + V3]
setnames(sim, old = paste0("V", 1:4), new = c("sim_[0,15)", "sim_[15,65)", "sim_[65,99]", "sim_[0,99]"))

# Combine into operating dataframe for model
# Duplicate 1950-1999 for burn in
prelimdemog <- dupfinal[Year %between% c(1950, 1999)][, Year := rep(1900:1949, each = 2)]
prelimdbirths <- dbirths[Year %between% c(1950, 1999)][, Year := rep(1900:1949, each = 2)]
prelimsim <- sim[Year %between% c(1950, 1999)][, Year := rep(1900:1949, each = 2)]

# Stack
manual_demog <- rbind(prelimdemog, dupfinal)[, Index := 1:402]
manual_births <- rbind(prelimdbirths, dbirths)[, Index := 1:402]
manual_sim <- rbind(prelimsim, sim)[, Index := 1:402]

# Merge births and mort/ao
manual_demog <- merge(merge(manual_demog, manual_births, by = c("Index", "Year")), manual_sim, by = c("Index", "Year"))

# Annual total pop
fwrite(manual_demog[Index %in% seq(1, 400, 2), .(Year, `sim_[0,15)`, `sim_[15,65)`, `sim_[65,99]`, `sim_[0,99]`)], file = "data/proc_data/3ac_sim_total.csv")

# Write out
fwrite(final, file = "data/proc_data/3ac_demog.csv")
fwrite(manual_demog, file = "data/proc_data/3ac_manual.csv")
