library(data.table)
library(Rcpp)
sourceCpp("src/screen_SA.cpp")

# # Read initial seed from command line
args <- commandArgs(trailingOnly = F)
file_num <-  args[length(args)]
file_num <- as.numeric(unlist(strsplit(file_num,'-')))

seed_num<-file_num[1]
prog_speed<-file_num[2]
iter_num<-file_num[3]
print(seed_num)

seed <- as.data.frame(fread("calibrated_seeds.csv"))[seed_num,1]
print(seed)

# split data into cancer/non-cancer
full_data <- fread(paste0("Results/cancer/screening_data51_", seed, ".csv"))
# changing to run age 55-79 for the trial design
#full_data <- full_data[V2 > 69 & V2 < 75]
full_data <- full_data[V2 > 44 & V2 < 80]
#full_data <- full_data[V2 > 54 & V2 < 80]

infection_data <- as.matrix(full_data[V3 == 0])

cancer_data <- as.matrix(full_data[V3 == 1])

# stage simulation data
ajcc7 = as.matrix(fread("Data/AJCC7.by.age.csv"))

stage_progression <- as.matrix(fread("Data/stage.progression.csv"))

weights <- as.matrix(fread("Data/longitudinal_weights.csv"))
pop_size <- as.matrix(fread("Data/pop_size.csv"))
weights <- t(t(weights) %*% diag(as.vector(pop_size))) / 100000

mortality_data <- as.matrix(fread("Data/mortality32.csv"))

tnm <- fread("Data/AJCC7 TNM proportions_numeric.csv")
tnm <- split(tnm, by = "AJCC7\ stage", keep.by = F)
tnm <- lapply(tnm, function(x) as.matrix(x))

screen_scenarios <- fread("Data/screening_scenarios.csv")
screen_scenarios[, scenario := 1:nrow(screen_scenarios)]
screen_scenarios <- split(screen_scenarios, by = "scenario", keep.by = F)
screen_scenarios <- lapply(screen_scenarios, function(x) na.omit(as.numeric(x)))

cumu_survival <- fread("Data/cumulative.cause.specific.survival.by.smoking.csv")
cumu_survival <- split(cumu_survival, by = "V1", keep.by = F)
cumu_survival <- lapply(cumu_survival, function(x)
  lapply(split(x, by = "V2", keep.by = F), as.matrix))


no_infection_mortality <-
  as.matrix(fread(paste0("Results/cancer/mortality_denominator_51_", seed, "_1.csv")))

system.time(cpp_sim(ajcc7,
     stage_progression,
     cancer_data,
     infection_data,
     weights,
     mortality_data,
     tnm,
     prog_speed,                       # Progression speed (1-4)
     seed,
     screen_scenarios,
     iter_num,                       # number of iterations
     no_infection_mortality,
     cumu_survival,
     0,                           # sensitivity_1: redundant
     0,                            # sensitivity_2: change the sensitivity of ultrasound by a fixed percentage for all cancers
     1))                          # sensitivity_3: allow men to die earlier with screening than without
