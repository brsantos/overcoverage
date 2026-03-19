library(overcoverage2)

# Provide input paths via environment variables to avoid hard-coding confidential data.
data_sameyear <- Sys.getenv("OC2_DATA_SAMEYEAR")
data_sameyear_v2 <- Sys.getenv("OC2_DATA_SAMEYEAR_V2")

if (data_sameyear == "" || data_sameyear_v2 == "") {
  stop("Set OC2_DATA_SAMEYEAR and OC2_DATA_SAMEYEAR_V2 to the input .rds paths.")
}

data_model <- readRDS(data_sameyear)
data_model <- data_model[data_model$year != 2021, ]

data_model_v2 <- readRDS(data_sameyear_v2)
data_model_v2 <- data_model_v2[data_model_v2$year == 2021, ]

data_all <- rbind(data_model, data_model_v2)

data_all <- data_all[, c(
  "id", "year", "rtb", "age", "cob", "female", "death",
  "firstimmig", "immig", "emig", "reimmig",
  "married", "divorced", "amf", "studies",
  "intmove", "child", "swecit", "group_countries",
  "age_cat", "time_sweden", "time_sweden_cat",
  "faminc_b", "indinc_b", "aldpens", "forvers", "socink"
)]

year_beginning <- 2002
final_year <- 2022

# Apply consistency checks
checked <- oc2_check_register_data(
  data_all,
  year_beginning = year_beginning,
  final_year = final_year
)
data_all <- checked$data

# Create binary registers
data_all <- oc2_prepare_register_data(
  data_all,
  register_cols = c("married", "divorced", "amf", "studies", "intmove",
                    "child", "pension", "job", "social", "faminc_b"),
  binary_rules = list(
    pension = list(source = "aldpens", threshold = 0),
    job = list(source = "forvers", threshold = 0),
    social = list(source = "socink", threshold = 0)
  )
)

# Age dummy covariates
data_all$age1 <- ifelse(data_all$age_cat == "[18, 35]", 1, 0)
data_all$age2 <- ifelse(data_all$age_cat == "[36, 60]", 1, 0)
data_all$age3 <- ifelse(data_all$age_cat == "[60, ...]", 1, 0)

ages <- data_all[, c("id", "year", "firstimmig", "age", "age1", "age2", "age3")]
ages <- ages[order(ages$id, ages$year), ]
ages$lastyear <- ave(ages$year, ages$id, FUN = function(x) {
  if (max(x) == final_year) final_year else max(x)
})

unique_id <- unique(ages$id)
age_covariate <- array(data = NA, dim = c(length(unique_id), 7, final_year - year_beginning))
age_covariate[, 1, ] <- unique_id
for (yr in 1:20) {
  ages_year <- ages[ages$year == 2002 + yr, ]
  j <- match(ages_year$id, age_covariate[, 1, 1])
  age_covariate[j, 2, ] <- ages_year$firstimmig
  age_covariate[j, 3, ] <- ages_year$lastyear
  age_covariate[j, 4, yr] <- ages_year$age
  age_covariate[j, 5, yr] <- ages_year$age1
  age_covariate[j, 6, yr] <- ages_year$age2
  age_covariate[j, 7, yr] <- ages_year$age3
}

for (b in 1:20) {
  j <- which(is.na(age_covariate[, 4, b]) & (age_covariate[, 2, b] < (2002 + b)))
  age_covariate[j, 4, b] <- age_covariate[j, 4, b - 1] + 1
  x <- which(age_covariate[j, 4, b] <= 35)
  if (length(x) != 0) {
    if (length(j) == 1) {
      age_covariate[j, , ][5, b] <- 1
      age_covariate[j, , ][6:7, b] <- 0
    } else {
      age_covariate[j, , ][x, 5, b] <- 1
      age_covariate[j, , ][x, 6:7, b] <- 0
    }
  }
  y <- which(age_covariate[j, 4, b] > 35 & age_covariate[j, 4, b] <= 60)
  if (length(y) != 0) {
    if (length(j) == 1) {
      age_covariate[j, , ][5, b] <- 0
      age_covariate[j, , ][6, b] <- 1
      age_covariate[j, , ][7, b] <- 0
    } else {
      age_covariate[j, , ][y, 5, b] <- 0
      age_covariate[j, , ][y, 6, b] <- 1
      age_covariate[j, , ][y, 7, b] <- 0
    }
  }
  z <- which(age_covariate[j, 4, b] > 60)
  if (length(z) != 0) {
    if (length(j) == 1) {
      age_covariate[j, , ][5:6, b] <- 0
      age_covariate[j, , ][7, b] <- 1
    } else {
      age_covariate[j, , ][z, 5:6, b] <- 0
      age_covariate[j, , ][z, 7, b] <- 1
    }
  }
}

age_covariate <- age_covariate[, c(-1, -2, -3, -4), ]
age_covariate <- age_covariate[, c(-1), ]

# Time in Sweden dummy covariates
data_all$tis1 <- ifelse(data_all$time_sweden_cat == "[0, 5]", 1, 0)
data_all$tis2 <- ifelse(data_all$time_sweden_cat == "[6, 10]", 1, 0)
data_all$tis3 <- ifelse(data_all$time_sweden_cat == "10+", 1, 0)

sweden_time <- data_all[, c("id", "year", "firstimmig", "time_sweden", "tis1", "tis2", "tis3")]
sweden_time <- sweden_time[order(sweden_time$id, sweden_time$year), ]
sweden_time$lastyear <- ave(sweden_time$year, sweden_time$id, FUN = function(x) {
  if (max(x) == final_year) final_year else max(x)
})

unique_id <- unique(sweden_time$id)
tis_covariate <- array(data = NA, dim = c(length(unique_id), 7, 20))
tis_covariate[, 1, ] <- unique_id
for (yr in 1:20) {
  tis_year <- sweden_time[sweden_time$year == 2002 + yr, ]
  j <- match(tis_year$id, tis_covariate[, 1, 1])
  tis_covariate[j, 2, ] <- tis_year$firstimmig
  tis_covariate[j, 3, ] <- tis_year$lastyear
  tis_covariate[j, 4, yr] <- tis_year$time_sweden
  tis_covariate[j, 5, yr] <- tis_year$tis1
  tis_covariate[j, 6, yr] <- tis_year$tis2
  tis_covariate[j, 7, yr] <- tis_year$tis3
}

for (b in 1:20) {
  j <- which(is.na(tis_covariate[, 4, b]) & (tis_covariate[, 2, b] < (2002 + b)))
  tis_covariate[j, 4, b] <- tis_covariate[j, 4, b - 1] + 1
  x <- which(tis_covariate[j, 4, b] == 0)
  if (length(x) != 0) {
    if (length(j) == 1) {
      tis_covariate[j, , ][5, b] <- 1
      tis_covariate[j, , ][6:7, b] <- 0
    } else {
      tis_covariate[j, , ][x, 5, b] <- 1
      tis_covariate[j, , ][x, 6:7, b] <- 0
    }
  }
  y <- which(tis_covariate[j, 4, b] > 0 & tis_covariate[j, 4, b] <= 5)
  if (length(y) != 0) {
    if (length(j) == 1) {
      tis_covariate[j, , ][5, b] <- 0
      tis_covariate[j, , ][6, b] <- 1
      tis_covariate[j, , ][7, b] <- 0
    } else {
      tis_covariate[j, , ][y, 5, b] <- 0
      tis_covariate[j, , ][y, 6, b] <- 1
      tis_covariate[j, , ][y, 7, b] <- 0
    }
  }
  z <- which(tis_covariate[j, 4, b] > 5)
  if (length(z) != 0) {
    if (length(j) == 1) {
      tis_covariate[j, , ][5:6, b] <- 0
      tis_covariate[j, , ][7, b] <- 1
    } else {
      tis_covariate[j, , ][z, 5:6, b] <- 0
      tis_covariate[j, , ][z, 7, b] <- 1
    }
  }
}

tis_covariate <- tis_covariate[, c(-1, -2, -3, -4), ]
tis_covariate <- tis_covariate[, c(-1), ]

# Time independent covariates (sex + country of birth groups)
data_all$sex <- ifelse(data_all$female == "female", 1, 0)
data_all$cob1 <- ifelse(data_all$group_countries == "Denmark/Norway", 1, 0)
data_all$cob2 <- ifelse(data_all$group_countries == "Eastern Europe", 1, 0)
data_all$cob3 <- ifelse(data_all$group_countries == "Iceland/Finland", 1, 0)
data_all$cob4 <- ifelse(data_all$group_countries == "MENA", 1, 0)
data_all$cob5 <- ifelse(data_all$group_countries == "USA/Canada/Oceania", 1, 0)
data_all$cob6 <- ifelse(data_all$group_countries == "Western Europe", 1, 0)
data_all$cob7 <- ifelse(data_all$group_countries == "World", 1, 0)

covariates_to_use <- data_all[, c("id", "sex", "cob1", "cob2", "cob3", "cob4", "cob5", "cob6", "cob7")]
covariates_to_use <- covariates_to_use[!duplicated(covariates_to_use), ]

covariates <- cbind(rep(1, nrow(covariates_to_use)), covariates_to_use[, 2:9])
all_covariates <- covariates[, c(-3)]

# Observation combinations
register_cols <- c("married", "divorced", "amf", "studies", "intmove",
                   "child", "pension", "job", "social", "faminc_b")

X <- expand.grid(rep(list(c(1, 0)), length(register_cols) + 1))
X <- as.matrix(X)

n_age <- 2
dummies <- matrix(0, nrow = nrow(X), ncol = n_age)
for (i in 1:n_age) {
  mat <- matrix(0, nrow = nrow(X), ncol = n_age)
  mat[, i] <- 1
  dummies <- rbind(dummies, mat)
}
X_extended <- X[rep(1:nrow(X), n_age + 1), ]
X <- cbind(X_extended, dummies)

n_tis <- 2
dummies <- matrix(0, nrow = nrow(X), ncol = n_tis)
for (i in 1:n_tis) {
  mat <- matrix(0, nrow = nrow(X), ncol = n_tis)
  mat[, i] <- 1
  dummies <- rbind(dummies, mat)
}
X_extended <- X[rep(1:nrow(X), n_tis + 1), ]
X <- cbind(X_extended, dummies)

X <- as.data.frame(X)
colnames(X) <- c("married", "divorced", "amf", "studies", "intmove",
                 "child", "pension", "job", "social", "faminc_b",
                 "sex", "age2", "age3", "tis2", "tis3")
X$row_id <- seq_len(nrow(X))

X0 <- X
X0$reimmig <- 0
X0$y <- X0$row_id
X1 <- X
X1$reimmig <- 1
X1$y <- X1$row_id + (nrow(X) + 2)

X_all <- rbind(X0, X1)

data_all <- merge(
  data_all,
  X_all,
  by = c("married", "divorced", "amf", "studies", "intmove",
         "child", "pension", "job", "social", "faminc_b",
         "sex", "age2", "age3", "tis2", "tis3", "reimmig"),
  all.x = TRUE
)

data_all$y[data_all$death == 1] <- nrow(X) + 1
data_all$y[data_all$emig == 1] <- nrow(X) + 2

# Wide format
data_to_use <- reshape(
  data_all[, c("id", "year", "y")],
  idvar = "id",
  timevar = "year",
  direction = "wide"
)

# Replace NAs after observations with "unobserved"
removeNA <- aggregate(year ~ id, data_all, function(x) min(x))
names(removeNA) <- c("id", "first_year")
data_to_use <- merge(data_to_use, removeNA, by = "id", all.x = TRUE)

unobs <- which(apply(X[, 1:10], 1, function(r) all(r == 0)))
cov_comb <- X[unobs, 11:15]

combins <- matrix(NA, nrow = nrow(data_to_use), ncol = 20)
for (j in 1:nrow(data_to_use)) {
  for (t in (data_to_use$first_year[j] - 2002):20) {
    cov <- c(all_covariates[j, 2], age_covariate[j, , t], tis_covariate[j, , t])
    combins[j, t] <- which(apply(cov_comb, 1, function(row) all(row == cov)))
  }
}

prov_y <- as.matrix(data_to_use)
for (i in 1:nrow(data_to_use)) {
  unobs_idx <- which(is.na(prov_y[i, (data_to_use$first_year[i] - year_beginning + 2):(ncol(prov_y) - 1)]))
  unobs_idx <- unobs_idx + data_to_use$first_year[i] - year_beginning + 1
  prov_y[i, unobs_idx] <- combins[i, unobs_idx - 1] * 2^10
}

y_matrix <- as.matrix(prov_y[, 2:(ncol(prov_y) - 1)])

out_dir <- Sys.getenv("OC2_PREP_OUTPUT", unset = "inst/tests/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(y_matrix, file.path(out_dir, "y.rds"))
saveRDS(all_covariates, file.path(out_dir, "all_covariates.rds"))
saveRDS(age_covariate, file.path(out_dir, "age_covariate.rds"))
saveRDS(tis_covariate, file.path(out_dir, "tis_covariate.rds"))
saveRDS(combins, file.path(out_dir, "combins.rds"))
