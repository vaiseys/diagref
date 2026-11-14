# Build the `gss2024` example dataset from the full GSS 2024 file.
#
# Input:  data-raw/gss2024.rds, the GSS 2024 cross-section (3,309 respondents,
#         813 variables, haven-labelled).
# Output: data/gss2024.rda, a small subset for examples and teaching.
#
# Run from the package root:  source("data-raw/gss2024.R")

full <- readRDS("data-raw/gss2024.rds")
num  <- function(v) as.numeric(haven::zap_labels(full[[v]]))
int  <- function(v) as.integer(num(v))

madeg <- int("madeg")
padeg <- int("padeg")

gss2024 <- tibble::tibble(
  id       = int("id"),
  # Dyadic variables (0-4, same coding)
  degree   = int("degree"),
  madeg    = madeg,
  padeg    = padeg,
  pardeg   = pmax(madeg, padeg, na.rm = TRUE),
  # Covariates
  female   = as.integer(int("sex") == 2L),
  age      = int("age"),
  race     = int("race"),
  # Outcomes
  polviews = int("polviews"),
  partyid  = ifelse(int("partyid") == 7L, NA_integer_, int("partyid")),
  fefam    = int("fefam"),
  attend   = int("attend"),
  happy    = int("happy"),
  health   = int("health"),
  childs   = int("childs"),
  wordsum  = int("wordsum"),
  coninc   = num("coninc"),
  prestg10 = int("prestg10"),
  # Survey weight
  wtssps   = num("wtssps")
)

# Keep respondents with both dyadic variables observed
gss2024 <- gss2024[!is.na(gss2024$degree) & !is.na(gss2024$pardeg), ]

usethis::use_data(gss2024, overwrite = TRUE, compress = "xz")
