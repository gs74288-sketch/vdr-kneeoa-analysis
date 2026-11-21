# NHANES III: E-value, non-linearity, and adiposity stratification
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(splines))
options(survey.lonely.psu = "adjust")

csv <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/clinical_nhanes3_20261011/nhanes3_knee_vitD_design.csv"
df0 <- read.csv(csv, stringsAsFactors = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
df <- df0
df$vitD <- num(df$vitD_nmol)
df$kll <- num(df$KL_left); df$klr <- num(df$KL_right)
df$age <- num(df$age); df$bmi <- num(df$bmi); df$month <- num(df$exam_month)
df$wt <- num(df$wt); df$psu <- as.integer(df$psu); df$stratum <- as.integer(df$stratum)
df <- df[!is.na(df$vitD) & df$kll %in% 0:4 & df$klr %in% 0:4, ]
df$bmi[df$bmi < 12 | df$bmi > 70] <- NA
df$klmax <- pmax(df$kll, df$klr)
df$roa <- as.integer(df$klmax >= 2)
df$female <- as.integer(df$sex == 2)
df$smoke <- ifelse(df$smoke100 == 1, 1L, ifelse(df$smoke100 == 2, 0L, NA_integer_))
df$season <- ifelse(df$month %in% c(12,1,2), "winter",
             ifelse(df$month %in% c(3,4,5), "spring",
             ifelse(df$month %in% c(6,7,8), "summer", "autumn")))
df$vitD25 <- df$vitD/25
df$race <- factor(df$race); df$season <- factor(df$season)

d <- df[complete.cases(df[, c("roa","vitD25","age","female","race","bmi","season","smoke","wt","psu","stratum")]), ]
cat("analysis N:", nrow(d), "\n")

des <- svydesign(ids = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = d)
p0 <- as.numeric(coef(svymean(~roa, des)))
cat(sprintf("weighted prevalence: %.3f | design df: %.0f\n", p0, degf(des)))

fit1 <- function(f, design = des) {
  m <- svyglm(f, design = design, family = quasibinomial())
  b <- unname(coef(m)["vitD25"]); se <- unname(SE(m)["vitD25"])
  dfr <- max(1, degf(design)); tq <- qt(0.975, dfr)
  c(OR = exp(b), lo = exp(b - tq*se), hi = exp(b + tq*se),
    p = 2*pt(-abs(b/se), dfr), b = b, se = se)
}

## ---------- 1. E-value ----------
evalue_or <- function(or, p0) {
  rr <- if (or < 1) or / (1 - p0 + p0*or) else or
  if (rr < 1) { rr <- 1/rr }
  rr + sqrt(rr * (rr - 1))
}
crude <- fit1(roa ~ vitD25)
adj <- fit1(roa ~ vitD25 + age + female + race + bmi + season + smoke)
cat("\n=== E-value ===\n")
cat(sprintf("crude     OR %.3f (%.3f-%.3f)\n", crude["OR"], crude["lo"], crude["hi"]))
cat(sprintf("adjusted  OR %.3f (%.3f-%.3f)\n", adj["OR"], adj["lo"], adj["hi"]))
cat(sprintf("E-value, crude point estimate      : %.2f\n", evalue_or(crude["OR"], p0)))
cat(sprintf("E-value, crude CI bound near null  : %.2f\n", evalue_or(crude["hi"], p0)))
cat(sprintf("E-value, adjusted point estimate   : %.2f\n", evalue_or(adj["OR"], p0)))
cat(sprintf("(OR->RR conversion used p0 = %.3f)\n", p0))

## ---------- 2. non-linearity ----------
cat("\n=== non-linearity (natural cubic spline, 3 df) ===\n")
m_lin <- svyglm(roa ~ vitD25 + age + female + race + bmi + season + smoke, design = des, family = quasibinomial())
m_spl <- svyglm(roa ~ ns(vitD, df = 3) + age + female + race + bmi + season + smoke, design = des, family = quasibinomial())
a <- anova(m_lin, m_spl, method = "Wald")
print(a)
cat(sprintf("test for non-linearity: F = %.3f, df = %d, P = %.3f\n",
            a$F[2], a$ddf[2] - a$ddf[1], a$`Pr(>F)`[2]))

## ---------- 3. adiposity stratification ----------
cat("\n=== by adiposity (fully adjusted, BMI removed from the model) ===\n")
for (g in list(list("BMI < 25", 0, 25), list("BMI 25-29.9", 25, 30), list("BMI >= 30", 30, 200))) {
  sub <- d[d$bmi >= g[[2]] & d$bmi < g[[3]], ]
  if (nrow(sub) < 100) next
  ds <- svydesign(ids = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = sub)
  m <- try(svyglm(roa ~ vitD25 + age + female + race + season + smoke, design = ds,
                  family = quasibinomial()), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-11s n=%4d  model failed\n", g[[1]], nrow(sub))); next }
  b <- unname(coef(m)["vitD25"]); se <- unname(SE(m)["vitD25"])
  dfr <- max(1, degf(ds)); tq <- qt(0.975, dfr)
  p <- 2*pt(-abs(b/se), dfr)
  cat(sprintf("  %-11s n=%4d  OR %6.3f (%.3f-%.3f)  P = %.3f\n",
              g[[1]], nrow(sub), exp(b), exp(b - tq*se), exp(b + tq*se), p))
}
