# Load packages for mixed-effects models and spline functions
library(lme4)
library(splines)
library(optimx)

# input dataset df: 
# data frame with columns pseudonymised patient episode identifier,
# time from blood culture collection to clinical parameter measurements,
# clinical parameter values (CRP, WBC, HR, RR, Temp), age, sex,
# comorbidities (Elixhauser, Charlson), source of infection, community-onset status,
# immunosuppression status, baseline antimicrobial susceptibility, and pathogen groups.


# Box-Cox transformation for CRP and WBC
BCtrans <- MASS::boxcox(lm(df$CRP ~ 1))
lambda <- BCtrans$x[which.max(BCtrans$y)]
df$CRP <- ifelse(lambda == 0, log(df$CRP), (df$CRP^lambda - 1) / lambda)

BCtrans <- MASS::boxcox(lm(df$WBC ~ 1))
lambda <- BCtrans$x[which.max(BCtrans$y)]
df$WBC <- ifelse(lambda == 0, log(df$WBC), (df$WBC^lambda - 1) / lambda)

# Define a natural spline for time using identified knots and boundary knots
TimeSplineKnots <- quantile(df$time, c(0.20, 0.40, 0.60, 0.80))
TimeSplineBoundary <- quantile(df$time, c(0.01, 0.99))
TimeSpline <- ns(df$time, knots = TimeSplineKnots, Boundary.knots = TimeSplineBoundary)

# Define a natural spline for age using identified knots and boundary knots
AgeSplineKnots <- quantile(df$Age, c(0.50))
AgeSplineBoundary <- quantile(df$Age, c(0.05, 0.95))
AgeSpline <- ns(df$Age, knots = AgeSplineKnots, Boundary.knots = AgeSplineBoundary)

# Fit linear mixed-effects models with the defined natural splines
# C-reactive protein
# Main model
RE_CRP <- lmer(
  CRP ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID), # Random intercept and slope for time
  data = df,
  REML = TRUE)

# For source of infection
RE_CRP_Source <- lmer(
  CRP ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For baseline antimicrobial susceptibility
RE_CRP_Suscept <- lmer(
  CRP ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression + BaseSuscept) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# White blood cell count
# Main model
RE_WBC <- lmer(
  WBC ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For source of infection
RE_WBC_Source <- lmer(
  WBC ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For baseline antimicrobial susceptibility
RE_WBC_Suscept <- lmer(
  WBC ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression + BaseSuscept) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# Heart rate
# Main model
RE_HR <- lmer(
  HR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For source of infection
RE_HR_Source <- lmer(
  HR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For baseline antimicrobial susceptibility
RE_HR_Suscept <- lmer(
  HR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression + BaseSuscept) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# Respiratory rate
# Main model
RE_RR <- lmer(
  RR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For source of infection
RE_RR_Source <- lmer(
  RR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For baseline antimicrobial susceptibility
RE_RR_Suscept <- lmer(
  RR ~ Elixhauser + Charlson + Age + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression + BaseSuscept) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE
# Temperature
# Main model
RE_Temp <- lmer(
  Temp ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For source of infection
RE_Temp_Source <- lmer(
  Temp ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (Source + CommunityOnset + ImmunoSuppression) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE)

# For baseline antimicrobial susceptibility
RE_Temp_Suscept <- lmer(
  Temp ~ Elixhauser + Charlson + AgeSpline + Sex +
    TimeSpline * (BugGroup + Source + CommunityOnset + ImmunoSuppression + BaseSuscept) +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE,
  control = lmerControl(
    check.nobs.vs.nRE = "ignore", 
    optimizer = "optimx", 
    optCtrl = list(method = "nlminb")
  )