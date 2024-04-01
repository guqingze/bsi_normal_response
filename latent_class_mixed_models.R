# Load the necessary packages for latent class mixed models and spline functions
library(lcmm)
library(splines)

# Define the natural spline basis for 'time'
# using the 25th, 50th, and 75th percentiles as knots and
# 1st and 99th percentiles as boundary knots
timeSplineKnots <- quantile(df$time, c(0.25, 0.5, 0.75))
timeSplineBoundary <- quantile(df$time, c(0.01, 0.99))
timeSpline <- ns(df$time, knots = timeSplineKnots, Boundary.knots = timeSplineBoundary)

# Renaming columns for interpretability
timeSplineColNames <- paste0("T", 1:ncol(timeSpline))
colnames(timeSpline) <- timeSplineColNames

# Binding the natural spline terms to the original dataframe
df <- cbind(df, timeSpline)

# Fitting a series of latent class mixed models with varying number of classes (ng)
# Beginning with a standard linear mixed model (LCMM_CRP_1)
LCMM_CRP_1 <- lcmm(
  CRP_trans ~ T1 + T2 + T3 + T4,
  random = ~ T1 + T2 + T3 + T4,
  subject = "EpisodeID",
  nproc = 20, # Number of cores to use for parallel processing
  data = df
)

# Fitting latent class mixed models with multiple classes
# using LCMM_CRP_1 as the starting values (B = LCMM_CRP_1)
LCMM_CRP_2 <- lcmm(
  CRP_trans ~ T1 + T2 + T3 + T4,
  mixture = ~ T1 + T2 + T3 + T4,
  random = ~ T1 + T2 + T3 + T4, maxiter = 10000,
  subject = "EpisodeID", ng = 2, B = LCMM_CRP_1,
  nproc = 20, # Number of cores to use for parallel processing
  data = df
)

# The process continues, increasing ng for each model to test for better fitting class numbers
# The `update` function is used to change the number of classes while keeping the rest of the model the same.
LCMM_CRP_3 <- update(LCMM_CRP_2, ng = 3)
LCMM_CRP_4 <- update(LCMM_CRP_3, ng = 4)
LCMM_CRP_5 <- update(LCMM_CRP_4, ng = 5)
LCMM_CRP_6 <- update(LCMM_CRP_5, ng = 6)