# Load packages for gamlss, mixed-effects models and spline functions
library(data.table)
library(gamlss)
library(lme4)
library(splines)

# df: dataset for standard responders
df_unique <- df[, .SD[1], by = EpisodeID]
df_unique <- df_unique[, .(EpisodeID, BugGroup, Source, Age, Sex,
                           CommunityOnset, ImmunoSuppression, Elixhauser, Charlson)]
# Bootstrap
set.seed(7)
df_resampled <- df_unique[sample(.N, 100000, replace = TRUE)]

# Replace duplicate EpisodeID with random unique numbers
duplicates <- duplicated(df_resampled$EpisodeID)
df_resampled$EpisodeID[duplicates] <- sample(900000:999999,
                                             size = sum(duplicates), replace = FALSE)

# Expand the number of rows per group and add random time
df_resampled <- df_resampled[rep(1:.N, each = 9)]
df_resampled[, time := runif(9, min = 0, max = 8), by = EpisodeID]

# Order by EpisodeID and time
df_resampled <- df_resampled[order(EpisodeID, time)]

# Fit linear mixed model with only fixed and random effects for time (natural cubic spline)
# Define the natural spline basis for 'time'
TimeSplineKnots <- quantile(df$time, c(0.20, 0.40, 0.60, 0.80))
TimeSplineBoundary <- quantile(df$time, c(0.01, 0.99))
TimeSpline <- ns(df$time, knots = TimeSplineKnots, Boundary.knots = TimeSplineBoundary)

RE_CRP <- lmer(
  CRP ~ TimeSpline +
    (TimeSpline | EpisodeID),
  data = df,
  REML = TRUE
)

# Simulation
set.seed(1)
sim <- simulate(RE_CRP, nsim = 1, re.form = NA, cond.sim = FALSE,
                     newdata = df_resampled, allow.new.levels = TRUE)
setDT(sim)
setnames(sim, old = c("sim_1"), new = c("sim"))
df_sim <- cbind(df_resampled, sim)
df_sim[, sim := ifelse(lambda == 0, exp(sim), (lambda * sim + 1)^(1 / lambda))] # invert the box-cox transformation

# Fit a GAMLSS model (LMS method) with the simulated data
log(nrow(df_sim))
LMS_CRP <- gamlss(sim ~ pb(time, method = "GAIC", k = log(nrow(df_sim))),
                        sigma.formula = ~ pb(time, method = "GAIC", k = log(nrow(df_sim))),
                        nu.formula = ~ pb(time, method = "GAIC", k = log(nrow(df_sim))),
                       family = BCTo, data = df_sim)
# Predict centiles
Centiles_CRP <- centiles.pred(LMS_CRP,
                              xname = "time",
                              xvalues = seq(0, 8, 0.1),
                              cent = c(5, 10, 25, 50, 75, 90, 95))
setDT(Centiles_CRP)

# Reshape data
Centiles_CRP <- melt(Centiles_CRP,
                     id.vars = c("x"),
                     measure.vars = c("5", "10", "25", "50", "75", "90", "95"),
                     variable.name = "Centiles",
                     value.name = "CRP")

Centiles_CRP$Centiles <- factor(Centiles_CRP$Centiles, levels = c("5", "10", "25", "50", "75", "90", "95"))

# Plot centile reference chart
p_centiles_crp <- 
ggplot(Centiles_CRP, aes(x = x, y = CRP)) +
  geom_line(aes(color = Centiles, linetype = Centiles), linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 8, 1)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 600, 50)) +
  scale_colour_manual(values = c("red", "black", "red", "black", "red", "black", "red")) +
  scale_linetype_manual(values = c("dotted", "dotdash", "dashed", "solid", "dashed", "dotdash", "dotted")) +
  labs(x = "Time (days after blood culture collection)",
       y = "C-Reactive Protein (mg/L)") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.text.align = 0,
        legend.background = element_rect(fill = "transparent")) +
  guides(color = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE))

print(p_centiles_crp)