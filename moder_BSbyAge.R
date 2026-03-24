# =========================================================
# Moderation analysis:
# Does the relation between brain scores and behavior differ by age group?
# Model for each behavioral variable:
#   behavior ~ brain_score * age_group
# =========================================================
setwd("~/OneDrive/Documents/Entropy_in_aging")
# ---------- packages ----------
library(tidyverse)
library(broom)
library(emmeans)
library(readxl)

# ---------- load data ----------
# Replace with your actual file path
df <- read_excel("EntropyAging_MultReg.xlsx")
names(df) <- gsub(" ", "_", names(df))
names(df) <- gsub("-", "_", names(df))
df$brainscore_blink <- as.numeric(gsub(",", ".", df$brainscore_blink))
df$Age_Group <- factor(df$Age_Group,
                       levels = c(-1, 1),
                       labels = c("young", "old"))

# ---------- define variables ----------
# Replace these names with your actual column names
brain_var <- "brainscore_blink"
group_var <- "Age_Group"

# List your 11 behavioral variables here
behavior_vars <- c(
  "d2_KL",
  "VLMT_Dg1_5",
  "VLMT_Dg1",
  "VLMT_Dg5",
  "VLMT_I",
  "VLMT_Dg7",
  "WMT_2",
  "Zahlen_ges",
  "Zahlen_vor",
  "Zahlen_rueck",
  "MWT_B"
)

# "VLMT W",
# "VLMT W-F",

# ---------- recode / prepare ----------
# Make sure age group is a factor and set reference level if needed
df[[group_var]] <- factor(df[[group_var]])

# Optional: set "younger" as reference if those are your labels
# df[[group_var]] <- relevel(df[[group_var]], ref = "younger")

# Optional but recommended: z-score brain score
df[[brain_var]] <- scale(df[[brain_var]], center = TRUE, scale = TRUE)[, 1]

# Optional: z-score all behavioral outcomes too
# This makes coefficients more comparable across measures
for (v in behavior_vars) {
  df[[v]] <- scale(df[[v]], center = TRUE, scale = TRUE)[, 1]
}

# ---------- function to run one moderation model ----------
run_moderation <- function(outcome, data, brain_var, group_var) {
  
  # formula: outcome ~ brain_score * age_group
  form <- as.formula(
    paste(outcome, "~", brain_var, "*", group_var)
  )
  
  model <- lm(form, data = data)
  
  # tidy coefficients
  coef_tab <- broom::tidy(model) %>%
    mutate(outcome = outcome)
  
  # model fit
  fit_tab <- broom::glance(model) %>%
    mutate(outcome = outcome)
  
  # simple slopes: effect of brain score within each age group
  slopes <- emtrends(model, specs = group_var, var = brain_var) %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(outcome = outcome)
  
  # estimated marginal means at +/- 1 SD for plotting interpretation
  # since brain score is z-scored, -1 and +1 are convenient values
  emm <- emmeans(model, specs = group_var, at = setNames(list(c(-1, 1)), brain_var)) %>%
    as.data.frame() %>%
    mutate(outcome = outcome)
  
  list(
    model = model,
    coef_tab = coef_tab,
    fit_tab = fit_tab,
    slopes = slopes,
    emm = emm
  )
}

# ---------- run all models ----------
results <- lapply(behavior_vars, run_moderation,
                  data = df,
                  brain_var = brain_var,
                  group_var = group_var)

names(results) <- behavior_vars

# ---------- combine outputs ----------
coef_all <- bind_rows(lapply(results, `[[`, "coef_tab"))
fit_all  <- bind_rows(lapply(results, `[[`, "fit_tab"))
slopes_all <- bind_rows(lapply(results, `[[`, "slopes"))
emm_all <- bind_rows(lapply(results, `[[`, "emm"))

# ---------- extract interaction terms ----------
interaction_results <- coef_all %>%
  filter(str_detect(term, paste0("^", brain_var, ":|:", brain_var, "$"))) %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(p_adj_fdr)

print(interaction_results)

# ---------- save results ----------
write.csv(coef_all, "moderation_all_coefficients.csv", row.names = FALSE)
write.csv(fit_all, "moderation_model_fit.csv", row.names = FALSE)
write.csv(slopes_all, "moderation_simple_slopes.csv", row.names = FALSE)
write.csv(emm_all, "moderation_emmeans.csv", row.names = FALSE)
write.csv(interaction_results, "moderation_interactions_only.csv", row.names = FALSE)

# ---------- show compact summary ----------
cat("\n============================\n")
cat("Interaction effects summary\n")
cat("============================\n")

interaction_results %>%
  select(outcome, term, estimate, std.error, statistic, p.value, p_adj_fdr) %>%
  print(n = Inf)

# =========================================================
# Optional: plot interaction for one selected behavioral measure
# =========================================================

plot_interaction <- function(outcome_name, data, brain_var, group_var) {
  graphics.off()
  
  p <- ggplot(
    data,
    aes(
      x = .data[[brain_var]],
      y = .data[[outcome_name]],
      color = .data[[group_var]]
    )
  ) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic()
  
  print(p)
  invisible(p)
}

# Example:
p <- plot_interaction("VLMT_Dg1_5", df, brain_var, group_var)


cor(df[df$Age_Group=="young", c("brainscore_blink", behavior_vars)])
cor(df[df$Age_Group=="old", c("brainscore_blink", behavior_vars)])

