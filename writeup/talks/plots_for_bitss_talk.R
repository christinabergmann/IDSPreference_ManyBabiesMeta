# Forest plots for BITSS talks 2021 (M. Lewis)
# https://docs.google.com/presentation/d/1CFBuytkxks9l4UytFjlqyWaXuNNY4H-LbcRyzlQTFX0/edit#slide=id.gb3889ff16e_0_59

library(tidyverse)

# Es estimates from dunst et al and manybabies 2020 (data from original_es)
A1 <- "dark"
A2 <- "light"
es_estimates <- data.frame(type = c("MA", "MLR", "MA", "MLR", "MA", "MLR", "MA", "MLR"),
                          source = c("original", "original", "re-analysis", "re-analysis",
                                     "naive", "naive", "moderated", "moderated"),
                          es_estimate = c(.67,.35, .70, .35, .68, .35, .61, .13),
                          es_upper = c(.76,.42, 1.02, .43, 1, .42, .96, .35 ),
                          es_lower = c(.57,.29, .37, .28, .35, .28, .26, -.08),
                          alpha = c(A1, A1,  A2, A2, A2, A2, A2, A2))

# original + naive + moderated
es_estimates %>%
  filter(source != "re-analysis") %>%
  mutate(type = fct_recode(type,"Meta-analysis" = "MA", "Multi-Lab Replication" = "MLR"),
         source = fct_relevel(source, "original", "naive", "moderated")) %>%
  ggplot(aes(x = es_estimate, y = fct_rev(type), color = type, shape = source)) +
  geom_pointrange(aes(xmin = es_lower, xmax = es_upper, alpha = alpha),
                  size = 1.1, position = position_dodge2(0.5, reverse = T)) +
  xlim(-.10, 1) +
  #scale_alpha_manual(guide = "none") +
  geom_vline(aes(xintercept = 0), linetype = 2)+
  xlab(expression(paste("Effect Size (Cohen's ", italic(d), ")"))) +
  scale_color_manual(values = c("blue", "red"), guide = "none") +
  #scale_shape_manual(name = "model type") +
  scale_alpha_manual(values = c(1, .6), guide = FALSE) +
  ylab("") +
  ggtitle("Aggregate Effect Size Estimates") +
  labs(shape = "Model ") +
  theme_classic(base_size = 16) +
  theme(# axis.text.y = element_text(angle = 90, hjust = .5),
         legend.position = "bottom")

# original only
es_estimates %>%
  filter(source != "re-analysis") %>%
  mutate(type = fct_recode(type,"Meta-analysis" = "MA", "Multi-Lab Replication" = "MLR"),
         source = fct_relevel(source, "original", "naive", "moderated")) %>%
  ggplot(aes(x = es_estimate, y = fct_rev(type), color = type, shape = source)) +
  geom_pointrange(aes(xmin = es_lower, xmax = es_upper, alpha = alpha),
                  size = 1.1, position = position_dodge2(0.5, reverse = T)) +
  xlim(-.10, 1) +
  #scale_alpha_manual(guide = "none") +
  geom_vline(aes(xintercept = 0), linetype = 2)+
  xlab(expression(paste("Effect Size (Cohen's ", italic(d), ")"))) +
  scale_color_manual(values = c("blue", "red"), guide = "none") +
  #scale_shape_manual(name = "model type") +
  scale_alpha_manual(values = c(1, 0), guide = FALSE) +
  ylab("") +
  ggtitle("Aggregate Effect Size Estimates\nfor IDS preference") +
  labs(shape = "Model ") +
  theme_classic(base_size = 16) +
  theme(# axis.text.y = element_text(angle = 90, hjust = .5),
    legend.position = "bottom")

# original + naive
es_estimates %>%
  mutate(es_lower = case_when(es_lower < 0 ~ .1,
                              TRUE ~ es_lower)) %>%
  filter(source != "re-analysis") %>%
  mutate(type = fct_recode(type,"Meta-analysis" = "MA", "Multi-Lab Replication" = "MLR"),
         source = fct_relevel(source, "original", "naive", "moderated")) %>%
  ggplot(aes(x = es_estimate, y = fct_rev(type), color = type, shape = source)) +
  geom_pointrange(aes(xmin = es_lower, xmax = es_upper, alpha = alpha),
                  size = 1.1, position = position_dodge2(0.5, reverse = T)) +
  xlim(-.10, 1) +
  #scale_alpha_manual(guide = "none") +
  geom_vline(aes(xintercept = 0), linetype = 2)+
  xlab(expression(paste("Effect Size (Cohen's ", italic(d), ")"))) +
  scale_color_manual(values = c("blue", "red"), guide = "none") +
  #scale_shape_manual(name = "model type") +
  scale_alpha_manual(values = c(1, .6), guide = FALSE) +
  ylab("") +
  ggtitle("Aggregate Effect Size Estimates\nfor IDS preference") +
  labs(shape = "Model ") +
  theme_classic(base_size = 16) +
  theme(# axis.text.y = element_text(angle = 90, hjust = .5),
    legend.position = "bottom")


# manybabies exclusion criteria plot (from table 6 in manybabies plot)

exclusion_criteria = data.frame(method = c("CF", "CF", "CF", "HPP", "HPP", "HPP"),
                                min_trials = c(2,4,8, 2,4,8),
                                estimate = c(.29, .34, .40, .51, .56, .63),
                                se = c(.06, .06, .06, .06, .06, .07))

exclusion_criteria_with_ci <- exclusion_criteria %>%
  mutate(ci_lower = estimate - (1.96*se),
         ci_upper = estimate + (1.96*se))

ggplot(exclusion_criteria_with_ci, aes(x = estimate, y = min_trials, color = method)) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 1.1, position = position_dodge2(0.5, reverse = T)) +
  geom_line() +
  ggtitle("MLR effect size varying inclusion criteria") +
  xlab(expression(paste("Effect Size (Cohen's ", italic(d), ")"))) +
  ylab("Minimum number of trials for inclusion") +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_classic(base_size = 16)
