library(here)
library(tidyverse)

data <- read_csv(here("data","prepped_with_augmented_ma_extended","mb_ma_combined_prepped_0.125.csv"))

ggplot(aes(y = method, x = mean_age, color = test_lang), data = data) +
  labs(y = "Method", x = "Infant Age in Days") +
  scale_colour_manual(breaks=c("b.nonnative","a.native","c.artificial"),labels=c("Non-native","Native","Other"),values = c("#046C9A", "#C93312", "#a6761d")) +
  ggtitle('Distribution of Data across Different Conditions') +
  geom_point(aes(color = test_lang), alpha = 0.5, size = 3, position = position_jitterdodge(jitter.height = 0.25,jitter.width=0.2,dodge.width=0.6), data = data) +
  facet_grid(~sourcePretty) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=18),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(shape = "Speech Type", color = "Test Language") +
  scale_y_discrete(breaks=c("a.cf","b.hpp","c.other"),labels=c("CF","HPP","Other"))+
  theme(legend.position=c(0.8,0.8),
        strip.text = element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12))
  

ggsave(here("results_from_R","results_with_augmented_ma_extended","DataDifferentConditions.pdf"),width=9,height=6)
ggsave(here("Overleaf","augmented_ma_extended","DataDifferentConditions.pdf"),width=9,height=6)
