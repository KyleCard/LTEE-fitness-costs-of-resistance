setwd("C:\\Users\\cardk\\Box Sync\\Manuscripts\\Fitness costs of ABR mutations in the LTEE\\Card et al. 2020 - Fitness costs of ABR mutations in the LTEE")

library(tidyverse)
library(agricolae)
library(ggpubr)
library(cowplot)

counts_df <- read_csv("competition_data.csv")
counts_df$paired_ID <- as.character(counts_df$paired_ID)
strain_MICs_df <- read_csv("TET_strains_MICs.csv") # MIC data from Card et al. 2019


## Data wrangling ##
fitness_df <- counts_df %>% 
  mutate(malth_comp1 = log(comp1_d3 * (100^3) / comp1_d0)) %>%
  mutate(malth_comp2 = log(comp2_d3 * (100^3) / comp2_d0)) %>% 
  mutate(relative_fitness = malth_comp1 / malth_comp2) %>% 
  select(!c(antibiotic, competitor_2, comp1_d0:comp2_d3, malth_comp1, malth_comp2))

CompareFitness <- function(dat) {
  fitness_vec <- c()
  
  for (i in unique(dat$block)) {
    block_frame <- dat %>% 
      filter(block == i)
    
    for (j in unique(block_frame$paired_ID)) {
      paired_frame <- block_frame %>% 
        filter(paired_ID == j)
      
      true_fitness <- paired_frame[1, 6] / paired_frame[2, 6]

      fitness_vec <- bind_rows(fitness_vec, true_fitness)
    }
  }
  return(fitness_vec)
}

# Normalizes relative fitness based upon common competitor values
normalized_fitness_col <- CompareFitness(fitness_df) %>% 
  drop_na %>% 
  transmute(ln_relative_fitness = log(relative_fitness)) # Log(e) transform the relative fitness values

fitness_df <- fitness_df %>%
  drop_na %>%
  bind_cols(., normalized_fitness_col) %>% 
  select(!paired_ID:relative_fitness) %>% 
  rename(strain = competitor_1)

fitness_df$background <- as.factor(fitness_df$background)
fitness_df$strain <- as.factor(fitness_df$strain)

# Create data frame with average fitness of each strain
avg_fitness_df <- fitness_df %>%
  group_by(strain, background) %>% 
  summarize(average = mean(ln_relative_fitness)) %>% 
  mutate(lower_CI = average - 0.094547) %>% 
  mutate(upper_CI = average + 0.094547)


## Analyses ##


# One-sample t-test of average fitness costs of resistant lines
t.test(avg_fitness_df$average, mu = 0, alternative = "less")


# ANOVA tests of fitness differences between strains (both including and excluding strains without identified mutations; see Card et al. 2021)
FitnessAnova <- function(dat, exclude = FALSE, tukey = FALSE) {
  if (exclude == TRUE) {
    exclude_df <- dat %>% 
      filter(!strain %in% c("KJC65", "KJC66"))
    
    exclude_aov <- aov(ln_relative_fitness ~ strain, data = exclude_df)
    
    if (tukey == FALSE) {
      return(summary(exclude_aov))
    }
    else {
      return(exclude_aov)
    }
  }
  else { 
    include_aov <- aov(ln_relative_fitness ~ strain, data = dat)
    
    return(summary(include_aov))
  }
}

FitnessAnova(fitness_df)
FitnessAnova(fitness_df, exclude = TRUE)

# More data wrangling - Merge data frame of average strain fitness with MIC values from Card et al. 2019
avg_fitness_combined_df <- strain_MICs_df %>% 
  select(!background) %>% 
  left_join(avg_fitness_df, strain_MICs_df, by = "strain") %>% 
  filter(!strain %in% c("KJC65", "KJC66")) %>% 
  mutate(MIC_parent = log2(MIC_parent)) %>% 
  mutate(MIC_daughter = log2(MIC_daughter)) %>% 
  mutate(fold_change = (MIC_daughter - MIC_parent))


# Correlations between strain fitness and MIC
cor.test(avg_fitness_combined_df$average, avg_fitness_combined_df$MIC_daughter, alternative = "two.sided")

# Correlations between strain fitness and level of resistance conferred by mutation
cor.test(avg_fitness_combined_df$average, avg_fitness_combined_df$fold_change, alternative = "two.sided")

# ANOVA to test for possible interaction between genetic background and MIC on average fitness costs
aov(average ~ background * MIC_daughter, data = avg_fitness_combined_df) %>% summary()

# ANOVA to test for main effect of genetic background on average fitness costs
aov(average ~ background, data = avg_fitness_combined_df) %>% summary()

# ANOVA of fitness differences between strains (EXCLUDING KJC73 and KJC80)
fitness_df %>% 
  filter(!strain %in% c("KJC65", "KJC66","KJC73", "KJC80")) %>% 
  aov(ln_relative_fitness ~ strain, data = .) %>% summary()

# Compare resistant lines with single mutations against lines with multiple mutations (i.e., hitchhiking hypothesis)
avg_fitness_mult_df <- avg_fitness_combined_df %>% 
  filter(strain %in% c("KJC61", "KJC74", "KJC75", "KJC80", "KJC81"))

avg_fitness_single_df <- anti_join(avg_fitness_combined_df, avg_fitness_mult_df)

mean(avg_fitness_single_df$average)
mean(avg_fitness_mult_df$average)

t.test(avg_fitness_single_df$average, avg_fitness_mult_df$average, alternative = "greater", )

# ANOVA of ancestor-derived mutants
fitness_df %>% 
  filter(background == "Ancestor") %>% 
  aov(ln_relative_fitness ~ strain, data = .) %>% summary()

# ANOVA of ancestor-derived mutants with SINGLE mutation
fitness_df %>% 
  filter(strain %in% c("KJC62", "KJC63")) %>% 
  aov(ln_relative_fitness ~ strain, data = .) %>% summary()

# Tukey test

tukey_test <- FitnessAnova(fitness_df, exclude = TRUE, tukey = TRUE) %>% 
  HSD.test(., trt = "strain")

tukey_groups <- data.frame(group = c("bc", "cd", "cd", "d", "d", "cd", "cd", "ab", "c", "c", "a", "cd", "cd", "cd"))
avg_fitness_combined_df <- bind_cols(avg_fitness_combined_df, tukey_groups)


# Visualizations of the data


fitness_versus_MIC_plot <- avg_fitness_combined_df %>% 
  ggplot(aes(x = MIC_daughter, y = average)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
    stat_cor(method = "pearson", cor.coef.name = "r", digits = 4, label.x = 2, label.y = -0.15, geom = "text", size = 5) +
    scale_x_continuous(limits = c(1, 3), breaks = c(1, 2, 3)) +
    labs(y = "Log"[e]~"relative fitness", x = "Log"[2]~"MIC") +
    theme_cowplot()

fitness_versus_foldMIC_plot <- avg_fitness_combined_df %>% 
  ggplot(aes(x = fold_change, y = average)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 1) +
    stat_cor(method = "pearson", cor.coef.name = "r", digits = 4, label.x = 1, label.y = -0.15, geom = "text", size = 5) +
    scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
    labs(y = "Log"[e]~"relative fitness", x = "Log"[2]~"fold-increase in MIC") +
    theme_cowplot() +
    theme(axis.title.y = element_text(color = "white"),
          axis.text.y = element_text(color = "white"),
          axis.ticks.y = element_blank())

corr_plot <- plot_grid(fitness_versus_MIC_plot, fitness_versus_foldMIC_plot,
                       labels = "AUTO")

# ggsave("corr_plot.pdf", corr_plot, path = "Figures", device = "pdf", width = 10, height = 6, units = "in")


avg_fitness_plot <- avg_fitness_combined_df %>% 
  ggplot(aes(x = strain, y = average)) +
    geom_point(aes(x = reorder(strain, average) , y = average), size = 3) +
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.3) +
    geom_text(aes(label = group), vjust = -7) +
    scale_y_continuous(limits = c(-0.5, 0.2), breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2)) +
    scale_x_discrete(labels = c("Ara+4-3", "Ara+5-2", "Ancestor-1", "Ara\u20135-2", "Ara\u20136-2", "Ara+5-3", "Ara\u20136-1", "Ancestor-3", "Ancestor-2",
                                "Ara\u20136-3", "Ara\u20135-3", "Ara+4-2", "Ancestor-4", "Ara+4-1")) +
    labs(y = expression("Log"[e]~"relative fitness")) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

avg_fitness_plot

# ggsave("avg_fitness_plot.pdf", avg_fitness_plot, path = "Figures", device = "pdf", width = 10, height = 6, units = "in")












