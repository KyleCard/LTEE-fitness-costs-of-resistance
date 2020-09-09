setwd("C:\\Users\\cardk\\Box Sync\\Manuscripts\\Fitness costs of ABR mutations in the LTEE\\Card et al. 2020 - Fitness costs of ABR mutations in the LTEE")

library(tidyverse)
library(agricolae)
library(cowplot)

counts_df <- read_csv("competition_data.csv")
counts_df$paired_ID <- as.character(counts_df$paired_ID)
strain_MICs <- read_csv("TET_strains_MICs.csv") # MIC data from Card et al. 2019


## Data wrangling ##
relative_fitness_df <- counts_df %>% 
  mutate(malth_comp1 = log(comp1_d3 * (100^3) / comp1_d0)) %>%
  mutate(malth_comp2 = log(comp2_d3 * (100^3) / comp2_d0)) %>% 
  mutate(relative_fitness = malth_comp1 / malth_comp2)

CompareFitness <- function(dat) {
  fitness_vec <- c()
  
  for (i in unique(dat$block)) {
    block_frame <- dat %>% 
      filter(block == i)
    
    for (j in unique(block_frame$paired_ID)) {
      paired_frame <- block_frame %>% 
        filter(paired_ID == j)
      
      true_fitness <- paired_frame[1, 16] / paired_frame[2, 16]

      fitness_vec <- bind_rows(fitness_vec, true_fitness)
    }
  }
  return(fitness_vec)
}


# Normalizes relative fitness based upon common competitor values
normalized_fitness_col <- CompareFitness(relative_fitness_df) %>% 
  drop_na %>% 
  transmute(ln_relative_fitness = log(relative_fitness)) # Log(e) transform the relative fitness values

subset_fitness_df <- relative_fitness_df %>% 
  drop_na

fitness_df <- bind_cols(subset_fitness_df, normalized_fitness_col)
fitness_df$background <- as.factor(fitness_df$background)
fitness_df$competitor_1 <- as.factor(fitness_df$competitor_1)


# Create data frame with average fitness of each strain
avg_fitness_strain <- fitness_df %>%
  group_by(competitor_1, background) %>% 
  summarize(average = mean(ln_relative_fitness)) %>% 
  mutate(lower_CI = average - 0.094547) %>% 
  mutate(upper_CI = average + 0.094547)


# One-sample t-test of average fitness costs of resistant mutants
ttest_avg_costs <- t.test(avg_fitness_strain$average, mu = 0, alternative = "less")
ttest_avg_costs

# ANOVA of fitness differences between strains (INCLUDING strains without identified mutations; see Card et al. 2020)
strains_lm <- lm(ln_relative_fitness ~ competitor_1, data = fitness_df)
strains_aov <- aov(strains_lm)
summary(strains_aov)


# ANOVA of fitness differences between strains (EXCLUDING strains without identified mutations)
fitness_df_ex <- fitness_df %>% 
  filter(!competitor_1 %in% c("KJC65", "KJC66"))

strains_lm_ex <- lm(ln_relative_fitness ~ competitor_1, data = fitness_df_ex)
strains_aov_ex <- aov(strains_lm_ex)
summary(strains_aov_ex)


# Merge data frame of average strain fitness with MIC values from Card et al. 2019
avg_fitness_MICs <- avg_fitness_strain %>% 
  rename(strain = competitor_1)

avg_fitness_MICs <- left_join(avg_fitness_MICs, strain_MICs, by = "strain") %>% 
  filter(!strain %in% c("KJC65", "KJC66"))

avg_fitness_MICs$MIC_parent <- log2(avg_fitness_MICs$MIC_parent)
avg_fitness_MICs$MIC_daughter <- log2(avg_fitness_MICs$MIC_daughter)

avg_fitness_MICs <- avg_fitness_MICs %>% mutate(fold_change = (MIC_daughter - MIC_parent))


# Correlations between strain fitness and MIC
correl_MIC <- cor.test(avg_fitness_MICs$average, avg_fitness_MICs$MIC_daughter, alternative = "two.sided")
correl_MIC


# Correlations between strain fitness and level of resistance conferred by mutation
correl_fold_change <- cor.test(avg_fitness_MICs$average, avg_fitness_MICs$fold_change, alternative = "two.sided")
correl_fold_change


# ANOVA of fitness differences between genetic backgrounds
background_lm <- lm(average ~ background.x, data = avg_fitness_MICs)
background_aov <- aov(background_lm)
summary(background_aov)


# ANOVA of fitness differences between strains (EXCLUDING KJC73 and KJC80)
excluding_lowest <- fitness_df %>% 
  filter(!competitor_1 %in% c("KJC65", "KJC66","KJC73", "KJC80"))

excluding_lowest_lm <- lm(ln_relative_fitness ~ competitor_1, data = excluding_lowest)
excluding_lowest_aov <- aov(excluding_lowest_lm)
summary(excluding_lowest_aov)


# Compare resistant lines with single mutations against lines with multiple mutations (i.e., hitchhiking hypothesis)
avg_fitness_mult <- avg_fitness_MICs %>% 
  filter(strain %in% c("KJC61", "KJC74", "KJC75", "KJC80", "KJC81"))

avg_fitness_single <- anti_join(avg_fitness_MICs, avg_fitness_mult)

mean_single <- mean(avg_fitness_single$average)
mean_single

mean_mult <- mean(avg_fitness_mult$average)
mean_mult

ttest_single_mult <- t.test(avg_fitness_single$average, avg_fitness_mult$average, alternative = "greater", )
ttest_single_mult


# ANOVA of ancestor-derived mutants
fitness_ancestors <- fitness_df %>% 
  filter(background == "Ancestor")

ancestors_lm <- lm(ln_relative_fitness ~ competitor_1, data = fitness_ancestors)
ancestors_aov <- aov(ancestors_lm)
summary(ancestors_aov)


# ANOVA of ancestor-derived mutants with SINGLE mutation
fitness_ancestors_single <- fitness_ancestors %>% 
  filter(competitor_1 %in% c("KJC62", "KJC63"))

ancestors_single_lm <- lm(ln_relative_fitness ~ competitor_1, data = fitness_ancestors_single)
ancestors_single_aov <- aov(ancestors_single_lm)
summary(ancestors_single_aov)


# Tukey test
tukey_test <- HSD.test(strains_aov_ex, trt = "competitor_1")
tukey_groups <- data.frame(group = c("bc", "cd", "cd", "d", "d", "cd", "cd", "ab", "c", "c", "a", "cd", "cd", "cd"))

avg_fitness_MICs <- bind_cols(avg_fitness_MICs, tukey_groups)


# Visualizations of the data
avg_rel_fitness_plot <- avg_fitness_MICs %>% 
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

avg_rel_fitness_plot

# ggsave("avg_fitness_strain.tif", avg_rel_fitness_plot, path = "Figures", device = "tiff", width = 10, height = 6, units = "in", dpi = 300, compression = "lzw")












