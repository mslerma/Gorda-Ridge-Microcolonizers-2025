# 5/7/25 
#### RDA Analysis ####
library(vegan)
library(phyloseq)
library(ggplot2)
subs_for_rda <- mc_tmp_18 %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "left") %>%
  filter(!(SAMPLE == "Background seawater" | SAMPLE == "Mt Edwards diffuse fluid"))

subs_for_rda <- subs_for_rda %>%
  mutate(Time = case_when(
    MC == "MC1" ~ 7,
    MC == "MC2" ~ 6,
    MC == "MC3" ~ 6,
    MC == "MC4" ~ 7,
    MC == "MC5" ~ 8,
    MC == "MC6" ~ 8
  )) %>%
  mutate(Temperature = case_when(
    MC == "MC1" ~ 8.5,
    MC == "MC2" ~ 40.5,
    MC == "MC3" ~ 16.58,
    MC == "MC4" ~ 8.23,
    MC == "MC5" ~ 9.29,
    MC == "MC6" ~ 7.18
  ))

comm_matrix <- subs_for_rda %>%
  ungroup() %>%
  select(SAMPLE, FeatureID, SEQ_AVG) %>%
  pivot_wider(names_from = FeatureID, values_from = SEQ_AVG, values_fill = 0) %>%
  column_to_rownames(var = "SAMPLE")

env_data <- subs_for_rda %>%
  ungroup() %>%
  select(SAMPLE, MC, SUBSTRATE, Time, Temperature) %>%
  distinct() %>%
  column_to_rownames(var = "SAMPLE")

##
samp_env_data <- subs_for_rda %>%
  ungroup() %>%
  select(SAMPLE, Time, Temperature) %>%
  distinct()

comm_hell <- decostand(comm_matrix, method = "hellinger")

rda_MCs_18S <- rda(comm_hell ~ Time + Temperature, data = samp_env_data)
summary(rda_MCs_18S)
anova.cca(rda_MCs_18S, by = "term")
RsquareAdj(rda_MCs_18S)
all_MCs_TnT_rda <- plot(rda_MCs_18S, main = "All MCs")
  
##

rda_all_18S <- rda(comm_hell ~ SUBSTRATE + Time + Temperature, data = env_data)
summary(rda_all_18S)
anova.cca(rda_all_18S, by = "term")
RsquareAdj(rda_all_18S)
all_MCs_STT_rda <- plot(rda_all_18S, main = "All MCs")

RsquareAdj(rda(comm_hell ~ Time, data = env_data))
RsquareAdj(rda(comm_hell ~ Temperature, data = env_data))
RsquareAdj(rda(comm_hell ~ SUBSTRATE, data = env_data))

rda_time_18S <- rda(comm_hell ~ Time, data = env_data)
summary(rda_time_18S)
anova.cca(rda_time_18S, by = "term")
RsquareAdj(rda_time_18S)
all_time_rda <- plot(rda_time_18S)

rda_temp_18S <- rda(comm_hell ~ Temperature, data = env_data)
summary(rda_temp_18S)
anova.cca(rda_temp_18S, by = "term")
RsquareAdj(rda_temp_18S)
all_temp_rda <- plot(rda_temp_18S)

#### Quartz ####
env_quartz <- env_data %>%
  filter(SUBSTRATE == "Quartz")
comm_quartz <- comm_hell[rownames(env_quartz), ]
rda_quartz_18S <- rda(comm_quartz ~ Time + Temperature, data = env_quartz)

summary(rda_quartz_18S)
anova.cca(rda_quartz_18S, by = "term")
RsquareAdj(rda_quartz_18S)
quartz_TnT_rda <- plot(rda_quartz_18S)

#### Riftia ####
env_riftia <- env_data %>%
  filter(SUBSTRATE == "Riftia")
comm_riftia <- comm_hell[rownames(env_riftia), ]
rda_riftia_18S <- rda(comm_riftia ~ Time + Temperature, data = env_riftia)

summary(rda_riftia_18S)
anova.cca(rda_riftia_18S, by = "term")
RsquareAdj(rda_riftia_18S)
riftia_TnT_rda <- plot(rda_riftia_18S)

#### Shell ####
env_shell <- env_data %>%
  filter(SUBSTRATE == "Shell")
comm_shell <- comm_hell[rownames(env_shell), ]
rda_shell_18S <- rda(comm_shell ~ Time + Temperature, data = env_shell)

summary(rda_shell_18S)
anova.cca(rda_shell_18S, by = "term")
RsquareAdj(rda_shell_18S)
shell_TnT_rda <- plot(rda_shell_18S)

#### MC 2&3 ####
env_MC23 <- env_data %>%
  filter(MC %in% c("MC2", "MC3"))

comm_MC23 <- comm_matrix[rownames(env_MC23), ]

rda_MC23_time <- rda(comm_MC23 ~ Time, data = env_MC23)
summary(rda_MC23_time)
plot(rda_MC23_time)
rda_MC23_temp <- rda(comm_MC23 ~ Temperature, data = env_MC23)
summary(rda_MC23_temp)
anova.cca(rda_MC23_temp, by = "term")
RsquareAdj(rda_MC23_temp)
plot(rda_MC23_temp)
rda_MC23_subs <- rda(comm_MC23 ~ SUBSTRATE, data = env_MC23)
summary(rda_MC23_subs)
anova.cca(rda_MC23_subs, by = "term")
plot(rda_MC23_subs)

rda_MC23_subs_temp <- rda(comm_MC23 ~ Temperature + SUBSTRATE, data = env_MC23)
summary(rda_MC23_subs_temp)
plot(rda_MC23_subs)
#### ####
matrix_rdasubs_wide <- subs_for_rda %>%
  select(SAMPLE, FeatureID, SEQ_AVG) %>%
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE, values_from = SEQ_AVG, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()
covariance_rdasubs <- as.matrix(matrix_rdasubs_wide) %*% t(matrix_rdasubs_wide)
det(covariance_rdasubs)
stand_matrix_rdasubs <- decostand(matrix_rdasubs_wide, method = "hellinger", MARGIN = 2)
clr
####
rda_plots <- list()
wrap_plots(rda_plots, ncol = 3)