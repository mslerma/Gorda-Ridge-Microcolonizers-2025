#### ASVs on Substrate and in Vent fluid ####

format_18S <- mc_18S_df %>%
  mutate(SAMPLE_TYPE = case_when(
    SAMPLE == "Background seawater" ~ "Background", 
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent", 
    TRUE ~ "Substrate")) %>%
  select(FeatureID, SEQ_AVG, SAMPLE_TYPE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE_TYPE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

format_18S_assigned <- format_18S %>%
  mutate(DISTRIBUTION = case_when(
    (is.na(Background) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent only", 
    (is.na(Background) & !(is.na(Substrate)) & !(is.na(Vent))) ~ "Vent & Substrate only", 
    (is.na(Background) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate only", 
    (!(is.na(Background)) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent & Background only",
    (!(is.na(Background)) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate & Background only", 
    (!(is.na(Background)) & is.na(Substrate) & is.na(Vent)) ~ "Background only", 
    (!(is.na(Background)) & !(is.na(Substrate)) & !(is.na(Vent)) ~ "All")
  ))

key_dist_18S <- format_18S_assigned %>%
  select(FeatureID, DISTRIBUTION) %>%
  distinct()

df_18S_all <- mc_18S_df %>%
  left_join(key_dist_18S)

df_18S_vs <- mc_18S_df %>%
  left_join(key_dist_18S) %>%
  filter(DISTRIBUTION == "Vent & Substrate only") %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right") %>%
  filter(Phylum != "Metazoa")

shared_MC2 <- df_18S_vs %>% 
  filter(MC %in% c("MC2", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()
shared_MC3 <- df_18S_vs %>%
  filter(MC %in% c("MC3", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()
shared_MC1 <- df_18S_vs %>%
  filter(MC %in% c("MC1", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()
shared_MC4 <- df_18S_vs %>%
  filter(MC %in% c("MC4", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()
shared_MC5 <- df_18S_vs %>%
  filter(MC %in% c("MC5", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()
shared_MC6 <- df_18S_vs %>%
  filter(MC %in% c("MC6", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%   
  ungroup()

across_all <- df_18S_vs %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 7) %>%
  ungroup()

df_18S_all_filt <- df_18S_all %>%
  filter(Phylum != "Metazoa")

unique_18S_vs <- df_18S_vs %>%
  group_by(FeatureID, Taxon) %>%
  summarise(ASV_OCCURANCE = n()) 

## Sequence averages > 99
df_18s_vs_filt <- df_18s_vs %>%
  filter(Phylum != "Metazoa") %>%
  filter(SEQ_AVG > 99)

vs_dist_plot <- df_18s_vs_filt %>%
  add_column(COUNT = 1) %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(SUM_PHYLUM= sum(COUNT)) %>% 
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)
length(unique(df_18s_vs_filt$FeatureID))

df_18s_vs_filt$ordered_SAMPLES <- factor(df_18s_vs_filt$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid")) 
vs_dist_plot$ordered_SAMPLES <- factor(vs_dist_plot$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid")) 

ggplot(vs_dist_plot, aes(x = ordered_SAMPLES, y = SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(cols = vars(SUBSTRATE), scales = "free", space = "free") +
  scale_fill_manual(values = c("#a6859f", "#d9bdc8", "#ffffff", "#aee2ff", "#8db7ff", "#6d80fa", 
                               "#8465ec", "#834dc4", "#7d2da0", "#4e187c")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 

major_unique_18S_vs <- df_18S_vs %>%
  filter(SEQ_AVG > 99) %>%
  group_by(FeatureID, Taxon) %>%
  summarise(ASV_OCCURANCE = n()) 

#### ASVs on Substrates only ####
mc_df_18s <- df_18s_all %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid") %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta") %>%
  filter(!is.na(Phylum)) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  add_column(PRESENCE = 1)

MC_only_18S <- mc_df_18s %>%
  filter(DISTRIBUTION == "Substrate only")

### ASVs unique to MCs ###
view(MC_only_18S)
asv_counts <- MC_only_18S %>%
  distinct(FeatureID, MC) %>%
  group_by(FeatureID) %>%
  summarise(mc_count = n())

unique_asvs <- asv_counts %>%
  filter(mc_count == 1) %>%
  pull(FeatureID)

df_unique <- MC_only_18S %>%
  filter(FeatureID %in% unique_asvs)

df_unique_by_mc <- df_unique %>%
  group_by(MC) %>%
  summarise(unique_ASVs = list(unique(FeatureID)))

key_ASV_taxons <- MC_only_18S %>%
  ungroup() %>%
  select(FeatureID, Supergroup, Phylum, Class, Order, Family, Genus, Species, SEQ_AVG) %>%
  distinct()

asv_mc_sets <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarise(mc_combo = paste(sort(unique(MC)), collapse = ","), .groups = "drop")

joined_ASVs_MC_only <- asv_mc_sets %>%
  left_join(key_ASV_taxons)

asv_subs_sets <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarise(subs_combo = paste(sort(unique(SUBSTRATE)), collapse = ","), .groups = "drop")

joined_ASVs_SUBS_only <- asv_subs_sets %>%
  left_join(key_ASV_taxons)

# MC only ASVs on quartz only
quartz_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Quartz")

# MC only ASVs on riftia only
riftia_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Riftia")

# MC only ASVs on shell only
shell_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Shell")

# MC only ASVs on all substrates 
all_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Quartz,Riftia,Shell")
length(unique(joined_ASVs_SUBS_only$FeatureID))

for_subs_table <- joined_ASVs_SUBS_only %>%
  mutate(SUBS_DIST = case_when(
    subs_combo == "Quartz,Riftia,Shell" ~ "Across all",
    subs_combo == "Quartz,Shell" ~ "Quartz and Shell",
    subs_combo == "Quartz,Riftia" ~ "Quartz and Riftia",
    subs_combo == "Riftia,Shell" ~ "Riftia and Shell",
    subs_combo == "Quartz" ~ "Quartz only",
    subs_combo == "Riftia" ~ "Riftia only",
    subs_combo == "Shell" ~ "Shell only")
  )

asv_counts_table <- MC_only_18S %>%
  group_by(FeatureID, SUBSTRATE) %>%
  summarise(PRESENCE = max(PRESENCE), .groups = "drop") %>%  
  filter(PRESENCE == 1) %>%
  group_by(FeatureID) %>%
  summarise(SUBS_DIST = paste(sort(unique(SUBSTRATE)), collapse = " and "), .groups = "drop") %>%
  mutate(SUBS_DIST = case_when(
    SUBS_DIST == "Quartz and Riftia and Shell" ~ "Across all",
    TRUE ~ SUBS_DIST
  )) %>%
  distinct(FeatureID, SUBS_DIST) %>%   
  count(SUBS_DIST, name = "Num_distinct_ASVs") %>%
  arrange(desc(Num_distinct_ASVs))

MC_only_18S_filt <- MC_only_18S %>%
  filter(SEQ_AVG > 99) 
length(unique(MC_only_18S_filt$FeatureID))

# ASV distribution
mc_dist_plot <- MC_only_18S_filt %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(SUM_PHYLUM= sum(PRESENCE)) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)
unique(mc_dist_plot$Phylum)

mc_dist_plot$ordered_SAMPLES <- factor(mc_dist_plot$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid")) 

ggplot(mc_dist_plot, aes(x = ordered_SAMPLES, y = SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(cols = vars(SUBSTRATE), scales = "free", space = "free") +
  scale_fill_manual(values = c("#fcef8d", "#ffb879", "#ea6262", "#cc425e",
                               "#a32858", "#751756", "#390947", "#611851", 
                               "#873555", "#a6555f", "#c97373", "#f2ae99", 
                               "#ffc3f2", "#ee8fcb", "#d46eb3", "#873e84")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 

## Sequence averages > 99
more_than_one_18S <- for_subs_table %>%
  filter(subs_combo != "Quartz" & subs_combo != "Riftia" & subs_combo != "Shell") %>%
  filter(SEQ_AVG > 99)

df_more_than_one <- more_than_one_18S %>%
  select(-SEQ_AVG) %>%
  distinct() 



