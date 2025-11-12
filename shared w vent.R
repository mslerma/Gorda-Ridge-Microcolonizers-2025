### look into unique vs shared taxa on MCs and Vent ####
format_18s_assigned
?n_distinct()
#### 16S ####
view(df_16s_all)
sep_16s_vs <- df_16s_all %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right")
  
shared_MC2_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC2", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

shared_MC3_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC3", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC1_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC1", "Mt Edwards diffuse fluid")) %>% 
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC4_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC4", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC5_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC5", "Mt Edwards diffuse fluid")) %>% 
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC6_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC6", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

#### 18S ####
df_18s_all_filt <- df_18s_all %>%
  filter(Phylum != "Metazoa")

df_18s_vs <- df_18s_vs %>%
  filter(Phylum != "Metazoa")

df_18s_vs_filt <- df_18s_vs %>%
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

MC_only_18S_filt <- MC_only_18S %>%
  filter(SEQ_AVG > 99) 
length(unique(MC_only_18S_filt$FeatureID))

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
#
unique(sep_18s_vs$SAMPLE)
sep_18s_vs <- df_18s_vs %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right") %>%
  filter(Phylum != "Metazoa")

no_met_18s_vs <- df_18s_all %>%
  filter(Phylum != "Metazoa")

# Find FeatureIDs that are present in both locations
shared_MC2 <- sep_18s_vs %>% 
  filter(MC %in% c("MC2", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC3 <- sep_18s_vs %>%
  filter(MC %in% c("MC3", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC1 <- sep_18s_vs %>%
  filter(MC %in% c("MC1", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

## MC1 has 553 unique ASVs, 705 total 
MC1 <- mc_df_18s %>%
  filter(MC == "MC1")
length(MC1$FeatureID)

shared_MC4 <- sep_18s_vs %>%
  filter(MC %in% c("MC4", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC5 <- sep_18s_vs %>%
  filter(MC %in% c("MC5", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC6 <- sep_18s_vs %>%
  filter(MC %in% c("MC6", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

unique(shared_MC6$Class)

across_all <- sep_18s_vs %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 7) %>%
  ungroup()


unique(across_all$FeatureID)
df_18s_all_nomet <- df_18s_all %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta")

sep_18s_vs <- df_18s_vs %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right") %>%
  filter(Phylum != "Metazoa")
unique(sep_18s_vs$SAMPLE)

#### ####
length(unique(sep_18s_vs$FeatureID))
length(unique(df_18s_vs$FeatureID))


MC2_MC3_asvs <- sep_18s_vs %>%
  filter(MC %in% c("MC2", "MC3"))

vent_asvs <- sep_18s_vs %>%
  filter(SAMPLE == "Mt Edwards diffuse fluid")

shared_asvs <- sep_18s_vs %>%
  filter(FeatureID %in% MC2_MC3_asvs$FeatureID & FeatureID %in% vent_asvs$FeatureID)
shared_asvs

vs_presence <- sep_18s_vs %>%
  filter(!is.na(Phylum)) %>%
  add_column(COUNT = 1)

ggplot(vs_presence, aes(x = SUBSTRATE, y = FeatureID, fill = COUNT)) +
  geom_tile(stat = "identity") + 
  theme_classic() +
  facet_grid(cols = vars(MC), rows = vars(Phylum), scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = -45, hjust = 0), strip.text.y = element_text(angle = 0)) +
  labs(x = "", fill="Presence")

#### diversity of taxa at vent and substrates only, compare with sub only 
#### 16S ####
unique(filt_df_16s_svs$Phylum)
svs_taxa <- df_16s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

filt_df_16s_svs <- df_16s_svs %>%
  filter(SEQ_AVG > 10)

filt_svs_taxa <- filt_df_16s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(svs_taxa, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()
ggplot(filt_svs_taxa, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()
 
wide_svs_16s <- df_16s_svs %>%
  group_by(DISTRIBUTION, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>%
  select(DISTRIBUTION, FeatureID, SUM) %>%
  pivot_wider(names_from = DISTRIBUTION, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_svs_16s_mat <- as.matrix(wide_svs_16s)
stand_wide_svs_16S_mat <- decostand(wide_svs_16s_mat, method = "hellinger", MARGIN = 2)
shannon_svs_16S <- diversity(wide_svs_16s_mat, index = "shannon", MARGIN = 2)
shannon_svs_16S
plot_shannon_svs_16S <- as.data.frame(shannon_svs_16S) %>%
  rownames_to_column(var = "DISTRIBUTION")
plot_shannon_svs_16S
ggplot(plot_shannon_svs_16S, aes(x= DISTRIBUTION, y= shannon_svs_16S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "", y= "16S Shannon diversity")

#### 18S ####

unique(filt_df_18s_svs$Phylum)
svs_taxa_18S <- df_18s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  filter(!(is.na(Phylum)))

filt_df_18s_svs <- df_18s_svs %>%
  filter(SEQ_AVG > 10) %>%
  filter(!(is.na(Phylum)))

filt_svs_taxa_18S <- filt_df_18s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(svs_taxa_18S, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()
ggplot(filt_svs_taxa_18S, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()

wide_svs_16s <- df_16s_svs %>%
  group_by(DISTRIBUTION, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>%
  select(DISTRIBUTION, FeatureID, SUM) %>%
  pivot_wider(names_from = DISTRIBUTION, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_svs_16s_mat <- as.matrix(wide_svs_16s)
stand_wide_svs_16S_mat <- decostand(wide_svs_16s_mat, method = "hellinger", MARGIN = 2)
shannon_svs_16S <- diversity(wide_svs_16s_mat, index = "shannon", MARGIN = 2)
shannon_svs_16S
plot_shannon_svs_16S <- as.data.frame(shannon_svs_16S) %>%
  rownames_to_column(var = "DISTRIBUTION")
plot_shannon_svs_16S
ggplot(plot_shannon_svs_16S, aes(x= DISTRIBUTION, y= shannon_svs_16S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "", y= "16S Shannon diversity")

unique(mc_df_18s$DISTRIBUTION)
MC_only_18S <- mc_df_18s %>%
  filter(DISTRIBUTION == "Substrate only")

MC_only_16S <- mc_df_16s %>%
  filter(DISTRIBUTION == "Substrate only")
length(unique(MC_only_18S$FeatureID))
length(unique(MC_only_16S$FeatureID))
MC_only_16S <- MC_only_16S %>%
  unite(SAMPLE, MC, SUBSTRATE, sep = "-", remove = FALSE)

plot_MC_only_18S <- MC_only_18S %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')
ggplot(plot_MC_only_18S, aes(x= SAMPLE, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45))

library(dplyr)

ASVs_18S_mult_sites <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarize(site_count = n_distinct(SAMPLE)) %>%
  filter(site_count > 1)
  
ASV_df_18S_mult_site <- MC_only_18S %>%
  filter(FeatureID %in% ASVs_18S_mult_sites$FeatureID)

ASVs_18S_df <- MC_only_18S %>%
  left_join(ASVs_18S_mult_sites) %>%
  filter(!(is.na(site_count)))
length(ASVs_18S_df)

plot_MC_ASVs_18S <- ASVs_18S_df %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') 
ggplot(plot_MC_ASVs_18S, aes(x= SAMPLE, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45))
unique(ASVs_18S_df$Class)

ASVs_16S_mult_sites <- MC_only_16S %>%
  group_by(FeatureID) %>%
  summarize(site_count = n_distinct(SAMPLE)) %>%
  filter(site_count > 1)

ASV_df_16S_mult_site <- MC_only_16S %>%
  filter(FeatureID %in% ASVs_16S_mult_sites$FeatureID)

ASVs_16S_df <- MC_only_16S %>%
  left_join(ASVs_16S_mult_sites) %>%
  filter(!(is.na(site_count)))
length(unique(ASVs_16S_df$FeatureID))


length(unique(sep_18s_vs$FeatureID))
major_unique_18S_vs <- df_18S_vs %>%
  filter(SEQ_AVG > 99) %>%
  group_by(FeatureID, Taxon) %>%
  summarise(ASV_OCCURANCE = n()) 

unique_18S_vs <- df_18S_vs %>%
  group_by(FeatureID, Taxon) %>%
  summarise(ASV_OCCURANCE = n()) 

distinct(FeatureID, MC) %>%
  group_by(FeatureID) %>%
  summarise(mc_count = n())
summarise(mc_count = n())
summarise(unique_ASVs = list(unique(FeatureID)))
mutate(SUBS_DIST = case_when(
  SUBS_DIST == "Quartz and Riftia and Shell" ~ "Across all",
  TRUE ~ SUBS_DIST
)) %>%
  