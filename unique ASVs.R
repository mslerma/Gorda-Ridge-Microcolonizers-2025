MC1_only <- MC_only_18S %>%
  filter(MC == "MC1")
MC6_only <- MC_only_18S %>%
  filter(MC == "MC6")
MC5_only <- MC_only_18S %>%
  filter(MC == "MC5")

filt_mc_18S_df <- mc_18S_df %>%
  ungroup() %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta")
unique(filt_mc_18S_df$FeatureID)
length(filt_mc_18S_df$FeatureID)
more_filt_mc_18S_df <- filt_mc_18S_df %>%
  filter(SEQ_AVG > 10)
length(MC_only_18S$FeatureID)


#### ASVs unique to MCs ####
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

##
asv_mc_counts <- MC_only_18S %>%
  distinct(FeatureID, MC) %>%
  group_by(FeatureID) %>%
  summarise(mc_count = n(), .groups = "drop")

df_labeled <- MC_only_18S %>%
  left_join(asv_mc_counts, by = "FeatureID") %>%
  mutate(asv_type = ifelse(mc_count == 1, "unique", "shared"))

##
total_mcs <- MC_only_18S %>% distinct(MC) %>% nrow()

df_labeled2 <- df_labeled %>%
  mutate(asv_type = case_when(
    mc_count == 1 ~ "unique",
    mc_count == total_mcs ~ "ubiquitous",
    TRUE ~ "shared"
  ))

##
asv_mc_sets <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarise(mc_combo = paste(sort(unique(MC)), collapse = ","), .groups = "drop")

df_tagged <- MC_only_18S %>%
  left_join(asv_mc_sets, by = "FeatureID")

df_tagged %>%
  distinct(FeatureID, mc_combo) %>%
  group_by(mc_combo) %>%
  summarise(ASVs = list(FeatureID), count = n(), .groups = "drop")

install.packages("UpSetR")
library(UpSetR)
# Create a binary presence/absence matrix
asv_matrix <- MC_only_18S %>%
  ungroup() %>%
  select(MC, FeatureID, PRESENCE) %>%
  group_by(FeatureID, MC) %>%
  summarise(PRESENCE = max(PRESENCE), .groups = "drop") %>%  # collapse duplicates
  pivot_wider(
    names_from = MC,
    values_from = PRESENCE,
    values_fill = list(PRESENCE = 0)   # fill absent with 0
  ) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

# convert your matrix to data.frame
asv_df <- as.data.frame(asv_matrix)

# run upset plot
upset(
  asv_df,
  sets = colnames(asv_df),
  order.by = "freq"
)
upset(
  asv_df,
  sets = colnames(asv_df),
  order.by = "freq",            # order intersections by frequency
  main.bar.color = "steelblue", # color for intersection bars
  sets.bar.color = "gray40",    # color for set size bars
  text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.2), # resize labels
  mainbar.y.label = "ASV Intersection Size",
  sets.x.label = "ASVs per MC"
)

? upset()

key_ASV_taxons <- MC_only_18S %>%
  ungroup() %>%
  select(FeatureID, Supergroup, Phylum, Class, Order, Family, Genus, Species, SEQ_AVG) %>%
  distinct()

joined_ASVs_MC_only <- asv_mc_sets %>%
  left_join(key_ASV_taxons)

#### ASVs unique to substrates ####
asv_matrix_sub <- MC_only_18S %>%
  ungroup() %>%
  select(SUBSTRATE, FeatureID, PRESENCE) %>%
  group_by(FeatureID, SUBSTRATE) %>%
  summarise(PRESENCE = max(PRESENCE), .groups = "drop") %>%  # collapse duplicates
  pivot_wider(
    names_from = SUBSTRATE,
    values_from = PRESENCE,
    values_fill = list(PRESENCE = 0)
  ) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

asv_df_sub <- as.data.frame(asv_matrix_sub)

upset(
  asv_df_sub,
  sets = colnames(asv_df_sub),      # Quartz, Riftia, Shell
  order.by = "freq",
  main.bar.color = "lightpink",
  sets.bar.color = "gray40",
  text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.2),
  mainbar.y.label = "ASV Intersection Size",
  sets.x.label = "ASVs per
  Substrate"
)

asv_subs_sets <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarise(subs_combo = paste(sort(unique(SUBSTRATE)), collapse = ","), .groups = "drop")

joined_ASVs_SUBS_only <- asv_subs_sets %>%
  left_join(key_ASV_taxons)

quartz_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Quartz")
riftia_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Riftia")
shell_only_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Shell")
 
all_SUB <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Quartz,Riftia,Shell")
length(unique(joined_ASVs_SUBS_only$FeatureID))

###
unique(joined_ASVs_SUBS_only$Phylum)
unique(joined_ASVs_MC_only$Phylum)
library(dplyr)

class_sub_summary <- MC_only_18S %>%
  group_by(SUBSTRATE, Class) %>%
  summarise(total = sum(PRESENCE), .groups = "drop") %>%
  group_by(SUBSTRATE) %>%
  mutate(rel_abund = total / sum(total))

phyl_sub_summary <- MC_only_18S %>%
  group_by(SUBSTRATE, Phylum) %>%
  summarise(total = sum(PRESENCE), .groups = "drop") %>%
  group_by(SUBSTRATE) %>%
  mutate(rel_abund = total / sum(total))

library(ggplot2)

ggplot(class_sub_summary, aes(x = SUBSTRATE, y = rel_abund, fill = Class)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Substrate", y = "Relative Abundance", fill = "Class") +
  theme_minimal(base_size = 14)

ggplot(phyl_sub_summary, aes(x = SUBSTRATE, y = rel_abund, fill = Phylum)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Substrate", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal(base_size = 14)

class_sub_summary <- MC_only_18S %>%
  group_by(SUBSTRATE, Class) %>%
  summarise(total = sum(PRESENCE), .groups = "drop")


# Identify rare classes
rare_classes <- class_sub_summary %>%
  group_by(Class) %>%
  summarise(global_total = sum(total)) %>%
  filter(global_total < 10) %>%
  pull(Class)

# Recode rare ones to "Other"
class_sub_summary <- class_sub_summary %>%
  mutate(Class = ifelse(Class %in% rare_classes, "Other", Class)) %>%
  group_by(SUBSTRATE, Class) %>%
  summarise(total = sum(total), .groups = "drop") %>%
  group_by(SUBSTRATE) %>%
  mutate(rel_abund = total / sum(total))

ggplot(class_sub_summary, aes(x = SUBSTRATE, y = rel_abund, fill = Class)) +
  geom_col(color = "black") +
  scale_fill_hue() +
  labs(x = "Substrate", y = "Relative Abundance", fill = "Class") +
  theme_minimal(base_size = 14)


###
install.packages("ggupset")
library(ggupset)
alv <- c("Alveolata-Ellobiopsidae", "Alveolata-Perkinsea", "Alveolata-Unknown", "Alveolata-Chrompodellids", "Alveolata-Apicomplexa")
all_taxa_order <- c("Alveolata-Ciliophora", "Alveolata-Dinoflagellata", "Protalveolata", "Other Alveolata", "Amoebozoa", "Apusozoa", "Excavata", "Hacrobia", "Archaeplastida", "Rhizaria", "Rhizaria-Radiolaria", "Rhizaria-Cercozoa", "Stramenopiles", "Stramenopiles-Opalozoa", "Stramenopiles-Sagenista", "Stramenopiles-Ochrophyta", "Opisthokonta", "Unknown Eukaryota")
all_taxa_color = c("#fa9fb5", "#c51b8a", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c", "#c7e9b4", "#1d91c0", "#deebf7", "#253494", "#9e9ac8", "#238b45", "#54278f", "#bdbdbd", "#252525", "#fa9fb5", "#c51b8a", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c", "#c7e9b4", "#1d91c0", "#253494", "#9e9ac8", "#238b45", "#54278f", "#bdbdbd", "#252525")

asv_overlap <- MC_only_18S %>%
  filter(Supergroup != "Opisthokonta") %>% 
  mutate(Supergroup = ifelse(is.na(Supergroup), "Unknown Eukaryota", Supergroup),
         Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
         Phylum = ifelse(Phylum == "Alveolata_X", "Ellobiopsidae", Phylum),
         Supergroup = ifelse(Supergroup == "Alveolata", paste(Supergroup, Phylum, sep = "-"), Supergroup)) %>% 
  mutate(SUPERGROUP = case_when(
    Supergroup %in% alv ~ "Other Alveolata",
    Supergroup == "Eukaryota_X" ~ "Unknown Eukaryota",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Ochrophyta" ~ "Stramenopiles-Ochrophyta",
    Phylum == "Opalozoa" ~ "Stramenopiles-Opalozoa",
    Phylum == "Sagenista" ~ "Stramenopiles-Sagenista",
    Supergroup == "Protalveolata" ~ "Other Alveolata",
    TRUE ~ Supergroup
  )) %>% 
  # Taxa to supergroup
  mutate(SupergroupPhylum = SUPERGROUP) %>%
  # Average across replicates
  group_by(FeatureID, SAMPLE, MC, SUBSTRATE, SupergroupPhylum) %>% 
  ungroup() %>% 
  unite(SAMPLE, MC, SUBSTRATE, sep = " ", remove = FALSE) %>% 
  group_by(FeatureID, SupergroupPhylum, SAMPLE) %>% 
  summarise(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>% 
  mutate(SUPERGROUP_ORDER = factor(SupergroupPhylum, levels = all_taxa_order)) %>%
  distinct(FeatureID, SUPERGROUP_ORDER, SUM, SAMPLE, .keep_all = TRUE) %>% 
  group_by(FeatureID, SUPERGROUP_ORDER) %>% 
  summarise(SAMPLE = list(SAMPLE)) %>% 
  ggplot(aes(x = SAMPLE)) +
  geom_hline(yintercept = 200, linetype = "dashed", alpha = 0.6) +
  geom_bar(color = NA, width = 0.7, aes(fill = SUPERGROUP_ORDER)) +
  scale_x_upset(order_by = "freq", n_intersections = 25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Shared ASVs") +
  theme_linedraw() +
  theme(axis.text.y = element_text(color="black", size=8, face = "bold"),
        axis.text.x = element_text(color="black", size=8, face = "bold"),
        axis.title = element_text(color="black", size=8, face = "bold"),
        legend.text = element_text(color = "black", size = 6),
        plot.margin = margin(1, 1, 1, 5, "cm"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = all_taxa_color) +
  theme_combmatrix(
    combmatrix.panel.striped_background.color.two = "#878686",
    combmatrix.panel.point.size = 3.5
  )

asv_overlap <- MC_only_18S %>%
  filter(Supergroup != "Opisthokonta") %>% 
  mutate(Supergroup = ifelse(is.na(Supergroup), "Unknown Eukaryota", Supergroup),
         Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
         Phylum = ifelse(Phylum == "Alveolata_X", "Ellobiopsidae", Phylum),
         Supergroup = ifelse(Supergroup == "Alveolata", paste(Supergroup, Phylum, sep = "-"), Supergroup)) %>% 
  mutate(SUPERGROUP = case_when(
    Supergroup %in% alv ~ "Other Alveolata",
    Supergroup == "Eukaryota_X" ~ "Unknown Eukaryota",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Ochrophyta" ~ "Stramenopiles-Ochrophyta",
    Phylum == "Opalozoa" ~ "Stramenopiles-Opalozoa",
    Phylum == "Sagenista" ~ "Stramenopiles-Sagenista",
    Supergroup == "Protalveolata" ~ "Other Alveolata",
    TRUE ~ Supergroup
  )) %>% 
  # Taxa to supergroup
  mutate(SupergroupPhylum = SUPERGROUP) %>%
  # Average across replicates
  group_by(FeatureID, SAMPLE, MC, SUBSTRATE, SupergroupPhylum) %>% 
  ungroup() %>% 
  unite(SAMPLE, MC, SUBSTRATE, sep = " ", remove = FALSE) %>% 
  group_by(FeatureID, SupergroupPhylum, MC) %>% 
  summarise(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>% 
  mutate(SUPERGROUP_ORDER = factor(SupergroupPhylum, levels = all_taxa_order)) %>%
  distinct(FeatureID, SUPERGROUP_ORDER, SUM, MC, .keep_all = TRUE) %>% 
  group_by(FeatureID, SUPERGROUP_ORDER) %>% 
  summarise(MC = list(MC)) %>% 
  ggplot(aes(x = MC)) +
  geom_hline(yintercept = 200, linetype = "dashed", alpha = 0.6) +
  geom_bar(color = NA, width = 0.7, aes(fill = SUPERGROUP_ORDER)) +
  scale_x_upset(order_by = "freq", n_intersections = 25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Shared ASVs") +
  theme_linedraw() +
  theme(axis.text.y = element_text(color="black", size=8, face = "bold"),
        axis.text.x = element_text(color="black", size=8, face = "bold"),
        axis.title = element_text(color="black", size=8, face = "bold"),
        legend.text = element_text(color = "black", size = 6),
        plot.margin = margin(1, 1, 1, 5, "cm"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = all_taxa_color) +
  theme_combmatrix(
    combmatrix.panel.striped_background.color.two = "#878686",
    combmatrix.panel.point.size = 3.5
  )

asv_overlap

asv_overlap_2 <- MC_only_18S %>% 
  filter(SEQ_AVG > 10) %>%
  filter(Domain == "Eukaryota") %>% #select eukaryotes only
  filter(Supergroup != "Opisthokonta") %>% # remove multicellular metazoa
  mutate(Supergroup = ifelse(is.na(Supergroup), "Unknown Eukaryota", Supergroup),
         Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
         Phylum = ifelse(Phylum == "Alveolata_X", "Ellobiopsidae", Phylum),
         Supergroup = ifelse(Supergroup == "Alveolata", paste(Supergroup, Phylum, sep = "-"), Supergroup)) %>% 
  mutate(SUPERGROUP = case_when(
    Supergroup %in% alv ~ "Other Alveolata",
    Supergroup == "Eukaryota_X" ~ "Unknown Eukaryota",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Ochrophyta" ~ "Stramenopiles-Ochrophyta",
    Phylum == "Opalozoa" ~ "Stramenopiles-Opalozoa",
    Phylum == "Sagenista" ~ "Stramenopiles-Sagenista",
    TRUE ~ Supergroup
  )) %>% 
  # Taxa to supergroup
  mutate(SupergroupPhylum = SUPERGROUP) %>% #add modified "supergroup-phylum category"
  # Average across replicates
  group_by(FeatureID, SAMPLE, MC, SUBSTRATE, SupergroupPhylum, Taxon) %>% 
  summarise(AVG = mean(SEQ_AVG)) %>% 
  ungroup() %>% 
  group_by(SupergroupPhylum, Taxon, SAMPLE) %>% 
  summarise(SUM = sum(AVG)) %>%
  ungroup() %>%
  distinct(Taxon, SupergroupPhylum, SUM, SAMPLE, .keep_all = TRUE) %>% 
  group_by(SupergroupPhylum, Taxon) %>% 
  summarise(SAMPLE = list(SAMPLE)) %>% 
  ggplot(aes(x = SAMPLE)) +
  geom_bar(color = "black", width = 0.7, aes(fill = SupergroupPhylum)) +
  scale_x_upset(n_intersections = 25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Shared at taxonomic level") +
  theme_linedraw() +
  theme(axis.text.y = element_text(color="black", size=8, face = "bold"),
        axis.text.x = element_text(color="black", size=8, face = "bold"),
        axis.title = element_text(color="black", size=8, face = "bold"),
        legend.text = element_text(color = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 5, "cm")) + 
  scale_fill_manual(values = all_taxa_color)

asv_overlap_2 <- MC_only_18S %>% 
  filter(SEQ_AVG > 10) %>%
  filter(Domain == "Eukaryota") %>% #select eukaryotes only
  filter(Supergroup != "Opisthokonta") %>% # remove multicellular metazoa
  mutate(Supergroup = ifelse(is.na(Supergroup), "Unknown Eukaryota", Supergroup),
         Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
         Phylum = ifelse(Phylum == "Alveolata_X", "Ellobiopsidae", Phylum),
         Supergroup = ifelse(Supergroup == "Alveolata", paste(Supergroup, Phylum, sep = "-"), Supergroup)) %>% 
  mutate(SUPERGROUP = case_when(
    Supergroup %in% alv ~ "Other Alveolata",
    Supergroup == "Eukaryota_X" ~ "Unknown Eukaryota",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Ochrophyta" ~ "Stramenopiles-Ochrophyta",
    Phylum == "Opalozoa" ~ "Stramenopiles-Opalozoa",
    Phylum == "Sagenista" ~ "Stramenopiles-Sagenista",
    TRUE ~ Supergroup
  )) %>% 
  # Taxa to supergroup
  mutate(SupergroupPhylum = SUPERGROUP) %>% #add modified "supergroup-phylum category"
  # Average across replicates
  group_by(FeatureID, SAMPLE, MC, SUBSTRATE, SupergroupPhylum, Taxon) %>% 
  summarise(AVG = mean(SEQ_AVG)) %>% 
  ungroup() %>% 
  group_by(SupergroupPhylum, Taxon, MC) %>% 
  summarise(SUM = sum(AVG)) %>%
  ungroup() %>%
  distinct(Taxon, SupergroupPhylum, SUM, MC, .keep_all = TRUE) %>% 
  group_by(SupergroupPhylum, Taxon) %>% 
  summarise(MC = list(MC)) %>% 
  ggplot(aes(x = MC)) +
  geom_bar(color = "black", width = 0.7, aes(fill = SupergroupPhylum)) +
  scale_x_upset(n_intersections = 25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Shared at taxonomic level") +
  theme_linedraw() +
  theme(axis.text.y = element_text(color="black", size=8, face = "bold"),
        axis.text.x = element_text(color="black", size=8, face = "bold"),
        axis.title = element_text(color="black", size=8, face = "bold"),
        legend.text = element_text(color = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 5, "cm")) + 
  scale_fill_manual(values = all_taxa_color)

asv_overlap_2 <- MC_only_18S %>% 
  filter(SEQ_AVG > 10) %>%
  filter(Domain == "Eukaryota") %>% #select eukaryotes only
  filter(Supergroup != "Opisthokonta") %>% # remove multicellular metazoa
  mutate(Supergroup = ifelse(is.na(Supergroup), "Unknown Eukaryota", Supergroup),
         Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
         Phylum = ifelse(Phylum == "Alveolata_X", "Ellobiopsidae", Phylum),
         Supergroup = ifelse(Supergroup == "Alveolata", paste(Supergroup, Phylum, sep = "-"), Supergroup)) %>% 
  mutate(SUPERGROUP = case_when(
    Supergroup %in% alv ~ "Other Alveolata",
    Supergroup == "Eukaryota_X" ~ "Unknown Eukaryota",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Ochrophyta" ~ "Stramenopiles-Ochrophyta",
    Phylum == "Opalozoa" ~ "Stramenopiles-Opalozoa",
    Phylum == "Sagenista" ~ "Stramenopiles-Sagenista",
    TRUE ~ Supergroup
  )) %>% 
  # Taxa to supergroup
  mutate(SupergroupPhylum = SUPERGROUP) %>% #add modified "supergroup-phylum category"
  # Average across replicates
  group_by(FeatureID, SAMPLE, MC, SUBSTRATE, SupergroupPhylum, Taxon) %>% 
  summarise(AVG = mean(SEQ_AVG)) %>% 
  ungroup() %>% 
  group_by(SupergroupPhylum, Taxon, SUBSTRATE) %>% 
  summarise(SUM = sum(AVG)) %>%
  ungroup() %>%
  distinct(Taxon, SupergroupPhylum, SUM, SUBSTRATE, .keep_all = TRUE) %>% 
  group_by(SupergroupPhylum, Taxon) %>% 
  summarise(SUBSTRATE = list(SUBSTRATE)) %>% 
  ggplot(aes(x = SUBSTRATE)) +
  geom_bar(color = "black", width = 0.7, aes(fill = SupergroupPhylum)) +
  scale_x_upset(n_intersections = 25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Shared at taxonomic level") +
  theme_linedraw() +
  theme(axis.text.y = element_text(color="black", size=8, face = "bold"),
        axis.text.x = element_text(color="black", size=8, face = "bold"),
        axis.title = element_text(color="black", size=8, face = "bold"),
        legend.text = element_text(color = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 5, "cm")) + 
  scale_fill_manual(values = all_taxa_color)
asv_overlap_2

shared_across <- df_18s_s_v_vs %>% 
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta") %>%
  ungroup() %>%
  filter(SEQ_AVG > 1) %>%
  select(FeatureID, SAMPLE, DISTRIBUTION, COUNT) %>% 
  select(FeatureID, DISTRIBUTION, COUNT) %>% 
  pivot_wider(names_from = DISTRIBUTION, values_from = COUNT, values_fn = sum) 
  select(FeatureID, DISTRIBUTION)

  
unique(for_subs_table$SUBS_DIST)
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

library(dplyr)

asv_counts_table <- MC_only_18S %>%
  group_by(FeatureID, SUBSTRATE) %>%
  summarise(PRESENCE = max(PRESENCE), .groups = "drop") %>%  # one row per ASV–substrate
  filter(PRESENCE == 1) %>%
  group_by(FeatureID) %>%
  summarise(SUBS_DIST = paste(sort(unique(SUBSTRATE)), collapse = " and "), .groups = "drop") %>%
  mutate(SUBS_DIST = case_when(
    SUBS_DIST == "Quartz and Riftia and Shell" ~ "Across all",
    TRUE ~ SUBS_DIST
  )) %>%
  distinct(FeatureID, SUBS_DIST) %>%   # ensure unique pairs
  count(SUBS_DIST, name = "Num_distinct_ASVs") %>%
  arrange(desc(Num_distinct_ASVs))

asv_counts_table$SUBS_DIST <- factor(
  asv_counts_table$SUBS_DIST,
  levels = c("Across all",
             "Quartz",
             "Riftia",
             "Shell",
             "Quartz and Riftia",
             "Quartz and Shell",
             "Riftia and Shell")
)
length(unique(for_subs_table$FeatureID))

subs_dist_filt <- for_subs_table %>%
  filter(SEQ_AVG > 99)

more_than_one_18S <- for_subs_table %>%
  filter(subs_combo != "Quartz" & subs_combo != "Riftia" & subs_combo != "Shell") %>%
  filter(SEQ_AVG > 99)
length(unique(more_than_one_18S$FeatureID))
unique(more_than_one_18S$Phylum)

df_more_than_one <- more_than_one_18S %>%
  select(-SEQ_AVG) %>%
  distinct() 

asv_counts_table <- asv_counts_table %>%
  arrange(SUBS_DIST)
asv_counts_table <- for_subs_table %>%
  count(SUBS_DIST, name = "ASVs") %>%
  arrange(desc(ASVs))

SUBS_TABLE <- for_subs_table %>%
  bind_rows(summary_16s, summary_18s) %>%
  select(Domain, Total_sequences, Total_ASVs, Total_Unassigned)

#### ####
ASVs_MC_18s <- df_18s_all %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!(is.na(SUBSTRATE)))
length(ASVs_MC_18s$FeatureID)
length(unique(ASVs_MC_18s$FeatureID))

# ASVs on substrate only (MC only - no B/V)
df_18s_subs <- mc_tmp_18 %>%
  left_join(key_dist_18) %>%
  filter(DISTRIBUTION == "Substrate only")
length(df_18s_subs$FeatureID)
length(unique(df_18s_subs$FeatureID))

# MC only ASVs on more than one substrate
asvs_subs_only_18S <- df_18s_subs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)

format_18s_subonly <- asvs_subs_only_18S %>%
  select(FeatureID, SEQ_AVG, SUBSTRATE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SUBSTRATE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

format_18s_subonly_assigned <- format_18s_subonly %>%
  mutate(Distribution = case_when(
    (is.na(Riftia) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz only", 
    (is.na(Riftia) & !(is.na(Shell)) & !(is.na(Quartz))) ~ "Quartz & Shell only", 
    (is.na(Riftia) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell only", 
    (!(is.na(Riftia)) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz & Riftia only",
    (!(is.na(Riftia)) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell & Riftia only", 
    (!(is.na(Riftia)) & is.na(Shell) & is.na(Quartz)) ~ "Riftia only", 
    (!(is.na(Riftia)) & !(is.na(Shell)) & !(is.na(Quartz)) ~ "All")
  ))

key_sub_dist_18 <- format_18s_subonly_assigned %>%
  select(FeatureID, Distribution) %>%
  distinct()
unique(key_sub_dist_18$Distribution)

more_than_one <- key_sub_dist_18 %>%
  filter(!(Distribution == "Quartz only" | Distribution == "Riftia only" | Distribution == "Shell only"))
length(more_than_one$FeatureID)

# MC only ASVs on shell only
shell_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Shell only")
length(shell_only_18s$FeatureID)

# MC only ASVs on riftia only
riftia_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Riftia only")
length(riftia_only_18s$FeatureID)

# MC only ASVs on quartz only
quartz_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Quartz only")
length(quartz_only_18s$FeatureID)

# MC only ASVs on all substrates 
onall_sub_18s <- key_sub_dist_18 %>%
  filter(Distribution == "All")
length(unique(onall_sub_18s$FeatureID))