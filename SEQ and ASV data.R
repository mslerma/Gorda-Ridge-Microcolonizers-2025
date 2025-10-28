load(file = "C:/Users/Madeleine/OneDrive - Texas A&M University/Research lab/microcolonizers/GR_Microcolonizer_data.RData")

df_18S_data <- asv_wtax_18 %>%
  select(-SAMPLE) |> 
  filter(SAMPLETYPE == "Microcolonizer" | VENT == "Mt Edwards" | VENT == "Deep seawater") |> 
  filter(Domain == "Eukaryota") |> 
  mutate(SUBSTRATE = str_replace_all(Substrate, "z2", "z")) |> 
  mutate(SAMPLE = case_when(
    SAMPLETYPE == "Microcolonizer" ~ paste(MC, SUBSTRATE, sep = "-"),
    VENT == "Deep seawater" ~ "Background seawater",
    VENT == "Mt Edwards" ~ "Mt Edwards diffuse fluid"
  )) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(SEQUENCE_COUNT > 0) %>%
  select(-c("VENT", "COORDINATES", "SITE", "SAMPLEID", "SAMPLETYPE", "DEPTH", "temp", "pH", "percseawater", "mg", "h2", "h2s", "ch4", "ProkConc", "dataset", "Substrate")) %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta" & Phylum != "NA")


mc_18S_df_filt <- mc_18S_df %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta" & Phylum != "NA")
unique(mc_18S_df_filt$Phylum)  
df_18s_all_filt <- df_18s_all %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta" & Phylum != "NA")
unique(df_18s_all_filt$Phylum)  
tot_ASVs_all_samp <- unique(df_18S_data$FeatureID)

df_18S_MC_data <- df_18S_data %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid")
tot_ASVs_MC <- unique(df_18S_MC_data$FeatureID)

ASVs_MC_only <- df_18s_all_filt %>%
filter(DISTRIBUTION == "Substrate only")
length(unique(ASVs_MC_only$FeatureID))

asv_subs_sets <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarise(subs_combo = paste(sort(unique(SUBSTRATE)), collapse = ","), .groups = "drop")

key_ASV_taxons <- MC_only_18S %>%
  ungroup() %>%
  select(FeatureID, Supergroup, Phylum, Class, Order, Family, Genus, Species, SEQ_AVG) %>%
  distinct()

joined_ASVs_SUBS_only <- asv_subs_sets %>%
  left_join(key_ASV_taxons)

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

###SEQs
# total
total_seq_all_MCs <- sum(df_18S_data$SEQUENCE_COUNT)
 
# MC1
MC1_18S_data <- df_18S_data %>%
  filter(MC == "MC1")
total_seq_MC1 <- sum(MC1_18S_data$SEQUENCE_COUNT)

# MC2
MC2_18S_data <- df_18S_data %>%
  filter(MC == "MC2")
total_seq_MC2 <- sum(MC2_18S_data$SEQUENCE_COUNT)

# MC3
MC3_18S_data <- df_18S_data %>%
  filter(MC == "MC3")
total_seq_MC3 <- sum(MC3_18S_data$SEQUENCE_COUNT)

# MC4
MC4_18S_data <- df_18S_data %>%
  filter(MC == "MC4")
total_seq_MC4 <- sum(MC4_18S_data$SEQUENCE_COUNT)

# MC5
MC5_18S_data <- df_18S_data %>%
  filter(MC == "MC5")
total_seq_MC5 <- sum(MC5_18S_data$SEQUENCE_COUNT)

# MC6
MC6_18S_data <- df_18S_data %>%
  filter(MC == "MC6")
total_seq_MC6 <- sum(MC6_18S_data$SEQUENCE_COUNT)
unique(joined_ASVs_SUBS_only$subs_combo)
mult_subs_ASVs <- joined_ASVs_SUBS_only %>%
  filter(subs_combo != "Quartz" & subs_combo != "Riftia" & subs_combo != "Shell")
length(unique(mult_subs_ASVs$FeatureID))

ASVs_shell <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Shell")
length(unique(ASVs_shell$FeatureID))

ASVs_riftia <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Riftia")
length(unique(ASVs_riftia$FeatureID))

ASVs_quartz <- joined_ASVs_SUBS_only %>%
  filter(subs_combo == "Quartz")
length(unique(ASVs_quartz$FeatureID))

###ASVs
# total
total_ASVs_18S <- length(unique(df_18S_data$FeatureID))

# MC1
total_ASVs_MC1 <- length(unique(MC1_18S_data$FeatureID))

# MC2
total_ASVs_MC2 <- length(unique(MC2_18S_data$FeatureID))

# MC3
total_ASVs_MC3 <- length(unique(MC3_18S_data$FeatureID))

# MC4
total_ASVs_MC4 <- length(unique(MC4_18S_data$FeatureID))

# MC5
total_ASVs_MC5 <- length(unique(MC5_18S_data$FeatureID))

# MC6
total_ASVs_MC6 <- length(unique(MC6_18S_data$FeatureID))

##Table
summary_SEQs <- df_18S_data %>%
  summarise(Total_sequences = total_seq_all_MCs,
            MC1 = total_seq_MC1,
            MC2 = total_seq_MC2,
            MC3 = total_seq_MC3,
            MC4 = total_seq_MC4,
            MC5 = total_seq_MC5, 
            MC6 = total_seq_MC6)
summary_ASVs <- df_18S_data %>%
  summarise(Total_ASVs = total_ASVs_18S, 
            MC1 = total_ASVs_MC1,
            MC2 = total_ASVs_MC2,
            MC3 = total_ASVs_MC3,
            MC4 = total_ASVs_MC4,
            MC5 = total_ASVs_MC5, 
            MC6 = total_ASVs_MC6)


summary_table <- bind_rows(summary_SEQs, summary_ASVs) 
  select(-c())
rows(summary_16s, summary_18s) %>%
  select(Domain, Total_sequences, Total_ASVs, Total_Unassigned)

