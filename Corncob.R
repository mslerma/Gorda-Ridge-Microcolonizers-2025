#### Indicator species #### 

#### Corncob ####
# 4/17/25

## recruitment vs nonrecruitment (enrichment vs not)
phy <- requireNamespace("phyloseq", quietly = TRUE) == TRUE
remotes::install_github("statdivlab/corncob")

library(corncob)
library(phyloseq)
library(tidyr)
library(dplyr)

unique(sep_18s_vs$SAMPLE)


#####
# Make ASV table (wide format, samples in rows, ASVs in cols)
cc_asv_table <- sep_18s_vs %>%
  ungroup() %>%
  select(SAMPLE, FeatureID, SEQ_AVG) %>%
  pivot_wider(names_from = FeatureID, values_from = SEQ_AVG, values_fill = 0) %>%
  column_to_rownames("SAMPLE")

# Sample metadata (one row per sample)
sep_18s_vs <- sep_18s_vs %>%
  mutate(Time = case_when(
    MC == "MC1" ~ 7,
    MC == "MC2" ~ 6,
    MC == "MC3" ~ 6,
    MC == "MC4" ~ 7,
    MC == "MC5" ~ 8,
    MC == "MC6" ~ 8, 
    MC == "Mt Edwards diffuse fluid" ~ 0
  )) %>%
  mutate(Temperature = case_when(
    MC == "MC1" ~ 8.5,
    MC == "MC2" ~ 40.5,
    MC == "MC3" ~ 16.58,
    MC == "MC4" ~ 8.23,
    MC == "MC5" ~ 9.29,
    MC == "MC6" ~ 7.18, 
    MC == "Mt Edwards diffuse fluid" ~ 13
  )) %>%
  mutate(SUBSTRATE = ifelse(is.na(SUBSTRATE), "None", SUBSTRATE))

sep_18s_vs <- sep_18s_vs %>%
  mutate(SAMPLETYPE_BIN = case_when(
    MC == "MC1" ~ "non-vent",
    MC == "MC2" ~ "non-vent",
    MC == "MC3" ~ "non-vent",
    MC == "MC4" ~ "non-vent",
    MC == "MC5" ~ "non-vent",
    MC == "MC6" ~ "non-vent",
    MC == "Mt Edwards diffuse fluid" ~ "vent"
  ))
  
cc_env_data <- sep_18s_vs %>%
  ungroup() %>%
  select(SAMPLE, MC, SUBSTRATE, SAMPLETYPE_BIN) %>%
  distinct()

# Convert to phyloseq objects
cc_otu <- otu_table(as.matrix(cc_asv_table[ , -1]), taxa_are_rows = FALSE)
rownames(cc_otu) <- cc_asv_table$SAMPLE

cc_sdata <- sample_data(cc_env_data)
rownames(cc_sdata) <- cc_env_data$SAMPLE

####
cc_tax_table <- sep_18s_vs %>%
  ungroup() %>%
  select(FeatureID, Domain, Supergroup, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct()
view(cc_tax_table)
# Convert tibble -> dataframe
cc_tax_table <- as.data.frame(cc_tax_table)

# Set FeatureID as rownames
rownames(cc_tax_table) <- cc_tax_table$FeatureID

# Drop FeatureID column now that it’s rownames
cc_tax_table <- cc_tax_table[ , -which(names(cc_tax_table) == "FeatureID")]

# Convert to tax_table object
cc_tax <- tax_table(as.matrix(cc_tax_table))
view(cc_tax)
view(cc_otu)
cc_otu <- otu_table(data.matrix(cc_asv_table[ , ]), taxa_are_rows = FALSE)
rownames(cc_otu) <- cc_asv_table$SAMPLE
view(cc_sdata)
cc_phyloseq <- phyloseq(cc_otu, cc_tax, cc_sdata)

####
unique(sample_data(cc_phyloseq)$SAMPLETYPE_BIN)
grmc_super <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>%
  tax_glom("Supergroup")

grmc_phy <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Phylum")

grmc_class <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Class")

grmc_order <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Order")

grmc_fam <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Family")

grmc_genera <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Genus")

grmc_spp <- cc_phyloseq %>%
  phyloseq::subset_samples(SAMPLETYPE_BIN %in% c("non-vent", "vent")) %>%
  phyloseq::subset_taxa(Domain == "Eukaryota") %>% 
  tax_glom("Species")

sample_data(grmc_order)

####
# compare between MCs, between substrates
corncob_grmc <- function(df_in){
    ## SAMPLETYPE_BIN specifically compares vent to non-vent.
  da_analysis_output <- differentialTest(formula = ~ SAMPLETYPE_BIN,
                                         phi.formula = ~ SAMPLETYPE_BIN,
                                         formula_null = ~ 1,
                                         phi.formula_null = ~ SAMPLETYPE_BIN,
                                         test = "Wald", boot = FALSE,
                                         data = df_in,
                                         fdr_cutoff = 0.05)
  # da_analysis_output <- differentialTest(formula = ~ SAMPLE,
  #                               phi.formula = ~ SAMPLE,
  #                               formula_null = ~ 1,
  #                               phi.formula_null = ~ SAMPLE,
  #                               test = "Wald", boot = FALSE,
  #                               data = df_in,
  #                               fdr_cutoff = 0.05)
  #
  # da_analysis_output <- differentialTest(formula = ~ MC + SUBSTRATE,
  #                               phi.formula = ~ MC + SUBSTRATE,
  #                               formula_null = ~ 1,
  #                               phi.formula_null = ~ MC + SUBSTRATE,
  #                               test = "Wald", boot = FALSE,
  #                               data = df_in,
  #                               fdr_cutoff = 0.05)
  
  list_ofsig <- as.character(da_analysis_output$significant_taxa)
  total_number <- length(list_ofsig)[1]
  #
  sig_taxa_names <- as.data.frame(tax_table(df_in)) %>% 
    rownames_to_column(var = "FEATURE") %>% 
    filter(FEATURE %in% list_ofsig) %>% 
    rownames_to_column(var = "NUMBER")
  #
  for(var in 1:total_number){
    out_0 <- data.frame(da_analysis_output$significant_models[[var]]$coefficients) %>% 
      add_column(NUMBER = as.character(var))
    cat("extracted # ",var, "\n")
    if (!exists("extracted_coef")){
      extracted_coef <- out_0 # create the final count table
    } else {
      extracted_coef <- rbind(extracted_coef, out_0)
    }
    rm(out_0) # remove excess df
  }
  output_full <- extracted_coef %>% 
    rownames_to_column(var = "variable") %>%
    filter(grepl("mu.", variable)) %>% 
    left_join(sig_taxa_names, by = c("NUMBER" = "NUMBER")) %>% 
    mutate(VARIABLE = str_remove(variable, "[:digit:]+")) %>% select(-variable) #%>% 
  # pivot_wider(names_from = VARIABLE, names_glue = "{VARIABLE}_{.value}", values_from = c(Estimate, Std..Error, t.value, Pr...t..))
  rm(extracted_coef)
  return(output_full)}

#####
cc_phyloseq
sample_names(cc_phyloseq) %>% head()
taxa_names(cc_phyloseq) %>% head()
sample_variables(cc_phyloseq)


# Fit beta-binomial regression on a specific taxon (e.g., genus, ASV)

"2133fc07124882f822e01d15c48278fb" %in% taxa_names(cc_phyloseq)
otu_table(cc_phyloseq)[ , "2133fc07124882f822e01d15c48278fb"]
dim(otu_table(cc_phyloseq))
ntaxa(cc_phyloseq)   # number of taxa
nsamples(cc_phyloseq)
rownames(cc_otu)
subset_taxa(cc_phyloseq, taxa_names(cc_phyloseq) == "2133fc07124882f822e01d15c48278fb")
sum(otu_table(cc_phyloseq)["2133fc07124882f822e01d15c48278fb", ])
sum(as(otu_table(cc_phyloseq), "matrix")["2133fc07124882f822e01d15c48278fb", ])
sum(as(otu_table(cc_phyloseq), "matrix")[, "2133fc07124882f822e01d15c48278fb"])
taxa_are_rows(cc_phyloseq)
target_taxon <- "2133fc07124882f822e01d15c48278fb"
phylo_sub <- subset_taxa(cc_phyloseq, taxa_names(cc_phyloseq) == target_taxon)
bb_mod <- bbdml(
  formula = Abundance ~ SUBSTRATE + Time + Temperature,
  phi.formula = ~ SUBSTRATE + Time + Temperature,
  data = phylo_sub,
  taxa = target_taxon
)
str(sample_data(cc_phyloseq))


bb_mod <- bbdml(
  formula = Abundance ~ SUBSTRATE + Time + Temperature,
  phi.formula = ~ SUBSTRATE + Time + Temperature,
  data = cc_phyloseq,
  taxa = "2133fc07124882f822e01d15c48278fb")

summary(bb_mod)

da_test <- differentialTest(
  formula = ~ SUBSTRATE + Time + Temperature,
  phi.formula = ~ SUBSTRATE + Time + Temperature,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  boot = FALSE,
  data = cc_phyloseq,
  fdr_cutoff = 0.05
)

da_test$significant_taxa


da_substrate <- differentialTest(
  formula = ~ SUBSTRATE,
  phi.formula = ~ SUBSTRATE,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  boot = FALSE,
  data = cc_phyloseq,
  fdr_cutoff = 0.05
)

plot(da_substrate)   # visualize significant taxa
####
library(phyloseq)
table(sample_data(cc_phyloseq)$SUBSTRATE)
table(sample_data(cc_phyloseq)$Time)         # if Time is categorical
summary(sample_data(cc_phyloseq)$Temperature) # if numeric


####
vs_18s_cc <- sep_18s_vs %>% ungroup %>% 
  filter(SEQ_AVG > 0) %>%
  mutate(SAMPLE = case_when(`SAMPLE` == "Mt Edwards diffuse fluid" ~ "Vent", .default = `SAMPLE`)) %>%
  select(FeatureID, SAMPLE, SEQ_AVG) %>%
  pivot_wider(names_from = SAMPLE, values_from = SEQ_AVG, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

######
samplenames_order <- c("MC1-Quartz", "MC1-Riftia", "MC1-Shell",
                       "MC2-Quartz", "MC3-Riftia", "MC3-Shell",               
                       "MC4-Quartz", "MC4-Shell", "MC5-Quartz",              
                       "MC5-Riftia", "MC5-Shell", "MC6-Quartz",              
                       "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid")
samplenames_color <- c("#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93", "#2a2b47")

names(samplenames_color) <- samplenames_order
#####
