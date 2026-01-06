########                       MARCO-BOLO Data Analysis Challenge                       ########
########                     12S, 16S, COI fish metabarcoding data                      ########

########         The influence of bioinformatic choices on fish diversity patterns      ########
##########         script by Emilie Boulanger on R version 4.2.3 (2023-03-15)         ##########

## load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(hillR)

## load data

# if working on csv file
challengedata_12s16sCOI <- read.csv(file = "../data/challengedata_12s16sCOI.csv", row.names = 1) %>% 
  mutate(readAbundance = as.numeric(readAbundance)) %>% 
  rename(submitterID = submitterID_anon) %>% 
  filter(!(pipeline %in% c("Anacapa40", "Anacapa50", "Anacapa60", "Anacapa80", "Anacapa90"))) # Only keep thresholds 95 and 100 for the Anacapa pipeline

# if working on rds file
challengedata_12s16sCOI <- readRDS(file = "../data/challengedata_12s16sCOI.rds") %>% 
  mutate(readAbundance = as.numeric(readAbundance)) %>% 
  rename(submitterID = submitterID_anon) %>% 
  filter(!(pipeline %in% c("Anacapa40", "Anacapa50", "Anacapa60", "Anacapa80", "Anacapa90"))) # Only keep thresholds 95 and 100 for the Anacapa pipeline

# aquarium species lists
aquarium_residents <- read.csv("../data/aquarium_12s16sCOI_residents.csv")
aquarium_feed      <- read.csv("../data/aquarium_12s16sCOI_feed.csv")
aquarium_species <- c(aquarium_residents$Species, aquarium_feed$Potential_Contaminant_Species)

# metadata file built as I go
metadata_12s16sCOI <- read.csv("../data/metadata_12s16sCOI.csv") 

#### Calculate read abundances and summary metrics ####
##### Overall #####

# Read Abundance
total_readAb <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "species") %>%  # select one taxonomic level to count every sequence only once
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_read_total = sum(readAbundance), .groups = "drop") %>% 
  mutate(n_read_total = as.numeric(n_read_total))

sample_readAb <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "species") %>%  # select one taxonomic level to count every sequence only once
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_read = sum(readAbundance), .groups = "drop") %>% 
  mutate(n_read = as.numeric(n_read)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_read,
    names_prefix = "n_read_")

# Number ASVs or OTUs
total_seq <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "species") %>%  
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_seqID_total = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_total = as.numeric(n_seqID_total))

sample_seq <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "species") %>%  
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_seqID = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID = as.numeric(n_seqID)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_seqID,
    names_prefix = "n_seqID_")

results <- total_readAb %>% 
  full_join(sample_readAb, by = c("submitterID", "marker", "referenceDB", "pipeline")) %>% 
  full_join(total_seq, by = c("submitterID", "marker", "referenceDB", "pipeline")) %>% 
  full_join(sample_seq, by = c("submitterID", "marker", "referenceDB", "pipeline"))
# add metadata
results <- results %>% 
  mutate(submissionID = paste(submitterID, marker, referenceDB, pipeline, sep = "_")) %>% 
  left_join(metadata_12s16sCOI)

# Some quick boxplots
ggplot(data = results, aes(x= submitterID, y = n_read_total)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_grid(marker ~.) +
  ggtitle("Total read abundance")

ggplot(data = results, aes(x= submitterID, y = n_seqID_total)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_grid(marker ~.) +
  ggtitle("Total number of ASVs or OTUs identified")

##### Target ASV/OTU ####
# get the total read abundance looking only at fish identifications
# so get a list of all ASVs assigned to class c("Actinopteri", "Chondrichthyes", "Elasmobranchii", "Teleostei"))
# or by filtering for all Chordata, but then also get Mammalia etc. 
challengedata_12s16sCOI %>% 
  pivot_wider(names_from = taxonomicLevel, values_from = assignment) %>%
  filter(phylum == "Chordata") %>% pull(class) %>% unique() #contains Actinopteri, Chondrichthyes, Elasmobranchii, Teleostei, NA, Mammalia, Holostei, Ascidiacea 
challengedata_12s16sCOI %>% 
  pivot_wider(names_from = taxonomicLevel, values_from = assignment) %>%
  filter(class %in% c("Mammalia", "Holostei", "Ascidiacea")) %>% pull(species) %>% unique() #some weird IDs. Keep out of target here.

target_ASV <- challengedata_12s16sCOI %>% 
  filter(assignment %in% c("Actinopteri", "Chondrichthyes", "Elasmobranchii", "Teleostei")) %>% 
  pull(uniqueID) #get a list of all the ASVs assigned to fish, sharks and rays

# now get read abundances and number of ASVs for these target ASVs ----
total_readAb_target <- challengedata_12s16sCOI %>% 
  filter(uniqueID %in% target_ASV) %>% 
  filter(taxonomicLevel == "species") %>%  # select one taxonomic level to count every sequence only once
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_read_total_target = sum(readAbundance), .groups = "drop") %>% 
  mutate(n_read_total_target = as.numeric(n_read_total_target))

sample_readAb_target <- challengedata_12s16sCOI %>% 
  filter(uniqueID %in% target_ASV) %>% 
  filter(taxonomicLevel == "species") %>%  # select one taxonomic level to count every sequence only once
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_read_target = sum(readAbundance), .groups = "drop") %>% 
  mutate(n_read_target = as.numeric(n_read_target)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_read_target,
    names_prefix = "n_read_target_")

total_seq_target <- challengedata_12s16sCOI %>% 
  filter(uniqueID %in% target_ASV) %>% 
  filter(taxonomicLevel == "species") %>%  
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_seqID_total_target = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_total_target = as.numeric(n_seqID_total_target))

sample_seq_target <- challengedata_12s16sCOI %>% 
  filter(uniqueID %in% target_ASV) %>% 
  filter(taxonomicLevel == "species") %>%  
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_seqID_target = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_target = as.numeric(n_seqID_target)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_seqID_target,
    names_prefix = "n_seqID_target_")

#add to results
results <- results %>% 
  left_join(total_readAb_target, by = c("submitterID", "marker", "referenceDB", "pipeline")) %>% 
  left_join(sample_readAb_target, by = c("submitterID", "marker", "referenceDB", "pipeline")) %>% 
  left_join(total_seq_target, by = c("submitterID", "marker", "referenceDB", "pipeline")) %>% 
  left_join(sample_seq_target, by = c("submitterID", "marker", "referenceDB", "pipeline"))


#boxplots
ggplot(data = results, aes(x= submitterID, y = n_read_total_target)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_grid(marker ~.) +
  ggtitle("Total read abundance for target ASV/OTU")

ggplot(data = results, aes(x= submitterID, y = n_seqID_total_target)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_grid(marker ~.) +
  ggtitle("Total number of ASVs or OTUs identified for target ASV/OTU")

ggplot(data = total_readAb_target, aes(x= submitterID, y = n_read_total_target)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_grid(marker ~.) +
  ggtitle("Total read abundance target ASVs & OTUs")

##### By Assignment level #####
# calculate number of ASV/OTU assigned to species/genus/family level
# and append to (previous) results #EDIT: NEED TO FIX THIS

total_seq_species <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "species") %>%  
  filter(readAbundance > 0) %>% 
  filter(!is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_seqID_species_total = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_species_total = as.numeric(n_seqID_species_total)) 

total_seq_qenus <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "genus") %>%  
  filter(readAbundance > 0) %>% 
  filter(!is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_seqID_genus_total = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_genus_total = as.numeric(n_seqID_genus_total)) 

total_seq_family <- challengedata_12s16sCOI %>% 
  filter(taxonomicLevel == "family") %>%  
  filter(readAbundance > 0) %>% 
  filter(!is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline) %>% 
  summarise(n_seqID_family_total = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_family_total = as.numeric(n_seqID_family_total)) %>% 
  left_join(results, ., by = c("submitterID", "marker", "referenceDB", "pipeline"))

colnames(results)

# same but per sample
# number of ASV/OTU assigned to species level per sample
sample_seq_species <- challengedata_12s16sCOI %>%
  filter(taxonomicLevel == "species") %>%  
  filter(readAbundance > 0, !is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_seqID_species_sample = n_distinct(seqID), .groups = "drop") %>%
  mutate(n_seqID_species_sample = as.numeric(n_seqID_species_sample)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_seqID_species_sample,
    names_prefix = "n_seqID_species_")


# number of species per sample
sample_specnumber <- challengedata_12s16sCOI %>%
  filter(taxonomicLevel == "species") %>%  
  filter(readAbundance > 0, !is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_species_sample = n_distinct(assignment), .groups = "drop") %>%
  mutate(n_species_sample = as.numeric(n_species_sample)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_species_sample,
    names_prefix = "n_species_") %>% 
  left_join(results, ., by = c("submitterID", "marker", "referenceDB", "pipeline"))


# number of target species per sample
sample_target_specnumber <- challengedata_12s16sCOI %>%
  filter(uniqueID %in% target_ASV) %>% 
  filter(taxonomicLevel == "species") %>%  
  filter(readAbundance > 0, !is.na(assignment)) %>% 
  group_by(submitterID, marker, referenceDB, pipeline, sampleID) %>% 
  summarise(n_targetspecies_sample = n_distinct(assignment), .groups = "drop") %>%
  mutate(n_targetspecies_sample = as.numeric(n_targetspecies_sample)) %>% 
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline),
    names_from = sampleID,
    values_from = n_targetspecies_sample,
    names_prefix = "n_targetspecies_") %>% 
  left_join(results, ., by = c("submitterID", "marker", "referenceDB", "pipeline"))

#### Diversity numbers: Species richness, Shannon, Simpson, Hill numbers ----

all_samples <- challengedata_12s16sCOI %>%
  select(submissionID, submitterID, marker, referenceDB, pipeline, sampleID) %>%
  distinct()

# species level
diversity_numbers_species <- challengedata_12s16sCOI %>%
  filter(uniqueID %in% target_ASV,
         taxonomicLevel == "species",
         !is.na(assignment)) %>%
  group_by(submitterID, marker, referenceDB, pipeline, sampleID, assignment) %>%
  # pool duplicate species assignments 
  summarise(abundance = sum(readAbundance), .groups = "drop") %>%
  # widen to produce a species x site matrix for diversity calculations
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline, sampleID),
    names_from = assignment,
    values_from = abundance,
    values_fill = 0) %>%
  ungroup() %>%     # now keep working on the table as a single group
  {                 # create a code block within the %>% 
    # extract numeric abundance matrix that both vegan and hillR functions can use directly
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    # compute diversity metrics directly on the matrix
    bind_cols(                                                         # merge the sample metadata and diversity metrics into a single data.frame
      select(., submitterID, marker, referenceDB, pipeline, sampleID), # keep sample metadata to combine with results later
      data.frame(                                                      # calculate the diversity metrics with vegan and hillR
        specnumber = specnumber(abundance_mat),
        shannon  = diversity(abundance_mat, index = "shannon"),
        simpson  = diversity(abundance_mat, index = "simpson"),
        hillq0   = hill_taxa(abundance_mat, q = 0),
        hillq1   = hill_taxa(abundance_mat, q = 1),
        hillq15  = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% # close code block
  left_join(all_samples, ., 
            by = c("submitterID", "marker", "referenceDB", "pipeline", "sampleID"))

# genus level
diversity_numbers_genus <- challengedata_12s16sCOI %>%
  filter(uniqueID %in% target_ASV,
         taxonomicLevel == "genus",
         !is.na(assignment)) %>%
  group_by(submitterID, marker, referenceDB, pipeline, sampleID, assignment) %>%
  summarise(abundance = sum(readAbundance), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline, sampleID),
    names_from = assignment,
    values_from = abundance,
    values_fill = 0) %>%
  ungroup() %>%     
  {               
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    bind_cols(                                                         
      select(., submitterID, marker, referenceDB, pipeline, sampleID), 
      data.frame(                                                      
        specnumber = specnumber(abundance_mat),
        shannon  = diversity(abundance_mat, index = "shannon"),
        simpson  = diversity(abundance_mat, index = "simpson"),
        hillq0   = hill_taxa(abundance_mat, q = 0),
        hillq1   = hill_taxa(abundance_mat, q = 1),
        hillq15  = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% 
  left_join(all_samples, ., 
            by = c("submitterID", "marker", "referenceDB", "pipeline", "sampleID"))

# ASV or OTU level
diversity_numbers_unit <- challengedata_12s16sCOI %>%
  filter(uniqueID %in% target_ASV,
         taxonomicLevel == "species") %>% # keep unassigned sequences
  group_by(submitterID, marker, referenceDB, pipeline, sampleID, seqID) %>%
  summarise(abundance = sum(readAbundance), .groups = "drop") %>% # probably unnecessary but keep anyway
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline, sampleID),
    names_from = seqID,                     # warning, different entries will have the same seqID. As long as I work per-sample, this is not a problem. 
    values_from = abundance,
    values_fill = 0) %>%
  ungroup() %>%     
  {               
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    bind_cols(                                                         
      select(., submitterID, marker, referenceDB, pipeline, sampleID), 
      data.frame(                                                      
        specnumber = specnumber(abundance_mat),
        shannon  = diversity(abundance_mat, index = "shannon"),
        simpson  = diversity(abundance_mat, index = "simpson"),
        hillq0   = hill_taxa(abundance_mat, q = 0),
        hillq1   = hill_taxa(abundance_mat, q = 1),
        hillq15  = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% 
  left_join(all_samples, ., 
            by = c("submitterID", "marker", "referenceDB", "pipeline", "sampleID"))

# ASV or OTU level but without filtering for Actino & Chondri, but ALL ASV&OTU
diversity_numbers_unit_noTarget <- challengedata_12s16sCOI %>%
  filter(#uniqueID %in% target_ASV,
    taxonomicLevel == "species") %>% # keep unassigned sequences
  group_by(submitterID, marker, referenceDB, pipeline, sampleID, seqID) %>%
  summarise(abundance = sum(readAbundance), .groups = "drop") %>% # probably unnecessary but keep anyway
  pivot_wider(
    id_cols = c(submitterID, marker, referenceDB, pipeline, sampleID),
    names_from = seqID,                     # warning, different entries will have the same seqID. As long as I work per-sample, this is not a problem. 
    values_from = abundance,
    values_fill = 0) %>%
  ungroup() %>%     
  {               
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    bind_cols(                                                         
      select(., submitterID, marker, referenceDB, pipeline, sampleID), 
      data.frame(                                                      
        specnumber = specnumber(abundance_mat),
        shannon  = diversity(abundance_mat, index = "shannon"),
        simpson  = diversity(abundance_mat, index = "simpson"),
        hillq0   = hill_taxa(abundance_mat, q = 0),
        hillq1   = hill_taxa(abundance_mat, q = 1),
        hillq15  = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% 
  left_join(all_samples, ., 
            by = c("submitterID", "marker", "referenceDB", "pipeline", "sampleID"))

# compare
plot(diversity_numbers_species$specnumber ~ diversity_numbers_species$hillq0)
plot(diversity_numbers_species$specnumber ~ diversity_numbers_genus$specnumber)
plot(diversity_numbers_species$hillq1     ~ diversity_numbers_genus$hillq1)
plot(diversity_numbers_species$hillq15    ~ diversity_numbers_genus$hillq15)

plot(diversity_numbers_unit$specnumber ~ diversity_numbers_unit_noTarget$specnumber)
plot(diversity_numbers_unit$hillq1     ~ diversity_numbers_unit_noTarget$hillq1)
plot(diversity_numbers_unit$hillq15    ~ diversity_numbers_unit_noTarget$hillq15)



# export
#write.csv(diversity_numbers_species,       file = "../results/,row.names = FALSE")
#write.csv(diversity_numbers_genus,         file = "../results/,row.names = FALSE)
#write.csv(diversity_numbers_unit,          file = "../results/,row.names = FALSE)
#write.csv(diversity_numbers_unit_noTarget, file = "../results/,row.names = FALSE)

