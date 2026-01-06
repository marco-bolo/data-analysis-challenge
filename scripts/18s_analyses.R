########                       MARCO-BOLO Data Analysis Challenge                       ########
########                         18s protist metabarcoding data                         ########

########         The influence of bioinformatic choices on time series patterns         ########
##########         script by Emilie Boulanger on R version 4.2.3 (2023-03-15)         ##########

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

library(vegan)
library(hillR)
library(mgcv)

library(ggplot2)
library(patchwork)
library(ggrepel)
library(viridis)

# import data 

# if working on csv file
challengedata_18s <- read.csv(file = "../data/challengedata_18s.csv", row.names = 1) %>% 
  mutate(submitterID = str_split_fixed(.$submissionID, '_', 4)[,1]) %>% 
  mutate(marker =  str_split_fixed(.$submissionID, '_', 4)[,2]) %>% 
  mutate(referenceDB =  str_split_fixed(.$submissionID, '_', 4)[,3]) %>% 
  mutate(pipeline =  str_split_fixed(.$submissionID, '_', 4)[,4]) %>% 
  select(submissionID, submitterID, marker, referenceDB, pipeline, everything())

# if working on rds file
challengedata_18s <- readRDS(file = "../data/challengedata_18s.rds") %>% 
  mutate(submitterID = str_split_fixed(.$submissionID, '_', 4)[,1]) %>% 
  mutate(marker =  str_split_fixed(.$submissionID, '_', 4)[,2]) %>% 
  mutate(referenceDB =  str_split_fixed(.$submissionID, '_', 4)[,3]) %>% 
  mutate(pipeline =  str_split_fixed(.$submissionID, '_', 4)[,4]) %>% 
  select(submissionID, submitterID, marker, referenceDB, pipeline, everything())

metadata_astan <- read.csv("../data/metadata_astan.csv")

microscopy_for_challengedata <- read.csv("../data/microscopy_harmonised.csv")

metadata_microscopy<- read.csv("../data/metadata_microscopy.csv")

# merge the metabarcoding and microscopy data
challengedata_18s_merged <- full_join(challengedata_18s , microscopy_for_challengedata)
metadata_18s_samples     <- full_join(metadata_astan, metadata_microscopy)

# set global parameters 
# define an ordination plotting theme
theme_ordination <- function() {
  list(theme_classic() + # base
         theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
               axis.title = element_text(size = 12, face = "bold"),
               axis.text = element_text(size = 10),
               axis.line = element_line(color = "black"),
               panel.grid = element_blank()),
       # dashed zero lines
       geom_hline(yintercept = 0, linetype = "dashed", color = "grey50"),
       geom_vline(xintercept = 0, linetype = "dashed", color = "grey50"))}

# Seasonal diversity fluctuations ----
## calculate diversity indices ----

# work at the "species" level
# first produce a species x site matrix for diversity calculation
diversity_df_18s <- challengedata_18s_merged %>%
  filter(!is.na(species) & species != "") %>% # some empty species strings otherwise cause problems
  pivot_longer(cols = starts_with("RA"), names_to = "sampleID", values_to = "readAbundance") %>%
  mutate(readAbundance = ifelse(is.na(readAbundance), 0, readAbundance)) %>%  # replace NA readAbundance by zero but keep all samples
  group_by(submissionID, submitterID, marker, referenceDB, pipeline, sampleID, species) %>%
  # pool duplicate species assignments 
  summarise(abundance = sum(readAbundance), .groups = "drop") %>%
  # widen to produce a species x site matrix for diversity calculations
  pivot_wider(
    id_cols = c(submissionID, submitterID, marker, referenceDB, pipeline, sampleID),
    names_from = species,
    values_from = abundance,
    values_fill = 0) %>%
  ungroup() 
colnames(diversity_df_18s)

# measure diversity numbers on the diversity matrix

## separateley for metabarcoding results
diversity_numbers_metab <- diversity_df_18s %>% 
  # only work on 18s data to not skew results because different samples
  filter(marker == "18S",
         !str_detect(sampleID, "^RAMORPHO")) %>% 
  {                 # create a code block within the %>% 
    # extract numeric abundance matrix that both vegan and hillR functions can use directly
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    # compute diversity metrics directly on the matrix
    bind_cols(                                                         # merge the sample metadata and diversity metrics into a single data.frame
      select(., submissionID, submitterID, marker, referenceDB, pipeline, sampleID), # keep sample metadata to combine with results later
      data.frame(                                                      # calculate the diversity metrics with vegan and hillR
        specnumber = specnumber(abundance_mat),
        shannon    = diversity(abundance_mat, index = "shannon"),
        simpson    = diversity(abundance_mat, index = "simpson"),
        hillq0     = hill_taxa(abundance_mat, q = 0),
        hillq1     = hill_taxa(abundance_mat, q = 1),
        hillq15    = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% # close code block
  left_join( ., metadata_18s_samples, by = "sampleID")

## separately for microscopy results
diversity_numbers_micros <- diversity_df_18s %>% 
  # only work on microscopy data this time
  filter(marker == "morphology",
         str_detect(sampleID, "^RAMORPHO")) %>% 
  {                 # create a code block within the %>% 
    # extract numeric abundance matrix that both vegan and hillR functions can use directly
    abundance_mat <- as.matrix(select(., where(is.numeric)))
    # compute diversity metrics directly on the matrix
    bind_cols(                                                         # merge the sample metadata and diversity metrics into a single data.frame
      select(., submissionID, submitterID, marker, referenceDB, pipeline, sampleID), # keep sample metadata to combine with results later
      data.frame(                                                      # calculate the diversity metrics with vegan and hillR
        specnumber = specnumber(abundance_mat),
        shannon    = diversity(abundance_mat, index = "shannon"),
        simpson    = diversity(abundance_mat, index = "simpson"),
        hillq0     = hill_taxa(abundance_mat, q = 0),
        hillq1     = hill_taxa(abundance_mat, q = 1),
        hillq15    = hill_taxa(abundance_mat, q = 1.5)
      ))} %>% # close code block
  left_join( ., metadata_18s_samples, by = "sampleID") 


# merge both into an overall dataset
colnames(diversity_numbers_metab)
colnames(diversity_numbers_micros)
diversity_numbers_merged <- full_join(diversity_numbers_metab, 
                                      diversity_numbers_micros)

## Boxplots to explore results ----
# boxplot are shown per month, pooling all years, and separated by pipeline
summary(diversity_numbers_merged$shannon)
summary(diversity_numbers_merged$specnumber)
summary(diversity_numbers_merged$hillq15)

pBOX_shannon <- diversity_numbers_merged %>% 
  ggplot(aes(x=factor(month), y=shannon, fill=season)) +
  geom_boxplot() + #need factor levels for month though
  scale_fill_manual(values=c("#0000ff", "#00ff02", "#ffa602","#a52a2b")) +
  labs(x = "Month", y = "Shannon Diversity") +
  facet_wrap(.~pipeline) +
  theme_classic()
pBOX_shannon

pBOX_richness <- diversity_numbers_merged %>% 
  ggplot(aes(x=factor(month), y=specnumber, fill=season)) +
  geom_boxplot() + #need factor levels for month though
  scale_fill_manual(values=c("#0000ff", "#00ff02", "#ffa602","#a52a2b")) +
  labs(x = "Month", y = "Species/Taxon richness") +
  facet_wrap(.~pipeline) +
  theme_classic()
pBOX_richness

pBOX_hillq15 <- diversity_numbers_merged %>% 
  ggplot(aes(x=factor(month), y=hillq15, fill=season)) +
  geom_boxplot() + #need factor levels for month though
  scale_fill_manual(values=c("#0000ff", "#00ff02", "#ffa602","#a52a2b")) +
  labs(x = "Month", y = "Hill number q = 1.5") +
  facet_wrap(.~pipeline) +
  theme_classic()
pBOX_hillq15

filter(diversity_numbers_merged, hillq15 == "Inf") # a few entries have no species assigned and thus infinite values for hill q=1.5

## Generalized Additive Models (GAM) of diversity indices in time ----
###  Compute GAM across indices, years, months and pipelines ----

pipeline_gam_shannon <- diversity_numbers_merged %>%
  group_by(pipeline) %>%
  do({fit <- gam(shannon ~ s(month, bs = "cc", k = 12), data = .)
  pred <- predict(fit, newdata = data.frame(month = 1:12), se.fit = TRUE)
  tibble(month = 1:12, fit = pred$fit, se = pred$se.fit,
         lwr = pred$fit - 1.96 * pred$se.fit, 
         upr = pred$fit + 1.96 * pred$se.fit)})

pipeline_gam_richness <- diversity_numbers_merged %>%
  group_by(pipeline) %>%
  do({fit <- gam(specnumber ~ s(month, bs = "cc", k = 12), data = .)
  pred <- predict(fit, newdata = data.frame(month = 1:12), se.fit = TRUE)
  tibble(month = 1:12, fit = pred$fit, se = pred$se.fit,
         lwr = pred$fit - 1.96 * pred$se.fit, 
         upr = pred$fit + 1.96 * pred$se.fit)})

pipeline_gam_hill15 <- diversity_numbers_merged %>%
  filter(pipeline != "SLIM") %>% # leave out because Inf values cause problems
  group_by(pipeline) %>%
  do({fit <- gam(hillq15 ~ s(month, bs = "cc", k = 12), data = .)
  pred <- predict(fit, newdata = data.frame(month = 1:12), se.fit = TRUE)
  tibble(month = 1:12, fit = pred$fit, se = pred$se.fit,
         lwr = pred$fit - 1.96 * pred$se.fit, 
         upr = pred$fit + 1.96 * pred$se.fit)})

### Plot the GAM outputs by defining labels for each ----
# too many labels to make sense of the color legend
# instead, add labels to the end of each line
# and order by last value

#### Shannon diversity ----
# order each pipeline by the last fitted (Shannon) GAM value:
pipeline_order_shannon <- pipeline_gam_shannon %>%
  group_by(pipeline) %>%
  summarize(last_fit = fit[which.max(month)]) %>%
  arrange(last_fit) %>%
  pull(pipeline)

# apply the order to the dataset
pipeline_gam_shannon <- pipeline_gam_shannon %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_shannon))

diversity_numbers_merged <- diversity_numbers_merged %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_shannon))

# define label locations
end_labels_shannon <- pipeline_gam_shannon %>%
  group_by(pipeline) %>%
  slice_max(month, n = 1) %>%   # last point in the line
  ungroup()

# set custom colour scale withe only mid to dark shades and set as palette
n_cols <- length(unique(diversity_numbers_merged$pipeline))
dark_mako <- viridis(n = n_cols, option = "mako",
                     begin = 0,      # start of palette (dark)
                     end = 0.8)       # cutoff before it gets too light
palette_shannon <- dark_mako
names(palette_shannon) <- pipeline_order_shannon

# Override microscopy color
palette_shannon["microscopy"] <- "red"  # or any other color you like

# ggplot
pGAM_shannon <- ggplot(diversity_numbers_merged, aes(x = month, y = shannon, color = pipeline)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pipeline_gam_shannon,
            aes(x = month, y = fit, color = pipeline),
            linewidth = 1) +
  geom_ribbon(data = pipeline_gam_shannon,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              aes(x = month, ymin = lwr, ymax = upr, fill = pipeline),
              alpha = 0.15, inherit.aes = FALSE) +
  geom_label_repel(data = end_labels_shannon,
                   aes(x = month, y = fit, label = pipeline, color = pipeline),
                   hjust = 0, size = 3, show.legend = FALSE,
                   direction = "y", nudge_x = 0.5, #fontface = "bold",
                   segment.color = "grey50") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.05, 0.2))) +
  scale_color_manual(values = palette_shannon) +
  scale_fill_manual(values  = palette_shannon) +
  labs(x = "Month", y = "Shannon Diversity") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("ASTAN 18S Time Series \nGAM fit + SE across pipelines and microscopy data")
pGAM_shannon
#ggsave("../figures/18s/plankton_timeseries_gam_shannon.png", width = 8.79, height = 7.5)

#### Species richness ----
# order each pipeline by the last fitted (Richness) GAM value:
pipeline_order_richness <- pipeline_gam_richness %>%
  group_by(pipeline) %>%
  summarize(last_fit = fit[which.max(month)]) %>%
  arrange(last_fit) %>%
  pull(pipeline)

# apply the order to the dataset
pipeline_gam_richness <- pipeline_gam_richness %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_richness))

diversity_numbers_merged <- diversity_numbers_merged %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_richness))

# define label locations
end_labels_richness <- pipeline_gam_richness %>%
  group_by(pipeline) %>%
  slice_max(month, n = 1) %>%   # last point in the line
  ungroup()

pGAM_richness <- ggplot(diversity_numbers_merged, aes(x = month, y = specnumber, color = pipeline)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pipeline_gam_richness,
            aes(x = month, y = fit, color = pipeline),
            linewidth = 1) +
  geom_ribbon(data = pipeline_gam_richness,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              aes(x = month, ymin = lwr, ymax = upr, fill = pipeline),
              alpha = 0.15, inherit.aes = FALSE) +
  geom_label_repel(data = end_labels_richness,
                   aes(x = month, y = fit, label = pipeline, color = pipeline),
                   hjust = 0, size = 3, show.legend = FALSE,
                   direction = "y", nudge_x = 0.5, #fontface = "bold",
                   segment.color = "grey50") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.05, 0.2))) +
  scale_color_manual(values = palette_shannon) + # use the same colours across plots
  scale_fill_manual(values  = palette_shannon) +
  labs(x = "Month", y = "Species richness") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("ASTAN 18S Time Series \nGAM fit + SE across pipelines and microscopy data")
pGAM_richness
#ggsave("../figures/18s/plankton_timeseries_gam_richness.png", width = 8.79, height = 7.5)

#### Hill number q = 1.5 ----
# order each pipeline by the last fitted (Richness) GAM value:
pipeline_order_hill15 <- pipeline_gam_hill15 %>%
  group_by(pipeline) %>%
  summarize(last_fit = fit[which.max(month)]) %>%
  arrange(last_fit) %>%
  pull(pipeline)

# apply the order to the dataset
pipeline_gam_hill15 <- pipeline_gam_hill15 %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_hill15))

diversity_numbers_merged <- diversity_numbers_merged %>%
  mutate(pipeline = factor(pipeline, levels = pipeline_order_hill15))

# define label locations
end_labels_hill15 <- pipeline_gam_hill15 %>%
  group_by(pipeline) %>%
  slice_max(month, n = 1) %>%   # last point in the line
  ungroup()

pGAM_hill15 <- diversity_numbers_merged %>% 
  filter(pipeline != "SLIM") %>% 
  ggplot(., aes(x = month, y = hillq15, color = pipeline)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pipeline_gam_hill15,
            aes(x = month, y = fit, color = pipeline),
            linewidth = 1) +
  geom_ribbon(data = pipeline_gam_hill15,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              aes(x = month, ymin = lwr, ymax = upr, fill = pipeline),
              alpha = 0.15, inherit.aes = FALSE) +
  geom_label_repel(data = end_labels_hill15,
                   aes(x = month, y = fit, label = pipeline, color = pipeline),
                   hjust = 0, size = 3, show.legend = FALSE,
                   direction = "y", nudge_x = 0.5, #fontface = "bold",
                   segment.color = "grey50") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.05, 0.2))) +
  scale_color_manual(values = palette_shannon) + # use the same colours across plots
  scale_fill_manual(values  = palette_shannon) +
  labs(x = "Month", y = "Hill numbers q = 1.5") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("ASTAN 18S Time Series \nGAM fit + SE across pipelines and microscopy data")
pGAM_hill15
#ggsave("../figures/18s/plankton_timeseries_gam_hillq15.png", width = 8.79, height = 7.5)


# Community composition ----

## NMDS plots across pipelines ----

# extract the matrix of read or microscopy abundances
abundance_mat <- as.matrix(select(diversity_df_18s, where(is.numeric)))
# remove empty rows
empty_rows <- rowSums(abundance_mat) == 0
abundance_mat <- abundance_mat[!empty_rows, ]
# apply hellinger transformation on relative abundances for Bray-Curtis
relabundance_hellinger <- decostand(abundance_mat, method = "hellinger")
# transform to presence/absence data for Jaccard
pa_mat <- abundance_mat # work on matrix without NAs
pa_mat[pa_mat>0] <- 1

### apply NMDS 

# with Bray-Curtis dissimilarity
nmds_bray <- metaMDS(relabundance_hellinger, distance = "bray", k = 2, trymax = 100)  # this will take some time, approx. 60-80 minutes
# temporary export

# with Jaccard distance
nmds_jacc <- metaMDS(pa_mat, distance = "jaccard", k = 2, trymax = 100)

# check stress
nmds_bray$stress #  0.1198466
nmds_jacc$stress #  0.07030733

# export nmds files
#save(nmds_bray, nmds_jacc, file = "../data/intermediate_nmds_18s.RData")

# Extract plotting coordinates
nmds <- nmds_bray               # here define nmds_jacc for the jaccard nmds plots. be mindful to also change the plotting lables and figure export names to jaccard instead of bray
nmds_df <- data.frame(diversity_df_18s[!empty_rows,1:6], NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2]) %>% 
  left_join(metadata_18s_samples)

# define label locations
center_labels <- nmds_df %>% 
  group_by(submitterID, pipeline) %>% 
  summarise(NMDS1 = mean(NMDS1, na.rm = T),
            NMDS2 = mean (NMDS2, na.rm = T)) 

# plot 
nmds_pipeline <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = pipeline)) +
  geom_point(size = 3) +
  scale_color_manual(values = palette_shannon, guide = guide_legend(ncol = 1)) +
  geom_label_repel(data = center_labels,
                   aes(x = NMDS1, y = NMDS2, label = pipeline, color = pipeline),
                   hjust = 1, size = 3, show.legend = FALSE,
                   max.overlaps = 30,
                   direction = "both", nudge_x = 0.2, nudge_y = 0.2, #fontface = "bold",
                   segment.color = "grey50") +
  annotate(geom = "text", label = paste0("Stress = ", round(nmds$stress,5)),
           x = -2.9, y = -3.5, hjust = 1, vjust = 1) +    # for Bray-Curtis
 #          x = -7.5, y = -13.5, hjust = 1, vjust = 1) +    # for Jaccard
  theme_ordination() +
  theme(legend.position = "none") +
  labs(x = "NMDS1", y = "NMDS2")

nmds_season <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = season)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Winter" = "blue", "Spring" = "green","Summer" = "orange","Autumn" = "brown")) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2")

nmds_year <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = year)) +
  geom_point(size = 3) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2")

nmds_pipeline / (nmds_season + nmds_year) +
  plot_annotation(title = "NMDS on Hellinger-transformed relative abundances & Bray-Curtis distance")
#  plot_annotation(title = "NMDS on presence-absence data & Jaccard distance")
#ggsave(file = "../figures/18s/plankton_nmds_bray.png", width = 12, height = 12)
#ggsave(file = "../figures/18s/plankton_nmds_jacc.png", width = 12, height = 12)


## Separate NMDS by pipeline ----
# it's hard to distinguish finer patterns for each of the pipeline submissions
# here, run separate NMDS analyses for each pipeline and plot them together

# first, split by submissionID
split_list <- diversity_df_18s %>%
  group_split(submissionID)

### Bray-Curtis dissimilarity ----
# use purrr::map to iterate the nmds code over each submissionID
# while keeping submission metadata

nmds_list <- map(split_list, function(df) {
  abund <- df %>% select(where(is.numeric))
  meta  <- df %>% select(submissionID, submitterID, marker, referenceDB, pipeline, sampleID)
  # compute mask removing empty rows once
  mask <- rowSums(abund, na.rm = TRUE) > 0
  abund <- abund[mask, ]
  meta  <- meta[mask, ]
  # apply hellinger transformation on relative abundances
  abund_hell <- decostand(abund, method = "hellinger")
  # bray-curtis nmds
  nmds_split  <- metaMDS(abund_hell, distance = "bray", k = 2, trymax = 100)
  data.frame(meta,NMDS1 = nmds_split$points[,1],NMDS2 = nmds_split$points[,2])
})

# combine results into one df and add metadata
nmds_split_bray <- bind_rows(nmds_list) %>% 
  left_join(metadata_18s_samples, by = "sampleID")

# plot seasonal and yearly patterns

nmds_split_bray_season <- ggplot(nmds_split_bray, aes(NMDS1, NMDS2, color = season)) +
  geom_point(size = 3) +
  facet_wrap(.~pipeline, scales = "free") +
  stat_ellipse(type = "norm", level = 0.95) +
  scale_color_manual(values = c("Winter" = "blue", "Spring" = "green","Summer" = "orange","Autumn" = "brown")) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2") 
nmds_split_bray_season
#ggsave(file = "../figures/18s/plankton_nmds_bray-season.png", width = 12, height = 10)

nmds_split_bray_year <- ggplot(nmds_split_bray, aes(NMDS1, NMDS2, color = year)) +
  geom_point(size = 3) +
  facet_wrap(.~pipeline, scales = "free") +
  scale_color_gradient(breaks = unique(nmds_split_bray$year)) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2", color = "Year") 
nmds_split_bray_year
#ggsave(file = "../figures/18s/plankton_nmds_bray-year.png", width = 12, height = 10)

### Jaccard dissimilarity ----

# use the same split list, but change to p/a & jaccard in the map

nmds_list_jacc <- map(split_list, function(df) {
  abund <- df %>% select(where(is.numeric))
  meta  <- df %>% select(submissionID, submitterID, marker, referenceDB, pipeline, sampleID)
  # compute mask removing empty rows once
  mask <- rowSums(abund, na.rm = TRUE) > 0
  abund <- abund[mask, ]
  meta  <- meta[mask, ]
  # presence/absence
  pa_mat <- abund # work on matrix without NAs
  pa_mat[pa_mat>0] <- 1
  # jaccard nmds
  nmds_split_jacc  <- metaMDS(pa_mat, distance = "jaccard", k = 2, trymax = 100)
  data.frame(meta,NMDS1 = nmds_split_jacc$points[,1],NMDS2 = nmds_split_jacc$points[,2])
})

# combine results into one df and add metadata
nmds_split_jacc <- bind_rows(nmds_list_jacc) %>% 
  left_join(metadata_18s_samples, by = "sampleID")


nmds_jacc_season <- ggplot(nmds_split_jacc, aes(NMDS1, NMDS2, color = season)) +
  geom_point(size = 3) +
  facet_wrap(.~pipeline, scales = "free") +
  stat_ellipse(type = "norm", level = 0.95) +
  scale_color_manual(values = c("Winter" = "blue", "Spring" = "green","Summer" = "orange","Autumn" = "brown")) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS on presence-absence matrix with Jaccard distance") 
nmds_jacc_season
#ggsave(file = "../figures/18s/plankton_nmds_jacc-season.png", width = 12, height = 10)


nmds_jacc_year <- ggplot(nmds_split_jacc, aes(NMDS1, NMDS2, color = year)) +
  geom_point(size = 3) +
  facet_wrap(.~pipeline, scales = "free") +
  scale_color_gradient(breaks = unique(nmds_split_jacc$year)) +
  theme_ordination() +
  labs(x = "NMDS1", y = "NMDS2", color = "Year", title = "NMDS on presence-absence matrix with Jaccard distance") 
nmds_jacc_year
#ggsave(file = "../figures/18s/plankton_nmds_jacc-year.png", width = 12, height = 10)


