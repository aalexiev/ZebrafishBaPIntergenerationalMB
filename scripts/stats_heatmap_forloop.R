# Make a figure showing all the stats in a heatmap

## calculate all alpha and beta div metrics that I don't already have


########## read in taxonomic data phyloseq object and separate into otu table and metadata ###########

library(phyloseq)

taxa_dm <- readRDS(paste0(input_dir, "ps.qc.rds"))
tax_tab <- as.data.frame(taxa_dm@otu_table)
rownames(tax_tab) <- substr(rownames(tax_tab), 1, 8) # clean up row names

tax_meta <- as.data.frame(taxa_dm@sam_data)
tax_meta$Morpholino <- as.factor(tax_meta$Morph)
tax_meta$Exposure <- as.factor(tax_meta$Exposure)
tax_meta$Sex <- as.factor(tax_meta$Sex)
rownames(tax_meta) <- NULL # remove the repeated rownames since we have a column to work with
tax_meta$SampleID <- substr(tax_meta$SampleID, 1, 8) # clean up row names

# make a presence/absence data frame
tax_cov <- tax_tab
tax_cov[tax_cov > 0] <- 1
tax_cov[tax_cov <= 0] <- 0

# calculate all metrics for taxa
tax_richness <- rowSums(tax_cov)
tax_shannon <- diversity(tax_tab, "shannon", base = exp(1))
tax_simpson <- diversity(tax_tab, "simpson")
tax_dm <- vegdist(tax_tab, method = "bray")
tax_sor <- vegdist(tax_tab, binary = TRUE)

# calculate sorenson metric for pathways
pw_sor <- vegdist(pathabund_tab_ALL, binary = TRUE)

#######################


######### add alpha div metrics to metadata files ##########

# filter to only important variables we are testing and relevel so baseline category is controls
meta_modALL <- metaALL %>%
  dplyr::select(c("SampleID", "Exposure", "Morpholino", "Generation", "Sex"))

# add taxonomic alpha diversity to the metadata file
alphatax_meta <- data.frame(tax_richness, tax_shannon, tax_simpson) %>%
  rownames_to_column("SampleID") %>%
  left_join(tax_meta, by = "SampleID") %>%
  mutate(tax_richness = (tax_richness - min(tax_richness))/(max(tax_richness) - min(tax_richness)),
         tax_shannon = (tax_shannon - min(tax_shannon))/(max(tax_shannon) - min(tax_shannon)),
         tax_simpson = (tax_simpson - min(tax_simpson))/(max(tax_simpson) - min(tax_simpson)),
         Exposure = relevel(Exposure, "DMSO"), Morpholino = relevel(Morpholino, "CoMo")) # normalize alpha diversity metrics

# add pathway alpha diveristy to the metadata file
alphapw_meta <- data.frame(pw_richness, pw_shannon, pw_simpson) %>%
  rownames_to_column("SampleID") %>%
  right_join(meta_modALL, by = "SampleID")%>%
  mutate(pw_richness = (pw_richness - min(pw_richness))/(max(pw_richness) - min(pw_richness)),
         pw_shannon = (pw_shannon - min(pw_shannon))/(max(pw_shannon) - min(pw_shannon)),
         pw_simpson = (pw_simpson - min(pw_simpson))/(max(pw_simpson) - min(pw_simpson)),
         Exposure = relevel(Exposure, "DMSO"), Morpholino = relevel(Morpholino, "CoMo")) # normalize alpha diversity metrics

#######################


########## prep for loops ##########

# make a list of all variables
covs <- c("Exposure", "Morpholino", "Generation", "Sex")

#######################

######### for loops to calculate alpha diversity stats ############

# make an empty output file
stat_hm_out <- c()

# taxonomy
alphatax <- c("tax_richness", "tax_shannon", "tax_simpson")
for (i in alphatax) {
  for (y in covs) {
    formula <- paste0(i, "~", y) # define the formula as each alpha div metric ~ covariate
    model <- summary(lm(formula = formula, data = alphatax_meta)) # run lm and summarize
    coeffs <- as.data.frame(model$coefficients) %>% # take model coefficients
      rownames_to_column("vars") %>% # rownames are variables, move to column called vars
      filter(vars != "(Intercept)") # remove (Intercept) row
    df <- data.frame(div = "alpha", data = "taxa", i, y, coeffs$vars, coeffs$Estimate, coeffs$`Pr(>|t|)`) # make a data frame of the relevant info for graphing
    stat_hm_out <- rbind(stat_hm_out, df) # add it to the output data frame
  }
}

# pathways
alphapws <- c("pw_richness", "pw_shannon", "pw_simpson")
for (i in alphapws) {
  for (y in covs) {
    formula <- paste0(i, "~", y)
    model <- summary(lm(formula = formula, data = alphapw_meta))
    coeffs <- as.data.frame(model$coefficients) %>%
      rownames_to_column("vars") %>%
      filter(vars != "(Intercept)")
    df <- data.frame(div = "alpha", data = "pathways", i, y, coeffs$vars, coeffs$Estimate, coeffs$`Pr(>|t|)`)
    stat_hm_out <- rbind(stat_hm_out, df)
  }
}
View(stat_hm_out)

#######################

########### graphing alpha div stats heatmap #############

# do some clean up and renaming of columns and values for easier graphing
alphastat_hm_outclean <- stat_hm_out %>%
  mutate(coeffs.vars = replace(coeffs.vars, coeffs.vars == "ExposureBaP", "Exposure"),
         coeffs.vars = replace(coeffs.vars, coeffs.vars == "MorpholinoAHR2Mo", "Morpholino"),
         coeffs.vars = replace(coeffs.vars, coeffs.vars == "SexM", "Sex"),
         i = replace(i, i == "tax_richness", "Taxonomic Richness"),
         i = replace(i, i == "tax_shannon", "Taxonomic Shannon diversity"),
         i = replace(i, i == "tax_simpson", "Taxonomic Simpson Index"),
         i = replace(i, i == "pw_richness", "Pathway Richness"),
         i = replace(i, i == "pw_shannon", "Pathway Shannon diversity"),
         i = replace(i, i == "pw_simpson", "Pathway Simpson Index"),
         coeffs..Pr...t... = replace(coeffs..Pr...t..., coeffs..Pr...t... <= 0.05, "*"),
         coeffs..Pr...t... = replace(coeffs..Pr...t..., coeffs..Pr...t... > 0.05, "n.s."))

alphastat_hm_outclean$coeffs.vars <- factor(alphastat_hm_outclean$coeffs.vars, 
                                       levels = c("Sex",
                                                  "Exposure",
                                                  "Morpholino",
                                                  "GenerationF1",
                                                  "GenerationF2"))

# heatmap ggplot
heat_alphastats <- ggplot(data = alphastat_hm_outclean, aes(x = coeffs.vars, y = i)) +
  geom_tile(aes(fill = coeffs.Estimate), colour = "white", linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "tan", mid = "white", high = "darkgreen", 
                       midpoint = 0, na.value = "white")  + 
  labs(y = "", x = "", fill = "Linear model\ncoefficient estimate") +
  geom_text(aes(label = coeffs..Pr...t...)) + 
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Alpha diversity metrics")
heat_alphastats

#######################



######## beta diversity for loops #############

# make a list of beta diversity metrics
beta_tax <- c("tax_dm", "tax_sor")
beta_pw <- c("pw_dm", "pw_sor")

# make an empty output file
stat_hm_out <- c()

# taxonomy
for (i in beta_tax) {
  for (y in covs) {
    formula <- as.formula(paste0(i, " ~ ", y)) # define the formula as each alpha div metric ~ covariate
    model <- adonis2(formula = formula, data = alphatax_meta)
    coeffs <- data.frame(vars = rownames(model), model$R2, model$`Pr(>F)`) %>%
      dplyr::filter(vars == paste(y))
    df <- data.frame(div = "beta", data = "taxa", i = i, y = y, coeffs$vars, coeffs$`model.R2`, coeffs$`model..Pr..F..`)
    stat_hm_out <- rbind(stat_hm_out, df)
  }
}
  
# pathways
for (i in beta_pw) {
  for (y in covs) {
    formula <- as.formula(paste0(i, " ~ ", y))
    model <- adonis2(formula = formula, data = alphapw_meta)
    coeffs <- data.frame(vars = rownames(model), model$R2, model$`Pr(>F)`) %>%
      dplyr::filter(vars == paste(y))
    df <- data.frame(div = "beta", data = "pathways", i = i, y = y, coeffs$vars, coeffs$`model.R2`, coeffs$`model..Pr..F..`)
    stat_hm_out <- rbind(stat_hm_out, df)
  }
}
# View(stat_hm_out)

#######################

########### graphing beta div stats heatmap #############

# do some clean up and renaming of columns and values for easier graphing
betastat_hm_outclean <- stat_hm_out %>%
  mutate(i = replace(i, i == "tax_dm", "Taxonomic Bray-Curtis dissimilarity"),
         i = replace(i, i == "tax_sor", "Taxonomic Sorenson metric"),
         i = replace(i, i == "pw_dm", "Pathway Bray-Curtis dissimilarity"),
         i = replace(i, i == "pw_sor", "Pathway Sorenson metric"),
         coeffs.model..Pr..F.. = replace(coeffs.model..Pr..F.., coeffs.model..Pr..F.. <= 0.05, "*"),
         coeffs.model..Pr..F.. = replace(coeffs.model..Pr..F.., coeffs.model..Pr..F.. > 0.05, "n.s."))

betastat_hm_outclean$coeffs.vars <- factor(betastat_hm_outclean$coeffs.vars, 
                                            levels = c("Sex",
                                                       "Exposure",
                                                       "Morpholino",
                                                       "Generation"))

# heatmap ggplot
heat_betastats <- ggplot(data = betastat_hm_outclean, aes(x = coeffs.vars, y = i)) +
  geom_tile(aes(fill = coeffs.model.R2), colour = "white", linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "tan", mid = "white", high = "darkgreen", 
                       midpoint = 0, na.value = "white")  + 
  labs(y = "", x = "", fill = "PERMANOVA\nR-squared") +
  geom_text(aes(label = coeffs.model..Pr..F..)) + 
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Beta diversity metrics")
heat_betastats

#######################

  
############ for loops to calculate beta disper stats ############
# make a list of the beta diversity results object names
beta_tax <- list(tax_dm, tax_sor)
beta_pw <- list(pw_dm, pw_sor)

# make an empty output file
stat_hm_out <- c()

# taxonomy
for (i in 1:length(beta_tax)) {
  for (y in 1:length(covs)) {
    # calculate beta dispersion
    dispers <- betadisper(beta_tax[[i]], alphatax_meta$Treatments)
    
    # clean up the dispersion file
    Disper1 <- data.frame(dispers$distances)
    colnames(Disper1) <- "dispers"
    Disper2 <- Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = alphatax_meta, 
                 by = "SampleID")
    
    # run model
    dist <- dispers$distances
    formula <- as.formula(paste("dist ~ ", covs[[y]])) # define the formula as each alpha div metric ~ covariate
    model <- summary(lm(formula = formula, data = Disper2))
    
    # clean up output and add to results file
    coeffs <- as.data.frame(model$coefficients) %>% # take model coefficients
      rownames_to_column("vars") %>% # rownames are variables, move to column called vars
      filter(vars != "(Intercept)") # remove (Intercept) row
    df <- data.frame(div = "betadisp", data = "taxa", i = paste(i, "tax"), covs[[y]], coeffs$vars, coeffs$Estimate, coeffs$`Pr(>|t|)`) # make a data frame of the relevant info for graphing
    stat_hm_out <- rbind(stat_hm_out, df) # add it to the output data frame
  }
}

# pathways
for (i in 1:length(beta_pw)) {
  for (y in 1:length(covs)) {
    # calculate beta dispersion
    dispers <- betadisper(beta_tax[[i]], alphatax_meta$Treatments)
    
    # clean up the dispersion file
    Disper1 <- data.frame(dispers$distances)
    colnames(Disper1) <- "dispers"
    Disper2 <- Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = alphatax_meta, 
                 by = "SampleID")
    
    # run model
    dist <- dispers$distances
    formula <- as.formula(paste("dist ~ ", covs[[y]])) # define the formula as each alpha div metric ~ covariate
    model <- summary(lm(formula = formula, data = Disper2))
    
    # clean up output and add to results file
    coeffs <- as.data.frame(model$coefficients) %>% # take model coefficients
      rownames_to_column("vars") %>% # rownames are variables, move to column called vars
      filter(vars != "(Intercept)") # remove (Intercept) row
    df <- data.frame(div = "betadisp", data = "pathways", i = paste(i, "pw"), covs[[y]], coeffs$vars, coeffs$Estimate, coeffs$`Pr(>|t|)`) # make a data frame of the relevant info for graphing
    stat_hm_out <- rbind(stat_hm_out, df) # add it to the output data frame
  }
}

View(stat_hm_out)

#######################
  
########### graphing stats heatmap #############

# do some clean up and renaming of columns and values for easier graphing
bdispstat_hm_outclean <- stat_hm_out %>%
  mutate(coeffs.vars = replace(coeffs.vars, coeffs.vars == "ExposureBaP", "Exposure"),
         coeffs.vars = replace(coeffs.vars, coeffs.vars == "MorpholinoAHR2Mo", "Morpholino"),
         coeffs.vars = replace(coeffs.vars, coeffs.vars == "SexM", "Sex"),
         i = replace(i, i == "1 tax", "Taxonomic Bray-Curtis dissimilarity"),
         i = replace(i, i == "2 tax", "Taxonomic Sorenson metric"),
         i = replace(i, i == "1 pw", "Pathway Bray-Curtis dissimilarity"),
         i = replace(i, i == "2 pw", "Pathway Sorenson metric"),
         coeffs..Pr...t... = replace(coeffs..Pr...t..., coeffs..Pr...t... <= 0.05, "*"),
         coeffs..Pr...t... = replace(coeffs..Pr...t..., coeffs..Pr...t... > 0.05, "n.s."))

bdispstat_hm_outclean$coeffs.vars <- factor(bdispstat_hm_outclean$coeffs.vars, 
                                            levels = c("Sex",
                                                       "Exposure",
                                                       "Morpholino",
                                                       "GenerationF1",
                                                       "GenerationF2"))

# heatmap ggplot
heat_bdispstats <- ggplot(data = bdispstat_hm_outclean, aes(x = coeffs.vars, y = i)) +
  geom_tile(aes(fill = coeffs.Estimate), colour = "white", linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "tan", mid = "white", high = "darkgreen", 
                       midpoint = 0, na.value = "white")  + 
  labs(y = "", x = "", fill = "Linear model\ncoefficient estimate") +
  geom_text(aes(label = coeffs..Pr...t...)) + 
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Beta dispersion")
heat_bdispstats

#######################

