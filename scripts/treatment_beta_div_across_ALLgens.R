# A script to produce a figure showing treatment beta diversity trends between generations
# author: Alexandra Alexiev

# Goal is to illustrate that the differences between generations in control fish microbiomes are lost when fish are treated

# Math, using point coordinates from ordination (per each treatment): 
# F0 = [ (distance from centroid F0 to each point in F0) / (median distance from centroid F0 to each point in F0) ] - 1
# F1 = [ (distance from centroid F0 to each point in F1) / (median distance from centroid F0 to each point in F0) ] - 1
# F2 = [ (distance from centroid F0 to each point in F2) / (median distance from centroid F0 to each point in F0) ] - 1

# write first with control treatment

# # extract F0 centroid
# cent_f0 <- t(as.data.frame(pw.disper[["centroids"]]["CoMo_DMSO_F0",1:2]))
# # calculate metric
# F0_metric <- pw_PC1 %>% # take all ordination coords
#   filter(treatgen == "CoMo_DMSO_F0") %>% # filter the points for focal treatment
#   dplyr::select(c("CAP1", "CAP2")) %>% # select only x and y coords for ordination
#   mutate(distF0 = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
#   mutate(distNorm = (distF0 / median(distF0)) - 1) # normalize to distance from F0 centroid to each F0 point (should be 0)
# # save F0 baseline dist for F1 and F2 calcs for this treatment
# distF0_CoMo_BaP <- median(F0_metric$distF0)
# F1_metric <- pw_PC1 %>% # take all ordination coords
#   filter(treatgen == "CoMo_DMSO_F1") %>% # filter the points for focal treatment
#   dplyr::select(c("CAP1", "CAP2")) %>% # select only x and y coords for ordination
#   mutate(distF1 = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
#   mutate(distNorm = (distF1 / distF0_CoMo_BaP) - 1) # normalize to distance from F0 centroid to each F0 point
# F2_metric <- pw_PC1 %>% # take all ordination coords
#   filter(treatgen == "CoMo_DMSO_F2") %>% # filter the points for focal treatment
#   dplyr::select(c("CAP1", "CAP2")) %>% # select only x and y coords for ordination
#   mutate(distF2 = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
#   mutate(distNorm = (distF2 / distF0_CoMo_BaP) - 1) # normalize to distance from F0 centroid to each F0 point

# for loop to calc metrics for each treatment and generation
treats <- c("CoMo_DMSO", "CoMo_BaP", "AHR2Mo_DMSO", "AHR2Mo_BaP")
pw_PC1 <- rownames_to_column(pw_PC1, "SampleID")
out <- c()

for (i in 1:length(treats)) {
  # extract F0 centroid
  cent_f0 <- t(as.data.frame(pw.disper[["centroids"]][paste0(treats[i],"_F0"),1:2]))
  
  # calculate metric
  F0_metric <- pw_PC1 %>% # take all ordination coords
    filter(treatgen == paste0(treats[i],"_F0")) %>% # filter the points for focal treatment
    dplyr::select(c("SampleID","CAP1", "CAP2")) %>% # select only x and y coords for ordination
    mutate(dist = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
    mutate(distNorm = (dist / median(dist)) - 1) # normalize to distance from F0 centroid to each F0 point (should be 0)
  
  # save F0 baseline dist for F1 and F2 calcs for this treatment
  distF0 <- median(F0_metric$dist)

  # calc F1 metrics
  F1_metric <- pw_PC1 %>% # take all ordination coords
    filter(treatgen == paste0(treats[i],"_F1")) %>% # filter the points for focal treatment
    dplyr::select(c("SampleID","CAP1", "CAP2")) %>% # select only x and y coords for ordination
    mutate(dist = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
    mutate(distNorm = (dist / distF0) - 1) # normalize to distance from F0 centroid to each F0 point
  
  # calc F2 metrics
  F2_metric <- pw_PC1 %>% # take all ordination coords
    filter(treatgen == paste0(treats[i],"_F2")) %>% # filter the points for focal treatment
    dplyr::select(c("SampleID","CAP1", "CAP2")) %>% # select only x and y coords for ordination
    mutate(dist = sqrt((cent_f0[1] - CAP1)^2 + (cent_f0[2] - CAP2)^2)) %>% # calc distance from F0 centroid to each point
    mutate(distNorm = (dist / distF0) - 1) # normalize to distance from F0 centroid to each F0 point
  
  # save metrics from each into output for graphing
  out[[i]] <- rbind(F0_metric, F1_metric, F2_metric)
  i <- i + 1

}

# clean up data frame
big_data <- do.call(rbind, out) %>%
  inner_join(pw_PC1, by = "SampleID") %>%
  dplyr::select(-CAP1.x, -CAP2.x) %>%
  rename(CAP1 = CAP1.y, CAP2 = CAP2.y) %>%
  group_by(Treatments, Generation) %>%
  summarise(median = median(distNorm), sd <- sd(distNorm)) %>%
  rename(sd = `sd <- sd(distNorm)`)

## the outcome should be a metric that represents the distance between generations normalized to F0's distance to itself.
# normalizing to F0 was just so that F0 could have a ~0 value and it would be easier to compare distance between each subsequent generation to 0.

