## Housekeeping

    # libraries
    library(dplyr)
    library(tidyr)
    library(plyr)
    library(ggplot2)
    library(MASS)
    library(EnvStats)
    library(ggpubr)
    library(stringr)
    # set seed so stats like permanova are reproducible
    set.seed(3)

    # to knit an Rmd file, you need explicit file paths typed out, so I created these objects as shortcuts so I could call the object for commonly used file paths instead of retyping the whole thing each time
    home_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/"
    input_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/input_files_pub/"
    output_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/output_files_pub/"

## Workflow summary (from Knecht et al., 2018 for adult behavior assays)

Clean up data: - remove dead and morphologically deficient fish, rows of
all 0’s and/or NA’s - truncate time period by removing first 5 min (when
fish are acclimating to new environment) and last portion of time
periods where some generations don’t have data compared to others (so
all generations have the same time endpoint) - change formatting of
generation and treatment columns to be simpler

Statistics

This part is confusing because the Knecht et al adult learning paper
does adult behavior assays like we did and claims to use a linear
regression approach that references Lisa’s paper, but the script Lisa
sent with graphs looks consistent with the other Knecht et al paper
(embryo and larval data) rather than the adult learning paper. That
approach in Knecht et al transgenerational paper uses non-adult (embryo
and larval) samples, just like Lisa’s paper, so probably not 100%
comparable to our samples. The curves we see here are also definitely
not linear (they are S-shaped) in CDF, although they are normal and also
could fit a gamma distribution. The image that Lisa last sent me to try
to recreate was also from the Transgenerational paper, not the adult
learning one. So below, I’ve done the transgenerational paper approach
but I’m not sure if that was correct and need clarification.

Knecht et al., transgenerational paper approach: - test assumptions of
statistical model are met: double truncated gamma distribution,
non-normal - identify statistical outliers - graph eCDF of each test
group (treatment per generation) - stats significance determined with
(independent samples) Kolmogorov-Smirnov test comparing exposed to
control - stats were computed as time-averaged differences between
individual time series and corresponding mean control response over the
truncated time interval

Knecht et al., adult learning paper approach: - test if normal and
linear - ANOVA with Tukey’s HSD used to compare between regression lines
(intercepts and slopes for each parameter)

## Read in and clean up data

-   Lisa already removed dead or morphologically deficient fish

## Freeswim

-   measured in distance (cm)

<!-- -->

    # behavior data: freeswim
    freeswim <- read.csv(file = paste0(input_dir, "freeswim_allCLEAN.csv"), 
                               header = TRUE)
    freeswim <-  mutate(freeswim, "group" = paste0(freeswim$generation, "_", freeswim$treatment))
    freeswim$MeanPerFish <- apply(freeswim[,4:(length(freeswim)-1)], 1, mean, na.rm = TRUE)
    which(freeswim$MeanPerFish == 0) # none
    ## integer(0)

    freeswim <- subset(freeswim, select = -c(T0:T10, T20:T28)) # truncate time period

    freeswim_long <- freeswim %>%
      mutate("FishNumber" = 1:nrow(freeswim)) %>% # add fish number
      # dplyr::filter(MeanPerFish > 0) %>% # filter out fish with a row mean of zero
      gather(timepoint, Distance, T11:T19, factor_key = T) # from wide to long data type

    ggplot(data = freeswim_long, 
           aes(x = timepoint, y = Distance, group = FishNumber)) +
      geom_line(aes(color = generation), alpha = 0.3) +
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/data cleanup freeswim-1.png" width="98%" height="98%" />

-   take out outliers per group using IQR method (this is code from
    Ebony Stretch that I modified into dplyr code)

<!-- -->

    # calculating the quartiles for each group
    freeswim_quartiles <- freeswim_long %>%
      group_by(group) %>%
      group_modify(~ {
         quantile(.x$MeanPerFish, probs = c(0.25, 0.75)) %>%
         tibble::enframe(name = "prob", value = "quantile") %>%
         spread(prob, quantile)
      }) %>%
      mutate(IQR = `75%` - `25%`,
             lower = `25%` - (1.5 * IQR),
             upper = `75%` + (1.5 * IQR)) %>%
      dplyr::select(c(group, lower, upper)) %>%
      full_join(freeswim_long, by = "group", multiple = "all") %>%
      ungroup() %>%
      as.data.frame()

    # subset 
    freeswim_long_outlierd <- subset(freeswim_quartiles, freeswim_quartiles$MeanPerFish > freeswim_quartiles$lower & freeswim_quartiles$MeanPerFish < freeswim_quartiles$upper)
    str(freeswim_long)
    ## 'data.frame':    1143 obs. of  8 variables:
    ##  $ generation : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment  : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID         : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ group      : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ MeanPerFish: num  1.951 0.691 0.569 1.714 1.866 ...
    ##  $ FishNumber : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint  : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Distance   : num  0.871 0.587 0.535 0.535 0.626 0.698 0.574 0.287 0.544 0.489 ...

    str(freeswim_long_outlierd)
    ## 'data.frame':    1107 obs. of  10 variables:
    ##  $ group      : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ lower      : num  -1.15 -1.15 -1.15 -1.15 -1.15 ...
    ##  $ upper      : num  3.43 3.43 3.43 3.43 3.43 ...
    ##  $ generation : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment  : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID         : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ MeanPerFish: num  1.951 0.691 0.569 1.714 1.866 ...
    ##  $ FishNumber : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint  : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Distance   : num  0.871 0.587 0.535 0.535 0.626 0.698 0.574 0.287 0.544 0.489 ...

    # some were removed so use the outlier removed file moving forward

## Shoaling: nearest-neighbor distance

-   The nearest-neighbor distance is the average of the nearest fish to
    each of the four fish in the shoal

<!-- -->

    # behavior data: shoaling nnd
    shoal_nnd <- read.csv(file = paste0(input_dir, "shoalingnnd_allCLEAN.csv"),
                               header = TRUE)
    shoal_nnd <-  mutate(shoal_nnd, "group" = paste0(shoal_nnd$generation, "_", shoal_nnd$treatment))
    shoal_nnd$MeanPerFish <- apply(shoal_nnd[,4:(length(shoal_nnd)-1)], 1, mean, na.rm = TRUE)
    which(shoal_nnd$MeanPerFish == 0) # none
    ## integer(0)

    shoal_nnd <- subset(shoal_nnd, select = -c(T1:T10, T20:T23)) # truncate time period

    # double check with line graphs of each fish's movement
    shoal_nnd_long <- shoal_nnd %>%
      mutate("FishNumber" = 1:nrow(shoal_nnd)) %>% # add fish number
      gather(timepoint, NearestNeighborDist, T11:T19, factor_key = T) # from wide to long data type

    ggplot(data = shoal_nnd_long, 
           aes(x = timepoint, y = NearestNeighborDist, group = FishNumber)) +
      geom_line(aes(color = generation), alpha = 0.3) +
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/data cleanup shoaling nnd-1.png" width="98%" height="98%" />

    # weird low and high set of fish dependent on generation, but Lisa said this could be biologically expected based on previous studies they've done at SARL

-   take out outliers per group using IQR method (this is code from
    Ebony Stretch that I modified into dplyr code)

<!-- -->

    # calculating the quartiles for each group
    shoal_nnd_quartiles <- shoal_nnd_long %>%
      group_by(group) %>%
      group_modify(~ {
         quantile(.x$MeanPerFish, probs = c(0.25, 0.75)) %>%
         tibble::enframe(name = "prob", value = "quantile") %>%
         spread(prob, quantile)
      }) %>%
      mutate(IQR = `75%` - `25%`,
             lower = `25%` - (1.5 * IQR),
             upper = `75%` + (1.5 * IQR)) %>%
      dplyr::select(c(group, lower, upper)) %>%
      full_join(shoal_nnd_long, by = "group", multiple = "all") %>%
      ungroup() %>%
      as.data.frame()

    # subset 
    shoal_nnd_long_outlierd <- subset(shoal_nnd_quartiles, shoal_nnd_quartiles$MeanPerFish > shoal_nnd_quartiles$lower & shoal_nnd_quartiles$MeanPerFish < shoal_nnd_quartiles$upper)
    str(shoal_nnd_long)
    ## 'data.frame':    1242 obs. of  8 variables:
    ##  $ generation         : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment          : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID                 : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ group              : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ MeanPerFish        : num  0.407 0.402 0.324 0.321 0.364 ...
    ##  $ FishNumber         : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint          : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ NearestNeighborDist: num  0.539 0.488 0.292 0.292 0.341 0.416 0.317 0.317 0.286 0.348 ...

    str(shoal_nnd_long_outlierd)
    ## 'data.frame':    1188 obs. of  10 variables:
    ##  $ group              : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ lower              : num  0.204 0.204 0.204 0.204 0.204 ...
    ##  $ upper              : num  0.475 0.475 0.475 0.475 0.475 ...
    ##  $ generation         : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment          : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID                 : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ MeanPerFish        : num  0.407 0.402 0.324 0.321 0.364 ...
    ##  $ FishNumber         : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint          : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ NearestNeighborDist: num  0.539 0.488 0.292 0.292 0.341 0.416 0.317 0.317 0.286 0.348 ...

    # some were removed so use the outlier removed file moving forward

## Shoaling: inter-individual distance

-   the inter-individual distance is the average distance between all
    four fish to each of the remaining 3 fish within the shoal

<!-- -->

    # behavior data: shoaling iid
    shoal_iid <- read.csv(file = paste0(input_dir, "shoalingiid_allCLEAN.csv"),
                               header = TRUE)
    shoal_iid <-  mutate(shoal_iid, "group" = paste0(shoal_iid$generation, "_", shoal_iid$treatment))
    shoal_iid$MeanPerFish <- apply(shoal_iid[,4:(length(shoal_iid)-1)], 1, mean, na.rm = TRUE)
    which(shoal_iid$MeanPerFish == 0) # no zero's in this data
    ## integer(0)

    shoal_iid <- subset(shoal_iid, select = -c(T0:T10, T20:T28)) # truncate time period

    # double check with line graphs of each fish's movement
    shoal_iid_long <- shoal_iid %>%
      mutate("FishNumber" = 1:nrow(shoal_iid)) %>% # add fish number
      gather(timepoint, InterIndivDist, T11:T19, factor_key = T) # from wide to long data type

    ggplot(data = shoal_iid_long, 
           aes(x = timepoint, y = InterIndivDist, group = FishNumber)) +
      geom_line(aes(color = generation), alpha = 0.3) +
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/data cleanup shoaling iid-1.png" width="98%" height="98%" />

    # same situation as nnd, there is a low and high set of fish

-   take out outliers per group using IQR method (this is code from
    Ebony Stretch that I modified into dplyr code)

<!-- -->

    # calculating the quartiles for each group
    shoal_iid_quartiles <- shoal_iid_long %>%
      group_by(group) %>%
      group_modify(~ {
         quantile(.x$MeanPerFish, probs = c(0.25, 0.75)) %>%
         tibble::enframe(name = "prob", value = "quantile") %>%
         spread(prob, quantile)
      }) %>%
      mutate(IQR = `75%` - `25%`,
             lower = `25%` - (1.5 * IQR),
             upper = `75%` + (1.5 * IQR)) %>%
      dplyr::select(c(group, lower, upper)) %>%
      full_join(shoal_iid_long, by = "group", multiple = "all") %>%
      ungroup() %>%
      as.data.frame()

    # subset 
    shoal_iid_long_outlierd <- subset(shoal_iid_quartiles, shoal_iid_quartiles$MeanPerFish > shoal_iid_quartiles$lower & shoal_iid_quartiles$MeanPerFish < shoal_iid_quartiles$upper)
    str(shoal_iid_long)
    ## 'data.frame':    999 obs. of  8 variables:
    ##  $ generation    : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment     : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID            : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank02_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank02_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ group         : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ MeanPerFish   : num  0.691 0.569 0.722 0.64 0.556 ...
    ##  $ FishNumber    : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint     : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ InterIndivDist: num  0.587 0.535 0.698 0.574 0.544 0.489 0.533 0.455 0.445 0.739 ...

    str(shoal_iid_long_outlierd)
    ## 'data.frame':    963 obs. of  10 variables:
    ##  $ group         : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ lower         : num  0.352 0.352 0.352 0.352 0.352 ...
    ##  $ upper         : num  0.895 0.895 0.895 0.895 0.895 0.895 0.895 0.895 0.895 0.895 ...
    ##  $ generation    : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment     : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID            : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank02_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank02_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ MeanPerFish   : num  0.691 0.569 0.722 0.64 0.556 ...
    ##  $ FishNumber    : int  1 2 3 4 5 6 1 2 3 4 ...
    ##  $ timepoint     : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 2 2 2 2 ...
    ##  $ InterIndivDist: num  0.587 0.535 0.698 0.574 0.544 0.489 0.683 0.499 0.699 0.531 ...

    # some were removed so use the outlier removed file moving forward

## Shoaling: speed

-   average speed of all four fish within the shoal

<!-- -->

    # behavior data: shoaling speed
    shoal_speed <- read.csv(file = paste0(input_dir, "shoalingspeed_allCLEAN.csv"),
                               header = TRUE)
    shoal_speed <- mutate(shoal_speed, "group" = paste0(shoal_speed$generation, "_", shoal_speed$treatment))
    shoal_speed$MeanPerFish <- apply(shoal_speed[,4:(length(shoal_speed)-1)], 1, mean, na.rm = TRUE)
    which(shoal_speed$MeanPerFish == 0) # no zero's in this data
    ## integer(0)

    shoal_speed <- subset(shoal_speed, select = -c(T0:T10, T20:T28)) # truncate time period

    # double check with line graphs of each fish's movement
    shoal_speed_long <- shoal_speed %>%
      mutate("FishNumber" = 1:nrow(shoal_speed)) %>% # add fish number
      gather(timepoint, speed, T11:T19, factor_key = T) # from wide to long data type

    ggplot(data = shoal_speed_long, 
           aes(x = timepoint, y = speed, group = FishNumber)) +
      geom_line(aes(color = generation), alpha = 0.3) +
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/data cleanup shoaling speed-1.png" width="98%" height="98%" />


    # same pattern as above, high and low grouping of fish
    # also more noise with these very high peaks on a few fish that might be outliers

-   take out outliers per group using IQR method (this is code from
    Ebony Stretch that I modified into dplyr code)

<!-- -->

    # calculating the quartiles for each group
    shoal_speed_quartiles <- shoal_speed_long %>%
      group_by(group) %>%
      group_modify(~ {
         quantile(.x$MeanPerFish, probs = c(0.25, 0.75)) %>%
         tibble::enframe(name = "prob", value = "quantile") %>%
         spread(prob, quantile)
      }) %>%
      mutate(IQR = `75%` - `25%`,
             lower = `25%` - (1.5 * IQR),
             upper = `75%` + (1.5 * IQR)) %>%
      dplyr::select(c(group, lower, upper)) %>%
      full_join(shoal_speed_long, by = "group", multiple = "all") %>%
      ungroup() %>%
      as.data.frame()

    # subset 
    shoal_speed_long_outlierd <- subset(shoal_speed_quartiles, shoal_speed_quartiles$MeanPerFish > shoal_speed_quartiles$lower & shoal_speed_quartiles$MeanPerFish < shoal_speed_quartiles$upper)
    str(shoal_speed_long)
    ## 'data.frame':    1152 obs. of  8 variables:
    ##  $ generation : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment  : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID         : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ group      : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ MeanPerFish: num  0.856 1.132 0.715 0.71 0.645 ...
    ##  $ FishNumber : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint  : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ speed      : num  0.948 0.997 0.678 0.678 0.829 ...

    str(shoal_speed_long_outlierd)
    ## 'data.frame':    1134 obs. of  10 variables:
    ##  $ group      : chr  "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" "F0_AHRMo_BaP" ...
    ##  $ lower      : num  0.244 0.244 0.244 0.244 0.244 ...
    ##  $ upper      : num  1.35 1.35 1.35 1.35 1.35 ...
    ##  $ generation : chr  "F0" "F0" "F0" "F0" ...
    ##  $ treatment  : chr  "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" "AHRMo_BaP" ...
    ##  $ ID         : chr  "Tank01_AHR 5uM BaP_Video Camera 5 5-17-2018 11-20-17 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-22-2018 9-56-12 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" "Tank01_AHR 5uM BaP_Video Camera 5 5-25-2018 10-01-29 AM 1" ...
    ##  $ MeanPerFish: num  0.856 1.132 0.715 0.71 0.645 ...
    ##  $ FishNumber : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ timepoint  : Factor w/ 9 levels "T11","T12","T13",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ speed      : num  0.948 0.997 0.678 0.678 0.829 ...

    # some were removed so use the outlier removed file moving forward

## Statistical analysis: two sample Kolmogorov-Smirnov test

-   Statistical significance determined with two-sample
    Kolmogorov-Smirnov test comparing experimental treatments
    (AHRMo\_BaP, AHRMo\_DMSO, CoMo \_BaP) to control ( CoMo\_DMSO)
-   Null hypothesis: two distributions were drawn from the same
    continuous distribution
-   p-value &lt;= 0.05 means reject this null and therefore
    distributions are from different sources (so experimental behavior
    is atypical compared to control behavior measure)
-   using “two-sided” as the alternative because we don’t actually know
    if treatment distributions will be to the left or right of the
    control (less or greater); anything not like the control would be
    considered atypical behavior
-   stats were computed as time-averaged differences between individual
    time series and corresponding mean control response over the
    specified time interval

## Freeswim

-   test assumptions of statistical model are met: non-normal data, f(x)
    is non-decreasing and right-continuous

<!-- -->

    # check normality assumption/what type of distribution
    # histogram (for all distances pooled across time points)
    mu_freeswim <- ddply(freeswim_long_outlierd, .(generation, treatment), summarise, 
                         grp.mean = mean(Distance, na.rm = TRUE))
    mu_freeswim
    ##    generation  treatment  grp.mean
    ## 1          F0  AHRMo_BaP 0.6104222
    ## 2          F0 AHRMo_DMSO 0.6160370
    ## 3          F0  CoMo _BaP 0.6389074
    ## 4          F0  CoMo_DMSO 0.5839074
    ## 5          F1  AHRMo_BaP 4.4779778
    ## 6          F1 AHRMo_DMSO 4.1567222
    ## 7          F1  CoMo _BaP 4.2964167
    ## 8          F1  CoMo_DMSO 4.5420185
    ## 9          F2  AHRMo_BaP 4.5727350
    ## 10         F2 AHRMo_DMSO 4.2414444
    ## 11         F2  CoMo _BaP 4.3777619
    ## 12         F2  CoMo_DMSO 4.2987778


    ggplot(data = freeswim_long_outlierd, 
           aes(x = Distance, color = generation, fill = generation)) +
      geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.2) +
      geom_vline(data = mu_freeswim, aes(xintercept = grp.mean, color = generation),
               linetype = "dashed") +
      geom_density(alpha = 0.1) + 
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/test assumptions-1.png" width="98%" height="98%" />

-   Two sample K-S test comparing experimental treatments to control
    group where p-value &lt;= 0.05 means the two data sets likely came
    from different distributions

<!-- -->

    # summary list for full stats for easy reference later
    ks_freeswimsmry <- c()

    ## F0 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_freeswimsmry$F0_BaPeffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_freeswimsmry$F0_AhR2effect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_AHRMo_DMSO"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_freeswimsmry$F0_intereffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F0_AHRMo_BaP"]) # from different distributions


    ## F1 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_freeswimsmry$F1_BaPeffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_freeswimsmry$F1_AhR2effect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_freeswimsmry$F1_intereffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F1_AHRMo_BaP"]) # not from different distributions


    ## F2 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_freeswimsmry$F2_BaPeffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_freeswimsmry$F2_AhR2effect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_AHRMo_DMSO"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_freeswimsmry$F2_intereffect <- ks.test(freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       freeswim_long_outlierd$Distance[freeswim_long_outlierd$group %in% "F2_AHRMo_BaP"]) # from different distributions

    ## make summary data frame of just the useful stats for publishing
    # list of test stats
    ksfreeswim_teststats <- c(ks_freeswimsmry$F0_BaPeffect$statistic, 
                              ks_freeswimsmry$F0_AhR2effect$statistic,
                              ks_freeswimsmry$F0_intereffect$statistic,
                              ks_freeswimsmry$F1_BaPeffect$statistic, 
                              ks_freeswimsmry$F1_AhR2effect$statistic,
                              ks_freeswimsmry$F1_intereffect$statistic, 
                              ks_freeswimsmry$F2_BaPeffect$statistic, 
                              ks_freeswimsmry$F2_AhR2effect$statistic,
                              ks_freeswimsmry$F2_intereffect$statistic)
    # list of p-values
    ksfreeswim_pvals <- c(ks_freeswimsmry$F0_BaPeffect$p.value, 
                              ks_freeswimsmry$F0_AhR2effect$p.value,
                              ks_freeswimsmry$F0_intereffect$p.value,
                              ks_freeswimsmry$F1_BaPeffect$p.value, 
                              ks_freeswimsmry$F1_AhR2effect$p.value,
                              ks_freeswimsmry$F1_intereffect$p.value, 
                              ks_freeswimsmry$F2_BaPeffect$p.value, 
                              ks_freeswimsmry$F2_AhR2effect$p.value,
                              ks_freeswimsmry$F2_intereffect$p.value)

    # data frame of all summary stats of KS test for quick, easy reference and publishing
    ks_freeswimsmry_short <- data.frame("Generation" = c(rep("F0", 3), rep("F1", 3), rep("F2", 3)),
                                        "Experimental_treatment_group" = c(rep(c("AhR2Mo - / BaP +", "AhR2Mo + / BaP -", "AhR2Mo + / BaP +"), 3)), 
                                        "Test_statistic" = ksfreeswim_teststats,
                                        "P_values" = ksfreeswim_pvals)
    ks_freeswimsmry_short$signif <- apply(ks_freeswimsmry_short, 1, 
                                   function(x) {ifelse(x[["P_values"]] <= 0.05, "*", "ns")})
    ks_freeswimsmry_short
    ##   Generation Experimental_treatment_group Test_statistic     P_values signif
    ## 1         F0             AhR2Mo - / BaP +      0.2407407 0.0857760115     ns
    ## 2         F0             AhR2Mo + / BaP -      0.1728395 0.2593598184     ns
    ## 3         F0             AhR2Mo + / BaP +      0.2296296 0.0460684996      *
    ## 4         F1             AhR2Mo - / BaP +      0.1759259 0.0706872848     ns
    ## 5         F1             AhR2Mo + / BaP -      0.3333333 0.0005537942      *
    ## 6         F1             AhR2Mo + / BaP +      0.1703704 0.0994281426     ns
    ## 7         F2             AhR2Mo - / BaP +      0.1111111 0.4694024858     ns
    ## 8         F2             AhR2Mo + / BaP -      0.1239316 0.3542931756     ns
    ## 9         F2             AhR2Mo + / BaP +      0.2336182 0.0043525418      *

-   graphing the eCDF by generation with lines as treatments

<!-- -->

    freeswim_cdf_plot <- ggplot(freeswim_long_outlierd, aes(x = Distance, color = treatment)) +
      stat_ecdf(geom = "step") +
      theme_classic() +
      facet_wrap("generation") +
      labs(x = "Freeswim Distance", y = "Time (min)", color = "Treatment") +
      scale_color_brewer(palette = "PuOr", 
                         name = "Treatments", 
                        labels = c("AhR2Mo + / BaP +",
                                     "AhR2Mo + / BaP -",
                                     "AhR2Mo - / BaP +",
                                     "AhR2Mo - / BaP -"))
    freeswim_cdf_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/CDF graph for freeswim data-1.png" width="98%" height="98%" />

## Shoaling: nearest-neighbor distance

-   The nearest-neighbor distance is the average of the nearest fish to
    each of the four fish in the shoal

<!-- -->

    # check normality assumption/what type of distribution
    # histogram (for all distances pooled across time points)
    mu_shoal_nnd <- ddply(shoal_nnd_long_outlierd, .(generation, treatment), summarise, 
                         grp.mean = mean(NearestNeighborDist, na.rm = TRUE))
    mu_shoal_nnd
    ##    generation  treatment  grp.mean
    ## 1          F0  AHRMo_BaP 0.3341204
    ## 2          F0 AHRMo_DMSO 0.3051235
    ## 3          F0   CoMo_BaP 0.3560617
    ## 4          F0  CoMo_DMSO 0.3181852
    ## 5          F1  AHRMo_BaP 2.5429333
    ## 6          F1 AHRMo_DMSO 2.3259074
    ## 7          F1   CoMo_BaP 2.4380833
    ## 8          F1  CoMo_DMSO 2.5283519
    ## 9          F2  AHRMo_BaP 2.5105726
    ## 10         F2 AHRMo_DMSO 2.3617949
    ## 11         F2   CoMo_BaP 2.4573810
    ## 12         F2  CoMo_DMSO 2.5055299


    ggplot(data = shoal_nnd_long_outlierd, 
           aes(x = NearestNeighborDist, color = generation, fill = generation)) +
      geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2) +
      geom_vline(data = mu_shoal_nnd, aes(xintercept = grp.mean, color = generation),
               linetype = "dashed") +
      geom_density(alpha = 0.1) + 
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/test assumptions and check distribution shoaling nnd-1.png" width="98%" height="98%" />

-   Two sample K-S test comparing experimental treatments to control
    group where p-value &lt;= 0.05 means the two data sets likely came
    from different distributions

<!-- -->

    # summary list for full stats for easy reference later
    ks_shoal_nndsmry <- c()

    ## F0 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_nndsmry$F0_BaPeffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_CoMo_BaP"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_nndsmry$F0_AhR2effect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_AHRMo_DMSO"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_nndsmry$F0_intereffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F0_AHRMo_BaP"]) # not from different distributions


    ## F1 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_nndsmry$F1_BaPeffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_CoMo_BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_nndsmry$F1_AhR2effect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_nndsmry$F1_intereffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F1_AHRMo_BaP"]) # from different distributions


    ## F2 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_nndsmry$F2_BaPeffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_CoMo_BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_nndsmry$F2_AhR2effect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_nndsmry$F2_intereffect <- ks.test(shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_nnd_long_outlierd$NearestNeighborDist[shoal_nnd_long_outlierd$group %in% "F2_AHRMo_BaP"]) # from different distributions

    ## make summary data frame of just the useful stats for publishing
    # list of test stats
    ksshoal_nnd_teststats <- c(ks_shoal_nndsmry$F0_BaPeffect$statistic, 
                              ks_shoal_nndsmry$F0_AhR2effect$statistic,
                              ks_shoal_nndsmry$F0_intereffect$statistic,
                              ks_shoal_nndsmry$F1_BaPeffect$statistic, 
                              ks_shoal_nndsmry$F1_AhR2effect$statistic,
                              ks_shoal_nndsmry$F1_intereffect$statistic, 
                              ks_shoal_nndsmry$F2_BaPeffect$statistic, 
                              ks_shoal_nndsmry$F2_AhR2effect$statistic,
                              ks_shoal_nndsmry$F2_intereffect$statistic)
    # list of p-values
    ksshoal_nnd_pvals <- c(ks_shoal_nndsmry$F0_BaPeffect$p.value, 
                              ks_shoal_nndsmry$F0_AhR2effect$p.value,
                              ks_shoal_nndsmry$F0_intereffect$p.value,
                              ks_shoal_nndsmry$F1_BaPeffect$p.value, 
                              ks_shoal_nndsmry$F1_AhR2effect$p.value,
                              ks_shoal_nndsmry$F1_intereffect$p.value, 
                              ks_shoal_nndsmry$F2_BaPeffect$p.value, 
                              ks_shoal_nndsmry$F2_AhR2effect$p.value,
                              ks_shoal_nndsmry$F2_intereffect$p.value)

    # data frame of all summary stats of KS test for quick, easy reference and publishing
    ks_shoal_nndsmry_short <- data.frame("Generation" = c(rep("F0", 3), rep("F1", 3), rep("F2", 3)),
                                        "Experimental_treatment_group" = c(rep(c("AhR2Mo - / BaP +", "AhR2Mo + / BaP -", "AhR2Mo + / BaP +"), 3)), 
                                        "Test_statistic" = ksshoal_nnd_teststats,
                                        "P_values" = ksshoal_nnd_pvals)
    ks_shoal_nndsmry_short$signif <- apply(ks_shoal_nndsmry_short, 1, 
                                   function(x) {ifelse(x[["P_values"]] <= 0.05, "*", "ns")})
    ks_shoal_nndsmry_short
    ##   Generation Experimental_treatment_group Test_statistic     P_values signif
    ## 1         F0             AhR2Mo - / BaP +     0.22222222 0.0306982334      *
    ## 2         F0             AhR2Mo + / BaP -     0.27160494 0.0038011043      *
    ## 3         F0             AhR2Mo + / BaP +     0.12962963 0.3612485741     ns
    ## 4         F1             AhR2Mo - / BaP +     0.14814815 0.1867397197     ns
    ## 5         F1             AhR2Mo + / BaP -     0.32407407 0.0008351264      *
    ## 6         F1             AhR2Mo + / BaP +     0.18888889 0.0497374855      *
    ## 7         F2             AhR2Mo - / BaP +     0.08730159 0.7442496196     ns
    ## 8         F2             AhR2Mo + / BaP -     0.23076923 0.0039358750      *
    ## 9         F2             AhR2Mo + / BaP +     0.17948718 0.0461390384      *

-   graphing the eCDF by generation with lines as treatments

<!-- -->

    shoal_nnd_cdf_plot <- ggplot(shoal_nnd_long_outlierd, aes(x = NearestNeighborDist, 
                                                              color = treatment)) +
      stat_ecdf(geom = "step") +
      theme_classic() +
      facet_wrap("generation") +
      labs(x = "Shoaling Nearest Neighbor Distance", 
           y = "Time (min)", 
           color = "Treatment") +
      scale_color_brewer(palette = "PuOr", 
                         name = "Treatments", 
                        labels = c("AhR2Mo + / BaP +",
                                     "AhR2Mo + / BaP -",
                                     "AhR2Mo - / BaP +",
                                     "AhR2Mo - / BaP -"))
    shoal_nnd_cdf_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/CDF graph for shoaling nnd data-1.png" width="98%" height="98%" />

## Shoaling: inter-individual distance

-   the inter-individual distance is the average distance between all
    four fish to each of the remaining 3 fish within the shoal

<!-- -->

    # check normality assumption/what type of distribution
    # histogram (for all distances pooled across time points)
    mu_shoal_iid <- ddply(shoal_iid_long_outlierd, .(generation, treatment), summarise, 
                         grp.mean = mean(InterIndivDist, na.rm = TRUE))
    mu_shoal_iid
    ##    generation  treatment  grp.mean
    ## 1          F0  AHRMo_BaP 0.6037963
    ## 2          F0 AHRMo_DMSO 0.5177037
    ## 3          F0  CoMo _BaP 0.6407037
    ## 4          F0  CoMo_DMSO 0.6401481
    ## 5          F1  AHRMo_BaP 4.4779778
    ## 6          F1 AHRMo_DMSO 4.1567222
    ## 7          F1  CoMo _BaP 4.2964167
    ## 8          F1  CoMo_DMSO 4.5420185
    ## 9          F2  AHRMo_BaP 4.5727350
    ## 10         F2 AHRMo_DMSO 4.2414444
    ## 11         F2  CoMo _BaP 4.3777619
    ## 12         F2  CoMo_DMSO 4.2987778


    ggplot(data = shoal_iid_long_outlierd, 
           aes(x = InterIndivDist, color = generation, fill = generation)) +
      geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2) +
      geom_vline(data = mu_shoal_iid, aes(xintercept = grp.mean, color = generation),
               linetype = "dashed") +
      geom_density(alpha = 0.1) + 
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/test assumptions and check distribution shoaling iid-1.png" width="98%" height="98%" />

-   Two sample K-S test comparing experimental treatments to control
    group where p-value &lt;= 0.05 means the two data sets likely came
    from different distributions

<!-- -->

    # summary list for full stats for easy reference later
    ks_shoal_iidsmry <- c()

    ## F0 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_iidsmry$F0_BaPeffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_iidsmry$F0_AhR2effect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_iidsmry$F0_intereffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F0_AHRMo_BaP"]) # not from different distributions


    ## F1 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_iidsmry$F1_BaPeffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_iidsmry$F1_AhR2effect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_iidsmry$F1_intereffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F1_AHRMo_BaP"]) # not from different distributions


    ## F2 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_iidsmry$F2_BaPeffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_CoMo _BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_iidsmry$F2_AhR2effect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_AHRMo_DMSO"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_iidsmry$F2_intereffect <- ks.test(shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_iid_long_outlierd$InterIndivDist[shoal_iid_long_outlierd$group %in% "F2_AHRMo_BaP"]) # from different distributions

    ## make summary data frame of just the useful stats for publishing
    # list of test stats
    ksshoal_iid_teststats <- c(ks_shoal_iidsmry$F0_BaPeffect$statistic, 
                              ks_shoal_iidsmry$F0_AhR2effect$statistic,
                              ks_shoal_iidsmry$F0_intereffect$statistic,
                              ks_shoal_iidsmry$F1_BaPeffect$statistic, 
                              ks_shoal_iidsmry$F1_AhR2effect$statistic,
                              ks_shoal_iidsmry$F1_intereffect$statistic, 
                              ks_shoal_iidsmry$F2_BaPeffect$statistic, 
                              ks_shoal_iidsmry$F2_AhR2effect$statistic,
                              ks_shoal_iidsmry$F2_intereffect$statistic)
    # list of p-values
    ksshoal_iid_pvals <- c(ks_shoal_iidsmry$F0_BaPeffect$p.value, 
                              ks_shoal_iidsmry$F0_AhR2effect$p.value,
                              ks_shoal_iidsmry$F0_intereffect$p.value,
                              ks_shoal_iidsmry$F1_BaPeffect$p.value, 
                              ks_shoal_iidsmry$F1_AhR2effect$p.value,
                              ks_shoal_iidsmry$F1_intereffect$p.value, 
                              ks_shoal_iidsmry$F2_BaPeffect$p.value, 
                              ks_shoal_iidsmry$F2_AhR2effect$p.value,
                              ks_shoal_iidsmry$F2_intereffect$p.value)

    # data frame of all summary stats of KS test for quick, easy reference and publishing
    ks_shoal_iidsmry_short <- data.frame("Generation" = c(rep("F0", 3), rep("F1", 3), rep("F2", 3)),
                                        "Experimental_treatment_group" = c(rep(c("AhR2Mo - / BaP +", "AhR2Mo + / BaP -", "AhR2Mo + / BaP +"), 3)), 
                                        "Test_statistic" = ksshoal_iid_teststats,
                                        "P_values" = ksshoal_iid_pvals)
    ks_shoal_iidsmry_short$signif <- apply(ks_shoal_iidsmry_short, 1, 
                                   function(x) {ifelse(x[["P_values"]] <= 0.05, "*", "ns")})
    ks_shoal_iidsmry_short
    ##   Generation Experimental_treatment_group Test_statistic     P_values signif
    ## 1         F0             AhR2Mo - / BaP +      0.3333333 0.0995624511     ns
    ## 2         F0             AhR2Mo + / BaP -      0.4074074 0.0206394820      *
    ## 3         F0             AhR2Mo + / BaP +      0.2592593 0.1684132191     ns
    ## 4         F1             AhR2Mo - / BaP +      0.1759259 0.0706872848     ns
    ## 5         F1             AhR2Mo + / BaP -      0.3333333 0.0005537942      *
    ## 6         F1             AhR2Mo + / BaP +      0.1703704 0.0994281426     ns
    ## 7         F2             AhR2Mo - / BaP +      0.1111111 0.4694024858     ns
    ## 8         F2             AhR2Mo + / BaP -      0.1239316 0.3542931756     ns
    ## 9         F2             AhR2Mo + / BaP +      0.2336182 0.0043525418      *

-   graphing the eCDF by generation with lines as treatments

<!-- -->

    shoal_iid_cdf_plot <- ggplot(shoal_iid_long_outlierd, aes(x = InterIndivDist, 
                                                              color = treatment)) +
      stat_ecdf(geom = "step") +
      theme_classic() +
      facet_wrap("generation") +
      labs(x = "Shoaling Inter-individual Distance", 
           y = "Time (min)", 
           color = "Treatment") +
      scale_color_brewer(palette = "PuOr", 
                         name = "Treatments", 
                        labels = c("AhR2Mo + / BaP +",
                                     "AhR2Mo + / BaP -",
                                     "AhR2Mo - / BaP +",
                                     "AhR2Mo - / BaP -"))
    shoal_iid_cdf_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/CDF graph for shoaling iid data-1.png" width="98%" height="98%" />

## Shoaling: speed

-   average speed of all four fish within the shoal

-   test assumptions of statistical model are met: non-normal data, f(x)
    is non-decreasing and right-continuous

<!-- -->

    # check normality assumption/what type of distribution
    # histogram (for all distances pooled across time points)
    mu_shoal_speed <- ddply(shoal_speed_long_outlierd, .(generation, treatment), summarise, 
                         grp.mean = mean(speed, na.rm = TRUE))
    mu_shoal_speed
    ##    generation  treatment  grp.mean
    ## 1          F0  AHRMo_BaP 0.8052929
    ## 2          F0 AHRMo_DMSO 0.6423457
    ## 3          F0   CoMo_BaP 0.5979111
    ## 4          F0  CoMo_DMSO 0.5952407
    ## 5          F1  AHRMo_BaP 4.3618333
    ## 6          F1 AHRMo_DMSO 4.7823175
    ## 7          F1   CoMo_BaP 4.6058426
    ## 8          F1  CoMo_DMSO 4.2983889
    ## 9          F2  AHRMo_BaP 5.0081453
    ## 10         F2 AHRMo_DMSO 5.0220000
    ## 11         F2   CoMo_BaP 4.4509444
    ## 12         F2  CoMo_DMSO 4.3162650


    ggplot(data = shoal_speed_long_outlierd, 
           aes(x = speed, color = generation, fill = generation)) +
      geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2) +
      geom_vline(data = mu_shoal_speed, aes(xintercept = grp.mean, color = generation),
               linetype = "dashed") +
      geom_density(alpha = 0.1) + 
      facet_wrap(~ treatment, ncol = 2)

<img src="zfBaP_behavior_files/figure-markdown_strict/test assumptions and check distribution shoaling speed-1.png" width="98%" height="98%" />

-   Two sample K-S test comparing experimental treatments to control
    group where p-value &lt;= 0.05 means the two data sets likely came
    from different distributions

<!-- -->

    # summary list for full stats for easy reference later
    ks_shoal_speedsmry <- c()

    ## F0 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_speedsmry$F0_BaPeffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_CoMo_BaP"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_speedsmry$F0_AhR2effect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_AHRMo_DMSO"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_speedsmry$F0_intereffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F0_AHRMo_BaP"]) # from different distributions


    ## F1 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_speedsmry$F1_BaPeffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_CoMo_BaP"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_speedsmry$F1_AhR2effect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_speedsmry$F1_intereffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F1_AHRMo_BaP"]) # not from different distributions


    ## F2 generation tests

    # K-S test comparing CoMo_DMSO to CoMo_BaP
    ks_shoal_speedsmry$F2_BaPeffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_CoMo_BaP"]) # not from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_DMSO
    ks_shoal_speedsmry$F2_AhR2effect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_AHRMo_DMSO"]) # from different distributions

    # K-S test comparing CoMo_DMSO to AhR2_BaP
    ks_shoal_speedsmry$F2_intereffect <- ks.test(shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_CoMo_DMSO"],
                                       shoal_speed_long_outlierd$speed[shoal_speed_long_outlierd$group %in% "F2_AHRMo_BaP"]) # from different distributions

    ## make summary data frame of just the useful stats for publishing
    # list of test stats
    ksshoal_speed_teststats <- c(ks_shoal_speedsmry$F0_BaPeffect$statistic, 
                              ks_shoal_speedsmry$F0_AhR2effect$statistic,
                              ks_shoal_speedsmry$F0_intereffect$statistic,
                              ks_shoal_speedsmry$F1_BaPeffect$statistic, 
                              ks_shoal_speedsmry$F1_AhR2effect$statistic,
                              ks_shoal_speedsmry$F1_intereffect$statistic, 
                              ks_shoal_speedsmry$F2_BaPeffect$statistic, 
                              ks_shoal_speedsmry$F2_AhR2effect$statistic,
                              ks_shoal_speedsmry$F2_intereffect$statistic)
    # list of p-values
    ksshoal_speed_pvals <- c(ks_shoal_speedsmry$F0_BaPeffect$p.value, 
                              ks_shoal_speedsmry$F0_AhR2effect$p.value,
                              ks_shoal_speedsmry$F0_intereffect$p.value,
                              ks_shoal_speedsmry$F1_BaPeffect$p.value, 
                              ks_shoal_speedsmry$F1_AhR2effect$p.value,
                              ks_shoal_speedsmry$F1_intereffect$p.value, 
                              ks_shoal_speedsmry$F2_BaPeffect$p.value, 
                              ks_shoal_speedsmry$F2_AhR2effect$p.value,
                              ks_shoal_speedsmry$F2_intereffect$p.value)

    # data frame of all summary stats of KS test for quick, easy reference and publishing
    ks_shoal_speedsmry_short <- data.frame("Generation" = c(rep("F0", 3), rep("F1", 3), rep("F2", 3)),
                                        "Experimental_treatment_group" = c(rep(c("AhR2Mo - / BaP +", "AhR2Mo + / BaP -", "AhR2Mo + / BaP +"), 3)), 
                                        "Test_statistic" = ksshoal_speed_teststats,
                                        "P_values" = ksshoal_speed_pvals)
    ks_shoal_speedsmry_short$signif <- apply(ks_shoal_speedsmry_short, 1, 
                                   function(x) {ifelse(x[["P_values"]] <= 0.05, "*", "ns")})
    ks_shoal_speedsmry_short
    ##   Generation Experimental_treatment_group Test_statistic     P_values signif
    ## 1         F0             AhR2Mo - / BaP +     0.29259259 2.304482e-02     ns
    ## 2         F0             AhR2Mo + / BaP -     0.17283951 2.511765e-01     ns
    ## 3         F0             AhR2Mo + / BaP +     0.38383838 3.641625e-05     ns
    ## 4         F1             AhR2Mo - / BaP +     0.21296296 1.492077e-02     ns
    ## 5         F1             AhR2Mo + / BaP -     0.39814815 3.676915e-06     ns
    ## 6         F1             AhR2Mo + / BaP +     0.10925926 5.565010e-01     ns
    ## 7         F2             AhR2Mo - / BaP +     0.08058608 8.256580e-01     ns
    ## 8         F2             AhR2Mo + / BaP -     0.31562882 1.125985e-05     ns
    ## 9         F2             AhR2Mo + / BaP +     0.32478632 8.730105e-06     ns

    ## NOTE: for some reason the above signif column code not working just on this behavior measure??

-   graphing the eCDF by generation with lines as treatments

<!-- -->

    shoal_speed_cdf_plot <- ggplot(shoal_speed_long_outlierd, aes(x = speed, color = treatment)) +
      stat_ecdf(geom = "step") +
      theme_classic() +
      facet_wrap("generation") +
      labs(x = "Shoaling Speed", y = "Time (min)", color = "Treatment") +
      scale_color_brewer(palette = "PuOr", 
                         name = "Treatments", 
                        labels = c("AhR2Mo + / BaP +",
                                     "AhR2Mo + / BaP -",
                                     "AhR2Mo - / BaP +",
                                     "AhR2Mo - / BaP -"))
    shoal_speed_cdf_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/CDF graph for shoaling speed data-1.png" width="98%" height="98%" />

## summary graph of all behavior measures

    behavior_cdfplots <- ggarrange(freeswim_cdf_plot, shoal_nnd_cdf_plot, 
                                   shoal_iid_cdf_plot, shoal_speed_cdf_plot,
                                   common.legend = TRUE,
                                   legend = "bottom",
                                   ncol = 1, nrow = 4)
    behavior_cdfplots

<img src="zfBaP_behavior_files/figure-markdown_strict/all behavior measures graph per generation-1.png" width="98%" height="98%" />


    # ggsave(paste0(output_dir, "behavior_cdfplots.tiff"),
    #        width = 7,
    #        height = 8,
    #        dpi = 300)

## Calculate and graph AUC for each metric

# Freeswim

    ggplot(data = freeswim_long_outlierd, 
           aes(x = timepoint, y = Distance, group = treatment)) +
      geom_point(aes(color = treatment)) +
      geom_smooth(formula = y ~ x, method = "loess", se = TRUE,
                  linetype = "solid", aes(color = treatment)) +
      facet_wrap(~ generation, ncol = 1)

<img src="zfBaP_behavior_files/figure-markdown_strict/graph raw Freeswim data with best fit line-1.png" width="98%" height="98%" />

    library(DescTools)

    freeswim_cdf <- lapply(split(freeswim_long_outlierd$Distance, 
                           freeswim_long_outlierd$FishNumber), ecdf)
    freeswim_cdf_y <- sapply(freeswim_cdf, 
                             function(e) e(freeswim_long_outlierd$Distance))

    freeswim_long_AUC <- as.data.frame(freeswim_cdf_y) %>%
      mutate(x = 1:nrow(as.data.frame(freeswim_cdf_y))) %>%
      gather(key = "FishNumber", value = "CDF_y", 1:123) %>%
      mutate(FishNumber = as.numeric(FishNumber)) %>%
      inner_join(freeswim_long_outlierd, by = "FishNumber", multiple = "all") %>%
      group_by(group, FishNumber) %>%
      dplyr::summarise(AUC = AUC(x, CDF_y, na.rm = TRUE)) # default method is trapezoid
    freeswim_long_AUC
    ## # A tibble: 123 × 3
    ## # Groups:   group [12]
    ##    group        FishNumber   AUC
    ##    <chr>             <dbl> <dbl>
    ##  1 F0_AHRMo_BaP          1  894.
    ##  2 F0_AHRMo_BaP          2  943.
    ##  3 F0_AHRMo_BaP          3  964.
    ##  4 F0_AHRMo_BaP          4  964.
    ##  5 F0_AHRMo_BaP          5  900.
    ##  6 F0_AHRMo_BaP          6  892.
    ##  7 F0_AHRMo_BaP          7  960.
    ##  8 F0_AHRMo_BaP          8 1072.
    ##  9 F0_AHRMo_BaP          9  960.
    ## 10 F0_AHRMo_BaP         10 1053.
    ## # ℹ 113 more rows

    # linear model
    freeswim_wide_outlierd <- spread(freeswim_long_outlierd, timepoint, Distance)
    freeswim_lm_AUC <- freeswim_long_AUC %>%
      inner_join(freeswim_wide_outlierd, by = "FishNumber") %>%
      mutate(treatment = str_replace(treatment, "CoMo _BaP", "CoMo_BaP"), 
             treats = treatment) %>%
      extract(treats, into = c("Morpholino", "Exposure"), "(.*)_([^_]+)$") %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"),
             Morpholino = relevel(factor(Morpholino), "CoMo"),
             Treatments = relevel(factor(treatment), "CoMo_DMSO"))

    # now filter and run linear models per generation
    freeswim_lm_AUCF0 <- dplyr::filter(freeswim_lm_AUC, generation == "F0")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = freeswim_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = freeswim_lm_AUCF0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -114.228  -49.843    0.378   37.423  111.600 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   986.37      27.56  35.794   <2e-16 ***
    ## ExposureBaP                   -34.00      38.97  -0.872    0.391    
    ## MorpholinoAHRMo               -13.81      35.58  -0.388    0.701    
    ## ExposureBaP:MorpholinoAHRMo    21.50      49.81   0.432    0.669    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 67.5 on 27 degrees of freedom
    ## Multiple R-squared:  0.03371,    Adjusted R-squared:  -0.07365 
    ## F-statistic: 0.314 on 3 and 27 DF,  p-value: 0.8151

    summary(lm(AUC ~ Treatments, data = freeswim_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = freeswim_lm_AUCF0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -114.228  -49.843    0.378   37.423  111.600 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            986.37      27.56  35.794   <2e-16 ***
    ## TreatmentsAHRMo_BaP    -26.30      34.86  -0.755    0.457    
    ## TreatmentsAHRMo_DMSO   -13.81      35.58  -0.388    0.701    
    ## TreatmentsCoMo_BaP     -34.00      38.97  -0.872    0.391    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 67.5 on 27 degrees of freedom
    ## Multiple R-squared:  0.03371,    Adjusted R-squared:  -0.07365 
    ## F-statistic: 0.314 on 3 and 27 DF,  p-value: 0.8151


    freeswim_lm_AUCF1 <- dplyr::filter(freeswim_lm_AUC, generation == "F1")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = freeswim_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = freeswim_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -281.47 -155.46  -15.15  121.91  408.40 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   355.44      57.81   6.148  4.4e-07 ***
    ## ExposureBaP                    86.59      81.76   1.059    0.297    
    ## MorpholinoAHRMo               153.06     100.13   1.529    0.135    
    ## ExposureBaP:MorpholinoAHRMo  -216.98     131.83  -1.646    0.108    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 200.3 on 36 degrees of freedom
    ## Multiple R-squared:  0.07462,    Adjusted R-squared:  -0.002496 
    ## F-statistic: 0.9676 on 3 and 36 DF,  p-value: 0.4186

    summary(lm(AUC ~ Treatments, data = freeswim_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = freeswim_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -281.47 -155.46  -15.15  121.91  408.40 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            355.44      57.81   6.148  4.4e-07 ***
    ## TreatmentsAHRMo_BaP     22.66      85.75   0.264    0.793    
    ## TreatmentsAHRMo_DMSO   153.06     100.13   1.529    0.135    
    ## TreatmentsCoMo_BaP      86.59      81.76   1.059    0.297    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 200.3 on 36 degrees of freedom
    ## Multiple R-squared:  0.07462,    Adjusted R-squared:  -0.002496 
    ## F-statistic: 0.9676 on 3 and 36 DF,  p-value: 0.4186


    freeswim_lm_AUCF2 <- dplyr::filter(freeswim_lm_AUC, generation == "F2")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = freeswim_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = freeswim_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -334.99 -107.70   -5.39  133.60  314.29 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   446.32      46.18   9.664 7.66e-13 ***
    ## ExposureBaP                   -30.22      62.94  -0.480    0.633    
    ## MorpholinoAHRMo                21.09      64.05   0.329    0.743    
    ## ExposureBaP:MorpholinoAHRMo   -91.95      88.88  -1.035    0.306    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 160 on 48 degrees of freedom
    ## Multiple R-squared:  0.08234,    Adjusted R-squared:  0.02499 
    ## F-statistic: 1.436 on 3 and 48 DF,  p-value: 0.244

    summary(lm(AUC ~ Treatments, data = freeswim_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = freeswim_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -334.99 -107.70   -5.39  133.60  314.29 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            446.32      46.18   9.664 7.66e-13 ***
    ## TreatmentsAHRMo_BaP   -101.08      64.05  -1.578    0.121    
    ## TreatmentsAHRMo_DMSO    21.09      64.05   0.329    0.743    
    ## TreatmentsCoMo_BaP     -30.22      62.94  -0.480    0.633    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 160 on 48 degrees of freedom
    ## Multiple R-squared:  0.08234,    Adjusted R-squared:  0.02499 
    ## F-statistic: 1.436 on 3 and 48 DF,  p-value: 0.244


    # AUC doesn't vary with any of our treatment variables. boo.

    freeswim_lm_AUC$Treatments <- factor(x = freeswim_lm_AUC$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHRMo_DMSO",
                                            "AHRMo_BaP"))

    freeswimAUC_plot <- ggplot(data = freeswim_lm_AUC,
           aes(x = Treatments, y = AUC)) +
      geom_boxplot(aes(fill = Treatments), outlier.shape = NA) +
      theme_classic() +
      labs(title="Free swim Distance", 
           x = "", y = "AUC (cm)", fill = "Treatment") +
      theme(text = element_text(size = 20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuOr", direction = -1,
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      facet_wrap(~ generation, ncol = 3, scales = "free_y")
    freeswimAUC_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Freeswim boxplots-1.png" width="98%" height="98%" />

# Shoaling nnd

    ggplot(data = shoal_nnd_long_outlierd, 
           aes(x = timepoint, y = NearestNeighborDist, group = treatment)) +
      geom_point(aes(color = treatment)) +
      geom_smooth(formula = y ~ x, method = "loess", se = TRUE,
                  linetype = "solid", aes(color = treatment)) +
      facet_wrap(~ generation, ncol = 1)

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling nnd raw data with best fit line-1.png" width="98%" height="98%" />

    library(DescTools)

    shoal_nnd_cdf <- lapply(split(shoal_nnd_long_outlierd$NearestNeighborDist, 
                           shoal_nnd_long_outlierd$FishNumber), ecdf)
    shoal_nnd_cdf_y <- sapply(shoal_nnd_cdf, 
                             function(e) e(shoal_nnd_long_outlierd$NearestNeighborDist))

    shoal_nnd_long_AUC <- as.data.frame(shoal_nnd_cdf_y) %>%
      mutate(x = 1:nrow(as.data.frame(shoal_nnd_cdf_y))) %>%
      gather(key = "FishNumber", value = "CDF_y", 1:132) %>%
      mutate(FishNumber = as.numeric(FishNumber)) %>%
      inner_join(shoal_nnd_long_outlierd, by = "FishNumber", multiple = "all") %>%
      group_by(group, FishNumber) %>%
      dplyr::summarise(AUC = AUC(x, CDF_y, na.rm = TRUE)) # default method is trapezoid
    shoal_nnd_long_AUC
    ## # A tibble: 132 × 3
    ## # Groups:   group [12]
    ##    group        FishNumber   AUC
    ##    <chr>             <dbl> <dbl>
    ##  1 F0_AHRMo_BaP          1  886.
    ##  2 F0_AHRMo_BaP          2  907.
    ##  3 F0_AHRMo_BaP          3 1022.
    ##  4 F0_AHRMo_BaP          4 1022.
    ##  5 F0_AHRMo_BaP          5  930.
    ##  6 F0_AHRMo_BaP          6  883.
    ##  7 F0_AHRMo_BaP          7 1017.
    ##  8 F0_AHRMo_BaP          8 1017.
    ##  9 F0_AHRMo_BaP          9 1130.
    ## 10 F0_AHRMo_BaP         10  979.
    ## # ℹ 122 more rows

    # linear model
    shoal_nnd_wide_outlierd <- spread(shoal_nnd_long_outlierd, timepoint, NearestNeighborDist)
    shoal_nnd_lm_AUC <- shoal_nnd_long_AUC %>%
      inner_join(shoal_nnd_wide_outlierd, by = "FishNumber") %>%
      mutate(treatment = as.factor(treatment), 
             treats = treatment) %>%
      extract(treats, into = c("Morpholino", "Exposure"), "(.*)_([^_]+)$") %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"),
             Morpholino = relevel(factor(Morpholino), "CoMo"),
             Treatments = relevel(factor(treatment), "CoMo_DMSO"))

    # now filter and run linear models per generation
    shoal_nnd_lm_AUCF0 <- dplyr::filter(shoal_nnd_lm_AUC, generation == "F0")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_nnd_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_nnd_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -122.24  -67.26   11.99   25.43  161.49 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 1020.074     27.080  37.668   <2e-16 ***
    ## ExposureBaP                  -45.012     38.297  -1.175    0.248    
    ## MorpholinoAHRMo               34.611     38.297   0.904    0.372    
    ## ExposureBaP:MorpholinoAHRMo   -4.326     52.441  -0.082    0.935    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 81.24 on 35 degrees of freedom
    ## Multiple R-squared:  0.1141, Adjusted R-squared:  0.03814 
    ## F-statistic: 1.502 on 3 and 35 DF,  p-value: 0.231

    summary(lm(AUC ~ Treatments, data = shoal_nnd_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_nnd_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -122.24  -67.26   11.99   25.43  161.49 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           1020.07      27.08  37.668   <2e-16 ***
    ## TreatmentsAHRMo_BaP    -14.73      35.82  -0.411    0.684    
    ## TreatmentsAHRMo_DMSO    34.61      38.30   0.904    0.372    
    ## TreatmentsCoMo_BaP     -45.01      38.30  -1.175    0.248    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 81.24 on 35 degrees of freedom
    ## Multiple R-squared:  0.1141, Adjusted R-squared:  0.03814 
    ## F-statistic: 1.502 on 3 and 35 DF,  p-value: 0.231


    shoal_nnd_lm_AUCF1 <- dplyr::filter(shoal_nnd_lm_AUC, generation == "F1")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_nnd_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_nnd_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -344.36 -168.43  -36.53  147.49  424.04 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   370.69      61.68   6.010 6.74e-07 ***
    ## ExposureBaP                    62.55      87.23   0.717    0.478    
    ## MorpholinoAHRMo               159.28     106.83   1.491    0.145    
    ## ExposureBaP:MorpholinoAHRMo  -219.96     140.65  -1.564    0.127    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 213.7 on 36 degrees of freedom
    ## Multiple R-squared:  0.07058,    Adjusted R-squared:  -0.00687 
    ## F-statistic: 0.9113 on 3 and 36 DF,  p-value: 0.4452

    summary(lm(AUC ~ Treatments, data = shoal_nnd_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_nnd_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -344.36 -168.43  -36.53  147.49  424.04 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           370.694     61.680   6.010 6.74e-07 ***
    ## TreatmentsAHRMo_BaP     1.872     91.487   0.020    0.984    
    ## TreatmentsAHRMo_DMSO  159.278    106.834   1.491    0.145    
    ## TreatmentsCoMo_BaP     62.551     87.229   0.717    0.478    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 213.7 on 36 degrees of freedom
    ## Multiple R-squared:  0.07058,    Adjusted R-squared:  -0.00687 
    ## F-statistic: 0.9113 on 3 and 36 DF,  p-value: 0.4452


    shoal_nnd_lm_AUCF2 <- dplyr::filter(shoal_nnd_lm_AUC, generation == "F2")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_nnd_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_nnd_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -317.81  -77.54   28.43  117.08  304.74 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                  410.590     45.262   9.071 4.62e-12 ***
    ## ExposureBaP                    6.394     62.857   0.102    0.919    
    ## MorpholinoAHRMo               81.038     64.010   1.266    0.211    
    ## ExposureBaP:MorpholinoAHRMo -122.993     89.712  -1.371    0.177    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 163.2 on 49 degrees of freedom
    ## Multiple R-squared:  0.06695,    Adjusted R-squared:  0.009827 
    ## F-statistic: 1.172 on 3 and 49 DF,  p-value: 0.33

    summary(lm(AUC ~ Treatments, data = shoal_nnd_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_nnd_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -317.81  -77.54   28.43  117.08  304.74 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           410.590     45.262   9.071 4.62e-12 ***
    ## TreatmentsAHRMo_BaP   -35.560     64.010  -0.556    0.581    
    ## TreatmentsAHRMo_DMSO   81.038     64.010   1.266    0.211    
    ## TreatmentsCoMo_BaP      6.394     62.857   0.102    0.919    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 163.2 on 49 degrees of freedom
    ## Multiple R-squared:  0.06695,    Adjusted R-squared:  0.009827 
    ## F-statistic: 1.172 on 3 and 49 DF,  p-value: 0.33


    # no significance

    shoal_nnd_lm_AUC$Treatments <- factor(x = shoal_nnd_lm_AUC$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHRMo_DMSO",
                                            "AHRMo_BaP"))

    shoal_nndAUC_plot <- ggplot(data = shoal_nnd_lm_AUC,
           aes(x = Treatments, y = AUC)) +
      geom_boxplot(aes(fill = Treatments), outlier.shape = NA) +
      theme_classic() +
      labs(title="Shoaling nnd", 
           x = "", y = "AUC (cm)", fill = "Treatment") +
      theme(text = element_text(size = 20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuOr", direction = -1,
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      facet_wrap(~ generation, ncol = 3, scales = "free_y")
    shoal_nndAUC_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling nnd boxplots-1.png" width="98%" height="98%" />

# Shoaling iid

    ggplot(data = shoal_iid_long_outlierd, 
           aes(x = timepoint, y = InterIndivDist, group = treatment)) +
      geom_point(aes(color = treatment)) +
      geom_smooth(formula = y ~ x, method = "loess", se = TRUE,
                  linetype = "solid", aes(color = treatment)) +
      facet_wrap(~ generation, ncol = 1)

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling iid raw data with best fit line-1.png" width="98%" height="98%" />

    library(DescTools)

    shoal_iid_cdf <- lapply(split(shoal_iid_long_outlierd$InterIndivDist, 
                           shoal_iid_long_outlierd$FishNumber), ecdf)
    shoal_iid_cdf_y <- sapply(shoal_iid_cdf, 
                             function(e) e(shoal_iid_long_outlierd$InterIndivDist))

    shoal_iid_long_AUC <- as.data.frame(shoal_iid_cdf_y) %>%
      mutate(x = 1:nrow(as.data.frame(shoal_iid_cdf_y))) %>%
      gather(key = "FishNumber", value = "CDF_y", 1:107) %>%
      mutate(FishNumber = as.numeric(FishNumber)) %>%
      inner_join(shoal_iid_long_outlierd, by = "FishNumber", multiple = "all") %>%
      group_by(group, FishNumber) %>%
      dplyr::summarise(AUC = AUC(x, CDF_y, na.rm = TRUE)) # default method is trapezoid
    shoal_iid_long_AUC
    ## # A tibble: 107 × 3
    ## # Groups:   group [12]
    ##    group         FishNumber   AUC
    ##    <chr>              <dbl> <dbl>
    ##  1 F0_AHRMo_BaP           1  884.
    ##  2 F0_AHRMo_BaP           2  893.
    ##  3 F0_AHRMo_BaP           3  854.
    ##  4 F0_AHRMo_BaP           4  892.
    ##  5 F0_AHRMo_BaP           5  894.
    ##  6 F0_AHRMo_BaP           6  938 
    ##  7 F0_AHRMo_DMSO          7  913.
    ##  8 F0_AHRMo_DMSO          8  921.
    ##  9 F0_AHRMo_DMSO          9  921.
    ## 10 F0_CoMo _BaP          10  858.
    ## # ℹ 97 more rows

    # linear model
    shoal_iid_wide_outlierd <- spread(shoal_iid_long_outlierd, timepoint, InterIndivDist)
    shoal_iid_lm_AUC <- shoal_iid_long_AUC %>%
      inner_join(shoal_iid_wide_outlierd, by = "FishNumber") %>%
      mutate(treatment = str_replace(treatment, "CoMo _BaP", "CoMo_BaP"), 
             treats = treatment) %>%
      extract(treats, into = c("Morpholino", "Exposure"), "(.*)_([^_]+)$") %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"),
             Morpholino = relevel(factor(Morpholino), "CoMo"),
             Treatments = relevel(factor(treatment), "CoMo_DMSO"))

    # now filter and run linear models per generation
    shoal_iid_lm_AUCF0 <- dplyr::filter(shoal_iid_lm_AUC, generation == "F0")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_iid_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_iid_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -40.519 -15.370   0.556   2.500  62.926 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   892.74      18.01  49.569 2.77e-14 ***
    ## ExposureBaP                   -12.20      25.47  -0.479    0.641    
    ## MorpholinoAHRMo                25.70      25.47   1.009    0.335    
    ## ExposureBaP:MorpholinoAHRMo   -13.69      33.69  -0.406    0.692    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 31.19 on 11 degrees of freedom
    ## Multiple R-squared:  0.1785, Adjusted R-squared:  -0.04556 
    ## F-statistic: 0.7966 on 3 and 11 DF,  p-value: 0.521

    summary(lm(AUC ~ Treatments, data = shoal_iid_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_iid_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -40.519 -15.370   0.556   2.500  62.926 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          892.7407    18.0101  49.569 2.77e-14 ***
    ## TreatmentsAHRMo_BaP   -0.1852    22.0577  -0.008    0.993    
    ## TreatmentsAHRMo_DMSO  25.7037    25.4701   1.009    0.335    
    ## TreatmentsCoMo_BaP   -12.2037    25.4701  -0.479    0.641    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 31.19 on 11 degrees of freedom
    ## Multiple R-squared:  0.1785, Adjusted R-squared:  -0.04556 
    ## F-statistic: 0.7966 on 3 and 11 DF,  p-value: 0.521


    shoal_iid_lm_AUCF1 <- dplyr::filter(shoal_iid_lm_AUC, generation == "F1")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_iid_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_iid_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -281.47 -155.46  -15.15  121.91  408.40 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   355.44      57.81   6.148  4.4e-07 ***
    ## ExposureBaP                    86.59      81.76   1.059    0.297    
    ## MorpholinoAHRMo               153.06     100.13   1.529    0.135    
    ## ExposureBaP:MorpholinoAHRMo  -216.98     131.83  -1.646    0.108    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 200.3 on 36 degrees of freedom
    ## Multiple R-squared:  0.07462,    Adjusted R-squared:  -0.002496 
    ## F-statistic: 0.9676 on 3 and 36 DF,  p-value: 0.4186

    summary(lm(AUC ~ Treatments, data = shoal_iid_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_iid_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -281.47 -155.46  -15.15  121.91  408.40 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            355.44      57.81   6.148  4.4e-07 ***
    ## TreatmentsAHRMo_BaP     22.66      85.75   0.264    0.793    
    ## TreatmentsAHRMo_DMSO   153.06     100.13   1.529    0.135    
    ## TreatmentsCoMo_BaP      86.59      81.76   1.059    0.297    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 200.3 on 36 degrees of freedom
    ## Multiple R-squared:  0.07462,    Adjusted R-squared:  -0.002496 
    ## F-statistic: 0.9676 on 3 and 36 DF,  p-value: 0.4186


    shoal_iid_lm_AUCF2 <- dplyr::filter(shoal_iid_lm_AUC, generation == "F2")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_iid_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_iid_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -334.99 -107.70   -5.39  133.60  314.29 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   446.32      46.18   9.664 7.66e-13 ***
    ## ExposureBaP                   -30.22      62.94  -0.480    0.633    
    ## MorpholinoAHRMo                21.09      64.05   0.329    0.743    
    ## ExposureBaP:MorpholinoAHRMo   -91.95      88.88  -1.035    0.306    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 160 on 48 degrees of freedom
    ## Multiple R-squared:  0.08234,    Adjusted R-squared:  0.02499 
    ## F-statistic: 1.436 on 3 and 48 DF,  p-value: 0.244

    summary(lm(AUC ~ Treatments, data = shoal_iid_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_iid_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -334.99 -107.70   -5.39  133.60  314.29 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            446.32      46.18   9.664 7.66e-13 ***
    ## TreatmentsAHRMo_BaP   -101.08      64.05  -1.578    0.121    
    ## TreatmentsAHRMo_DMSO    21.09      64.05   0.329    0.743    
    ## TreatmentsCoMo_BaP     -30.22      62.94  -0.480    0.633    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 160 on 48 degrees of freedom
    ## Multiple R-squared:  0.08234,    Adjusted R-squared:  0.02499 
    ## F-statistic: 1.436 on 3 and 48 DF,  p-value: 0.244


    # no significance

    shoal_iid_lm_AUC$Treatments <- factor(x = shoal_iid_lm_AUC$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHRMo_DMSO",
                                            "AHRMo_BaP"))

    shoal_iidAUC_plot <- ggplot(data = shoal_iid_lm_AUC,
           aes(x = Treatments, y = AUC)) +
      geom_boxplot(aes(fill = Treatments), outlier.shape = NA) +
      theme_classic() +
      labs(title="Shoaling iid", 
           x = "", y = "AUC (cm)", fill = "Treatment") +
      theme(text = element_text(size = 20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuOr", direction = -1,
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      facet_wrap(~ generation, ncol = 3, scales = "free_y")
    shoal_iidAUC_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling iid boxplots-1.png" width="98%" height="98%" />

# Shoaling speed

    ggplot(data = shoal_speed_long_outlierd, 
           aes(x = timepoint, y = speed, group = treatment)) +
      geom_point(aes(color = treatment)) +
      geom_smooth(formula = y ~ x, method = "loess", se = TRUE,
                  linetype = "solid", aes(color = treatment)) +
      facet_wrap(~ generation, ncol = 1)

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling speed raw data with best fit line-1.png" width="98%" height="98%" />

    library(DescTools)

    shoal_speed_cdf <- lapply(split(shoal_speed_long_outlierd$speed, 
                           shoal_speed_long_outlierd$FishNumber), ecdf)
    shoal_speed_cdf_y <- sapply(shoal_speed_cdf, 
                             function(e) e(shoal_speed_long_outlierd$speed))

    shoal_speed_long_AUC <- as.data.frame(shoal_speed_cdf_y) %>%
      mutate(x = 1:nrow(as.data.frame(shoal_speed_cdf_y))) %>%
      gather(key = "FishNumber", value = "CDF_y", 1:126) %>%
      mutate(FishNumber = as.numeric(FishNumber)) %>%
      inner_join(shoal_speed_long_outlierd, by = "FishNumber", multiple = "all") %>%
      group_by(group, FishNumber) %>%
      dplyr::summarise(AUC = AUC(x, CDF_y, na.rm = TRUE)) # default method is trapezoid
    shoal_speed_long_AUC
    ## # A tibble: 126 × 3
    ## # Groups:   group [12]
    ##    group        FishNumber   AUC
    ##    <chr>             <dbl> <dbl>
    ##  1 F0_AHRMo_BaP          1  906 
    ##  2 F0_AHRMo_BaP          2  895.
    ##  3 F0_AHRMo_BaP          3  965 
    ##  4 F0_AHRMo_BaP          4  965 
    ##  5 F0_AHRMo_BaP          5 1007.
    ##  6 F0_AHRMo_BaP          6  864.
    ##  7 F0_AHRMo_BaP          7  963.
    ##  8 F0_AHRMo_BaP          8  963.
    ##  9 F0_AHRMo_BaP          9 1071.
    ## 10 F0_AHRMo_BaP         10  883 
    ## # ℹ 116 more rows

    # linear model
    shoal_speed_wide_outlierd <- spread(shoal_speed_long_outlierd, timepoint, speed)
    shoal_speed_lm_AUC <- shoal_speed_long_AUC %>%
      inner_join(shoal_speed_wide_outlierd, by = "FishNumber") %>%
      mutate(treatment = as.factor(treatment), 
             treats = treatment) %>%
      extract(treats, into = c("Morpholino", "Exposure"), "(.*)_([^_]+)$") %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"),
             Morpholino = relevel(factor(Morpholino), "CoMo"),
             Treatments = relevel(factor(treatment), "CoMo_DMSO"))

    # now filter and run linear models per generation
    shoal_speed_lm_AUCF0 <- dplyr::filter(shoal_speed_lm_AUC, generation == "F0")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_speed_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_speed_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -87.818 -41.733   7.067  37.492 119.348 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 1012.815     22.617  44.781   <2e-16 ***
    ## ExposureBaP                   28.896     33.547   0.861   0.3966    
    ## MorpholinoAHRMo               -5.006     29.199  -0.171   0.8651    
    ## ExposureBaP:MorpholinoAHRMo  -84.720     41.778  -2.028   0.0525 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 55.4 on 27 degrees of freedom
    ## Multiple R-squared:  0.2948, Adjusted R-squared:  0.2164 
    ## F-statistic: 3.761 on 3 and 27 DF,  p-value: 0.02234

    summary(lm(AUC ~ Treatments, data = shoal_speed_lm_AUCF0))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_speed_lm_AUCF0)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -87.818 -41.733   7.067  37.492 119.348 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          1012.815     22.617  44.781   <2e-16 ***
    ## TreatmentsAHRMo_BaP   -60.830     28.117  -2.163   0.0395 *  
    ## TreatmentsAHRMo_DMSO   -5.006     29.199  -0.171   0.8651    
    ## TreatmentsCoMo_BaP     28.896     33.547   0.861   0.3966    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 55.4 on 27 degrees of freedom
    ## Multiple R-squared:  0.2948, Adjusted R-squared:  0.2164 
    ## F-statistic: 3.761 on 3 and 27 DF,  p-value: 0.02234


    shoal_speed_lm_AUCF1 <- dplyr::filter(shoal_speed_lm_AUC, generation == "F1")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_speed_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_speed_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -377.47 -166.29  -12.58  166.88  311.06 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   517.69      58.08   8.914 9.57e-11 ***
    ## ExposureBaP                   -85.98      82.13  -1.047   0.3020    
    ## MorpholinoAHRMo              -170.45      95.68  -1.781   0.0831 .  
    ## ExposureBaP:MorpholinoAHRMo   225.03     128.75   1.748   0.0888 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 201.2 on 37 degrees of freedom
    ## Multiple R-squared:  0.08812,    Adjusted R-squared:  0.01418 
    ## F-statistic: 1.192 on 3 and 37 DF,  p-value: 0.3262

    summary(lm(AUC ~ Treatments, data = shoal_speed_lm_AUCF1))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_speed_lm_AUCF1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -377.47 -166.29  -12.58  166.88  311.06 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            517.69      58.08   8.914 9.57e-11 ***
    ## TreatmentsAHRMo_BaP    -31.40      86.14  -0.365   0.7175    
    ## TreatmentsAHRMo_DMSO  -170.45      95.68  -1.781   0.0831 .  
    ## TreatmentsCoMo_BaP     -85.98      82.13  -1.047   0.3020    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 201.2 on 37 degrees of freedom
    ## Multiple R-squared:  0.08812,    Adjusted R-squared:  0.01418 
    ## F-statistic: 1.192 on 3 and 37 DF,  p-value: 0.3262


    shoal_speed_lm_AUCF2 <- dplyr::filter(shoal_speed_lm_AUC, generation == "F2")
    summary(lm(AUC ~ Exposure + Morpholino + Exposure:Morpholino, data = shoal_speed_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = shoal_speed_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -425.40 -169.67   -2.19  156.93  388.46 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   501.92      56.62   8.864  7.9e-12 ***
    ## ExposureBaP                   -15.74      78.64  -0.200   0.8422    
    ## MorpholinoAHRMo              -182.72      78.64  -2.324   0.0243 *  
    ## ExposureBaP:MorpholinoAHRMo    19.47     111.21   0.175   0.8617    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 204.2 on 50 degrees of freedom
    ## Multiple R-squared:  0.1626, Adjusted R-squared:  0.1124 
    ## F-statistic: 3.236 on 3 and 50 DF,  p-value: 0.02983

    summary(lm(AUC ~ Treatments, data = shoal_speed_lm_AUCF2))
    ## 
    ## Call:
    ## lm(formula = AUC ~ Treatments, data = shoal_speed_lm_AUCF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -425.40 -169.67   -2.19  156.93  388.46 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            501.92      56.62   8.864  7.9e-12 ***
    ## TreatmentsAHRMo_BaP   -178.98      80.08  -2.235   0.0299 *  
    ## TreatmentsAHRMo_DMSO  -182.72      78.64  -2.324   0.0243 *  
    ## TreatmentsCoMo_BaP     -15.74      78.64  -0.200   0.8422    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 204.2 on 50 degrees of freedom
    ## Multiple R-squared:  0.1626, Adjusted R-squared:  0.1124 
    ## F-statistic: 3.236 on 3 and 50 DF,  p-value: 0.02983


    # some significant patterns here

    shoal_speed_lm_AUC$Treatments <- factor(x = shoal_speed_lm_AUC$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHRMo_DMSO",
                                            "AHRMo_BaP"))

    shoal_speedAUC_plot <- ggplot(data = shoal_speed_lm_AUC,
           aes(x = Treatments, y = AUC)) +
      geom_boxplot(aes(fill = Treatments), outlier.shape = NA) +
      theme_classic() +
      labs(title="Shoaling speed", 
           x = "", y = "AUC (cm/min)", fill = "Treatment") +
      theme(text = element_text(size = 20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuOr", direction = -1,
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      facet_wrap(~ generation, ncol = 3, scales = "free_y")
    shoal_speedAUC_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/graph Shoaling speed boxplots-1.png" width="98%" height="98%" />

# combine graphs into multi-panel

    # combine into panels
    AUC_plot <- ggarrange(freeswimAUC_plot, shoal_nndAUC_plot, 
                          shoal_iidAUC_plot, shoal_speedAUC_plot, 
                          ncol = 1, nrow = 4,
                          align = "v", 
                          common.legend = TRUE, legend = "right")
    AUC_plot

<img src="zfBaP_behavior_files/figure-markdown_strict/could use this graph to show the KS test patterns-1.png" width="98%" height="98%" />
