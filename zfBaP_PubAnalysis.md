## Housekeeping

    # libraries
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(gridExtra)
    library(knitr)
    library(stringr)
    library(vegan)
    library(MASS)
    # set seed so stats like permanova are reproducible
    set.seed(3)

    # to knit an Rmd file, you need explicit file paths typed out, so I created these objects as shortcuts so I could call the object for commonly used file paths instead of retyping the whole thing each time
    home_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/"
    input_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/input_files_pub/"
    output_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/02_finalAnalysis/output_files_pub/"

    ## whole dataset (all generations) files

    # presence/absence of pathways in whole dataset
    pathcov_tab_ALL <- read.csv(file = paste0(input_dir, "pathcov_tab_ALL.txt"), 
                               sep = "\t",
                               header = TRUE)

    # abundance of pathways in whole dataset
    pathabund_tab_ALL <- read.csv(file = paste0(input_dir, "pathabund_tab_ALL.txt"), 
                                 sep = "\t",
                                 header = TRUE)

    # metadata for whole dataset
    metaALL <- read.csv(file = paste0(input_dir, "meta_tab_ALL.txt"),
                       sep = "\t",
                       header = TRUE) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), Morpholino = relevel(factor(Morpholino), "CoMo"),
             Generation = factor(Generation),
             Sex = factor(Sex), # relevel and factor so these are accurate
             BMI = Weight_g / (Length_mm/1000)^2) # add BMI calculation
    # reorder factors for treatment
    metaALL$Treatments <- factor(x = metaALL$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHR2Mo_DMSO",
                                            "AHR2Mo_BaP"))

    ## F0 files

    # presence/absence of pathways in F0
    pathcov_tab_F0 <- read.csv(file = paste0(input_dir, "pathcov_tab_F0.txt"), 
                               sep = "\t",
                               header = TRUE)

    # abundance of pathways in F0
    pathabund_tab_F0 <- read.csv(file = paste0(input_dir, "pathabund_tab_F0.txt"), 
                                 sep = "\t",
                                 header = TRUE)

    # metadata for F0
    metaF0 <- read.csv(file = paste0(input_dir, "meta_tab_F0.txt"),
                       sep = "\t",
                       header = TRUE) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), Morpholino = relevel(factor(Morpholino), "CoMo"),
             Generation = factor(Generation),
             Sex = factor(Sex), # relevel and factor so these are accurate
             BMI = Weight_g / (Length_mm/1000)^2) # add BMI calculation

    metaF0$Treatments <- factor(x = metaF0$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHR2Mo_DMSO",
                                            "AHR2Mo_BaP"))

    ## F1 files

    # presence/absence of pathways in F1
    pathcov_tab_F1 <- read.csv(file = paste0(input_dir, "pathcov_tab_F1.txt"), 
                               sep = "\t",
                               header = TRUE)

    # abundance of pathways in F1
    pathabund_tab_F1 <- read.csv(file = paste0(input_dir, "pathabund_tab_F1.txt"), 
                                 sep = "\t",
                                 header = TRUE)

    # metadata for F1
    metaF1 <- read.csv(file = paste0(input_dir, "meta_tab_F1.txt"),
                       sep = "\t",
                       header = TRUE) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), Morpholino = relevel(factor(Morpholino), "CoMo"),
             Generation = factor(Generation),
             Sex = factor(Sex), # relevel and factor so these are accurate
             BMI = Weight_g / (Length_mm/1000)^2) # add BMI calculation

    metaF1$Treatments <- factor(x = metaF1$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHR2Mo_DMSO",
                                            "AHR2Mo_BaP"))

    ## F2 files

    # presence/absence of pathways in F2
    pathcov_tab_F2 <- read.csv(file = paste0(input_dir, "pathcov_tab_F2.txt"), 
                               sep = "\t",
                               header = TRUE)

    # abundance of pathways in F2
    pathabund_tab_F2 <- read.csv(file = paste0(input_dir, "pathabund_tab_F2.txt"), 
                                 sep = "\t",
                                 header = TRUE)

    # metadata for F2
    metaF2 <- read.csv(file = paste0(input_dir, "meta_tab_F2.txt"),
                       sep = "\t",
                       header = TRUE) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), Morpholino = relevel(factor(Morpholino), "CoMo"),
             Generation = factor(Generation),
             Sex = factor(Sex), # relevel and factor so these are accurate
             BMI = Weight_g / (Length_mm/1000)^2) # add BMI calculation

    metaF2$Treatments <- factor(x = metaF2$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHR2Mo_DMSO",
                                            "AHR2Mo_BaP"))

## Alpha diversity stats and plots

# Whole dataset

Calculate alpha diversity metrics (richness, shannon, simpson) on whole
dataset (this means all generations worth of data). This modeling
approach uses “Generation” as a term, along with other covariates like
Sex and Treatment, in the linear models of different alpha diversity
metrics, with the input being the whole dataset. AIC is used to
determine the optimal model terms.

    pw_richness <- rowSums(pathcov_tab_ALL)
    pw_shannon <- diversity(pathabund_tab_ALL, "shannon", base = exp(1))
    pw_simpson <- diversity(pathabund_tab_ALL, "simpson")
    meta_modALL <- metaALL %>%
      dplyr::select(c("SampleID", "Exposure", "Morpholino", "Generation", "Sex", "Treatments", "BMI")) %>% # filter to only important variables we are testing
      mutate(Treatments = factor(Treatments))
    alphadiv <- data.frame(pw_richness, pw_shannon, pw_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(meta_modALL, by = "SampleID")
    alphadiv
    ##     SampleID pw_richness pw_shannon pw_simpson Exposure Morpholino Generation
    ## 1   BaP_0001          34   4.550223  0.9885438      BaP     AHR2Mo         F2
    ## 2   BaP_0002          42   4.502771  0.9877448      BaP     AHR2Mo         F2
    ## 3   BaP_0003          23   4.504218  0.9881009     DMSO     AHR2Mo         F0
    ## 4   BaP_0004          42   4.545643  0.9886415     DMSO     AHR2Mo         F0
    ## 5   BaP_0005          26   4.520909  0.9881479     DMSO     AHR2Mo         F0
    ## 6   BaP_0006          25   4.353807  0.9844620     DMSO     AHR2Mo         F0
    ## 7   BaP_0007          33   4.439526  0.9869343     DMSO     AHR2Mo         F0
    ## 8   BaP_0008          31   4.522775  0.9880482     DMSO     AHR2Mo         F0
    ## 9   BaP_0009          49   4.445836  0.9875717     DMSO       CoMo         F1
    ## 10  BaP_0010          50   4.493771  0.9875042     DMSO       CoMo         F1
    ## 11  BaP_0011          34   4.438030  0.9866638     DMSO       CoMo         F1
    ## 12  BaP_0012          34   4.547843  0.9885828     DMSO       CoMo         F1
    ## 13  BaP_0013          24   4.425747  0.9868110     DMSO       CoMo         F1
    ## 14  BaP_0014          46   4.445500  0.9868567     DMSO       CoMo         F1
    ## 15  BaP_0015          21   4.497418  0.9879694     DMSO       CoMo         F1
    ## 16  BaP_0016          34   4.573694  0.9883337     DMSO       CoMo         F1
    ## 17  BaP_0017          35   4.490649  0.9877335      BaP       CoMo         F1
    ## 18  BaP_0018          24   4.516835  0.9882476      BaP       CoMo         F1
    ## 19  BaP_0019          24   4.332636  0.9860365     DMSO     AHR2Mo         F1
    ## 20  BaP_0020          39   4.497835  0.9878455     DMSO     AHR2Mo         F1
    ## 21  BaP_0021          30   4.494431  0.9877907     DMSO     AHR2Mo         F1
    ## 22  BaP_0022          39   4.426520  0.9869192     DMSO     AHR2Mo         F1
    ## 23  BaP_0023          40   4.436349  0.9869666     DMSO     AHR2Mo         F1
    ## 24  BaP_0024          32   4.541916  0.9886768     DMSO     AHR2Mo         F1
    ## 25  BaP_0025          38   4.462421  0.9871796     DMSO       CoMo         F2
    ## 26  BaP_0026          38   4.393515  0.9865654     DMSO       CoMo         F2
    ## 27  BaP_0027          38   4.325567  0.9856425     DMSO       CoMo         F2
    ## 28  BaP_0028          35   4.385863  0.9863195     DMSO       CoMo         F2
    ## 29  BaP_0029          24   4.411667  0.9864888     DMSO       CoMo         F2
    ## 30  BaP_0030          42   4.405563  0.9867490     DMSO       CoMo         F2
    ## 31  BaP_0031          34   4.482217  0.9877155     DMSO       CoMo         F2
    ## 32  BaP_0032          42   4.427853  0.9870069     DMSO       CoMo         F2
    ## 33  BaP_0033          40   4.408737  0.9869127     DMSO       CoMo         F2
    ## 34  BaP_0034          33   4.495006  0.9878998     DMSO       CoMo         F2
    ## 35  BaP_0035          40   4.498149  0.9879687     DMSO       CoMo         F2
    ## 36  BaP_0036          27   4.554704  0.9888240     DMSO       CoMo         F2
    ## 37  BaP_0037          26   4.496289  0.9878816     DMSO     AHR2Mo         F0
    ## 38  BaP_0038          33   4.316099  0.9853482     DMSO     AHR2Mo         F0
    ## 39  BaP_0039          28   4.504699  0.9879460     DMSO     AHR2Mo         F0
    ## 40  BaP_0040          34   4.425233  0.9873059     DMSO     AHR2Mo         F0
    ## 41  BaP_0041          26   4.520501  0.9883291     DMSO     AHR2Mo         F0
    ## 42  BaP_0042          27   4.534710  0.9883182     DMSO     AHR2Mo         F0
    ## 43  BaP_0043          41   4.501306  0.9879252     DMSO     AHR2Mo         F0
    ## 44  BaP_0044          26   4.525067  0.9882462     DMSO     AHR2Mo         F0
    ## 45  BaP_0045          24   4.490652  0.9878204     DMSO       CoMo         F1
    ## 46  BaP_0046          27   4.514056  0.9880311     DMSO       CoMo         F1
    ## 47  BaP_0047          26   4.539308  0.9883887     DMSO       CoMo         F1
    ## 48  BaP_0048          25   4.530884  0.9882733     DMSO       CoMo         F1
    ## 49  BaP_0049          31   4.483722  0.9872031     DMSO       CoMo         F1
    ## 50  BaP_0050          24   4.462053  0.9870267     DMSO       CoMo         F1
    ## 51  BaP_0051          29   4.522418  0.9880867     DMSO       CoMo         F1
    ## 52  BaP_0052          41   4.452999  0.9874027     DMSO       CoMo         F1
    ## 53  BaP_0053          30   4.548316  0.9885620     DMSO     AHR2Mo         F1
    ## 54  BaP_0054          26   4.328906  0.9852158     DMSO     AHR2Mo         F1
    ## 55  BaP_0055          32   4.428729  0.9868857     DMSO     AHR2Mo         F1
    ## 56  BaP_0056          32   4.517757  0.9884123     DMSO     AHR2Mo         F1
    ## 57  BaP_0057          24   4.418993  0.9867686     DMSO     AHR2Mo         F1
    ## 58  BaP_0058          33   4.481396  0.9877116     DMSO     AHR2Mo         F1
    ## 59  BaP_0059          25   4.502718  0.9880775     DMSO     AHR2Mo         F1
    ## 60  BaP_0060          38   4.417238  0.9866347     DMSO     AHR2Mo         F1
    ## 61  BaP_0061          29   4.518762  0.9883991     DMSO       CoMo         F2
    ## 62  BaP_0062          39   4.388660  0.9862537     DMSO       CoMo         F2
    ## 63  BaP_0063          37   4.374832  0.9864720     DMSO       CoMo         F2
    ## 64  BaP_0064          39   4.395850  0.9864343     DMSO       CoMo         F2
    ## 65  BaP_0065          35   4.368612  0.9862677     DMSO       CoMo         F2
    ## 66  BaP_0066          33   4.377396  0.9860513     DMSO       CoMo         F2
    ## 67  BaP_0067          25   4.542407  0.9887608     DMSO       CoMo         F2
    ## 68  BaP_0068          38   4.394275  0.9863954     DMSO       CoMo         F2
    ## 69  BaP_0069          35   4.389691  0.9863123     DMSO       CoMo         F2
    ## 70  BaP_0070          30   4.540525  0.9886508     DMSO       CoMo         F2
    ## 71  BaP_0071          26   4.361893  0.9863425     DMSO       CoMo         F2
    ## 72  BaP_0072          33   4.434827  0.9869978      BaP       CoMo         F2
    ## 73  BaP_0073          39   4.328620  0.9855610     DMSO     AHR2Mo         F0
    ## 74  BaP_0074          39   4.230060  0.9841720     DMSO     AHR2Mo         F0
    ## 75  BaP_0075          38   4.373401  0.9860194      BaP     AHR2Mo         F0
    ## 76  BaP_0076          36   4.341884  0.9855663      BaP     AHR2Mo         F0
    ## 77  BaP_0077          34   4.329710  0.9855582      BaP     AHR2Mo         F0
    ## 78  BaP_0078          34   4.438752  0.9869682      BaP     AHR2Mo         F0
    ## 79  BaP_0079          31   4.378341  0.9861319      BaP     AHR2Mo         F0
    ## 80  BaP_0080          33   4.287327  0.9850763      BaP     AHR2Mo         F0
    ## 81  BaP_0081          39   4.294424  0.9852336     DMSO       CoMo         F1
    ## 82  BaP_0083          38   4.294128  0.9851863     DMSO       CoMo         F1
    ## 83  BaP_0084          29   4.215196  0.9838677     DMSO       CoMo         F1
    ## 84  BaP_0085          38   4.263761  0.9845001     DMSO       CoMo         F1
    ## 85  BaP_0086          35   4.508152  0.9878353     DMSO       CoMo         F1
    ## 86  BaP_0087          43   4.344175  0.9853647     DMSO       CoMo         F1
    ## 87  BaP_0088          47   4.256614  0.9843278     DMSO       CoMo         F1
    ## 88  BaP_0089          35   4.472215  0.9874612      BaP     AHR2Mo         F2
    ## 89  BaP_0090          40   4.468907  0.9870658      BaP     AHR2Mo         F2
    ## 90  BaP_0091          40   4.529356  0.9884797     DMSO     AHR2Mo         F1
    ## 91  BaP_0092          33   4.489661  0.9873577     DMSO     AHR2Mo         F1
    ## 92  BaP_0093          40   4.397413  0.9865212     DMSO     AHR2Mo         F1
    ## 93  BaP_0094          45   4.436673  0.9866033     DMSO     AHR2Mo         F1
    ## 94  BaP_0095          32   4.401144  0.9865592     DMSO     AHR2Mo         F1
    ## 95  BaP_0096          37   4.352675  0.9860732     DMSO     AHR2Mo         F1
    ## 96  BaP_0097          34   4.532889  0.9884338      BaP       CoMo         F2
    ## 97  BaP_0098          41   4.451153  0.9869110      BaP     AHR2Mo         F2
    ## 98  BaP_0099          32   4.442931  0.9871331      BaP     AHR2Mo         F1
    ## 99  BaP_0100          37   4.472857  0.9872593      BaP     AHR2Mo         F1
    ## 100 BaP_0101          38   4.437928  0.9870747      BaP     AHR2Mo         F1
    ## 101 BaP_0102          32   4.367063  0.9861272      BaP     AHR2Mo         F1
    ## 102 BaP_0103          38   4.338099  0.9863603      BaP     AHR2Mo         F1
    ## 103 BaP_0104          37   4.386232  0.9862100      BaP     AHR2Mo         F1
    ## 104 BaP_0105          24   4.395963  0.9867431      BaP       CoMo         F2
    ## 105 BaP_0106          42   4.465500  0.9871323      BaP       CoMo         F2
    ## 106 BaP_0107          30   4.428108  0.9869742      BaP       CoMo         F2
    ## 107 BaP_0108          38   4.444523  0.9867667      BaP       CoMo         F2
    ## 108 BaP_0109          26   4.541155  0.9884644      BaP     AHR2Mo         F0
    ## 109 BaP_0110          27   4.566991  0.9888320      BaP     AHR2Mo         F0
    ## 110 BaP_0111          34   4.573263  0.9886815      BaP     AHR2Mo         F0
    ## 111 BaP_0112          25   4.518203  0.9883509      BaP     AHR2Mo         F0
    ## 112 BaP_0113          27   4.518732  0.9881668      BaP     AHR2Mo         F0
    ## 113 BaP_0114          34   4.555842  0.9886733      BaP     AHR2Mo         F0
    ## 114 BaP_0115          24   4.540651  0.9886155      BaP     AHR2Mo         F0
    ## 115 BaP_0116          34   4.554625  0.9885314      BaP     AHR2Mo         F0
    ## 116 BaP_0117          28   4.545590  0.9884380     DMSO       CoMo         F1
    ## 117 BaP_0118          36   4.525876  0.9881573     DMSO       CoMo         F1
    ## 118 BaP_0119          32   4.548033  0.9886903     DMSO       CoMo         F1
    ## 119 BaP_0121          25   4.540934  0.9885798      BaP       CoMo         F1
    ## 120 BaP_0122          23   4.517661  0.9882370      BaP       CoMo         F1
    ## 121 BaP_0123          24   4.532811  0.9883483      BaP       CoMo         F1
    ## 122 BaP_0124          25   4.477051  0.9877192      BaP       CoMo         F1
    ## 123 BaP_0125          28   4.561343  0.9887294     DMSO     AHR2Mo         F1
    ## 124 BaP_0126          30   4.527927  0.9882307     DMSO     AHR2Mo         F1
    ## 125 BaP_0127          39   4.434399  0.9870181     DMSO     AHR2Mo         F1
    ## 126 BaP_0128          36   4.530216  0.9885140     DMSO     AHR2Mo         F1
    ## 127 BaP_0129          29   4.484446  0.9878066     DMSO     AHR2Mo         F1
    ## 128 BaP_0130          28   4.532226  0.9885346     DMSO     AHR2Mo         F1
    ## 129 BaP_0131          28   4.495650  0.9880066     DMSO     AHR2Mo         F1
    ## 130 BaP_0132          28   4.523465  0.9883102     DMSO     AHR2Mo         F1
    ## 131 BaP_0133          41   4.470522  0.9874441      BaP     AHR2Mo         F1
    ## 132 BaP_0134          36   4.524628  0.9883918      BaP     AHR2Mo         F1
    ## 133 BaP_0135          38   4.482283  0.9877355      BaP     AHR2Mo         F1
    ## 134 BaP_0136          38   4.513877  0.9882242      BaP     AHR2Mo         F1
    ## 135 BaP_0137          38   4.508693  0.9881372     DMSO       CoMo         F2
    ## 136 BaP_0138          33   4.528405  0.9884916     DMSO       CoMo         F2
    ## 137 BaP_0139          26   4.501309  0.9880717     DMSO       CoMo         F2
    ## 138 BaP_0140          27   4.526699  0.9884871     DMSO       CoMo         F2
    ## 139 BaP_0141          24   4.526364  0.9885085      BaP       CoMo         F2
    ## 140 BaP_0142          36   4.495184  0.9879463      BaP     AHR2Mo         F2
    ## 141 BaP_0143          35   4.477761  0.9876701      BaP       CoMo         F2
    ## 142 BaP_0144          39   4.538766  0.9882522      BaP       CoMo         F2
    ## 143 BaP_0433          31   4.521076  0.9884088      BaP     AHR2Mo         F0
    ## 144 BaP_0434          24   4.496731  0.9879677      BaP     AHR2Mo         F0
    ## 145 BaP_0435          42   4.508273  0.9881797      BaP     AHR2Mo         F0
    ## 146 BaP_0436          38   4.510948  0.9876939      BaP     AHR2Mo         F0
    ## 147 BaP_0437          37   4.515367  0.9881783      BaP     AHR2Mo         F0
    ## 148 BaP_0438          42   4.559864  0.9886423      BaP     AHR2Mo         F0
    ## 149 BaP_0439          39   4.571648  0.9886480      BaP     AHR2Mo         F0
    ## 150 BaP_0440          40   4.537667  0.9880528      BaP     AHR2Mo         F0
    ## 151 BaP_0441          43   4.517536  0.9882380      BaP     AHR2Mo         F2
    ## 152 BaP_0442          28   4.572098  0.9888556      BaP     AHR2Mo         F2
    ## 153 BaP_0443          38   4.494443  0.9877333      BaP       CoMo         F1
    ## 154 BaP_0444          28   4.447137  0.9871850      BaP       CoMo         F1
    ## 155 BaP_0445          24   4.517276  0.9880003      BaP       CoMo         F1
    ## 156 BaP_0446          32   4.583506  0.9889707      BaP       CoMo         F1
    ## 157 BaP_0447          25   4.561361  0.9887000      BaP       CoMo         F1
    ## 158 BaP_0448          37   4.579595  0.9887938      BaP       CoMo         F1
    ## 159 BaP_0450          41   4.484328  0.9877138     DMSO     AHR2Mo         F1
    ## 160 BaP_0451          45   4.456636  0.9872566     DMSO     AHR2Mo         F1
    ## 161 BaP_0452          36   4.528606  0.9882602     DMSO     AHR2Mo         F1
    ## 162 BaP_0453          28   4.517043  0.9882430     DMSO     AHR2Mo         F1
    ## 163 BaP_0454          30   4.575520  0.9889817     DMSO     AHR2Mo         F1
    ## 164 BaP_0455          38   4.512575  0.9882346      BaP     AHR2Mo         F1
    ## 165 BaP_0456          34   4.546822  0.9885168      BaP     AHR2Mo         F1
    ## 166 BaP_0457          27   4.461372  0.9873845      BaP       CoMo         F2
    ## 167 BaP_0458          38   4.565708  0.9888064     DMSO     AHR2Mo         F2
    ## 168 BaP_0459          42   4.502028  0.9880679     DMSO     AHR2Mo         F2
    ## 169 BaP_0460          41   4.516054  0.9882113     DMSO     AHR2Mo         F2
    ## 170 BaP_0461          39   4.530688  0.9884452     DMSO     AHR2Mo         F2
    ## 171 BaP_0462          35   4.586090  0.9890140     DMSO     AHR2Mo         F2
    ## 172 BaP_0463          30   4.531274  0.9884541     DMSO     AHR2Mo         F2
    ## 173 BaP_0464          34   4.533524  0.9882487     DMSO     AHR2Mo         F2
    ## 174 BaP_0465          30   4.517963  0.9881035     DMSO       CoMo         F2
    ## 175 BaP_0466          40   4.494190  0.9878048     DMSO       CoMo         F2
    ## 176 BaP_0467          25   4.525114  0.9884395     DMSO       CoMo         F2
    ## 177 BaP_0468          31   4.599095  0.9891943     DMSO       CoMo         F2
    ## 178 BaP_0541          26   4.457410  0.9870694      BaP     AHR2Mo         F0
    ## 179 BaP_0542          38   4.497898  0.9879342      BaP     AHR2Mo         F0
    ## 180 BaP_0543          38   4.537138  0.9885087      BaP     AHR2Mo         F0
    ## 181 BaP_0544          43   4.484337  0.9878241      BaP     AHR2Mo         F0
    ## 182 BaP_0545          27   4.484008  0.9877804      BaP     AHR2Mo         F0
    ## 183 BaP_0546          25   4.482312  0.9878118      BaP     AHR2Mo         F0
    ## 184 BaP_0547          30   4.573616  0.9888375      BaP     AHR2Mo         F0
    ## 185 BaP_0548          45   4.486966  0.9877766      BaP     AHR2Mo         F0
    ## 186 BaP_0549          38   4.536344  0.9885234      BaP       CoMo         F1
    ## 187 BaP_0550          36   4.588442  0.9890731      BaP       CoMo         F1
    ## 188 BaP_0551          40   4.525920  0.9884031      BaP       CoMo         F1
    ## 189 BaP_0552          44   4.448724  0.9871658      BaP       CoMo         F1
    ## 190 BaP_0553          35   4.560406  0.9887393      BaP       CoMo         F1
    ## 191 BaP_0554          36   4.536312  0.9883334      BaP       CoMo         F1
    ## 192 BaP_0555          25   4.528842  0.9885057      BaP       CoMo         F1
    ## 193 BaP_0556          25   4.525814  0.9884885      BaP       CoMo         F1
    ## 194 BaP_0557          25   4.501041  0.9880131      BaP     AHR2Mo         F1
    ## 195 BaP_0558          24   4.501790  0.9881099      BaP     AHR2Mo         F1
    ## 196 BaP_0559          27   4.561566  0.9887271      BaP     AHR2Mo         F1
    ## 197 BaP_0560          47   4.490269  0.9877926      BaP     AHR2Mo         F1
    ## 198 BaP_0561          28   4.576161  0.9888381      BaP     AHR2Mo         F1
    ## 199 BaP_0562          44   4.463849  0.9874716      BaP     AHR2Mo         F1
    ## 200 BaP_0563          24   4.509535  0.9882125      BaP     AHR2Mo         F1
    ## 201 BaP_0564          47   4.498697  0.9877141      BaP     AHR2Mo         F1
    ## 202 BaP_0565          30   4.548769  0.9885137     DMSO     AHR2Mo         F2
    ## 203 BaP_0566          39   4.543399  0.9885190     DMSO     AHR2Mo         F2
    ## 204 BaP_0567          24   4.520689  0.9883654     DMSO     AHR2Mo         F2
    ## 205 BaP_0568          41   4.481282  0.9876492     DMSO     AHR2Mo         F2
    ## 206 BaP_0569          38   4.573315  0.9887190     DMSO     AHR2Mo         F2
    ## 207 BaP_0570          49   4.371779  0.9859827     DMSO     AHR2Mo         F2
    ## 208 BaP_0571          26   4.518165  0.9883440     DMSO     AHR2Mo         F2
    ## 209 BaP_0572          43   4.421724  0.9868799     DMSO     AHR2Mo         F2
    ## 210 BaP_0573          34   4.529663  0.9884505      BaP       CoMo         F2
    ## 211 BaP_0574          42   4.481864  0.9877075      BaP       CoMo         F2
    ## 212 BaP_0575          39   4.542426  0.9886135      BaP       CoMo         F2
    ## 213 BaP_0576          42   4.489621  0.9878006      BaP       CoMo         F2
    ## 214 BaP_0649          25   4.514694  0.9882774      BaP     AHR2Mo         F2
    ## 215 BaP_0650          38   4.502033  0.9880862      BaP     AHR2Mo         F2
    ## 216 BaP_0651          41   4.530872  0.9883613      BaP     AHR2Mo         F0
    ## 217 BaP_0652          37   4.635700  0.9894872      BaP     AHR2Mo         F0
    ## 218 BaP_0653          38   4.511315  0.9880871      BaP     AHR2Mo         F0
    ## 219 BaP_0654          30   4.575460  0.9887534      BaP     AHR2Mo         F0
    ## 220 BaP_0655          40   4.538502  0.9883690      BaP     AHR2Mo         F0
    ## 221 BaP_0656          24   4.553642  0.9886642      BaP     AHR2Mo         F0
    ## 222 BaP_0657          42   4.571462  0.9887573      BaP       CoMo         F1
    ## 223 BaP_0658          38   4.683788  0.9896237      BaP       CoMo         F1
    ## 224 BaP_0659          26   4.529908  0.9882084      BaP       CoMo         F1
    ## 225 BaP_0660          26   4.553212  0.9886697      BaP       CoMo         F1
    ## 226 BaP_0661          43   4.483877  0.9875114      BaP       CoMo         F1
    ## 227 BaP_0662          47   4.514247  0.9881306      BaP       CoMo         F1
    ## 228 BaP_0663          41   4.446053  0.9870659      BaP       CoMo         F1
    ## 229 BaP_0664          39   4.603789  0.9890413      BaP       CoMo         F1
    ## 230 BaP_0665          42   4.443380  0.9875153      BaP     AHR2Mo         F1
    ## 231 BaP_0666          26   4.538411  0.9883892      BaP     AHR2Mo         F1
    ## 232 BaP_0667          48   4.359747  0.9861924      BaP     AHR2Mo         F1
    ## 233 BaP_0668          45   4.602605  0.9889341      BaP     AHR2Mo         F1
    ## 234 BaP_0669          44   4.494491  0.9879774      BaP     AHR2Mo         F1
    ## 235 BaP_0670          41   4.507746  0.9879028      BaP     AHR2Mo         F1
    ## 236 BaP_0671          39   4.580277  0.9889097      BaP     AHR2Mo         F1
    ## 237 BaP_0672          25   4.490591  0.9879534      BaP     AHR2Mo         F1
    ## 238 BaP_0673          48   4.498142  0.9880294      BaP       CoMo         F2
    ## 239 BaP_0674          30   4.549705  0.9886907      BaP       CoMo         F2
    ## 240 BaP_0675          40   4.557115  0.9886677      BaP       CoMo         F2
    ## 241 BaP_0676          37   4.552293  0.9886820      BaP       CoMo         F2
    ## 242 BaP_0677          32   4.478903  0.9878334      BaP       CoMo         F2
    ## 243 BaP_0678          40   4.578955  0.9888786      BaP       CoMo         F2
    ## 244 BaP_0679          29   4.505667  0.9880266      BaP       CoMo         F2
    ## 245 BaP_0680          28   4.577815  0.9888212      BaP       CoMo         F2
    ## 246 BaP_0681          43   4.594604  0.9888714      BaP       CoMo         F2
    ## 247 BaP_0683          38   4.587992  0.9888858      BaP       CoMo         F2
    ## 248 BaP_0757          33   4.592464  0.9890881     DMSO       CoMo         F1
    ## 249 BaP_0758          37   4.490714  0.9875079     DMSO       CoMo         F1
    ## 250 BaP_0759          47   4.503360  0.9881061     DMSO       CoMo         F1
    ## 251 BaP_0760          43   4.580983  0.9887079     DMSO       CoMo         F1
    ## 252 BaP_0761          36   4.604834  0.9890739     DMSO       CoMo         F1
    ## 253 BaP_0762          40   4.516106  0.9880403     DMSO       CoMo         F1
    ## 254 BaP_0763          30   4.561636  0.9886822     DMSO       CoMo         F1
    ## 255 BaP_0764          35   4.592782  0.9891507     DMSO       CoMo         F1
    ## 256 BaP_0765          31   4.521452  0.9883179      BaP       CoMo         F1
    ## 257 BaP_0766          43   4.585418  0.9887643      BaP       CoMo         F1
    ## 258 BaP_0767          46   4.475837  0.9877168      BaP       CoMo         F1
    ## 259 BaP_0768          30   4.593279  0.9891742      BaP       CoMo         F1
    ## 260 BaP_0769          25   4.544934  0.9884985      BaP       CoMo         F1
    ## 261 BaP_0770          50   4.487092  0.9877380      BaP       CoMo         F1
    ## 262 BaP_0771          24   4.504819  0.9878869      BaP       CoMo         F1
    ## 263 BaP_0772          40   4.783915  0.9906523      BaP       CoMo         F1
    ## 264 BaP_0773          29   4.547295  0.9883929      BaP     AHR2Mo         F1
    ## 265 BaP_0774          42   4.601739  0.9887933      BaP     AHR2Mo         F1
    ## 266 BaP_0775          35   4.531207  0.9883433      BaP     AHR2Mo         F1
    ## 267 BaP_0776          40   4.497533  0.9877301      BaP     AHR2Mo         F1
    ## 268 BaP_0777          42   4.554920  0.9885032      BaP     AHR2Mo         F1
    ## 269 BaP_0778          31   4.561511  0.9888869      BaP     AHR2Mo         F1
    ## 270 BaP_0779          38   4.581208  0.9886861      BaP     AHR2Mo         F1
    ## 271 BaP_0780          35   4.494911  0.9879318      BaP     AHR2Mo         F1
    ## 272 BaP_0781          33   4.501919  0.9881847      BaP       CoMo         F2
    ## 273 BaP_0782          33   4.532279  0.9885477      BaP       CoMo         F2
    ## 274 BaP_0783          23   4.414131  0.9856429      BaP       CoMo         F2
    ## 275 BaP_0784          43   4.494852  0.9878793      BaP       CoMo         F2
    ## 276 BaP_0785          42   4.573259  0.9887914      BaP       CoMo         F2
    ## 277 BaP_0786          33   4.456109  0.9873294      BaP       CoMo         F2
    ## 278 BaP_0787          40   4.608135  0.9890402      BaP       CoMo         F2
    ## 279 BaP_0788          36   4.558937  0.9886826      BaP       CoMo         F2
    ## 280 BaP_0789          40   4.501261  0.9880634      BaP       CoMo         F2
    ## 281 BaP_0790          26   4.514630  0.9881663      BaP       CoMo         F2
    ## 282 BaP_0791          46   4.632574  0.9893353      BaP       CoMo         F2
    ## 283 BaP_0792          41   4.529855  0.9883095      BaP       CoMo         F2
    ## 284 BaP_0865          37   4.430378  0.9870195      BaP     AHR2Mo         F2
    ## 285 BaP_0866          26   4.525773  0.9884628      BaP     AHR2Mo         F2
    ## 286 BaP_0867          42   4.537241  0.9885059     DMSO       CoMo         F0
    ## 287 BaP_0868          36   4.559420  0.9887496     DMSO       CoMo         F0
    ## 288 BaP_0869          40   4.563146  0.9884425     DMSO       CoMo         F0
    ## 289 BaP_0870          25   4.532752  0.9885157     DMSO       CoMo         F0
    ## 290 BaP_0871          25   4.585741  0.9888897     DMSO       CoMo         F0
    ## 291 BaP_0872          45   4.529618  0.9884092     DMSO       CoMo         F0
    ## 292 BaP_0873          23   4.536499  0.9885496     DMSO       CoMo         F0
    ## 293 BaP_0874          29   4.562748  0.9887055     DMSO       CoMo         F0
    ## 294 BaP_0875          49   4.528531  0.9882597     DMSO       CoMo         F0
    ## 295 BaP_0876          25   4.533867  0.9885804     DMSO       CoMo         F0
    ## 296 BaP_0877          33   4.573732  0.9887946     DMSO       CoMo         F0
    ## 297 BaP_0878          32   4.559254  0.9887499     DMSO       CoMo         F0
    ## 298 BaP_0879          38   4.540871  0.9885075      BaP       CoMo         F0
    ## 299 BaP_0880          26   4.613909  0.9889981      BaP       CoMo         F0
    ## 300 BaP_0881          28   4.513351  0.9882516      BaP       CoMo         F0
    ## 301 BaP_0882          27   4.570542  0.9887851      BaP       CoMo         F0
    ## 302 BaP_0883          26   4.587557  0.9887919      BaP       CoMo         F0
    ## 303 BaP_0884          36   4.517204  0.9881091      BaP       CoMo         F0
    ## 304 BaP_0885          41   4.545778  0.9882872      BaP       CoMo         F0
    ## 305 BaP_0886          25   4.550773  0.9886811      BaP       CoMo         F0
    ## 306 BaP_0887          26   4.565892  0.9886970      BaP       CoMo         F0
    ## 307 BaP_0888          26   4.536289  0.9885451      BaP       CoMo         F0
    ## 308 BaP_0889          27   4.499857  0.9879037     DMSO     AHR2Mo         F2
    ## 309 BaP_0890          27   4.550619  0.9885459     DMSO     AHR2Mo         F2
    ## 310 BaP_0891          42   4.596655  0.9888162     DMSO     AHR2Mo         F2
    ## 311 BaP_0893          40   4.436523  0.9871808     DMSO     AHR2Mo         F2
    ## 312 BaP_0895          26   4.533494  0.9884286     DMSO     AHR2Mo         F2
    ## 313 BaP_0897          26   4.520421  0.9881975     DMSO     AHR2Mo         F2
    ## 314 BaP_0899          25   4.519657  0.9882915     DMSO     AHR2Mo         F2
    ## 315 BaP_0937          24   4.522185  0.9882154     DMSO       CoMo         F0
    ## 316 BaP_0938          46   4.508775  0.9881818     DMSO       CoMo         F0
    ## 317 BaP_0939          28   4.534658  0.9885164     DMSO       CoMo         F0
    ## 318 BaP_0940          45   4.584213  0.9885080     DMSO       CoMo         F0
    ## 319 BaP_0941          48   4.403116  0.9866326     DMSO       CoMo         F0
    ## 320 BaP_0942          31   4.594294  0.9889972     DMSO       CoMo         F0
    ## 321 BaP_0943          23   4.383083  0.9859967     DMSO       CoMo         F0
    ## 322 BaP_0944          46   4.579474  0.9883979     DMSO       CoMo         F0
    ## 323 BaP_0945          24   4.493030  0.9879130      BaP       CoMo         F0
    ## 324 BaP_0946          42   4.524984  0.9879366      BaP       CoMo         F0
    ## 325 BaP_0947          24   4.520566  0.9882318      BaP       CoMo         F0
    ## 326 BaP_0948          26   4.509400  0.9879252      BaP       CoMo         F0
    ## 327 BaP_0949          36   4.510599  0.9881336      BaP       CoMo         F0
    ## 328 BaP_0950          26   4.548811  0.9885605      BaP       CoMo         F0
    ## 329 BaP_0951          37   4.536197  0.9885085      BaP       CoMo         F0
    ## 330 BaP_0952          43   4.514813  0.9875756      BaP       CoMo         F0
    ## 331 BaP_0953          37   4.551988  0.9883968      BaP       CoMo         F0
    ## 332 BaP_0954          29   4.511145  0.9881735      BaP       CoMo         F0
    ## 333 BaP_0955          25   4.573314  0.9887514      BaP       CoMo         F0
    ## 334 BaP_0956          35   4.494903  0.9877437      BaP       CoMo         F0
    ## 335 BaP_0957          28   4.498406  0.9878888     DMSO     AHR2Mo         F0
    ## 336 BaP_0958          26   4.412859  0.9866246     DMSO     AHR2Mo         F0
    ## 337 BaP_0959          28   4.573228  0.9886918     DMSO     AHR2Mo         F0
    ## 338 BaP_0960          45   4.542668  0.9884344     DMSO     AHR2Mo         F0
    ## 339 BaP_0961          32   4.529736  0.9880896     DMSO     AHR2Mo         F2
    ## 340 BaP_0962          24   4.486389  0.9876956     DMSO     AHR2Mo         F2
    ## 341 BaP_0963          24   4.464565  0.9873303     DMSO     AHR2Mo         F2
    ## 342 BaP_0964          41   4.492813  0.9876432     DMSO     AHR2Mo         F2
    ## 343 BaP_0965          25   4.536049  0.9884101     DMSO     AHR2Mo         F2
    ## 344 BaP_0966          26   4.529795  0.9884617     DMSO     AHR2Mo         F2
    ## 345 BaP_0967          48   4.673877  0.9895179     DMSO     AHR2Mo         F2
    ## 346 BaP_0968          27   4.538327  0.9884462     DMSO     AHR2Mo         F2
    ## 347 BaP_0969          36   4.471660  0.9875338      BaP     AHR2Mo         F2
    ## 348 BaP_0970          38   4.535278  0.9883476      BaP     AHR2Mo         F2
    ## 349 BaP_0971          53   4.470996  0.9875761      BaP     AHR2Mo         F2
    ## 350 BaP_0972          26   4.469410  0.9873433      BaP     AHR2Mo         F2
    ## 351 BaP_1009          29   4.464234  0.9874148     DMSO       CoMo         F0
    ## 352 BaP_1010          31   4.470033  0.9873883     DMSO       CoMo         F0
    ## 353 BaP_1011          27   4.452881  0.9873293     DMSO       CoMo         F0
    ## 354 BaP_1012          45   4.274714  0.9853417     DMSO       CoMo         F0
    ## 355 BaP_1013          25   4.446758  0.9871227     DMSO       CoMo         F0
    ## 356 BaP_1014          42   4.487749  0.9876448     DMSO       CoMo         F0
    ## 357 BaP_1015          25   4.514678  0.9881170     DMSO       CoMo         F0
    ## 358 BaP_1016          42   4.612949  0.9892287     DMSO       CoMo         F0
    ## 359 BaP_1017          27   4.531667  0.9885300      BaP     AHR2Mo         F2
    ## 360 BaP_1018          49   4.425298  0.9872458      BaP     AHR2Mo         F2
    ## 361 BaP_1019          34   4.539563  0.9886260      BaP       CoMo         F0
    ## 362 BaP_1020          28   4.514816  0.9881500      BaP       CoMo         F0
    ## 363 BaP_1021          25   4.581829  0.9888790      BaP       CoMo         F0
    ## 364 BaP_1022          43   4.600920  0.9891380      BaP       CoMo         F0
    ## 365 BaP_1023          40   4.442561  0.9872198      BaP       CoMo         F0
    ## 366 BaP_1024          41   4.545772  0.9885967      BaP       CoMo         F0
    ## 367 BaP_1025          40   4.521984  0.9880028     DMSO     AHR2Mo         F0
    ## 368 BaP_1026          24   4.484919  0.9878194     DMSO     AHR2Mo         F0
    ## 369 BaP_1027          40   4.504881  0.9880215     DMSO     AHR2Mo         F0
    ## 370 BaP_1028          39   4.448848  0.9876286     DMSO     AHR2Mo         F0
    ## 371 BaP_1029          39   4.572114  0.9888226     DMSO     AHR2Mo         F0
    ## 372 BaP_1030          37   4.555932  0.9885582     DMSO     AHR2Mo         F0
    ## 373 BaP_1031          31   4.559720  0.9887723     DMSO     AHR2Mo         F0
    ## 374 BaP_1032          29   4.542443  0.9883773     DMSO     AHR2Mo         F0
    ## 375 BaP_1033          28   4.546102  0.9883440     DMSO     AHR2Mo         F2
    ## 376 BaP_1034          38   4.304793  0.9851062     DMSO     AHR2Mo         F2
    ## 377 BaP_1035          41   4.549272  0.9881433     DMSO     AHR2Mo         F2
    ## 378 BaP_1036          25   4.540316  0.9884947     DMSO     AHR2Mo         F2
    ## 379 BaP_1037          37   4.382060  0.9860371     DMSO     AHR2Mo         F2
    ## 380 BaP_1038          27   4.546698  0.9885625     DMSO     AHR2Mo         F2
    ## 381 BaP_1039          44   4.573911  0.9884341      BaP     AHR2Mo         F2
    ## 382 BaP_1040          32   4.460502  0.9875155      BaP     AHR2Mo         F2
    ## 383 BaP_1041          31   4.506099  0.9879788      BaP     AHR2Mo         F2
    ## 384 BaP_1042          24   4.563573  0.9886930      BaP     AHR2Mo         F2
    ## 385 BaP_1043          32   4.420825  0.9869515      BaP     AHR2Mo         F2
    ## 386 BaP_1044          41   4.540977  0.9880741      BaP     AHR2Mo         F2
    ## 387 BaP_1081          30   4.579584  0.9889271     DMSO       CoMo         F0
    ## 388 BaP_1082          35   4.614369  0.9892223     DMSO       CoMo         F0
    ## 389 BaP_1083          40   4.449053  0.9867586     DMSO       CoMo         F0
    ## 390 BaP_1084          45   4.629706  0.9891099     DMSO       CoMo         F0
    ## 391 BaP_1085          24   4.448315  0.9872674     DMSO       CoMo         F0
    ## 392 BaP_1086          36   4.624406  0.9891998     DMSO       CoMo         F0
    ## 393 BaP_1087          28   4.563440  0.9886356     DMSO       CoMo         F0
    ## 394 BaP_1088          39   4.609836  0.9890670     DMSO       CoMo         F0
    ## 395 BaP_1089          25   4.598048  0.9891037      BaP       CoMo         F0
    ## 396 BaP_1090          35   4.641423  0.9891380      BaP       CoMo         F0
    ## 397 BaP_1091          24   4.467356  0.9875381      BaP       CoMo         F0
    ## 398 BaP_1092          40   4.630763  0.9893792      BaP       CoMo         F0
    ## 399 BaP_1093          25   4.579059  0.9887814      BaP       CoMo         F0
    ## 400 BaP_1094          30   4.510950  0.9881143      BaP       CoMo         F0
    ## 401 BaP_1095          28   4.531201  0.9883345      BaP       CoMo         F0
    ## 402 BaP_1096          37   4.575171  0.9887221      BaP       CoMo         F0
    ## 403 BaP_1097          40   4.655449  0.9895118     DMSO     AHR2Mo         F0
    ## 404 BaP_1098          34   4.568771  0.9886379     DMSO     AHR2Mo         F0
    ## 405 BaP_1099          45   4.616147  0.9891246     DMSO     AHR2Mo         F0
    ## 406 BaP_1100          36   4.577251  0.9886513     DMSO     AHR2Mo         F0
    ## 407 BaP_1101          41   4.599216  0.9888807     DMSO     AHR2Mo         F0
    ## 408 BaP_1102          34   4.594431  0.9886945     DMSO     AHR2Mo         F0
    ## 409 BaP_1103          29   4.650466  0.9894567     DMSO     AHR2Mo         F0
    ## 410 BaP_1104          30   4.425107  0.9865785     DMSO     AHR2Mo         F0
    ## 411 BaP_1105          28   4.591510  0.9889295      BaP     AHR2Mo         F2
    ## 412 BaP_1106          45   4.605341  0.9889344      BaP     AHR2Mo         F2
    ## 413 BaP_1107          26   4.534361  0.9882333      BaP     AHR2Mo         F2
    ## 414 BaP_1108          42   4.519372  0.9880783      BaP     AHR2Mo         F2
    ## 415 BaP_1109          28   4.540853  0.9883700      BaP     AHR2Mo         F2
    ## 416 BaP_1110          28   4.613805  0.9891114      BaP     AHR2Mo         F2
    ## 417 BaP_1111          26   4.478225  0.9876314      BaP     AHR2Mo         F2
    ## 418 BaP_1112          23   4.423231  0.9864562      BaP     AHR2Mo         F2
    ## 419 BaP_1113          26   4.561057  0.9884440      BaP     AHR2Mo         F2
    ## 420 BaP_1114          29   4.598979  0.9887733      BaP     AHR2Mo         F2
    ## 421 BaP_1115          40   4.589359  0.9887910      BaP     AHR2Mo         F2
    ## 422 BaP_1116          36   4.671099  0.9893012      BaP     AHR2Mo         F2
    ##     Sex  Treatments       BMI
    ## 1     F  AHR2Mo_BaP  764.4444
    ## 2     F  AHR2Mo_BaP  641.2742
    ## 3     F AHR2Mo_DMSO 1094.1828
    ## 4     F AHR2Mo_DMSO  620.8912
    ## 5     F AHR2Mo_DMSO 1069.3878
    ## 6     M AHR2Mo_DMSO  884.4953
    ## 7     F AHR2Mo_DMSO  897.9592
    ## 8     M AHR2Mo_DMSO  717.9931
    ## 9     F   CoMo_DMSO  671.4410
    ## 10    M   CoMo_DMSO  372.3974
    ## 11    F   CoMo_DMSO  490.8642
    ## 12    M   CoMo_DMSO  534.8007
    ## 13    F   CoMo_DMSO  498.6479
    ## 14    M   CoMo_DMSO  378.1163
    ## 15    F   CoMo_DMSO  636.5432
    ## 16    M   CoMo_DMSO  471.6049
    ## 17    F    CoMo_BaP  725.2569
    ## 18    M    CoMo_BaP  489.3750
    ## 19    F AHR2Mo_DMSO  534.0265
    ## 20    M AHR2Mo_DMSO  520.4082
    ## 21    F AHR2Mo_DMSO  650.2268
    ## 22    M AHR2Mo_DMSO  507.9221
    ## 23    F AHR2Mo_DMSO  695.9600
    ## 24    M AHR2Mo_DMSO  477.5000
    ## 25    F   CoMo_DMSO  927.4691
    ## 26    M   CoMo_DMSO  668.0428
    ## 27    F   CoMo_DMSO  604.6713
    ## 28    M   CoMo_DMSO  659.5253
    ## 29    F   CoMo_DMSO  818.5604
    ## 30    M   CoMo_DMSO  647.3829
    ## 31    F   CoMo_DMSO 1070.6696
    ## 32    M   CoMo_DMSO  634.7699
    ## 33    F   CoMo_DMSO  842.2206
    ## 34    M   CoMo_DMSO  682.6172
    ## 35    F   CoMo_DMSO  588.7500
    ## 36    M   CoMo_DMSO  534.6260
    ## 37    F AHR2Mo_DMSO  770.0312
    ## 38    M AHR2Mo_DMSO  604.0816
    ## 39    F AHR2Mo_DMSO  977.5087
    ## 40    F AHR2Mo_DMSO  674.7405
    ## 41    F AHR2Mo_DMSO  993.7500
    ## 42    M AHR2Mo_DMSO  579.5848
    ## 43    F AHR2Mo_DMSO  889.7959
    ## 44    M AHR2Mo_DMSO  694.4444
    ## 45    F   CoMo_DMSO  645.5577
    ## 46    M   CoMo_DMSO  473.7696
    ## 47    F   CoMo_DMSO  578.5124
    ## 48    M   CoMo_DMSO  529.8765
    ## 49    F   CoMo_DMSO  510.3306
    ## 50    M   CoMo_DMSO  489.5895
    ## 51    F   CoMo_DMSO  660.7407
    ## 52    M   CoMo_DMSO  463.7500
    ## 53    F AHR2Mo_DMSO  571.2992
    ## 54    M AHR2Mo_DMSO  531.8519
    ## 55    F AHR2Mo_DMSO  642.8000
    ## 56    M AHR2Mo_DMSO  394.4773
    ## 57    F AHR2Mo_DMSO  792.9804
    ## 58    M AHR2Mo_DMSO  473.2288
    ## 59    F AHR2Mo_DMSO  615.7845
    ## 60    M AHR2Mo_DMSO  468.2540
    ## 61    F   CoMo_DMSO  787.1094
    ## 62    M   CoMo_DMSO  709.0947
    ## 63    F   CoMo_DMSO 1053.1250
    ## 64    M   CoMo_DMSO  551.1111
    ## 65    F   CoMo_DMSO  949.5465
    ## 66    M   CoMo_DMSO  589.9654
    ## 67    F   CoMo_DMSO  605.5515
    ## 68    M   CoMo_DMSO  528.3556
    ## 69    F   CoMo_DMSO  949.2326
    ## 70    F   CoMo_DMSO  831.9559
    ## 71    F   CoMo_DMSO  852.4470
    ## 72    F    CoMo_BaP  727.5383
    ## 73    F AHR2Mo_DMSO  683.3910
    ## 74    M AHR2Mo_DMSO  558.8585
    ## 75    F  AHR2Mo_BaP  942.2936
    ## 76    M  AHR2Mo_BaP  673.8281
    ## 77    F  AHR2Mo_BaP 1175.8585
    ## 78    M  AHR2Mo_BaP  651.9743
    ## 79    F  AHR2Mo_BaP  895.0617
    ## 80    M  AHR2Mo_BaP  718.0021
    ## 81    F   CoMo_DMSO  631.6928
    ## 82    F   CoMo_DMSO  655.0475
    ## 83    M   CoMo_DMSO  468.2540
    ## 84    F   CoMo_DMSO  536.8000
    ## 85    M   CoMo_DMSO  560.6576
    ## 86    F   CoMo_DMSO  697.1336
    ## 87    M   CoMo_DMSO  500.2469
    ## 88    M  AHR2Mo_BaP  549.4644
    ## 89    M  AHR2Mo_BaP  568.6728
    ## 90    F AHR2Mo_DMSO  583.7037
    ## 91    M AHR2Mo_DMSO  461.6300
    ## 92    F AHR2Mo_DMSO  617.7686
    ## 93    M AHR2Mo_DMSO  524.7934
    ## 94    F AHR2Mo_DMSO  691.8715
    ## 95    M AHR2Mo_DMSO  449.5465
    ## 96    F    CoMo_BaP  854.7816
    ## 97    F  AHR2Mo_BaP  827.1209
    ## 98    F  AHR2Mo_BaP  721.7882
    ## 99    M  AHR2Mo_BaP  482.4509
    ## 100   F  AHR2Mo_BaP  620.1901
    ## 101   M  AHR2Mo_BaP  390.0000
    ## 102   F  AHR2Mo_BaP  551.1364
    ## 103   M  AHR2Mo_BaP  420.6250
    ## 104   F    CoMo_BaP  701.2418
    ## 105   M    CoMo_BaP  569.3918
    ## 106   M    CoMo_BaP  588.0140
    ## 107   F    CoMo_BaP  897.5069
    ## 108   F  AHR2Mo_BaP  748.4568
    ## 109   M  AHR2Mo_BaP  622.8374
    ## 110   F  AHR2Mo_BaP  663.5802
    ## 111   M  AHR2Mo_BaP  551.5088
    ## 112   F  AHR2Mo_BaP  648.7889
    ## 113   M  AHR2Mo_BaP  693.8776
    ## 114   F  AHR2Mo_BaP  688.7052
    ## 115   M  AHR2Mo_BaP  605.4688
    ## 116   F   CoMo_DMSO  601.1342
    ## 117   M   CoMo_DMSO  513.7500
    ## 118   F   CoMo_DMSO  613.8453
    ## 119   F    CoMo_BaP  680.0000
    ## 120   M    CoMo_BaP  467.6871
    ## 121   F    CoMo_BaP  699.4329
    ## 122   M    CoMo_BaP  300.0000
    ## 123   F AHR2Mo_DMSO  643.0503
    ## 124   M AHR2Mo_DMSO  451.6765
    ## 125   F AHR2Mo_DMSO  555.0023
    ## 126   M AHR2Mo_DMSO  485.0000
    ## 127   F AHR2Mo_DMSO  719.0083
    ## 128   M AHR2Mo_DMSO  522.9854
    ## 129   F AHR2Mo_DMSO  574.3802
    ## 130   M AHR2Mo_DMSO  495.3086
    ## 131   F  AHR2Mo_BaP  682.3347
    ## 132   M  AHR2Mo_BaP  542.5342
    ## 133   F  AHR2Mo_BaP  599.0123
    ## 134   M  AHR2Mo_BaP  456.2760
    ## 135   F   CoMo_DMSO  868.7700
    ## 136   M   CoMo_DMSO  577.4222
    ## 137   F   CoMo_DMSO 1242.2145
    ## 138   M   CoMo_DMSO  619.1467
    ## 139   F    CoMo_BaP  955.6250
    ## 140   M  AHR2Mo_BaP  669.5502
    ## 141   M    CoMo_BaP  513.8504
    ## 142   F    CoMo_BaP  772.2681
    ## 143   F  AHR2Mo_BaP  717.9931
    ## 144   M  AHR2Mo_BaP  634.7656
    ## 145   F  AHR2Mo_BaP  770.0312
    ## 146   M  AHR2Mo_BaP  596.8858
    ## 147   F  AHR2Mo_BaP  832.3424
    ## 148   M  AHR2Mo_BaP  685.7143
    ## 149   F  AHR2Mo_BaP  778.5467
    ## 150   M  AHR2Mo_BaP  709.3426
    ## 151   F  AHR2Mo_BaP  839.3352
    ## 152   F  AHR2Mo_BaP  709.3750
    ## 153   F    CoMo_BaP  557.6560
    ## 154   M    CoMo_BaP  468.1737
    ## 155   F    CoMo_BaP  587.1605
    ## 156   M    CoMo_BaP  514.3750
    ## 157   F    CoMo_BaP  947.8306
    ## 158   M    CoMo_BaP  546.6984
    ## 159   M AHR2Mo_DMSO  474.7174
    ## 160   F AHR2Mo_DMSO  734.1270
    ## 161   M AHR2Mo_DMSO  510.0054
    ## 162   F AHR2Mo_DMSO  501.8929
    ## 163   F AHR2Mo_DMSO  668.6420
    ## 164   F  AHR2Mo_BaP  602.9630
    ## 165   M  AHR2Mo_BaP  423.6669
    ## 166   M    CoMo_BaP  662.0753
    ## 167   F AHR2Mo_DMSO  889.4558
    ## 168   M AHR2Mo_DMSO  768.4666
    ## 169   F AHR2Mo_DMSO 1011.7729
    ## 170   M AHR2Mo_DMSO  623.0469
    ## 171   F AHR2Mo_DMSO  780.3688
    ## 172   M AHR2Mo_DMSO  687.3630
    ## 173   F AHR2Mo_DMSO 1048.9796
    ## 174   F   CoMo_DMSO  929.8037
    ## 175   M   CoMo_DMSO  654.6939
    ## 176   F   CoMo_DMSO 1035.3186
    ## 177   M   CoMo_DMSO  886.0482
    ## 178   F  AHR2Mo_BaP  749.1082
    ## 179   M  AHR2Mo_BaP  457.8564
    ## 180   F  AHR2Mo_BaP  848.1262
    ## 181   M  AHR2Mo_BaP  624.3496
    ## 182   M  AHR2Mo_BaP  596.8858
    ## 183   F  AHR2Mo_BaP  685.7143
    ## 184   F  AHR2Mo_BaP  749.2196
    ## 185   M  AHR2Mo_BaP  759.1837
    ## 186   F    CoMo_BaP  686.9835
    ## 187   M    CoMo_BaP  450.1385
    ## 188   M    CoMo_BaP  498.9669
    ## 189   M    CoMo_BaP  457.4830
    ## 190   F    CoMo_BaP  726.3705
    ## 191   M    CoMo_BaP  416.8750
    ## 192   F    CoMo_BaP  447.3140
    ## 193   M    CoMo_BaP  352.3997
    ## 194   F  AHR2Mo_BaP  741.0208
    ## 195   M  AHR2Mo_BaP  271.1634
    ## 196   F  AHR2Mo_BaP  540.0000
    ## 197   M  AHR2Mo_BaP  468.1250
    ## 198   F  AHR2Mo_BaP  504.3750
    ## 199   M  AHR2Mo_BaP  507.5000
    ## 200   F  AHR2Mo_BaP  649.3827
    ## 201   M  AHR2Mo_BaP  364.9584
    ## 202   M AHR2Mo_DMSO  755.5556
    ## 203   F AHR2Mo_DMSO  803.3241
    ## 204   M AHR2Mo_DMSO  612.7385
    ## 205   F AHR2Mo_DMSO  950.1385
    ## 206   M AHR2Mo_DMSO  641.8685
    ## 207   F AHR2Mo_DMSO 1121.1911
    ## 208   M AHR2Mo_DMSO  633.4964
    ## 209   F AHR2Mo_DMSO  980.9336
    ## 210   M    CoMo_BaP  529.4118
    ## 211   F    CoMo_BaP  976.2655
    ## 212   M    CoMo_BaP  609.4183
    ## 213   F    CoMo_BaP  683.4568
    ## 214   M  AHR2Mo_BaP  671.8646
    ## 215   M  AHR2Mo_BaP  558.0716
    ## 216   F  AHR2Mo_BaP  818.1154
    ## 217   M  AHR2Mo_BaP  587.6951
    ## 218   M  AHR2Mo_BaP  722.2222
    ## 219   M  AHR2Mo_BaP  596.8779
    ## 220   M  AHR2Mo_BaP  795.8478
    ## 221   M  AHR2Mo_BaP  645.1613
    ## 222   F    CoMo_BaP  587.2934
    ## 223   M    CoMo_BaP  480.3719
    ## 224   F    CoMo_BaP  591.2698
    ## 225   M    CoMo_BaP  444.3783
    ## 226   F    CoMo_BaP  546.2412
    ## 227   M    CoMo_BaP  524.7934
    ## 228   F    CoMo_BaP  658.7902
    ## 229   M    CoMo_BaP  487.8049
    ## 230   F  AHR2Mo_BaP  632.3617
    ## 231   M  AHR2Mo_BaP  430.0000
    ## 232   F  AHR2Mo_BaP  646.1777
    ## 233   M  AHR2Mo_BaP  437.0748
    ## 234   F  AHR2Mo_BaP  580.4989
    ## 235   M  AHR2Mo_BaP  465.9091
    ## 236   F  AHR2Mo_BaP  535.7143
    ## 237   M  AHR2Mo_BaP  427.9778
    ## 238   M    CoMo_BaP  570.6122
    ## 239   F    CoMo_BaP  967.8598
    ## 240   M    CoMo_BaP  605.7099
    ## 241   F    CoMo_BaP  804.4077
    ## 242   M    CoMo_BaP  654.7291
    ## 243   F    CoMo_BaP  820.4082
    ## 244   M    CoMo_BaP  616.5333
    ## 245   F    CoMo_BaP  869.3750
    ## 246   M    CoMo_BaP  593.1337
    ## 247   M    CoMo_BaP  659.3896
    ## 248   F   CoMo_DMSO  544.7846
    ## 249   M   CoMo_DMSO  527.9012
    ## 250   F   CoMo_DMSO  453.1445
    ## 251   M   CoMo_DMSO  261.2847
    ## 252   F   CoMo_DMSO  578.7654
    ## 253   M   CoMo_DMSO  417.7778
    ## 254   F   CoMo_DMSO  559.5238
    ## 255   M   CoMo_DMSO  461.4512
    ## 256   F    CoMo_BaP  563.5479
    ## 257   M    CoMo_BaP  578.7982
    ## 258   F    CoMo_BaP  609.9773
    ## 259   M    CoMo_BaP  488.0952
    ## 260   F    CoMo_BaP  676.5432
    ## 261   M    CoMo_BaP  457.5000
    ## 262   F    CoMo_BaP  628.6420
    ## 263   M    CoMo_BaP  414.3750
    ## 264   F  AHR2Mo_BaP  561.5705
    ## 265   M  AHR2Mo_BaP  466.8750
    ## 266   F  AHR2Mo_BaP  587.8685
    ## 267   M  AHR2Mo_BaP  443.1886
    ## 268   F  AHR2Mo_BaP  596.5909
    ## 269   M  AHR2Mo_BaP  571.9955
    ## 270   F  AHR2Mo_BaP  521.9037
    ## 271   M  AHR2Mo_BaP  472.0579
    ## 272   M    CoMo_BaP  567.4611
    ## 273   F    CoMo_BaP  798.2543
    ## 274   M    CoMo_BaP  770.4264
    ## 275   F    CoMo_BaP  936.8836
    ## 276   M    CoMo_BaP  600.4383
    ## 277   F    CoMo_BaP  831.7175
    ## 278   M    CoMo_BaP  601.8519
    ## 279   F    CoMo_BaP  708.5896
    ## 280   M    CoMo_BaP  695.3125
    ## 281   F    CoMo_BaP 1007.8895
    ## 282   M    CoMo_BaP  722.8724
    ## 283   F    CoMo_BaP  755.3035
    ## 284   F  AHR2Mo_BaP  681.6108
    ## 285   F  AHR2Mo_BaP  637.7383
    ## 286   F   CoMo_DMSO  993.9446
    ## 287   M   CoMo_DMSO  451.1719
    ## 288   F   CoMo_DMSO 1296.0761
    ## 289   M   CoMo_DMSO  674.0129
    ## 290   F   CoMo_DMSO 1029.9489
    ## 291   M   CoMo_DMSO  404.0404
    ## 292   F   CoMo_DMSO  824.4898
    ## 293   M   CoMo_DMSO  507.8125
    ## 294   F   CoMo_DMSO  787.1972
    ## 295   M   CoMo_DMSO  546.8750
    ## 296   F   CoMo_DMSO  734.6189
    ## 297   M   CoMo_DMSO  532.4074
    ## 298   F    CoMo_BaP  856.4815
    ## 299   F    CoMo_BaP  726.5306
    ## 300   F    CoMo_BaP  759.1837
    ## 301   M    CoMo_BaP  605.5363
    ## 302   F    CoMo_BaP  709.8765
    ## 303   M    CoMo_BaP 1020.7612
    ## 304   F    CoMo_BaP 1092.7456
    ## 305   M    CoMo_BaP  849.6094
    ## 306   F    CoMo_BaP  676.3788
    ## 307   M    CoMo_BaP  562.4543
    ## 308   M AHR2Mo_DMSO  501.9531
    ## 309   F AHR2Mo_DMSO  767.3130
    ## 310   M AHR2Mo_DMSO  709.2768
    ## 311   F AHR2Mo_DMSO  755.5402
    ## 312   M AHR2Mo_DMSO  613.8135
    ## 313   F AHR2Mo_DMSO 1015.3483
    ## 314   M AHR2Mo_DMSO  641.9753
    ## 315   F   CoMo_DMSO  697.1904
    ## 316   M   CoMo_DMSO  569.3297
    ## 317   F   CoMo_DMSO  830.4498
    ## 318   M   CoMo_DMSO  606.4209
    ## 319   F   CoMo_DMSO  743.9446
    ## 320   M   CoMo_DMSO  679.5225
    ## 321   F   CoMo_DMSO  712.8906
    ## 322   M   CoMo_DMSO  432.6531
    ## 323   F    CoMo_BaP  790.8429
    ## 324   M    CoMo_BaP  576.1719
    ## 325   F    CoMo_BaP  601.8519
    ## 326   M    CoMo_BaP  894.9011
    ## 327   F    CoMo_BaP  878.9063
    ## 328   M    CoMo_BaP  657.4394
    ## 329   F    CoMo_BaP  902.7778
    ## 330   M    CoMo_BaP  576.1719
    ## 331   F    CoMo_BaP 1009.3652
    ## 332   M    CoMo_BaP  556.6406
    ## 333   F    CoMo_BaP  957.3361
    ## 334   M    CoMo_BaP  830.0781
    ## 335   F AHR2Mo_DMSO  524.6914
    ## 336   M AHR2Mo_DMSO  460.1899
    ## 337   F AHR2Mo_DMSO  973.3701
    ## 338   M AHR2Mo_DMSO  437.0447
    ## 339   F AHR2Mo_DMSO  769.8882
    ## 340   M AHR2Mo_DMSO  615.5102
    ## 341   F AHR2Mo_DMSO  963.7188
    ## 342   M AHR2Mo_DMSO  664.4898
    ## 343   F AHR2Mo_DMSO  953.7500
    ## 344   M AHR2Mo_DMSO  649.4141
    ## 345   F AHR2Mo_DMSO 1033.5306
    ## 346   M AHR2Mo_DMSO  553.0649
    ## 347   F  AHR2Mo_BaP  689.4923
    ## 348   M  AHR2Mo_BaP  612.2449
    ## 349   F  AHR2Mo_BaP  886.4901
    ## 350   M  AHR2Mo_BaP  597.5510
    ## 351   F   CoMo_DMSO  622.8374
    ## 352   M   CoMo_DMSO  676.3788
    ## 353   M   CoMo_DMSO  615.2344
    ## 354   F   CoMo_DMSO  787.1972
    ## 355   F   CoMo_DMSO  900.2770
    ## 356   M   CoMo_DMSO  663.5802
    ## 357   F   CoMo_DMSO  910.4938
    ## 358   M   CoMo_DMSO  677.1861
    ## 359   M  AHR2Mo_BaP  708.5714
    ## 360   M  AHR2Mo_BaP  538.5802
    ## 361   F    CoMo_BaP  725.3086
    ## 362   M    CoMo_BaP  856.4014
    ## 363   F    CoMo_BaP  775.5102
    ## 364   M    CoMo_BaP  769.8962
    ## 365   F    CoMo_BaP  620.4082
    ## 366   M    CoMo_BaP  650.1096
    ## 367   F AHR2Mo_DMSO  942.9066
    ## 368   M AHR2Mo_DMSO  613.9438
    ## 369   F AHR2Mo_DMSO  640.4321
    ## 370   M AHR2Mo_DMSO  700.6920
    ## 371   F AHR2Mo_DMSO  783.6735
    ## 372   M AHR2Mo_DMSO  582.6397
    ## 373   F AHR2Mo_DMSO  824.0997
    ## 374   M AHR2Mo_DMSO  458.4775
    ## 375   F AHR2Mo_DMSO  733.7963
    ## 376   M AHR2Mo_DMSO  646.1938
    ## 377   F AHR2Mo_DMSO  779.7784
    ## 378   M AHR2Mo_DMSO  636.8410
    ## 379   F AHR2Mo_DMSO  925.7143
    ## 380   M AHR2Mo_DMSO  697.2318
    ## 381   F  AHR2Mo_BaP  789.6416
    ## 382   M  AHR2Mo_BaP  721.4533
    ## 383   F  AHR2Mo_BaP  799.1234
    ## 384   M  AHR2Mo_BaP  688.8889
    ## 385   F  AHR2Mo_BaP  907.8947
    ## 386   M  AHR2Mo_BaP  610.9964
    ## 387   F   CoMo_DMSO  832.6531
    ## 388   M   CoMo_DMSO  780.4370
    ## 389   F   CoMo_DMSO  759.6786
    ## 390   M   CoMo_DMSO  732.4219
    ## 391   F   CoMo_DMSO  634.7656
    ## 392   M   CoMo_DMSO  692.0415
    ## 393   F   CoMo_DMSO  918.2099
    ## 394   F   CoMo_DMSO  653.0612
    ## 395   F    CoMo_BaP  808.1633
    ## 396   M    CoMo_BaP  642.0927
    ## 397   F    CoMo_BaP  712.8906
    ## 398   M    CoMo_BaP  546.9388
    ## 399   F    CoMo_BaP  763.8889
    ## 400   M    CoMo_BaP  654.2969
    ## 401   F    CoMo_BaP  854.6384
    ## 402   M    CoMo_BaP  560.1469
    ## 403   F AHR2Mo_DMSO  748.4568
    ## 404   M AHR2Mo_DMSO  477.5023
    ## 405   F AHR2Mo_DMSO  761.2457
    ## 406   M AHR2Mo_DMSO  635.5004
    ## 407   F AHR2Mo_DMSO  840.8163
    ## 408   M AHR2Mo_DMSO  622.8374
    ## 409   M AHR2Mo_DMSO  613.9438
    ## 410   M AHR2Mo_DMSO  661.2245
    ## 411   F  AHR2Mo_BaP  949.9072
    ## 412   M  AHR2Mo_BaP  601.8519
    ## 413   F  AHR2Mo_BaP  556.1857
    ## 414   M  AHR2Mo_BaP  533.8776
    ## 415   F  AHR2Mo_BaP  973.5410
    ## 416   M  AHR2Mo_BaP  606.5306
    ## 417   F  AHR2Mo_BaP  951.2500
    ## 418   M  AHR2Mo_BaP  657.6075
    ## 419   F  AHR2Mo_BaP  803.0764
    ## 420   M  AHR2Mo_BaP  661.7969
    ## 421   F  AHR2Mo_BaP  803.9032
    ## 422   M  AHR2Mo_BaP  651.6620

What covariates predict microbiome richness?

    # step AIC to test which model
    testmod_rich <- lm(pw_richness ~ Generation + Exposure + Morpholino + Sex + BMI + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino, data = alphadiv)
    AIC_wholedata <- stepAIC(testmod_rich)
    ## Start:  AIC=1663.82
    ## pw_richness ~ Generation + Exposure + Morpholino + Sex + BMI + 
    ##     Generation:Exposure + Generation:Morpholino + Exposure:Morpholino
    ## 
    ##                         Df Sum of Sq   RSS    AIC
    ## - Generation:Morpholino  2    46.381 20601 1660.8
    ## - BMI                    1     0.562 20555 1661.8
    ## - Generation:Exposure    2   123.780 20678 1662.3
    ## <none>                               20554 1663.8
    ## - Sex                    1   109.139 20664 1664.0
    ## - Exposure:Morpholino    1   127.014 20681 1664.4
    ## 
    ## Step:  AIC=1660.77
    ## pw_richness ~ Generation + Exposure + Morpholino + Sex + BMI + 
    ##     Generation:Exposure + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq   RSS    AIC
    ## - BMI                  1      0.60 20601 1658.8
    ## - Generation:Exposure  2    126.43 20727 1659.3
    ## <none>                             20601 1660.8
    ## - Sex                  1    108.89 20710 1661.0
    ## - Exposure:Morpholino  1    125.93 20727 1661.3
    ## 
    ## Step:  AIC=1658.78
    ## pw_richness ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq   RSS    AIC
    ## - Generation:Exposure  2    126.09 20727 1657.4
    ## <none>                             20601 1658.8
    ## - Exposure:Morpholino  1    125.44 20727 1659.3
    ## - Sex                  1    171.61 20773 1660.3
    ## 
    ## Step:  AIC=1657.36
    ## pw_richness ~ Generation + Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq   RSS    AIC
    ## - Generation           2    142.85 20870 1656.2
    ## <none>                             20727 1657.4
    ## - Exposure:Morpholino  1    128.07 20856 1658.0
    ## - Sex                  1    169.72 20897 1658.8
    ## 
    ## Step:  AIC=1656.25
    ## pw_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq   RSS    AIC
    ## <none>                             20870 1656.2
    ## - Exposure:Morpholino  1    126.77 20997 1656.8
    ## - Sex                  1    169.66 21040 1657.7

    # extract lm that is best model for results
    lm_rich <- summary(AIC_wholedata)
    lm_rich
    ## 
    ## Call:
    ## lm(formula = pw_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino, 
    ##     data = alphadiv)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.766  -6.679   0.005   5.736  18.899 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   33.7658     0.7736  43.648   <2e-16 ***
    ## ExposureBaP                   -0.7708     0.9795  -0.787   0.4318    
    ## MorpholinoAHR2Mo              -1.0870     0.9861  -1.102   0.2709    
    ## SexM                           1.2692     0.6894   1.841   0.0663 .  
    ## ExposureBaP:MorpholinoAHR2Mo   2.1932     1.3780   1.592   0.1123    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.075 on 417 degrees of freedom
    ## Multiple R-squared:  0.01471,    Adjusted R-squared:  0.005263 
    ## F-statistic: 1.557 on 4 and 417 DF,  p-value: 0.185


    # best model is with Exposure + Morpholino + Sex + Exposure:Morpholino but Sex (and nothing else) is barely not significant

What covarites predict shannon richness?

    # step AIC to test which model
    testmod_shan <- lm(pw_shannon ~ Generation + Exposure + Morpholino + Sex + BMI + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino, data = alphadiv)
    AIC_wholedata <- stepAIC(testmod_shan)
    ## Start:  AIC=-2214.46
    ## pw_shannon ~ Generation + Exposure + Morpholino + Sex + BMI + 
    ##     Generation:Exposure + Generation:Morpholino + Exposure:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - BMI                    1  0.000927 2.0981 -2216.3
    ## - Sex                    1  0.003694 2.1009 -2215.7
    ## <none>                               2.0972 -2214.5
    ## - Generation:Exposure    2  0.025582 2.1228 -2213.3
    ## - Exposure:Morpholino    1  0.038638 2.1359 -2208.8
    ## - Generation:Morpholino  2  0.078949 2.1762 -2202.9
    ## 
    ## Step:  AIC=-2216.27
    ## pw_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## <none>                               2.0981 -2216.3
    ## - Sex                    1  0.011313 2.1094 -2216.0
    ## - Generation:Exposure    2  0.026263 2.1244 -2215.0
    ## - Exposure:Morpholino    1  0.037837 2.1360 -2210.7
    ## - Generation:Morpholino  2  0.078868 2.1770 -2204.7


    # best model is Generation + Exposure + Morpholino + Sex + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino
    # GenerationF1, GenerationF2, GenerationF1:ExposureBaP, GenerationF2:MorpholinoAHR2Mo, ExposureBaP:MorpholinoAHR2Mo significant
    # extract lm that is best model for results

    lm_shan <-summary(AIC_wholedata) 
    lm_shan
    ## 
    ## Call:
    ## lm(formula = pw_shannon ~ Generation + Exposure + Morpholino + 
    ##     Sex + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino, 
    ##     data = alphadiv)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.284508 -0.039117  0.008548  0.046281  0.240522 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    4.517379   0.011329 398.755  < 2e-16 ***
    ## GenerationF1                  -0.048434   0.014773  -3.279  0.00113 ** 
    ## GenerationF2                  -0.062332   0.015008  -4.153 3.99e-05 ***
    ## ExposureBaP                    0.026185   0.013797   1.898  0.05842 .  
    ## MorpholinoAHR2Mo              -0.013178   0.013797  -0.955  0.34010    
    ## SexM                           0.010367   0.006964   1.489  0.13734    
    ## GenerationF1:ExposureBaP       0.037897   0.016997   2.230  0.02632 *  
    ## GenerationF2:ExposureBaP       0.024775   0.017006   1.457  0.14592    
    ## GenerationF1:MorpholinoAHR2Mo  0.012699   0.016997   0.747  0.45543    
    ## GenerationF2:MorpholinoAHR2Mo  0.063407   0.017006   3.729  0.00022 ***
    ## ExposureBaP:MorpholinoAHR2Mo  -0.037911   0.013925  -2.722  0.00676 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07145 on 411 degrees of freedom
    ## Multiple R-squared:  0.1165, Adjusted R-squared:  0.095 
    ## F-statistic: 5.419 on 10 and 411 DF,  p-value: 1.477e-07

    # richness
    pw_richplot <- ggplot(alphadiv, aes(x = Generation, y = pw_richness, 
                                        fill = Treatments)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Pathway Richness")
    pw_richplot

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/pathway richness plot-1.png" width="98%" height="98%" />

    # shannon
    alphadiv$Treatments <- factor(alphadiv$Treatments, 
                                    levels = c("CoMo_DMSO",
                                               "CoMo_BaP", 
                                               "AHR2Mo_DMSO",
                                               "AHR2Mo_BaP"))
    pw_shanplotint <- ggplot(alphadiv, aes(x = Treatments, y = pw_shannon,
                                        fill = Generation)) + 
      geom_boxplot() +
      theme_classic() +
      labs(x = NULL, y = "Pathway Shannon Diversity") +
      scale_fill_brewer(palette = "PuBu") +
      scale_x_discrete(labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 25))

    pw_shanplotint

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/intergenerational shannon graph for manuscript-1.png" width="98%" height="98%" />

The following sections use a subset of the data that is only for each
generation of animals. It uses a similar modeling approach as above but
without Generation as a variable since one generation at a time is
evaluated.

# F0

    ## calculate alpha diversity metrics and put into a data frame
    pwF0_richness <- rowSums(pathcov_tab_F0)
    pwF0_shannon <- diversity(pathabund_tab_F0, "shannon", base = exp(1))
    pwF0_simpson <- diversity(pathabund_tab_F0, "simpson")
    meta_modF0 <- metaF0 %>%
      dplyr::select(c("SampleID", "Exposure", "Morpholino", "Sex", "BMI", "Treatments")) %>% # filter to only important variables we are testing
      mutate(Treatments = factor(Treatments))
    alphadiv_pathF0 <- data.frame(pwF0_richness, pwF0_shannon, pwF0_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(meta_modF0, by = "SampleID")
    alphadiv_pathF0
    ##     SampleID pwF0_richness pwF0_shannon pwF0_simpson Exposure Morpholino Sex
    ## 1   BaP_0003            23     4.504218    0.9881009     DMSO     AHR2Mo   F
    ## 2   BaP_0004            42     4.545643    0.9886415     DMSO     AHR2Mo   F
    ## 3   BaP_0005            26     4.520909    0.9881479     DMSO     AHR2Mo   F
    ## 4   BaP_0006            25     4.353807    0.9844620     DMSO     AHR2Mo   M
    ## 5   BaP_0007            33     4.439526    0.9869343     DMSO     AHR2Mo   F
    ## 6   BaP_0008            31     4.522775    0.9880482     DMSO     AHR2Mo   M
    ## 7   BaP_0037            26     4.496289    0.9878816     DMSO     AHR2Mo   F
    ## 8   BaP_0038            33     4.316099    0.9853482     DMSO     AHR2Mo   M
    ## 9   BaP_0039            28     4.504699    0.9879460     DMSO     AHR2Mo   F
    ## 10  BaP_0040            34     4.425233    0.9873059     DMSO     AHR2Mo   F
    ## 11  BaP_0041            26     4.520501    0.9883291     DMSO     AHR2Mo   F
    ## 12  BaP_0042            27     4.534710    0.9883182     DMSO     AHR2Mo   M
    ## 13  BaP_0043            41     4.501306    0.9879252     DMSO     AHR2Mo   F
    ## 14  BaP_0044            26     4.525067    0.9882462     DMSO     AHR2Mo   M
    ## 15  BaP_0073            39     4.328620    0.9855610     DMSO     AHR2Mo   F
    ## 16  BaP_0074            39     4.230060    0.9841720     DMSO     AHR2Mo   M
    ## 17  BaP_0075            38     4.373401    0.9860194      BaP     AHR2Mo   F
    ## 18  BaP_0076            36     4.341884    0.9855663      BaP     AHR2Mo   M
    ## 19  BaP_0077            34     4.329710    0.9855582      BaP     AHR2Mo   F
    ## 20  BaP_0078            34     4.438752    0.9869682      BaP     AHR2Mo   M
    ## 21  BaP_0079            31     4.378341    0.9861319      BaP     AHR2Mo   F
    ## 22  BaP_0080            33     4.287327    0.9850763      BaP     AHR2Mo   M
    ## 23  BaP_0109            26     4.541155    0.9884644      BaP     AHR2Mo   F
    ## 24  BaP_0110            27     4.566991    0.9888320      BaP     AHR2Mo   M
    ## 25  BaP_0111            34     4.573263    0.9886815      BaP     AHR2Mo   F
    ## 26  BaP_0112            25     4.518203    0.9883509      BaP     AHR2Mo   M
    ## 27  BaP_0113            27     4.518732    0.9881668      BaP     AHR2Mo   F
    ## 28  BaP_0114            34     4.555842    0.9886733      BaP     AHR2Mo   M
    ## 29  BaP_0115            24     4.540651    0.9886155      BaP     AHR2Mo   F
    ## 30  BaP_0116            34     4.554625    0.9885314      BaP     AHR2Mo   M
    ## 31  BaP_0433            31     4.521076    0.9884088      BaP     AHR2Mo   F
    ## 32  BaP_0434            24     4.496731    0.9879677      BaP     AHR2Mo   M
    ## 33  BaP_0435            42     4.508273    0.9881797      BaP     AHR2Mo   F
    ## 34  BaP_0436            38     4.510948    0.9876939      BaP     AHR2Mo   M
    ## 35  BaP_0437            37     4.515367    0.9881783      BaP     AHR2Mo   F
    ## 36  BaP_0438            42     4.559864    0.9886423      BaP     AHR2Mo   M
    ## 37  BaP_0439            39     4.571648    0.9886480      BaP     AHR2Mo   F
    ## 38  BaP_0440            40     4.537667    0.9880528      BaP     AHR2Mo   M
    ## 39  BaP_0541            26     4.457410    0.9870694      BaP     AHR2Mo   F
    ## 40  BaP_0542            38     4.497898    0.9879342      BaP     AHR2Mo   M
    ## 41  BaP_0543            38     4.537138    0.9885087      BaP     AHR2Mo   F
    ## 42  BaP_0544            43     4.484337    0.9878241      BaP     AHR2Mo   M
    ## 43  BaP_0545            27     4.484008    0.9877804      BaP     AHR2Mo   M
    ## 44  BaP_0546            25     4.482312    0.9878118      BaP     AHR2Mo   F
    ## 45  BaP_0547            30     4.573616    0.9888375      BaP     AHR2Mo   F
    ## 46  BaP_0548            45     4.486966    0.9877766      BaP     AHR2Mo   M
    ## 47  BaP_0651            41     4.530872    0.9883613      BaP     AHR2Mo   F
    ## 48  BaP_0652            37     4.635700    0.9894872      BaP     AHR2Mo   M
    ## 49  BaP_0653            38     4.511315    0.9880871      BaP     AHR2Mo   M
    ## 50  BaP_0654            30     4.575460    0.9887534      BaP     AHR2Mo   M
    ## 51  BaP_0655            40     4.538502    0.9883690      BaP     AHR2Mo   M
    ## 52  BaP_0656            24     4.553642    0.9886642      BaP     AHR2Mo   M
    ## 53  BaP_0867            42     4.537241    0.9885059     DMSO       CoMo   F
    ## 54  BaP_0868            36     4.559420    0.9887496     DMSO       CoMo   M
    ## 55  BaP_0869            40     4.563146    0.9884425     DMSO       CoMo   F
    ## 56  BaP_0870            25     4.532752    0.9885157     DMSO       CoMo   M
    ## 57  BaP_0871            25     4.585741    0.9888897     DMSO       CoMo   F
    ## 58  BaP_0872            45     4.529618    0.9884092     DMSO       CoMo   M
    ## 59  BaP_0873            23     4.536499    0.9885496     DMSO       CoMo   F
    ## 60  BaP_0874            29     4.562748    0.9887055     DMSO       CoMo   M
    ## 61  BaP_0875            49     4.528531    0.9882597     DMSO       CoMo   F
    ## 62  BaP_0876            25     4.533867    0.9885804     DMSO       CoMo   M
    ## 63  BaP_0877            33     4.573732    0.9887946     DMSO       CoMo   F
    ## 64  BaP_0878            32     4.559254    0.9887499     DMSO       CoMo   M
    ## 65  BaP_0879            38     4.540871    0.9885075      BaP       CoMo   F
    ## 66  BaP_0880            26     4.613909    0.9889981      BaP       CoMo   F
    ## 67  BaP_0881            28     4.513351    0.9882516      BaP       CoMo   F
    ## 68  BaP_0882            27     4.570542    0.9887851      BaP       CoMo   M
    ## 69  BaP_0883            26     4.587557    0.9887919      BaP       CoMo   F
    ## 70  BaP_0884            36     4.517204    0.9881091      BaP       CoMo   M
    ## 71  BaP_0885            41     4.545778    0.9882872      BaP       CoMo   F
    ## 72  BaP_0886            25     4.550773    0.9886811      BaP       CoMo   M
    ## 73  BaP_0887            26     4.565892    0.9886970      BaP       CoMo   F
    ## 74  BaP_0888            26     4.536289    0.9885451      BaP       CoMo   M
    ## 75  BaP_0937            24     4.522185    0.9882154     DMSO       CoMo   F
    ## 76  BaP_0938            46     4.508775    0.9881818     DMSO       CoMo   M
    ## 77  BaP_0939            28     4.534658    0.9885164     DMSO       CoMo   F
    ## 78  BaP_0940            45     4.584213    0.9885080     DMSO       CoMo   M
    ## 79  BaP_0941            48     4.403116    0.9866326     DMSO       CoMo   F
    ## 80  BaP_0942            31     4.594294    0.9889972     DMSO       CoMo   M
    ## 81  BaP_0943            23     4.383083    0.9859967     DMSO       CoMo   F
    ## 82  BaP_0944            46     4.579474    0.9883979     DMSO       CoMo   M
    ## 83  BaP_0945            24     4.493030    0.9879130      BaP       CoMo   F
    ## 84  BaP_0946            42     4.524984    0.9879366      BaP       CoMo   M
    ## 85  BaP_0947            24     4.520566    0.9882318      BaP       CoMo   F
    ## 86  BaP_0948            26     4.509400    0.9879252      BaP       CoMo   M
    ## 87  BaP_0949            36     4.510599    0.9881336      BaP       CoMo   F
    ## 88  BaP_0950            26     4.548811    0.9885605      BaP       CoMo   M
    ## 89  BaP_0951            37     4.536197    0.9885085      BaP       CoMo   F
    ## 90  BaP_0952            43     4.514813    0.9875756      BaP       CoMo   M
    ## 91  BaP_0953            37     4.551988    0.9883968      BaP       CoMo   F
    ## 92  BaP_0954            29     4.511145    0.9881735      BaP       CoMo   M
    ## 93  BaP_0955            25     4.573314    0.9887514      BaP       CoMo   F
    ## 94  BaP_0956            35     4.494903    0.9877437      BaP       CoMo   M
    ## 95  BaP_0957            28     4.498406    0.9878888     DMSO     AHR2Mo   F
    ## 96  BaP_0958            26     4.412859    0.9866246     DMSO     AHR2Mo   M
    ## 97  BaP_0959            28     4.573228    0.9886918     DMSO     AHR2Mo   F
    ## 98  BaP_0960            45     4.542668    0.9884344     DMSO     AHR2Mo   M
    ## 99  BaP_1009            29     4.464234    0.9874148     DMSO       CoMo   F
    ## 100 BaP_1010            31     4.470033    0.9873883     DMSO       CoMo   M
    ## 101 BaP_1011            27     4.452881    0.9873293     DMSO       CoMo   M
    ## 102 BaP_1012            45     4.274714    0.9853417     DMSO       CoMo   F
    ## 103 BaP_1013            25     4.446758    0.9871227     DMSO       CoMo   F
    ## 104 BaP_1014            42     4.487749    0.9876448     DMSO       CoMo   M
    ## 105 BaP_1015            25     4.514678    0.9881170     DMSO       CoMo   F
    ## 106 BaP_1016            42     4.612949    0.9892287     DMSO       CoMo   M
    ## 107 BaP_1019            34     4.539563    0.9886260      BaP       CoMo   F
    ## 108 BaP_1020            28     4.514816    0.9881500      BaP       CoMo   M
    ## 109 BaP_1021            25     4.581829    0.9888790      BaP       CoMo   F
    ## 110 BaP_1022            43     4.600920    0.9891380      BaP       CoMo   M
    ## 111 BaP_1023            40     4.442561    0.9872198      BaP       CoMo   F
    ## 112 BaP_1024            41     4.545772    0.9885967      BaP       CoMo   M
    ## 113 BaP_1025            40     4.521984    0.9880028     DMSO     AHR2Mo   F
    ## 114 BaP_1026            24     4.484919    0.9878194     DMSO     AHR2Mo   M
    ## 115 BaP_1027            40     4.504881    0.9880215     DMSO     AHR2Mo   F
    ## 116 BaP_1028            39     4.448848    0.9876286     DMSO     AHR2Mo   M
    ## 117 BaP_1029            39     4.572114    0.9888226     DMSO     AHR2Mo   F
    ## 118 BaP_1030            37     4.555932    0.9885582     DMSO     AHR2Mo   M
    ## 119 BaP_1031            31     4.559720    0.9887723     DMSO     AHR2Mo   F
    ## 120 BaP_1032            29     4.542443    0.9883773     DMSO     AHR2Mo   M
    ## 121 BaP_1081            30     4.579584    0.9889271     DMSO       CoMo   F
    ## 122 BaP_1082            35     4.614369    0.9892223     DMSO       CoMo   M
    ## 123 BaP_1083            40     4.449053    0.9867586     DMSO       CoMo   F
    ## 124 BaP_1084            45     4.629706    0.9891099     DMSO       CoMo   M
    ## 125 BaP_1085            24     4.448315    0.9872674     DMSO       CoMo   F
    ## 126 BaP_1086            36     4.624406    0.9891998     DMSO       CoMo   M
    ## 127 BaP_1087            28     4.563440    0.9886356     DMSO       CoMo   F
    ## 128 BaP_1088            39     4.609836    0.9890670     DMSO       CoMo   F
    ## 129 BaP_1089            25     4.598048    0.9891037      BaP       CoMo   F
    ## 130 BaP_1090            35     4.641423    0.9891380      BaP       CoMo   M
    ## 131 BaP_1091            24     4.467356    0.9875381      BaP       CoMo   F
    ## 132 BaP_1092            40     4.630763    0.9893792      BaP       CoMo   M
    ## 133 BaP_1093            25     4.579059    0.9887814      BaP       CoMo   F
    ## 134 BaP_1094            30     4.510950    0.9881143      BaP       CoMo   M
    ## 135 BaP_1095            28     4.531201    0.9883345      BaP       CoMo   F
    ## 136 BaP_1096            37     4.575171    0.9887221      BaP       CoMo   M
    ## 137 BaP_1097            40     4.655449    0.9895118     DMSO     AHR2Mo   F
    ## 138 BaP_1098            34     4.568771    0.9886379     DMSO     AHR2Mo   M
    ## 139 BaP_1099            45     4.616147    0.9891246     DMSO     AHR2Mo   F
    ## 140 BaP_1100            36     4.577251    0.9886513     DMSO     AHR2Mo   M
    ## 141 BaP_1101            41     4.599216    0.9888807     DMSO     AHR2Mo   F
    ## 142 BaP_1102            34     4.594431    0.9886945     DMSO     AHR2Mo   M
    ## 143 BaP_1103            29     4.650466    0.9894567     DMSO     AHR2Mo   M
    ## 144 BaP_1104            30     4.425107    0.9865785     DMSO     AHR2Mo   M
    ##           BMI  Treatments
    ## 1   1094.1828 AHR2Mo_DMSO
    ## 2    620.8912 AHR2Mo_DMSO
    ## 3   1069.3878 AHR2Mo_DMSO
    ## 4    884.4953 AHR2Mo_DMSO
    ## 5    897.9592 AHR2Mo_DMSO
    ## 6    717.9931 AHR2Mo_DMSO
    ## 7    770.0312 AHR2Mo_DMSO
    ## 8    604.0816 AHR2Mo_DMSO
    ## 9    977.5087 AHR2Mo_DMSO
    ## 10   674.7405 AHR2Mo_DMSO
    ## 11   993.7500 AHR2Mo_DMSO
    ## 12   579.5848 AHR2Mo_DMSO
    ## 13   889.7959 AHR2Mo_DMSO
    ## 14   694.4444 AHR2Mo_DMSO
    ## 15   683.3910 AHR2Mo_DMSO
    ## 16   558.8585 AHR2Mo_DMSO
    ## 17   942.2936  AHR2Mo_BaP
    ## 18   673.8281  AHR2Mo_BaP
    ## 19  1175.8585  AHR2Mo_BaP
    ## 20   651.9743  AHR2Mo_BaP
    ## 21   895.0617  AHR2Mo_BaP
    ## 22   718.0021  AHR2Mo_BaP
    ## 23   748.4568  AHR2Mo_BaP
    ## 24   622.8374  AHR2Mo_BaP
    ## 25   663.5802  AHR2Mo_BaP
    ## 26   551.5088  AHR2Mo_BaP
    ## 27   648.7889  AHR2Mo_BaP
    ## 28   693.8776  AHR2Mo_BaP
    ## 29   688.7052  AHR2Mo_BaP
    ## 30   605.4688  AHR2Mo_BaP
    ## 31   717.9931  AHR2Mo_BaP
    ## 32   634.7656  AHR2Mo_BaP
    ## 33   770.0312  AHR2Mo_BaP
    ## 34   596.8858  AHR2Mo_BaP
    ## 35   832.3424  AHR2Mo_BaP
    ## 36   685.7143  AHR2Mo_BaP
    ## 37   778.5467  AHR2Mo_BaP
    ## 38   709.3426  AHR2Mo_BaP
    ## 39   749.1082  AHR2Mo_BaP
    ## 40   457.8564  AHR2Mo_BaP
    ## 41   848.1262  AHR2Mo_BaP
    ## 42   624.3496  AHR2Mo_BaP
    ## 43   596.8858  AHR2Mo_BaP
    ## 44   685.7143  AHR2Mo_BaP
    ## 45   749.2196  AHR2Mo_BaP
    ## 46   759.1837  AHR2Mo_BaP
    ## 47   818.1154  AHR2Mo_BaP
    ## 48   587.6951  AHR2Mo_BaP
    ## 49   722.2222  AHR2Mo_BaP
    ## 50   596.8779  AHR2Mo_BaP
    ## 51   795.8478  AHR2Mo_BaP
    ## 52   645.1613  AHR2Mo_BaP
    ## 53   993.9446   CoMo_DMSO
    ## 54   451.1719   CoMo_DMSO
    ## 55  1296.0761   CoMo_DMSO
    ## 56   674.0129   CoMo_DMSO
    ## 57  1029.9489   CoMo_DMSO
    ## 58   404.0404   CoMo_DMSO
    ## 59   824.4898   CoMo_DMSO
    ## 60   507.8125   CoMo_DMSO
    ## 61   787.1972   CoMo_DMSO
    ## 62   546.8750   CoMo_DMSO
    ## 63   734.6189   CoMo_DMSO
    ## 64   532.4074   CoMo_DMSO
    ## 65   856.4815    CoMo_BaP
    ## 66   726.5306    CoMo_BaP
    ## 67   759.1837    CoMo_BaP
    ## 68   605.5363    CoMo_BaP
    ## 69   709.8765    CoMo_BaP
    ## 70  1020.7612    CoMo_BaP
    ## 71  1092.7456    CoMo_BaP
    ## 72   849.6094    CoMo_BaP
    ## 73   676.3788    CoMo_BaP
    ## 74   562.4543    CoMo_BaP
    ## 75   697.1904   CoMo_DMSO
    ## 76   569.3297   CoMo_DMSO
    ## 77   830.4498   CoMo_DMSO
    ## 78   606.4209   CoMo_DMSO
    ## 79   743.9446   CoMo_DMSO
    ## 80   679.5225   CoMo_DMSO
    ## 81   712.8906   CoMo_DMSO
    ## 82   432.6531   CoMo_DMSO
    ## 83   790.8429    CoMo_BaP
    ## 84   576.1719    CoMo_BaP
    ## 85   601.8519    CoMo_BaP
    ## 86   894.9011    CoMo_BaP
    ## 87   878.9063    CoMo_BaP
    ## 88   657.4394    CoMo_BaP
    ## 89   902.7778    CoMo_BaP
    ## 90   576.1719    CoMo_BaP
    ## 91  1009.3652    CoMo_BaP
    ## 92   556.6406    CoMo_BaP
    ## 93   957.3361    CoMo_BaP
    ## 94   830.0781    CoMo_BaP
    ## 95   524.6914 AHR2Mo_DMSO
    ## 96   460.1899 AHR2Mo_DMSO
    ## 97   973.3701 AHR2Mo_DMSO
    ## 98   437.0447 AHR2Mo_DMSO
    ## 99   622.8374   CoMo_DMSO
    ## 100  676.3788   CoMo_DMSO
    ## 101  615.2344   CoMo_DMSO
    ## 102  787.1972   CoMo_DMSO
    ## 103  900.2770   CoMo_DMSO
    ## 104  663.5802   CoMo_DMSO
    ## 105  910.4938   CoMo_DMSO
    ## 106  677.1861   CoMo_DMSO
    ## 107  725.3086    CoMo_BaP
    ## 108  856.4014    CoMo_BaP
    ## 109  775.5102    CoMo_BaP
    ## 110  769.8962    CoMo_BaP
    ## 111  620.4082    CoMo_BaP
    ## 112  650.1096    CoMo_BaP
    ## 113  942.9066 AHR2Mo_DMSO
    ## 114  613.9438 AHR2Mo_DMSO
    ## 115  640.4321 AHR2Mo_DMSO
    ## 116  700.6920 AHR2Mo_DMSO
    ## 117  783.6735 AHR2Mo_DMSO
    ## 118  582.6397 AHR2Mo_DMSO
    ## 119  824.0997 AHR2Mo_DMSO
    ## 120  458.4775 AHR2Mo_DMSO
    ## 121  832.6531   CoMo_DMSO
    ## 122  780.4370   CoMo_DMSO
    ## 123  759.6786   CoMo_DMSO
    ## 124  732.4219   CoMo_DMSO
    ## 125  634.7656   CoMo_DMSO
    ## 126  692.0415   CoMo_DMSO
    ## 127  918.2099   CoMo_DMSO
    ## 128  653.0612   CoMo_DMSO
    ## 129  808.1633    CoMo_BaP
    ## 130  642.0927    CoMo_BaP
    ## 131  712.8906    CoMo_BaP
    ## 132  546.9388    CoMo_BaP
    ## 133  763.8889    CoMo_BaP
    ## 134  654.2969    CoMo_BaP
    ## 135  854.6384    CoMo_BaP
    ## 136  560.1469    CoMo_BaP
    ## 137  748.4568 AHR2Mo_DMSO
    ## 138  477.5023 AHR2Mo_DMSO
    ## 139  761.2457 AHR2Mo_DMSO
    ## 140  635.5004 AHR2Mo_DMSO
    ## 141  840.8163 AHR2Mo_DMSO
    ## 142  622.8374 AHR2Mo_DMSO
    ## 143  613.9438 AHR2Mo_DMSO
    ## 144  661.2245 AHR2Mo_DMSO

What covariates predict F1 gut metagenome richness?

    # step AIC to test which model
    testmod_richF0 <- lm(pwF0_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino + Sex:BMI, 
                          data = alphadiv_pathF0)
    AIC_F0 <- stepAIC(testmod_richF0)
    ## Start:  AIC=568.35
    ## pwF0_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Exposure:Sex         1    43.295 6622.5 567.29
    ## - Sex:BMI              1    51.814 6631.0 567.48
    ## - Exposure:Morpholino  1    83.079 6662.3 568.15
    ## <none>                             6579.2 568.35
    ## - Morpholino:Sex       1   143.175 6722.4 569.45
    ## 
    ## Step:  AIC=567.29
    ## pwF0_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Morpholino:Sex + Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Sex:BMI              1    41.206 6663.7 566.18
    ## - Exposure:Morpholino  1    90.586 6713.1 567.25
    ## <none>                             6622.5 567.29
    ## - Morpholino:Sex       1   141.251 6763.7 568.33
    ## 
    ## Step:  AIC=566.18
    ## pwF0_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - BMI                  1     0.003 6663.7 564.18
    ## <none>                             6663.7 566.18
    ## - Exposure:Morpholino  1    97.314 6761.0 566.27
    ## - Morpholino:Sex       1   132.813 6796.5 567.03
    ## 
    ## Step:  AIC=564.18
    ## pwF0_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## <none>                             6663.7 564.18
    ## - Exposure:Morpholino  1    97.504 6761.2 564.28
    ## - Morpholino:Sex       1   133.066 6796.8 565.03


    # best model is Exposure + Morpholino + Sex + Exposure:Morpholino + Morpholino:Sex and Sex significant

    # save object
    lm_richF0 <- summary(AIC_F0)
    lm_richF0
    ## 
    ## Call:
    ## lm(formula = pwF0_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex, data = alphadiv_pathF0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -11.3007  -5.9006  -0.3007   6.3505  16.3216 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   32.6784     1.3933  23.453   <2e-16 ***
    ## ExposureBaP                   -2.7778     1.6379  -1.696   0.0921 .  
    ## MorpholinoAHR2Mo               0.5985     1.9712   0.304   0.7619    
    ## SexM                           3.6223     1.6404   2.208   0.0289 *  
    ## ExposureBaP:MorpholinoAHR2Mo   3.2972     2.3204   1.421   0.1576    
    ## MorpholinoAHR2Mo:SexM         -3.8556     2.3226  -1.660   0.0992 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.949 on 138 degrees of freedom
    ## Multiple R-squared:  0.05476,    Adjusted R-squared:  0.02051 
    ## F-statistic: 1.599 on 5 and 138 DF,  p-value: 0.1644


    # test just sex + treatment + interaction
    testmod_rich2F0 <- lm(pwF0_richness ~ Treatments + Sex + BMI + Sex:Treatments + Sex:BMI, 
                          data = alphadiv_pathF0)
    summary(stepAIC(testmod_rich2F0))
    ## Start:  AIC=569.64
    ## pwF0_richness ~ Treatments + Sex + BMI + Sex:Treatments + Sex:BMI
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - Treatments:Sex  3    216.93 6763.7 568.33
    ## - Sex:BMI         1     45.71 6592.5 568.64
    ## <none>                        6546.8 569.64
    ## 
    ## Step:  AIC=568.33
    ## pwF0_richness ~ Treatments + Sex + BMI + Sex:BMI
    ## 
    ##              Df Sum of Sq    RSS    AIC
    ## - Treatments  3   120.075 6883.8 564.86
    ## - Sex:BMI     1    32.767 6796.5 567.03
    ## <none>                    6763.7 568.33
    ## 
    ## Step:  AIC=564.86
    ## pwF0_richness ~ Sex + BMI + Sex:BMI
    ## 
    ##           Df Sum of Sq    RSS    AIC
    ## - Sex:BMI  1     58.11 6941.9 564.07
    ## <none>                 6883.8 564.86
    ## 
    ## Step:  AIC=564.07
    ## pwF0_richness ~ Sex + BMI
    ## 
    ##        Df Sum of Sq    RSS    AIC
    ## - BMI   1     0.316 6942.2 562.08
    ## - Sex   1    68.765 7010.7 563.49
    ## <none>              6941.9 564.07
    ## 
    ## Step:  AIC=562.08
    ## pwF0_richness ~ Sex
    ## 
    ##        Df Sum of Sq    RSS    AIC
    ## <none>              6942.2 562.08
    ## - Sex   1    107.52 7049.8 562.29
    ## 
    ## Call:
    ## lm(formula = pwF0_richness ~ Sex, data = alphadiv_pathF0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -10.0845  -6.3562  -0.0845   6.0976  16.6438 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  32.3562     0.8184  39.538   <2e-16 ***
    ## SexM          1.7283     1.1655   1.483     0.14    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.992 on 142 degrees of freedom
    ## Multiple R-squared:  0.01525,    Adjusted R-squared:  0.008316 
    ## F-statistic: 2.199 on 1 and 142 DF,  p-value: 0.1403


    # Sex best model but not significant

What covariates predict F0 gut metagenome shannon diversity?

    # step AIC to test which model
    testmod_shanF0 <- lm(pwF0_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + Sex:Morpholino + Exposure:Morpholino + Sex:BMI, 
                        data = alphadiv_pathF0)
    AIC_F0 <- stepAIC(testmod_shanF0)
    ## Start:  AIC=-740.68
    ## pwF0_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + 
    ##     Sex:Morpholino + Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Sex         1 0.0000875 0.74181 -742.66
    ## - Exposure:Morpholino  1 0.0040750 0.74580 -741.89
    ## - Sex:BMI              1 0.0060847 0.74781 -741.50
    ## <none>                             0.74172 -740.68
    ## - Morpholino:Sex       1 0.0169411 0.75866 -739.43
    ## 
    ## Step:  AIC=-742.66
    ## pwF0_shannon ~ Exposure + Morpholino + Sex + BMI + Morpholino:Sex + 
    ##     Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Morpholino  1 0.0040190 0.74583 -743.88
    ## - Sex:BMI              1 0.0059975 0.74781 -743.50
    ## <none>                             0.74181 -742.66
    ## - Morpholino:Sex       1 0.0169135 0.75872 -741.41
    ## 
    ## Step:  AIC=-743.88
    ## pwF0_shannon ~ Exposure + Morpholino + Sex + BMI + Morpholino:Sex + 
    ##     Sex:BMI
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## - Exposure        1 0.0043412 0.75017 -745.05
    ## - Sex:BMI         1 0.0055126 0.75134 -744.82
    ## <none>                        0.74583 -743.88
    ## - Morpholino:Sex  1 0.0175451 0.76337 -742.53
    ## 
    ## Step:  AIC=-745.05
    ## pwF0_shannon ~ Morpholino + Sex + BMI + Morpholino:Sex + Sex:BMI
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## - Sex:BMI         1 0.0038867 0.75406 -746.30
    ## <none>                        0.75017 -745.05
    ## - Morpholino:Sex  1 0.0166381 0.76681 -743.89
    ## 
    ## Step:  AIC=-746.3
    ## pwF0_shannon ~ Morpholino + Sex + BMI + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## - BMI             1 0.0018133 0.75587 -747.96
    ## <none>                        0.75406 -746.30
    ## - Morpholino:Sex  1 0.0158646 0.76992 -745.30
    ## 
    ## Step:  AIC=-747.96
    ## pwF0_shannon ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## <none>                        0.75587 -747.96
    ## - Morpholino:Sex  1  0.015482 0.77135 -747.04


    # best model is Morpholino + Sex + Morpholino:Sex
    # none significant but some close

    lm_shanF0 <- summary(AIC_F0)
    lm_shanF0
    ## 
    ## Call:
    ## lm(formula = pwF0_shannon ~ Morpholino + Sex + Morpholino:Sex, 
    ##     data = alphadiv_pathF0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.26786 -0.02970  0.01125  0.05047  0.15255 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            4.52135    0.01192 379.315   <2e-16 ***
    ## MorpholinoAHR2Mo      -0.01160    0.01721  -0.674   0.5014    
    ## SexM                   0.02969    0.01735   1.711   0.0892 .  
    ## MorpholinoAHR2Mo:SexM -0.04152    0.02452  -1.693   0.0926 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07348 on 140 degrees of freedom
    ## Multiple R-squared:  0.06726,    Adjusted R-squared:  0.04728 
    ## F-statistic: 3.365 on 3 and 140 DF,  p-value: 0.0205


    # test just sex + treatment + interaction
    testmod_shan2F0 <- lm(pwF0_shannon ~ Sex + Treatments + BMI + Sex:Treatments, 
                          data = alphadiv_pathF0)
    summary(stepAIC(testmod_shan2F0))
    ## Start:  AIC=-743.06
    ## pwF0_shannon ~ Sex + Treatments + BMI + Sex:Treatments
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## - BMI             1  0.001799 0.73134 -744.71
    ## <none>                        0.72954 -743.06
    ## - Sex:Treatments  3  0.034064 0.76360 -742.49
    ## 
    ## Step:  AIC=-744.71
    ## pwF0_shannon ~ Sex + Treatments + Sex:Treatments
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## <none>                        0.73134 -744.71
    ## - Sex:Treatments  3  0.034219 0.76556 -744.12
    ## 
    ## Call:
    ## lm(formula = pwF0_shannon ~ Sex + Treatments + Sex:Treatments, 
    ##     data = alphadiv_pathF0)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.257365 -0.031382  0.004406  0.047412  0.163042 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 4.500976   0.016823 267.543   <2e-16 ***
    ## SexM                        0.054113   0.024482   2.210   0.0288 *  
    ## TreatmentsCoMo_BaP          0.040743   0.023792   1.712   0.0891 .  
    ## TreatmentsAHR2Mo_DMSO       0.019450   0.023792   0.817   0.4151    
    ## TreatmentsAHR2Mo_BaP       -0.003916   0.024882  -0.157   0.8752    
    ## SexM:TreatmentsCoMo_BaP    -0.048851   0.034622  -1.411   0.1605    
    ## SexM:TreatmentsAHR2Mo_DMSO -0.087114   0.034622  -2.516   0.0130 *  
    ## SexM:TreatmentsAHR2Mo_BaP  -0.044340   0.034703  -1.278   0.2035    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07333 on 136 degrees of freedom
    ## Multiple R-squared:  0.09754,    Adjusted R-squared:  0.05109 
    ## F-statistic:   2.1 on 7 and 136 DF,  p-value: 0.04756


    # best model is Sex + Treatments + Sex:Treatments
    # none significant

    plot(lm(formula = pwF0_shannon ~ Sex + Treatments + Sex:Treatments, data = alphadiv_pathF0))

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F0-1.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F0-2.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F0-3.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F0-4.png" width="98%" height="98%" />

    # linear
    # normal for the most part, not as much in the lower quantile
    # homoskedastic
    # no outliers

    pwF0_richplot <- ggplot(alphadiv_pathF0, aes(x = Treatments, y = pwF0_richness)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Pathway Richness") +
      scale_y_continuous(limits = c(20, 55))

    pwF0_shanplot <- ggplot(alphadiv_pathF0, aes(x = Treatments, y = pwF0_shannon)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Shannon Diversity Metric") +
      scale_y_continuous(limits = c(4.2, 4.8))

    pwF0_alphaplot <- ggarrange(pwF0_richplot, pwF0_shanplot,
                            labels = c("A", "B"),
                            ncol = 1, nrow = 2)
    pwF0_alphaplot

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F0 alpha diversity plotting-1.png" width="98%" height="98%" />

# F1

    ## calculate alpha diversity metrics and put into a data frame
    pwF1_richness <- rowSums(pathcov_tab_F1)
    pwF1_shannon <- diversity(pathabund_tab_F1, "shannon", base = exp(1))
    pwF1_simpson <- diversity(pathabund_tab_F1, "simpson")
    meta_modF1 <- metaF1 %>%
      dplyr::select(c("SampleID", "Exposure", "Morpholino", "Sex", "BMI", "Treatments")) %>% # filter to only important variables we are testing
      mutate(Treatments = factor(Treatments))
    alphadiv_pathF1 <- data.frame(pwF1_richness, pwF1_shannon, pwF1_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(meta_modF1, by = "SampleID")
    alphadiv_pathF1
    ##     SampleID pwF1_richness pwF1_shannon pwF1_simpson Exposure Morpholino Sex
    ## 1   BaP_0009            49     4.445836    0.9875717     DMSO       CoMo   F
    ## 2   BaP_0010            50     4.493771    0.9875042     DMSO       CoMo   M
    ## 3   BaP_0011            34     4.438030    0.9866638     DMSO       CoMo   F
    ## 4   BaP_0012            34     4.547843    0.9885828     DMSO       CoMo   M
    ## 5   BaP_0013            24     4.425747    0.9868110     DMSO       CoMo   F
    ## 6   BaP_0014            46     4.445500    0.9868567     DMSO       CoMo   M
    ## 7   BaP_0015            21     4.497418    0.9879694     DMSO       CoMo   F
    ## 8   BaP_0016            34     4.573694    0.9883337     DMSO       CoMo   M
    ## 9   BaP_0017            35     4.490649    0.9877335      BaP       CoMo   F
    ## 10  BaP_0018            24     4.516835    0.9882476      BaP       CoMo   M
    ## 11  BaP_0019            24     4.332636    0.9860365     DMSO     AHR2Mo   F
    ## 12  BaP_0020            39     4.497835    0.9878455     DMSO     AHR2Mo   M
    ## 13  BaP_0021            30     4.494431    0.9877907     DMSO     AHR2Mo   F
    ## 14  BaP_0022            39     4.426520    0.9869192     DMSO     AHR2Mo   M
    ## 15  BaP_0023            40     4.436349    0.9869666     DMSO     AHR2Mo   F
    ## 16  BaP_0024            32     4.541916    0.9886768     DMSO     AHR2Mo   M
    ## 17  BaP_0045            24     4.490652    0.9878204     DMSO       CoMo   F
    ## 18  BaP_0046            27     4.514056    0.9880311     DMSO       CoMo   M
    ## 19  BaP_0047            26     4.539308    0.9883887     DMSO       CoMo   F
    ## 20  BaP_0048            25     4.530884    0.9882733     DMSO       CoMo   M
    ## 21  BaP_0049            31     4.483722    0.9872031     DMSO       CoMo   F
    ## 22  BaP_0050            24     4.462053    0.9870267     DMSO       CoMo   M
    ## 23  BaP_0051            29     4.522418    0.9880867     DMSO       CoMo   F
    ## 24  BaP_0052            41     4.452999    0.9874027     DMSO       CoMo   M
    ## 25  BaP_0053            30     4.548316    0.9885620     DMSO     AHR2Mo   F
    ## 26  BaP_0054            26     4.328906    0.9852158     DMSO     AHR2Mo   M
    ## 27  BaP_0055            32     4.428729    0.9868857     DMSO     AHR2Mo   F
    ## 28  BaP_0056            32     4.517757    0.9884123     DMSO     AHR2Mo   M
    ## 29  BaP_0057            24     4.418993    0.9867686     DMSO     AHR2Mo   F
    ## 30  BaP_0058            33     4.481396    0.9877116     DMSO     AHR2Mo   M
    ## 31  BaP_0059            25     4.502718    0.9880775     DMSO     AHR2Mo   F
    ## 32  BaP_0060            38     4.417238    0.9866347     DMSO     AHR2Mo   M
    ## 33  BaP_0081            39     4.294424    0.9852336     DMSO       CoMo   F
    ## 34  BaP_0083            38     4.294128    0.9851863     DMSO       CoMo   F
    ## 35  BaP_0084            29     4.215196    0.9838677     DMSO       CoMo   M
    ## 36  BaP_0085            38     4.263761    0.9845001     DMSO       CoMo   F
    ## 37  BaP_0086            35     4.508152    0.9878353     DMSO       CoMo   M
    ## 38  BaP_0087            43     4.344175    0.9853647     DMSO       CoMo   F
    ## 39  BaP_0088            47     4.256614    0.9843278     DMSO       CoMo   M
    ## 40  BaP_0091            40     4.529356    0.9884797     DMSO     AHR2Mo   F
    ## 41  BaP_0092            33     4.489661    0.9873577     DMSO     AHR2Mo   M
    ## 42  BaP_0093            40     4.397413    0.9865212     DMSO     AHR2Mo   F
    ## 43  BaP_0094            45     4.436673    0.9866033     DMSO     AHR2Mo   M
    ## 44  BaP_0095            32     4.401144    0.9865592     DMSO     AHR2Mo   F
    ## 45  BaP_0096            37     4.352675    0.9860732     DMSO     AHR2Mo   M
    ## 46  BaP_0099            32     4.442931    0.9871331      BaP     AHR2Mo   F
    ## 47  BaP_0100            37     4.472857    0.9872593      BaP     AHR2Mo   M
    ## 48  BaP_0101            38     4.437928    0.9870747      BaP     AHR2Mo   F
    ## 49  BaP_0102            32     4.367063    0.9861272      BaP     AHR2Mo   M
    ## 50  BaP_0103            38     4.338099    0.9863603      BaP     AHR2Mo   F
    ## 51  BaP_0104            37     4.386232    0.9862100      BaP     AHR2Mo   M
    ## 52  BaP_0117            28     4.545590    0.9884380     DMSO       CoMo   F
    ## 53  BaP_0118            36     4.525876    0.9881573     DMSO       CoMo   M
    ## 54  BaP_0119            32     4.548033    0.9886903     DMSO       CoMo   F
    ## 55  BaP_0121            25     4.540934    0.9885798      BaP       CoMo   F
    ## 56  BaP_0122            23     4.517661    0.9882370      BaP       CoMo   M
    ## 57  BaP_0123            24     4.532811    0.9883483      BaP       CoMo   F
    ## 58  BaP_0124            25     4.477051    0.9877192      BaP       CoMo   M
    ## 59  BaP_0125            28     4.561343    0.9887294     DMSO     AHR2Mo   F
    ## 60  BaP_0126            30     4.527927    0.9882307     DMSO     AHR2Mo   M
    ## 61  BaP_0127            39     4.434399    0.9870181     DMSO     AHR2Mo   F
    ## 62  BaP_0128            36     4.530216    0.9885140     DMSO     AHR2Mo   M
    ## 63  BaP_0129            29     4.484446    0.9878066     DMSO     AHR2Mo   F
    ## 64  BaP_0130            28     4.532226    0.9885346     DMSO     AHR2Mo   M
    ## 65  BaP_0131            28     4.495650    0.9880066     DMSO     AHR2Mo   F
    ## 66  BaP_0132            28     4.523465    0.9883102     DMSO     AHR2Mo   M
    ## 67  BaP_0133            41     4.470522    0.9874441      BaP     AHR2Mo   F
    ## 68  BaP_0134            36     4.524628    0.9883918      BaP     AHR2Mo   M
    ## 69  BaP_0135            38     4.482283    0.9877355      BaP     AHR2Mo   F
    ## 70  BaP_0136            38     4.513877    0.9882242      BaP     AHR2Mo   M
    ## 71  BaP_0443            38     4.494443    0.9877333      BaP       CoMo   F
    ## 72  BaP_0444            28     4.447137    0.9871850      BaP       CoMo   M
    ## 73  BaP_0445            24     4.517276    0.9880003      BaP       CoMo   F
    ## 74  BaP_0446            32     4.583506    0.9889707      BaP       CoMo   M
    ## 75  BaP_0447            25     4.561361    0.9887000      BaP       CoMo   F
    ## 76  BaP_0448            37     4.579595    0.9887938      BaP       CoMo   M
    ## 77  BaP_0450            41     4.484328    0.9877138     DMSO     AHR2Mo   M
    ## 78  BaP_0451            45     4.456636    0.9872566     DMSO     AHR2Mo   F
    ## 79  BaP_0452            36     4.528606    0.9882602     DMSO     AHR2Mo   M
    ## 80  BaP_0453            28     4.517043    0.9882430     DMSO     AHR2Mo   F
    ## 81  BaP_0454            30     4.575520    0.9889817     DMSO     AHR2Mo   F
    ## 82  BaP_0455            38     4.512575    0.9882346      BaP     AHR2Mo   F
    ## 83  BaP_0456            34     4.546822    0.9885168      BaP     AHR2Mo   M
    ## 84  BaP_0549            38     4.536344    0.9885234      BaP       CoMo   F
    ## 85  BaP_0550            36     4.588442    0.9890731      BaP       CoMo   M
    ## 86  BaP_0551            40     4.525920    0.9884031      BaP       CoMo   M
    ## 87  BaP_0552            44     4.448724    0.9871658      BaP       CoMo   M
    ## 88  BaP_0553            35     4.560406    0.9887393      BaP       CoMo   F
    ## 89  BaP_0554            36     4.536312    0.9883334      BaP       CoMo   M
    ## 90  BaP_0555            25     4.528842    0.9885057      BaP       CoMo   F
    ## 91  BaP_0556            25     4.525814    0.9884885      BaP       CoMo   M
    ## 92  BaP_0557            25     4.501041    0.9880131      BaP     AHR2Mo   F
    ## 93  BaP_0558            24     4.501790    0.9881099      BaP     AHR2Mo   M
    ## 94  BaP_0559            27     4.561566    0.9887271      BaP     AHR2Mo   F
    ## 95  BaP_0560            47     4.490269    0.9877926      BaP     AHR2Mo   M
    ## 96  BaP_0561            28     4.576161    0.9888381      BaP     AHR2Mo   F
    ## 97  BaP_0562            44     4.463849    0.9874716      BaP     AHR2Mo   M
    ## 98  BaP_0563            24     4.509535    0.9882125      BaP     AHR2Mo   F
    ## 99  BaP_0564            47     4.498697    0.9877141      BaP     AHR2Mo   M
    ## 100 BaP_0657            42     4.571462    0.9887573      BaP       CoMo   F
    ## 101 BaP_0658            38     4.683788    0.9896237      BaP       CoMo   M
    ## 102 BaP_0659            26     4.529908    0.9882084      BaP       CoMo   F
    ## 103 BaP_0660            26     4.553212    0.9886697      BaP       CoMo   M
    ## 104 BaP_0661            43     4.483877    0.9875114      BaP       CoMo   F
    ## 105 BaP_0662            47     4.514247    0.9881306      BaP       CoMo   M
    ## 106 BaP_0663            41     4.446053    0.9870659      BaP       CoMo   F
    ## 107 BaP_0664            39     4.603789    0.9890413      BaP       CoMo   M
    ## 108 BaP_0665            42     4.443380    0.9875153      BaP     AHR2Mo   F
    ## 109 BaP_0666            26     4.538411    0.9883892      BaP     AHR2Mo   M
    ## 110 BaP_0667            48     4.359747    0.9861924      BaP     AHR2Mo   F
    ## 111 BaP_0668            45     4.602605    0.9889341      BaP     AHR2Mo   M
    ## 112 BaP_0669            44     4.494491    0.9879774      BaP     AHR2Mo   F
    ## 113 BaP_0670            41     4.507746    0.9879028      BaP     AHR2Mo   M
    ## 114 BaP_0671            39     4.580277    0.9889097      BaP     AHR2Mo   F
    ## 115 BaP_0672            25     4.490591    0.9879534      BaP     AHR2Mo   M
    ## 116 BaP_0757            33     4.592464    0.9890881     DMSO       CoMo   F
    ## 117 BaP_0758            37     4.490714    0.9875079     DMSO       CoMo   M
    ## 118 BaP_0759            47     4.503360    0.9881061     DMSO       CoMo   F
    ## 119 BaP_0760            43     4.580983    0.9887079     DMSO       CoMo   M
    ## 120 BaP_0761            36     4.604834    0.9890739     DMSO       CoMo   F
    ## 121 BaP_0762            40     4.516106    0.9880403     DMSO       CoMo   M
    ## 122 BaP_0763            30     4.561636    0.9886822     DMSO       CoMo   F
    ## 123 BaP_0764            35     4.592782    0.9891507     DMSO       CoMo   M
    ## 124 BaP_0765            31     4.521452    0.9883179      BaP       CoMo   F
    ## 125 BaP_0766            43     4.585418    0.9887643      BaP       CoMo   M
    ## 126 BaP_0767            46     4.475837    0.9877168      BaP       CoMo   F
    ## 127 BaP_0768            30     4.593279    0.9891742      BaP       CoMo   M
    ## 128 BaP_0769            25     4.544934    0.9884985      BaP       CoMo   F
    ## 129 BaP_0770            50     4.487092    0.9877380      BaP       CoMo   M
    ## 130 BaP_0771            24     4.504819    0.9878869      BaP       CoMo   F
    ## 131 BaP_0772            40     4.783915    0.9906523      BaP       CoMo   M
    ## 132 BaP_0773            29     4.547295    0.9883929      BaP     AHR2Mo   F
    ## 133 BaP_0774            42     4.601739    0.9887933      BaP     AHR2Mo   M
    ## 134 BaP_0775            35     4.531207    0.9883433      BaP     AHR2Mo   F
    ## 135 BaP_0776            40     4.497533    0.9877301      BaP     AHR2Mo   M
    ## 136 BaP_0777            42     4.554920    0.9885032      BaP     AHR2Mo   F
    ## 137 BaP_0778            31     4.561511    0.9888869      BaP     AHR2Mo   M
    ## 138 BaP_0779            38     4.581208    0.9886861      BaP     AHR2Mo   F
    ## 139 BaP_0780            35     4.494911    0.9879318      BaP     AHR2Mo   M
    ##          BMI  Treatments
    ## 1   671.4410   CoMo_DMSO
    ## 2   372.3974   CoMo_DMSO
    ## 3   490.8642   CoMo_DMSO
    ## 4   534.8007   CoMo_DMSO
    ## 5   498.6479   CoMo_DMSO
    ## 6   378.1163   CoMo_DMSO
    ## 7   636.5432   CoMo_DMSO
    ## 8   471.6049   CoMo_DMSO
    ## 9   725.2569    CoMo_BaP
    ## 10  489.3750    CoMo_BaP
    ## 11  534.0265 AHR2Mo_DMSO
    ## 12  520.4082 AHR2Mo_DMSO
    ## 13  650.2268 AHR2Mo_DMSO
    ## 14  507.9221 AHR2Mo_DMSO
    ## 15  695.9600 AHR2Mo_DMSO
    ## 16  477.5000 AHR2Mo_DMSO
    ## 17  645.5577   CoMo_DMSO
    ## 18  473.7696   CoMo_DMSO
    ## 19  578.5124   CoMo_DMSO
    ## 20  529.8765   CoMo_DMSO
    ## 21  510.3306   CoMo_DMSO
    ## 22  489.5895   CoMo_DMSO
    ## 23  660.7407   CoMo_DMSO
    ## 24  463.7500   CoMo_DMSO
    ## 25  571.2992 AHR2Mo_DMSO
    ## 26  531.8519 AHR2Mo_DMSO
    ## 27  642.8000 AHR2Mo_DMSO
    ## 28  394.4773 AHR2Mo_DMSO
    ## 29  792.9804 AHR2Mo_DMSO
    ## 30  473.2288 AHR2Mo_DMSO
    ## 31  615.7845 AHR2Mo_DMSO
    ## 32  468.2540 AHR2Mo_DMSO
    ## 33  631.6928   CoMo_DMSO
    ## 34  655.0475   CoMo_DMSO
    ## 35  468.2540   CoMo_DMSO
    ## 36  536.8000   CoMo_DMSO
    ## 37  560.6576   CoMo_DMSO
    ## 38  697.1336   CoMo_DMSO
    ## 39  500.2469   CoMo_DMSO
    ## 40  583.7037 AHR2Mo_DMSO
    ## 41  461.6300 AHR2Mo_DMSO
    ## 42  617.7686 AHR2Mo_DMSO
    ## 43  524.7934 AHR2Mo_DMSO
    ## 44  691.8715 AHR2Mo_DMSO
    ## 45  449.5465 AHR2Mo_DMSO
    ## 46  721.7882  AHR2Mo_BaP
    ## 47  482.4509  AHR2Mo_BaP
    ## 48  620.1901  AHR2Mo_BaP
    ## 49  390.0000  AHR2Mo_BaP
    ## 50  551.1364  AHR2Mo_BaP
    ## 51  420.6250  AHR2Mo_BaP
    ## 52  601.1342   CoMo_DMSO
    ## 53  513.7500   CoMo_DMSO
    ## 54  613.8453   CoMo_DMSO
    ## 55  680.0000    CoMo_BaP
    ## 56  467.6871    CoMo_BaP
    ## 57  699.4329    CoMo_BaP
    ## 58  300.0000    CoMo_BaP
    ## 59  643.0503 AHR2Mo_DMSO
    ## 60  451.6765 AHR2Mo_DMSO
    ## 61  555.0023 AHR2Mo_DMSO
    ## 62  485.0000 AHR2Mo_DMSO
    ## 63  719.0083 AHR2Mo_DMSO
    ## 64  522.9854 AHR2Mo_DMSO
    ## 65  574.3802 AHR2Mo_DMSO
    ## 66  495.3086 AHR2Mo_DMSO
    ## 67  682.3347  AHR2Mo_BaP
    ## 68  542.5342  AHR2Mo_BaP
    ## 69  599.0123  AHR2Mo_BaP
    ## 70  456.2760  AHR2Mo_BaP
    ## 71  557.6560    CoMo_BaP
    ## 72  468.1737    CoMo_BaP
    ## 73  587.1605    CoMo_BaP
    ## 74  514.3750    CoMo_BaP
    ## 75  947.8306    CoMo_BaP
    ## 76  546.6984    CoMo_BaP
    ## 77  474.7174 AHR2Mo_DMSO
    ## 78  734.1270 AHR2Mo_DMSO
    ## 79  510.0054 AHR2Mo_DMSO
    ## 80  501.8929 AHR2Mo_DMSO
    ## 81  668.6420 AHR2Mo_DMSO
    ## 82  602.9630  AHR2Mo_BaP
    ## 83  423.6669  AHR2Mo_BaP
    ## 84  686.9835    CoMo_BaP
    ## 85  450.1385    CoMo_BaP
    ## 86  498.9669    CoMo_BaP
    ## 87  457.4830    CoMo_BaP
    ## 88  726.3705    CoMo_BaP
    ## 89  416.8750    CoMo_BaP
    ## 90  447.3140    CoMo_BaP
    ## 91  352.3997    CoMo_BaP
    ## 92  741.0208  AHR2Mo_BaP
    ## 93  271.1634  AHR2Mo_BaP
    ## 94  540.0000  AHR2Mo_BaP
    ## 95  468.1250  AHR2Mo_BaP
    ## 96  504.3750  AHR2Mo_BaP
    ## 97  507.5000  AHR2Mo_BaP
    ## 98  649.3827  AHR2Mo_BaP
    ## 99  364.9584  AHR2Mo_BaP
    ## 100 587.2934    CoMo_BaP
    ## 101 480.3719    CoMo_BaP
    ## 102 591.2698    CoMo_BaP
    ## 103 444.3783    CoMo_BaP
    ## 104 546.2412    CoMo_BaP
    ## 105 524.7934    CoMo_BaP
    ## 106 658.7902    CoMo_BaP
    ## 107 487.8049    CoMo_BaP
    ## 108 632.3617  AHR2Mo_BaP
    ## 109 430.0000  AHR2Mo_BaP
    ## 110 646.1777  AHR2Mo_BaP
    ## 111 437.0748  AHR2Mo_BaP
    ## 112 580.4989  AHR2Mo_BaP
    ## 113 465.9091  AHR2Mo_BaP
    ## 114 535.7143  AHR2Mo_BaP
    ## 115 427.9778  AHR2Mo_BaP
    ## 116 544.7846   CoMo_DMSO
    ## 117 527.9012   CoMo_DMSO
    ## 118 453.1445   CoMo_DMSO
    ## 119 261.2847   CoMo_DMSO
    ## 120 578.7654   CoMo_DMSO
    ## 121 417.7778   CoMo_DMSO
    ## 122 559.5238   CoMo_DMSO
    ## 123 461.4512   CoMo_DMSO
    ## 124 563.5479    CoMo_BaP
    ## 125 578.7982    CoMo_BaP
    ## 126 609.9773    CoMo_BaP
    ## 127 488.0952    CoMo_BaP
    ## 128 676.5432    CoMo_BaP
    ## 129 457.5000    CoMo_BaP
    ## 130 628.6420    CoMo_BaP
    ## 131 414.3750    CoMo_BaP
    ## 132 561.5705  AHR2Mo_BaP
    ## 133 466.8750  AHR2Mo_BaP
    ## 134 587.8685  AHR2Mo_BaP
    ## 135 443.1886  AHR2Mo_BaP
    ## 136 596.5909  AHR2Mo_BaP
    ## 137 571.9955  AHR2Mo_BaP
    ## 138 521.9037  AHR2Mo_BaP
    ## 139 472.0579  AHR2Mo_BaP

What covariates predict F1 gut metagenome richness?

    # step AIC to test which model
    testmod_richF1 <- lm(pwF1_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino + Sex:BMI, 
                          data = alphadiv_pathF1)
    AIC_F1 <- stepAIC(testmod_richF1)
    ## Start:  AIC=561.44
    ## pwF1_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Exposure:Sex         1     7.390 6940.8 559.59
    ## - Morpholino:Sex       1    12.468 6945.9 559.69
    ## - Sex:BMI              1    19.106 6952.5 559.82
    ## <none>                             6933.4 561.44
    ## - Exposure:Morpholino  1   158.662 7092.0 562.58
    ## 
    ## Step:  AIC=559.59
    ## pwF1_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Morpholino:Sex + Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Morpholino:Sex       1    12.600 6953.4 557.84
    ## - Sex:BMI              1    20.377 6961.2 557.99
    ## <none>                             6940.8 559.59
    ## - Exposure:Morpholino  1   160.712 7101.5 560.77
    ## 
    ## Step:  AIC=557.84
    ## pwF1_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Sex:BMI              1    20.084 6973.5 556.24
    ## <none>                             6953.4 557.84
    ## - Exposure:Morpholino  1   157.967 7111.3 558.96
    ## 
    ## Step:  AIC=556.24
    ## pwF1_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - BMI                  1     0.016 6973.5 554.24
    ## - Sex                  1    84.461 7057.9 555.91
    ## <none>                             6973.5 556.24
    ## - Exposure:Morpholino  1   160.876 7134.3 557.41
    ## 
    ## Step:  AIC=554.24
    ## pwF1_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## <none>                             6973.5 554.24
    ## - Exposure:Morpholino  1    167.87 7141.4 555.55
    ## - Sex                  1    177.04 7150.5 555.72


    # best model is Exposure + Morpholino + Sex + Exposure:Morpholino but not significant (has periods tho)

    # save object
    lm_richF1 <- summary(AIC_F1)
    lm_richF1
    ## 
    ## Call:
    ## lm(formula = pwF1_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino, 
    ##     data = alphadiv_pathF1)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.4352  -6.4063  -0.1759   5.2661  15.3220 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    33.790      1.365  24.757   <2e-16 ***
    ## ExposureBaP                    -1.371      1.727  -0.794   0.4286    
    ## MorpholinoAHR2Mo               -1.643      1.763  -0.932   0.3531    
    ## SexM                            2.259      1.225   1.844   0.0673 .  
    ## ExposureBaP:MorpholinoAHR2Mo    4.400      2.450   1.796   0.0747 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.214 on 134 degrees of freedom
    ## Multiple R-squared:  0.05157,    Adjusted R-squared:  0.02326 
    ## F-statistic: 1.822 on 4 and 134 DF,  p-value: 0.1283


    # test just sex + treatment + interaction
    testmod_rich2F1 <- lm(pwF1_richness ~ Treatments + Sex + BMI + Sex:Treatments, 
                          data = alphadiv_pathF1)
    summary(stepAIC(testmod_rich2F1))
    ## Start:  AIC=561.73
    ## pwF1_richness ~ Treatments + Sex + BMI + Sex:Treatments
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - Treatments:Sex  3   25.4960 6973.5 556.24
    ## - BMI             1    0.0029 6948.0 559.73
    ## <none>                        6948.0 561.73
    ## 
    ## Step:  AIC=556.24
    ## pwF1_richness ~ Treatments + Sex + BMI
    ## 
    ##              Df Sum of Sq    RSS    AIC
    ## - Treatments  3   198.064 7171.5 554.13
    ## - BMI         1     0.016 6973.5 554.24
    ## - Sex         1    84.461 7057.9 555.91
    ## <none>                    6973.5 556.24
    ## 
    ## Step:  AIC=554.13
    ## pwF1_richness ~ Sex + BMI
    ## 
    ##        Df Sum of Sq    RSS    AIC
    ## - BMI   1     7.113 7178.6 552.27
    ## - Sex   1    52.208 7223.7 553.14
    ## <none>              7171.5 554.13
    ## 
    ## Step:  AIC=552.27
    ## pwF1_richness ~ Sex
    ## 
    ##        Df Sum of Sq    RSS    AIC
    ## <none>              7178.6 552.27
    ## - Sex   1    174.02 7352.7 553.60
    ## 
    ## Call:
    ## lm(formula = pwF1_richness ~ Sex, data = alphadiv_pathF1)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.6522  -5.6522   0.3478   5.3478  15.5857 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  33.4143     0.8652  38.621   <2e-16 ***
    ## SexM          2.2379     1.2280   1.822   0.0706 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.239 on 137 degrees of freedom
    ## Multiple R-squared:  0.02367,    Adjusted R-squared:  0.01654 
    ## F-statistic: 3.321 on 1 and 137 DF,  p-value: 0.07057

    # sex best model but insigificant

What covariates predict F1 gut metagenome shannon diversity?

    # step AIC to test which model
    testmod_shanF1 <- lm(pwF1_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + Sex:Morpholino + Exposure:Morpholino + Sex:BMI, 
                        data = alphadiv_pathF1)
    AIC_F1 <- stepAIC(testmod_shanF1)
    ## Start:  AIC=-705.08
    ## pwF1_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + 
    ##     Sex:Morpholino + Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Sex         1 0.0009684 0.76624 -706.90
    ## - Morpholino:Sex       1 0.0034297 0.76870 -706.46
    ## - Sex:BMI              1 0.0035569 0.76883 -706.43
    ## <none>                             0.76527 -705.08
    ## - Exposure:Morpholino  1 0.0146252 0.77990 -704.45
    ## 
    ## Step:  AIC=-706.9
    ## pwF1_shannon ~ Exposure + Morpholino + Sex + BMI + Morpholino:Sex + 
    ##     Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Sex:BMI              1 0.0033805 0.76962 -708.29
    ## - Morpholino:Sex       1 0.0034048 0.76965 -708.29
    ## <none>                             0.76624 -706.90
    ## - Exposure:Morpholino  1 0.0148485 0.78109 -706.23
    ## 
    ## Step:  AIC=-708.29
    ## pwF1_shannon ~ Exposure + Morpholino + Sex + BMI + Morpholino:Sex + 
    ##     Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - BMI                  1 0.0030374 0.77266 -709.74
    ## - Morpholino:Sex       1 0.0033428 0.77296 -709.69
    ## <none>                             0.76962 -708.29
    ## - Exposure:Morpholino  1 0.0145118 0.78413 -707.69
    ## 
    ## Step:  AIC=-709.74
    ## pwF1_shannon ~ Exposure + Morpholino + Sex + Morpholino:Sex + 
    ##     Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Morpholino:Sex       1 0.0032702 0.77593 -711.16
    ## <none>                             0.77266 -709.74
    ## - Exposure:Morpholino  1 0.0123810 0.78504 -709.53
    ## 
    ## Step:  AIC=-711.16
    ## pwF1_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Sex                  1 0.0087807 0.78471 -711.59
    ## <none>                             0.77593 -711.16
    ## - Exposure:Morpholino  1 0.0128604 0.78879 -710.87
    ## 
    ## Step:  AIC=-711.59
    ## pwF1_shannon ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## <none>                             0.78471 -711.59
    ## - Exposure:Morpholino  1  0.013317 0.79803 -711.25

    lm_shanF1 <- summary(AIC_F1)
    lm_shanF1
    ## 
    ## Call:
    ## lm(formula = pwF1_shannon ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = alphadiv_pathF1)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.258414 -0.036404  0.008127  0.049275  0.245217 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   4.4736106  0.0130752 342.145  < 2e-16 ***
    ## ExposureBaP                   0.0650879  0.0182325   3.570 0.000495 ***
    ## MorpholinoAHR2Mo              0.0001005  0.0186307   0.005 0.995705    
    ## ExposureBaP:MorpholinoAHR2Mo -0.0391796  0.0258849  -1.514 0.132463    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07624 on 135 degrees of freedom
    ## Multiple R-squared:  0.1124, Adjusted R-squared:  0.09264 
    ## F-statistic: 5.696 on 3 and 135 DF,  p-value: 0.001057


    # best model is Exposure + Morpholino + Exposure:Morpholino
    # ExposureBaP significant

    # test just sex + treatment + interaction
    testmod_shan2F1 <- lm(pwF1_shannon ~ Sex + Treatments + BMI + Sex:Treatments, 
                          data = alphadiv_pathF1)
    summary(stepAIC(testmod_shan2F1))
    ## Start:  AIC=-704.51
    ## pwF1_shannon ~ Sex + Treatments + BMI + Sex:Treatments
    ## 
    ##                  Df Sum of Sq     RSS     AIC
    ## - Sex:Treatments  3 0.0045682 0.77296 -709.69
    ## - BMI             1 0.0024845 0.77088 -706.06
    ## <none>                        0.76840 -704.51
    ## 
    ## Step:  AIC=-709.69
    ## pwF1_shannon ~ Sex + Treatments + BMI
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - Sex         1  0.000593 0.77356 -711.58
    ## - BMI         1  0.002965 0.77593 -711.16
    ## <none>                    0.77296 -709.69
    ## - Treatments  3  0.098433 0.87140 -699.03
    ## 
    ## Step:  AIC=-711.58
    ## pwF1_shannon ~ Treatments + BMI
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - BMI         1  0.011152 0.78471 -711.59
    ## <none>                    0.77356 -711.58
    ## - Treatments  3  0.100170 0.87373 -700.66
    ## 
    ## Step:  AIC=-711.59
    ## pwF1_shannon ~ Treatments
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## <none>                    0.78471 -711.59
    ## - Treatments  3  0.099334 0.88404 -701.02
    ## 
    ## Call:
    ## lm(formula = pwF1_shannon ~ Treatments, data = alphadiv_pathF1)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.258414 -0.036404  0.008127  0.049275  0.245217 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           4.4736106  0.0130752 342.145  < 2e-16 ***
    ## TreatmentsCoMo_BaP    0.0650879  0.0182325   3.570 0.000495 ***
    ## TreatmentsAHR2Mo_DMSO 0.0001005  0.0186307   0.005 0.995705    
    ## TreatmentsAHR2Mo_BaP  0.0260088  0.0182325   1.427 0.156030    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07624 on 135 degrees of freedom
    ## Multiple R-squared:  0.1124, Adjusted R-squared:  0.09264 
    ## F-statistic: 5.696 on 3 and 135 DF,  p-value: 0.001057


    # best model is Treatments
    # TreatmentsCoMo_BaP significant

    plot(lm(formula = pwF1_shannon ~ Treatments, data = alphadiv_pathF1))

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F1-1.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F1-2.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F1-3.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F1-4.png" width="98%" height="98%" />

    # linear
    # normal
    # homoskedastic
    # no outliers

    ## combine, print to console, and save
    lm_alphasF1 <- cbind(lm_richF1, lm_shanF1)

    pwF1_richplot <- ggplot(alphadiv_pathF1, aes(x = Treatments, y = pwF1_richness)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Pathway Richness") + 
      scale_y_continuous(limits = c(20, 55))

    pwF1_shanplot <- ggplot(alphadiv_pathF1, aes(x = Treatments, y = pwF1_shannon)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Shannon Diversity Metric") +
      scale_y_continuous(limits = c(4.2, 4.8))

    pwF1_alphaplot <- ggarrange(pwF1_richplot, pwF1_shanplot,
                            labels = c("A", "B"),
                            ncol = 1, nrow = 2)
    pwF1_alphaplot

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F1 alpha diversity plotting-1.png" width="98%" height="98%" />

# F2

    ## calculate alpha diversity metrics and put into a data frame
    pwF2_richness <- rowSums(pathcov_tab_F2)
    pwF2_shannon <- diversity(pathabund_tab_F2, "shannon", base = exp(1))
    pwF2_simpson <- diversity(pathabund_tab_F2, "simpson")
    meta_modF2 <- metaF2 %>%
      dplyr::select(c("SampleID", "Exposure", "Morpholino", "Sex", "BMI", "Treatments")) %>% # filter to only important variables we are testing
      mutate(Treatments = factor(Treatments))
    alphadiv_pathF2 <- data.frame(pwF2_richness, pwF2_shannon, pwF2_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(meta_modF2, by = "SampleID")
    alphadiv_pathF2
    ##     SampleID pwF2_richness pwF2_shannon pwF2_simpson Exposure Morpholino Sex
    ## 1   BaP_0001            34     4.550223    0.9885438      BaP     AHR2Mo   F
    ## 2   BaP_0002            42     4.502771    0.9877448      BaP     AHR2Mo   F
    ## 3   BaP_0025            38     4.462421    0.9871796     DMSO       CoMo   F
    ## 4   BaP_0026            38     4.393515    0.9865654     DMSO       CoMo   M
    ## 5   BaP_0027            38     4.325567    0.9856425     DMSO       CoMo   F
    ## 6   BaP_0028            35     4.385863    0.9863195     DMSO       CoMo   M
    ## 7   BaP_0029            24     4.411667    0.9864888     DMSO       CoMo   F
    ## 8   BaP_0030            42     4.405563    0.9867490     DMSO       CoMo   M
    ## 9   BaP_0031            34     4.482217    0.9877155     DMSO       CoMo   F
    ## 10  BaP_0032            42     4.427853    0.9870069     DMSO       CoMo   M
    ## 11  BaP_0033            40     4.408737    0.9869127     DMSO       CoMo   F
    ## 12  BaP_0034            33     4.495006    0.9878998     DMSO       CoMo   M
    ## 13  BaP_0035            40     4.498149    0.9879687     DMSO       CoMo   F
    ## 14  BaP_0036            27     4.554704    0.9888240     DMSO       CoMo   M
    ## 15  BaP_0061            29     4.518762    0.9883991     DMSO       CoMo   F
    ## 16  BaP_0062            39     4.388660    0.9862537     DMSO       CoMo   M
    ## 17  BaP_0063            37     4.374832    0.9864720     DMSO       CoMo   F
    ## 18  BaP_0064            39     4.395850    0.9864343     DMSO       CoMo   M
    ## 19  BaP_0065            35     4.368612    0.9862677     DMSO       CoMo   F
    ## 20  BaP_0066            33     4.377396    0.9860513     DMSO       CoMo   M
    ## 21  BaP_0067            25     4.542407    0.9887608     DMSO       CoMo   F
    ## 22  BaP_0068            38     4.394275    0.9863954     DMSO       CoMo   M
    ## 23  BaP_0069            35     4.389691    0.9863123     DMSO       CoMo   F
    ## 24  BaP_0070            30     4.540525    0.9886508     DMSO       CoMo   F
    ## 25  BaP_0071            26     4.361893    0.9863425     DMSO       CoMo   F
    ## 26  BaP_0072            33     4.434827    0.9869978      BaP       CoMo   F
    ## 27  BaP_0089            35     4.472215    0.9874612      BaP     AHR2Mo   M
    ## 28  BaP_0090            40     4.468907    0.9870658      BaP     AHR2Mo   M
    ## 29  BaP_0097            34     4.532889    0.9884338      BaP       CoMo   F
    ## 30  BaP_0098            41     4.451153    0.9869110      BaP     AHR2Mo   F
    ## 31  BaP_0105            24     4.395963    0.9867431      BaP       CoMo   F
    ## 32  BaP_0106            42     4.465500    0.9871323      BaP       CoMo   M
    ## 33  BaP_0107            30     4.428108    0.9869742      BaP       CoMo   M
    ## 34  BaP_0108            38     4.444523    0.9867667      BaP       CoMo   F
    ## 35  BaP_0137            38     4.508693    0.9881372     DMSO       CoMo   F
    ## 36  BaP_0138            33     4.528405    0.9884916     DMSO       CoMo   M
    ## 37  BaP_0139            26     4.501309    0.9880717     DMSO       CoMo   F
    ## 38  BaP_0140            27     4.526699    0.9884871     DMSO       CoMo   M
    ## 39  BaP_0141            24     4.526364    0.9885085      BaP       CoMo   F
    ## 40  BaP_0142            36     4.495184    0.9879463      BaP     AHR2Mo   M
    ## 41  BaP_0143            35     4.477761    0.9876701      BaP       CoMo   M
    ## 42  BaP_0144            39     4.538766    0.9882522      BaP       CoMo   F
    ## 43  BaP_0441            43     4.517536    0.9882380      BaP     AHR2Mo   F
    ## 44  BaP_0442            28     4.572098    0.9888556      BaP     AHR2Mo   F
    ## 45  BaP_0457            27     4.461372    0.9873845      BaP       CoMo   M
    ## 46  BaP_0458            38     4.565708    0.9888064     DMSO     AHR2Mo   F
    ## 47  BaP_0459            42     4.502028    0.9880679     DMSO     AHR2Mo   M
    ## 48  BaP_0460            41     4.516054    0.9882113     DMSO     AHR2Mo   F
    ## 49  BaP_0461            39     4.530688    0.9884452     DMSO     AHR2Mo   M
    ## 50  BaP_0462            35     4.586090    0.9890140     DMSO     AHR2Mo   F
    ## 51  BaP_0463            30     4.531274    0.9884541     DMSO     AHR2Mo   M
    ## 52  BaP_0464            34     4.533524    0.9882487     DMSO     AHR2Mo   F
    ## 53  BaP_0465            30     4.517963    0.9881035     DMSO       CoMo   F
    ## 54  BaP_0466            40     4.494190    0.9878048     DMSO       CoMo   M
    ## 55  BaP_0467            25     4.525114    0.9884395     DMSO       CoMo   F
    ## 56  BaP_0468            31     4.599095    0.9891943     DMSO       CoMo   M
    ## 57  BaP_0565            30     4.548769    0.9885137     DMSO     AHR2Mo   M
    ## 58  BaP_0566            39     4.543399    0.9885190     DMSO     AHR2Mo   F
    ## 59  BaP_0567            24     4.520689    0.9883654     DMSO     AHR2Mo   M
    ## 60  BaP_0568            41     4.481282    0.9876492     DMSO     AHR2Mo   F
    ## 61  BaP_0569            38     4.573315    0.9887190     DMSO     AHR2Mo   M
    ## 62  BaP_0570            49     4.371779    0.9859827     DMSO     AHR2Mo   F
    ## 63  BaP_0571            26     4.518165    0.9883440     DMSO     AHR2Mo   M
    ## 64  BaP_0572            43     4.421724    0.9868799     DMSO     AHR2Mo   F
    ## 65  BaP_0573            34     4.529663    0.9884505      BaP       CoMo   M
    ## 66  BaP_0574            42     4.481864    0.9877075      BaP       CoMo   F
    ## 67  BaP_0575            39     4.542426    0.9886135      BaP       CoMo   M
    ## 68  BaP_0576            42     4.489621    0.9878006      BaP       CoMo   F
    ## 69  BaP_0649            25     4.514694    0.9882774      BaP     AHR2Mo   M
    ## 70  BaP_0650            38     4.502033    0.9880862      BaP     AHR2Mo   M
    ## 71  BaP_0673            48     4.498142    0.9880294      BaP       CoMo   M
    ## 72  BaP_0674            30     4.549705    0.9886907      BaP       CoMo   F
    ## 73  BaP_0675            40     4.557115    0.9886677      BaP       CoMo   M
    ## 74  BaP_0676            37     4.552293    0.9886820      BaP       CoMo   F
    ## 75  BaP_0677            32     4.478903    0.9878334      BaP       CoMo   M
    ## 76  BaP_0678            40     4.578955    0.9888786      BaP       CoMo   F
    ## 77  BaP_0679            29     4.505667    0.9880266      BaP       CoMo   M
    ## 78  BaP_0680            28     4.577815    0.9888212      BaP       CoMo   F
    ## 79  BaP_0681            43     4.594604    0.9888714      BaP       CoMo   M
    ## 80  BaP_0683            38     4.587992    0.9888858      BaP       CoMo   M
    ## 81  BaP_0781            33     4.501919    0.9881847      BaP       CoMo   M
    ## 82  BaP_0782            33     4.532279    0.9885477      BaP       CoMo   F
    ## 83  BaP_0783            23     4.414131    0.9856429      BaP       CoMo   M
    ## 84  BaP_0784            43     4.494852    0.9878793      BaP       CoMo   F
    ## 85  BaP_0785            42     4.573259    0.9887914      BaP       CoMo   M
    ## 86  BaP_0786            33     4.456109    0.9873294      BaP       CoMo   F
    ## 87  BaP_0787            40     4.608135    0.9890402      BaP       CoMo   M
    ## 88  BaP_0788            36     4.558937    0.9886826      BaP       CoMo   F
    ## 89  BaP_0789            40     4.501261    0.9880634      BaP       CoMo   M
    ## 90  BaP_0790            26     4.514630    0.9881663      BaP       CoMo   F
    ## 91  BaP_0791            46     4.632574    0.9893353      BaP       CoMo   M
    ## 92  BaP_0792            41     4.529855    0.9883095      BaP       CoMo   F
    ## 93  BaP_0865            37     4.430378    0.9870195      BaP     AHR2Mo   F
    ## 94  BaP_0866            26     4.525773    0.9884628      BaP     AHR2Mo   F
    ## 95  BaP_0889            27     4.499857    0.9879037     DMSO     AHR2Mo   M
    ## 96  BaP_0890            27     4.550619    0.9885459     DMSO     AHR2Mo   F
    ## 97  BaP_0891            42     4.596655    0.9888162     DMSO     AHR2Mo   M
    ## 98  BaP_0893            40     4.436523    0.9871808     DMSO     AHR2Mo   F
    ## 99  BaP_0895            26     4.533494    0.9884286     DMSO     AHR2Mo   M
    ## 100 BaP_0897            26     4.520421    0.9881975     DMSO     AHR2Mo   F
    ## 101 BaP_0899            25     4.519657    0.9882915     DMSO     AHR2Mo   M
    ## 102 BaP_0961            32     4.529736    0.9880896     DMSO     AHR2Mo   F
    ## 103 BaP_0962            24     4.486389    0.9876956     DMSO     AHR2Mo   M
    ## 104 BaP_0963            24     4.464565    0.9873303     DMSO     AHR2Mo   F
    ## 105 BaP_0964            41     4.492813    0.9876432     DMSO     AHR2Mo   M
    ## 106 BaP_0965            25     4.536049    0.9884101     DMSO     AHR2Mo   F
    ## 107 BaP_0966            26     4.529795    0.9884617     DMSO     AHR2Mo   M
    ## 108 BaP_0967            48     4.673877    0.9895179     DMSO     AHR2Mo   F
    ## 109 BaP_0968            27     4.538327    0.9884462     DMSO     AHR2Mo   M
    ## 110 BaP_0969            36     4.471660    0.9875338      BaP     AHR2Mo   F
    ## 111 BaP_0970            38     4.535278    0.9883476      BaP     AHR2Mo   M
    ## 112 BaP_0971            53     4.470996    0.9875761      BaP     AHR2Mo   F
    ## 113 BaP_0972            26     4.469410    0.9873433      BaP     AHR2Mo   M
    ## 114 BaP_1017            27     4.531667    0.9885300      BaP     AHR2Mo   M
    ## 115 BaP_1018            49     4.425298    0.9872458      BaP     AHR2Mo   M
    ## 116 BaP_1033            28     4.546102    0.9883440     DMSO     AHR2Mo   F
    ## 117 BaP_1034            38     4.304793    0.9851062     DMSO     AHR2Mo   M
    ## 118 BaP_1035            41     4.549272    0.9881433     DMSO     AHR2Mo   F
    ## 119 BaP_1036            25     4.540316    0.9884947     DMSO     AHR2Mo   M
    ## 120 BaP_1037            37     4.382060    0.9860371     DMSO     AHR2Mo   F
    ## 121 BaP_1038            27     4.546698    0.9885625     DMSO     AHR2Mo   M
    ## 122 BaP_1039            44     4.573911    0.9884341      BaP     AHR2Mo   F
    ## 123 BaP_1040            32     4.460502    0.9875155      BaP     AHR2Mo   M
    ## 124 BaP_1041            31     4.506099    0.9879788      BaP     AHR2Mo   F
    ## 125 BaP_1042            24     4.563573    0.9886930      BaP     AHR2Mo   M
    ## 126 BaP_1043            32     4.420825    0.9869515      BaP     AHR2Mo   F
    ## 127 BaP_1044            41     4.540977    0.9880741      BaP     AHR2Mo   M
    ## 128 BaP_1105            28     4.591510    0.9889295      BaP     AHR2Mo   F
    ## 129 BaP_1106            45     4.605341    0.9889344      BaP     AHR2Mo   M
    ## 130 BaP_1107            26     4.534361    0.9882333      BaP     AHR2Mo   F
    ## 131 BaP_1108            42     4.519372    0.9880783      BaP     AHR2Mo   M
    ## 132 BaP_1109            28     4.540853    0.9883700      BaP     AHR2Mo   F
    ## 133 BaP_1110            28     4.613805    0.9891114      BaP     AHR2Mo   M
    ## 134 BaP_1111            26     4.478225    0.9876314      BaP     AHR2Mo   F
    ## 135 BaP_1112            23     4.423231    0.9864562      BaP     AHR2Mo   M
    ## 136 BaP_1113            26     4.561057    0.9884440      BaP     AHR2Mo   F
    ## 137 BaP_1114            29     4.598979    0.9887733      BaP     AHR2Mo   M
    ## 138 BaP_1115            40     4.589359    0.9887910      BaP     AHR2Mo   F
    ## 139 BaP_1116            36     4.671099    0.9893012      BaP     AHR2Mo   M
    ##           BMI  Treatments
    ## 1    764.4444  AHR2Mo_BaP
    ## 2    641.2742  AHR2Mo_BaP
    ## 3    927.4691   CoMo_DMSO
    ## 4    668.0428   CoMo_DMSO
    ## 5    604.6713   CoMo_DMSO
    ## 6    659.5253   CoMo_DMSO
    ## 7    818.5604   CoMo_DMSO
    ## 8    647.3829   CoMo_DMSO
    ## 9   1070.6696   CoMo_DMSO
    ## 10   634.7699   CoMo_DMSO
    ## 11   842.2206   CoMo_DMSO
    ## 12   682.6172   CoMo_DMSO
    ## 13   588.7500   CoMo_DMSO
    ## 14   534.6260   CoMo_DMSO
    ## 15   787.1094   CoMo_DMSO
    ## 16   709.0947   CoMo_DMSO
    ## 17  1053.1250   CoMo_DMSO
    ## 18   551.1111   CoMo_DMSO
    ## 19   949.5465   CoMo_DMSO
    ## 20   589.9654   CoMo_DMSO
    ## 21   605.5515   CoMo_DMSO
    ## 22   528.3556   CoMo_DMSO
    ## 23   949.2326   CoMo_DMSO
    ## 24   831.9559   CoMo_DMSO
    ## 25   852.4470   CoMo_DMSO
    ## 26   727.5383    CoMo_BaP
    ## 27   549.4644  AHR2Mo_BaP
    ## 28   568.6728  AHR2Mo_BaP
    ## 29   854.7816    CoMo_BaP
    ## 30   827.1209  AHR2Mo_BaP
    ## 31   701.2418    CoMo_BaP
    ## 32   569.3918    CoMo_BaP
    ## 33   588.0140    CoMo_BaP
    ## 34   897.5069    CoMo_BaP
    ## 35   868.7700   CoMo_DMSO
    ## 36   577.4222   CoMo_DMSO
    ## 37  1242.2145   CoMo_DMSO
    ## 38   619.1467   CoMo_DMSO
    ## 39   955.6250    CoMo_BaP
    ## 40   669.5502  AHR2Mo_BaP
    ## 41   513.8504    CoMo_BaP
    ## 42   772.2681    CoMo_BaP
    ## 43   839.3352  AHR2Mo_BaP
    ## 44   709.3750  AHR2Mo_BaP
    ## 45   662.0753    CoMo_BaP
    ## 46   889.4558 AHR2Mo_DMSO
    ## 47   768.4666 AHR2Mo_DMSO
    ## 48  1011.7729 AHR2Mo_DMSO
    ## 49   623.0469 AHR2Mo_DMSO
    ## 50   780.3688 AHR2Mo_DMSO
    ## 51   687.3630 AHR2Mo_DMSO
    ## 52  1048.9796 AHR2Mo_DMSO
    ## 53   929.8037   CoMo_DMSO
    ## 54   654.6939   CoMo_DMSO
    ## 55  1035.3186   CoMo_DMSO
    ## 56   886.0482   CoMo_DMSO
    ## 57   755.5556 AHR2Mo_DMSO
    ## 58   803.3241 AHR2Mo_DMSO
    ## 59   612.7385 AHR2Mo_DMSO
    ## 60   950.1385 AHR2Mo_DMSO
    ## 61   641.8685 AHR2Mo_DMSO
    ## 62  1121.1911 AHR2Mo_DMSO
    ## 63   633.4964 AHR2Mo_DMSO
    ## 64   980.9336 AHR2Mo_DMSO
    ## 65   529.4118    CoMo_BaP
    ## 66   976.2655    CoMo_BaP
    ## 67   609.4183    CoMo_BaP
    ## 68   683.4568    CoMo_BaP
    ## 69   671.8646  AHR2Mo_BaP
    ## 70   558.0716  AHR2Mo_BaP
    ## 71   570.6122    CoMo_BaP
    ## 72   967.8598    CoMo_BaP
    ## 73   605.7099    CoMo_BaP
    ## 74   804.4077    CoMo_BaP
    ## 75   654.7291    CoMo_BaP
    ## 76   820.4082    CoMo_BaP
    ## 77   616.5333    CoMo_BaP
    ## 78   869.3750    CoMo_BaP
    ## 79   593.1337    CoMo_BaP
    ## 80   659.3896    CoMo_BaP
    ## 81   567.4611    CoMo_BaP
    ## 82   798.2543    CoMo_BaP
    ## 83   770.4264    CoMo_BaP
    ## 84   936.8836    CoMo_BaP
    ## 85   600.4383    CoMo_BaP
    ## 86   831.7175    CoMo_BaP
    ## 87   601.8519    CoMo_BaP
    ## 88   708.5896    CoMo_BaP
    ## 89   695.3125    CoMo_BaP
    ## 90  1007.8895    CoMo_BaP
    ## 91   722.8724    CoMo_BaP
    ## 92   755.3035    CoMo_BaP
    ## 93   681.6108  AHR2Mo_BaP
    ## 94   637.7383  AHR2Mo_BaP
    ## 95   501.9531 AHR2Mo_DMSO
    ## 96   767.3130 AHR2Mo_DMSO
    ## 97   709.2768 AHR2Mo_DMSO
    ## 98   755.5402 AHR2Mo_DMSO
    ## 99   613.8135 AHR2Mo_DMSO
    ## 100 1015.3483 AHR2Mo_DMSO
    ## 101  641.9753 AHR2Mo_DMSO
    ## 102  769.8882 AHR2Mo_DMSO
    ## 103  615.5102 AHR2Mo_DMSO
    ## 104  963.7188 AHR2Mo_DMSO
    ## 105  664.4898 AHR2Mo_DMSO
    ## 106  953.7500 AHR2Mo_DMSO
    ## 107  649.4141 AHR2Mo_DMSO
    ## 108 1033.5306 AHR2Mo_DMSO
    ## 109  553.0649 AHR2Mo_DMSO
    ## 110  689.4923  AHR2Mo_BaP
    ## 111  612.2449  AHR2Mo_BaP
    ## 112  886.4901  AHR2Mo_BaP
    ## 113  597.5510  AHR2Mo_BaP
    ## 114  708.5714  AHR2Mo_BaP
    ## 115  538.5802  AHR2Mo_BaP
    ## 116  733.7963 AHR2Mo_DMSO
    ## 117  646.1938 AHR2Mo_DMSO
    ## 118  779.7784 AHR2Mo_DMSO
    ## 119  636.8410 AHR2Mo_DMSO
    ## 120  925.7143 AHR2Mo_DMSO
    ## 121  697.2318 AHR2Mo_DMSO
    ## 122  789.6416  AHR2Mo_BaP
    ## 123  721.4533  AHR2Mo_BaP
    ## 124  799.1234  AHR2Mo_BaP
    ## 125  688.8889  AHR2Mo_BaP
    ## 126  907.8947  AHR2Mo_BaP
    ## 127  610.9964  AHR2Mo_BaP
    ## 128  949.9072  AHR2Mo_BaP
    ## 129  601.8519  AHR2Mo_BaP
    ## 130  556.1857  AHR2Mo_BaP
    ## 131  533.8776  AHR2Mo_BaP
    ## 132  973.5410  AHR2Mo_BaP
    ## 133  606.5306  AHR2Mo_BaP
    ## 134  951.2500  AHR2Mo_BaP
    ## 135  657.6075  AHR2Mo_BaP
    ## 136  803.0764  AHR2Mo_BaP
    ## 137  661.7969  AHR2Mo_BaP
    ## 138  803.9032  AHR2Mo_BaP
    ## 139  651.6620  AHR2Mo_BaP

What covariates predict F2 gut metagenome richness?

    # step AIC to test which model
    testmod_richF2 <- lm(pwF2_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino + Sex:BMI, 
                          data = alphadiv_pathF2)
    AIC_F2 <- stepAIC(testmod_richF2)
    ## Start:  AIC=546.26
    ## pwF2_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Exposure:Morpholino  1     7.884 6224.0 544.44
    ## - Exposure:Sex         1    20.093 6236.2 544.71
    ## - Sex:BMI              1    39.719 6255.8 545.14
    ## <none>                             6216.1 546.26
    ## - Morpholino:Sex       1   231.073 6447.2 549.33
    ## 
    ## Step:  AIC=544.44
    ## pwF2_richness ~ Exposure + Morpholino + Sex + BMI + Exposure:Sex + 
    ##     Morpholino:Sex + Sex:BMI
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - Exposure:Sex    1    19.010 6243.0 542.86
    ## - Sex:BMI         1    41.062 6265.1 543.35
    ## <none>                        6224.0 544.44
    ## - Morpholino:Sex  1   233.833 6457.8 547.56
    ## 
    ## Step:  AIC=542.86
    ## pwF2_richness ~ Exposure + Morpholino + Sex + BMI + Morpholino:Sex + 
    ##     Sex:BMI
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - Exposure        1    50.314 6293.3 541.98
    ## - Sex:BMI         1    55.072 6298.1 542.08
    ## <none>                        6243.0 542.86
    ## - Morpholino:Sex  1   239.678 6482.7 546.10
    ## 
    ## Step:  AIC=541.98
    ## pwF2_richness ~ Morpholino + Sex + BMI + Morpholino:Sex + Sex:BMI
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - Sex:BMI         1    56.205 6349.5 541.21
    ## <none>                        6293.3 541.98
    ## - Morpholino:Sex  1   242.803 6536.1 545.24
    ## 
    ## Step:  AIC=541.21
    ## pwF2_richness ~ Morpholino + Sex + BMI + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - BMI             1     5.799 6355.3 539.34
    ## <none>                        6349.5 541.21
    ## - Morpholino:Sex  1   248.573 6598.1 544.55
    ## 
    ## Step:  AIC=539.34
    ## pwF2_richness ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## <none>                        6355.3 539.34
    ## - Morpholino:Sex  1    252.31 6607.6 542.75

    # save object
    lm_richF2 <- summary(AIC_F2)
    lm_richF2
    ## 
    ## Call:
    ## lm(formula = pwF2_richness ~ Morpholino + Sex + Morpholino:Sex, 
    ##     data = alphadiv_pathF2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.1875  -6.5278   0.4857   5.4790  17.7500 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             33.514      1.160  28.898   <2e-16 ***
    ## MorpholinoAHR2Mo         1.736      1.629   1.066   0.2885    
    ## SexM                     2.673      1.678   1.593   0.1135    
    ## MorpholinoAHR2Mo:SexM   -5.395      2.331  -2.315   0.0221 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.861 on 135 degrees of freedom
    ## Multiple R-squared:  0.04236,    Adjusted R-squared:  0.02107 
    ## F-statistic:  1.99 on 3 and 135 DF,  p-value: 0.1184


    # best model is Morpholino + Sex + Morpholino:Sex with MorpholinoAHR2Mo:SexM significant

    # test just sex + treatment + interaction
    testmod_rich2F2 <- lm(pwF2_richness ~ Treatments + Sex + BMI + Sex:Treatments, 
                          data = alphadiv_pathF2)
    summary(stepAIC(testmod_rich2F2))
    ## Start:  AIC=545.51
    ## pwF2_richness ~ Treatments + Sex + BMI + Sex:Treatments
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## - BMI             1      5.62 6188.5 543.64
    ## <none>                        6182.9 545.51
    ## - Treatments:Sex  3    349.91 6532.8 547.17
    ## 
    ## Step:  AIC=543.64
    ## pwF2_richness ~ Treatments + Sex + Treatments:Sex
    ## 
    ##                  Df Sum of Sq    RSS    AIC
    ## <none>                        6188.5 543.64
    ## - Treatments:Sex  3    346.55 6535.0 545.21
    ## 
    ## Call:
    ## lm(formula = pwF2_richness ~ Treatments + Sex + Treatments:Sex, 
    ##     data = alphadiv_pathF2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -13.722  -6.028  -0.500   5.333  18.500 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                  32.353      1.667  19.408   <2e-16 ***
    ## TreatmentsCoMo_BaP            2.258      2.325   0.971   0.3331    
    ## TreatmentsAHR2Mo_DMSO         3.647      2.325   1.569   0.1191    
    ## TreatmentsAHR2Mo_BaP          2.147      2.325   0.924   0.3574    
    ## SexM                          3.147      2.481   1.269   0.2068    
    ## TreatmentsCoMo_BaP:SexM      -1.036      3.377  -0.307   0.7595    
    ## TreatmentsAHR2Mo_DMSO:SexM   -8.203      3.377  -2.429   0.0165 *  
    ## TreatmentsAHR2Mo_BaP:SexM    -3.536      3.377  -1.047   0.2970    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.873 on 131 degrees of freedom
    ## Multiple R-squared:  0.0675, Adjusted R-squared:  0.01767 
    ## F-statistic: 1.355 on 7 and 131 DF,  p-value: 0.23


    # Treatments + Sex + Treatments:Sex best model but none significant

What covariates predict F2 gut metagenome shannon diversity?

    # step AIC to test which model
    testmod_shanF2 <- lm(pwF2_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + Sex:Morpholino + Exposure:Morpholino + Sex:BMI, 
                        data = alphadiv_pathF2)
    AIC_F2 <- stepAIC(testmod_shanF2)
    ## Start:  AIC=-753.3
    ## pwF2_shannon ~ Exposure + Morpholino + Sex + BMI + Sex:Exposure + 
    ##     Sex:Morpholino + Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Morpholino:Sex       1 0.0000060 0.54094 -755.30
    ## - Exposure:Sex         1 0.0008326 0.54176 -755.09
    ## - Sex:BMI              1 0.0077761 0.54871 -753.32
    ## <none>                             0.54093 -753.30
    ## - Exposure:Morpholino  1 0.0262486 0.56718 -748.72
    ## 
    ## Step:  AIC=-755.3
    ## pwF2_shannon ~ Exposure + Morpholino + Sex + BMI + Exposure:Sex + 
    ##     Exposure:Morpholino + Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Sex         1 0.0008281 0.54176 -757.09
    ## - Sex:BMI              1 0.0077838 0.54872 -755.31
    ## <none>                             0.54094 -755.30
    ## - Exposure:Morpholino  1 0.0262483 0.56719 -750.71
    ## 
    ## Step:  AIC=-757.09
    ## pwF2_shannon ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino + 
    ##     Sex:BMI
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Sex:BMI              1 0.0070836 0.54885 -757.28
    ## <none>                             0.54176 -757.09
    ## - Exposure:Morpholino  1 0.0258895 0.56765 -752.60
    ## 
    ## Step:  AIC=-757.28
    ## pwF2_shannon ~ Exposure + Morpholino + Sex + BMI + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - BMI                  1 0.0008126 0.54966 -759.08
    ## - Sex                  1 0.0018874 0.55074 -758.81
    ## <none>                             0.54885 -757.28
    ## - Exposure:Morpholino  1 0.0251473 0.57400 -753.06
    ## 
    ## Step:  AIC=-759.08
    ## pwF2_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Sex                  1 0.0010941 0.55076 -760.80
    ## <none>                             0.54966 -759.08
    ## - Exposure:Morpholino  1 0.0261655 0.57583 -754.61
    ## 
    ## Step:  AIC=-760.8
    ## pwF2_shannon ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## <none>                             0.55076 -760.80
    ## - Exposure:Morpholino  1   0.02644 0.57720 -756.28

    lm_shanF2 <- summary(AIC_F2)
    lm_shanF2
    ## 
    ## Call:
    ## lm(formula = pwF2_shannon ~ Exposure + Morpholino + Exposure:Morpholino, 
    ##     data = alphadiv_pathF2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.209721 -0.048126  0.007401  0.041737  0.159363 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   4.45502    0.01147 388.345  < 2e-16 ***
    ## ExposureBaP                   0.06022    0.01565   3.848 0.000183 ***
    ## MorpholinoAHR2Mo              0.05949    0.01565   3.801 0.000217 ***
    ## ExposureBaP:MorpholinoAHR2Mo -0.05528    0.02172  -2.546 0.012027 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06387 on 135 degrees of freedom
    ## Multiple R-squared:  0.1422, Adjusted R-squared:  0.1231 
    ## F-statistic: 7.457 on 3 and 135 DF,  p-value: 0.0001172


    # best model is Exposure + Morpholino + Exposure:Morpholino
    # ExposureBaP, MorpholinoAHR2Mo, ExposureBaP:MorpholinoAHR2Mo significant

    # test just sex + treatment + interaction
    testmod_shan2F2 <- lm(pwF2_shannon ~ Sex + Treatments + BMI + Sex:Treatments, 
                          data = alphadiv_pathF2)
    summary(stepAIC(testmod_shan2F2))
    ## Start:  AIC=-751.38
    ## pwF2_shannon ~ Sex + Treatments + BMI + Sex:Treatments
    ## 
    ##                  Df  Sum of Sq     RSS     AIC
    ## - Sex:Treatments  3 0.00037043 0.54885 -757.28
    ## - BMI             1 0.00076458 0.54924 -753.18
    ## <none>                         0.54848 -751.38
    ## 
    ## Step:  AIC=-757.28
    ## pwF2_shannon ~ Sex + Treatments + BMI
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - BMI         1  0.000813 0.54966 -759.08
    ## - Sex         1  0.001887 0.55074 -758.81
    ## <none>                    0.54885 -757.28
    ## - Treatments  3  0.091115 0.63996 -741.93
    ## 
    ## Step:  AIC=-759.08
    ## pwF2_shannon ~ Sex + Treatments
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - Sex         1  0.001094 0.55076 -760.80
    ## <none>                    0.54966 -759.08
    ## - Treatments  3  0.090320 0.63998 -743.93
    ## 
    ## Step:  AIC=-760.8
    ## pwF2_shannon ~ Treatments
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## <none>                    0.55076 -760.80
    ## - Treatments  3  0.091268 0.64202 -745.49
    ## 
    ## Call:
    ## lm(formula = pwF2_shannon ~ Treatments, data = alphadiv_pathF2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.209721 -0.048126  0.007401  0.041737  0.159363 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            4.45502    0.01147 388.345  < 2e-16 ***
    ## TreatmentsCoMo_BaP     0.06022    0.01565   3.848 0.000183 ***
    ## TreatmentsAHR2Mo_DMSO  0.05949    0.01565   3.801 0.000217 ***
    ## TreatmentsAHR2Mo_BaP   0.06443    0.01565   4.117 6.64e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06387 on 135 degrees of freedom
    ## Multiple R-squared:  0.1422, Adjusted R-squared:  0.1231 
    ## F-statistic: 7.457 on 3 and 135 DF,  p-value: 0.0001172


    # best model is Treatments
    # TreatmentsCoMo_DMSO significant

    plot(lm(formula = pwF2_shannon ~ Treatments, data = alphadiv_pathF2))

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F2-1.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F2-2.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F2-3.png" width="98%" height="98%" /><img src="zfBaP_PubAnalysis_files/figure-markdown_strict/linear model shannon with covariates F2-4.png" width="98%" height="98%" />

    # linear
    # normal
    # homoskedastic
    # no outliers

    ## combine, print to console, and save
    lm_alphasF2 <- cbind(lm_richF2, lm_shanF2)

    pwF2_richplot <- ggplot(alphadiv_pathF2, aes(x = Treatments, y = pwF2_richness)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Pathway Richness") +
      scale_y_continuous(limits = c(20, 55))

    pwF2_shanplot <- ggplot(alphadiv_pathF2, aes(x = Treatments, y = pwF2_shannon)) + 
      geom_boxplot() + 
      theme_classic() +
      labs(x = NULL, y = "Shannon Diversity Metric") +
      scale_y_continuous(limits = c(4.2, 4.8))

    pwF2_alphaplot <- ggarrange(pwF2_richplot, pwF2_shanplot,
                            labels = c("A", "B"),
                            ncol = 1, nrow = 2)
    pwF2_alphaplot

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F2 alpha diversity plotting-1.png" width="98%" height="98%" />

## Beta diversity stats and plots

-   Used distance based RDA with bray curtis distance to create
    ordination.
-   Exploratory analysis revealed that bray curtis, jaccard, and
    sorenson data frame were correlated significantly (mantel test),
    meaning that the distance/dissimilarity measure does not change the
    information and ultimate interpretation downstream. We chose to
    proceed with bray curtis.
-   Ordistep used to choose the optimal model for dbRDA.
-   This optimal model was then used in the PERMANOVA to get
    significance.
-   Beta dispersion was also calculated and evaluated with linear
    regressions and graphed on a box and whisker plot.

# Whole data

Is beta diversity and dispersion of the gut metagenome associated with
generation, exposure, morpholino, or other variables?

    ## check model with all possible options (have checked sex with ordistep in the past and doesn't come out important to model)
    meta_modALL <- metaALL %>%
      dplyr::select(c("Exposure", "Morpholino", "Generation", "Sex", "BMI"))
    ALLmod0 <- capscale(pathabund_tab_ALL ~ 1, meta_modALL, distance = "bray")  # Model with intercept only
    ALLmod1 <- capscale(pathabund_tab_ALL ~ . + Exposure:Morpholino + Morpholino:Exposure + Generation:Exposure + Generation:Morpholino + Sex:BMI, meta_modALL, distance = "bray")  # Model with all explanatory variables
    pw_ordiALL <- ordistep(ALLmod0, scope = formula(ALLmod1)) # this determines what the best model is to run RDA on
    pw_ordiALL

    # best model is pathabund_tab_ALL ~ Generation + Exposure + Morpholino + Generation:Exposure + Exposure:Morpholino + Generation:Morpholino
    # model inertia was 5.4988; Inertia is scaled Chi-square
    # proportion variance explained is 0.2874

    # calculate RDA (distance based)
    pw_rda <- capscale(formula = pathabund_tab_ALL ~ Generation + Exposure + Morpholino + Generation:Exposure + Exposure:Morpholino + Generation:Morpholino, 
                       data = meta_modALL, distance = "bray")
    # using bray because euclidean, bray, and jaccard were all correlated via mantel test, and bray is consistent with other analyses in NMDS and other graphs in this analysis
    RsquareAdj(pw_rda)
    ## $r.squared
    ## [1] 0.2418108
    ## 
    ## $adj.r.squared
    ## [1] 0.2252484

    # 0.2718067 of variation explained by this model

    # and with generation as condition just to compare
    pw_rda2 <- capscale(formula = pathabund_tab_ALL ~ Condition(Generation) + Exposure + Morpholino + Exposure:Morpholino, 
                       data = meta_modALL, distance = "bray")

    # create files for graphing
    smry_rda <- summary(pw_rda)

    # this is the coordinates for the points and metadata
    pw_PC1  <- data.frame(smry_rda$sites[,1:2]) %>%  # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaALL, c(SampleID, Treatments, Generation, Exposure, Morpholino)), by = "SampleID") %>%
      mutate(treatgen = paste0(Treatments, "_", Generation)) %>%
      column_to_rownames("SampleID")
    pw_PC1$Treatments <- factor(x = pw_PC1$Treatments,
                                 levels = c("CoMo_DMSO",
                                            "CoMo_BaP",
                                            "AHR2Mo_DMSO",
                                            "AHR2Mo_BaP"))
    mylims <- range(with(pw_PC1, c(CAP1, CAP2)))

    # the biplot scores are the correlations between your environmental variables and axes
    pw_PC2  <- data.frame(smry_rda$biplot)

    # put plot together
    pwrda_plot <- ggplot(pw_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      stat_ellipse(aes(group = treatgen, 
                       color = Treatments), 
                   show.legend = T) +
      theme_classic() +
      labs(x = paste0("CAP1 (",round(100*smry_rda$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rda$cont$importance[2, "CAP2"], digits = 2),"%)")) +
      scale_color_brewer(palette = "PuOr", direction = -1,
                         labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 20)) +
      scale_shape_discrete(name = "Generation shape") +
      scale_linetype_discrete(name = "Ellipse line type") +
      facet_wrap("Generation")

    pwrda_plot

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/graph ordination-1.png" width="98%" height="98%" />

    pw_dm <- vegdist(pathabund_tab_ALL, method = "bray")
    PermExpandMod_pw <- adonis2(pw_dm ~ Generation + Exposure + Morpholino + Generation:Exposure + Exposure:Morpholino + Generation:Morpholino, data = metaALL)
    PermExpandMod_pw
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pw_dm ~ Generation + Exposure + Morpholino + Generation:Exposure + Exposure:Morpholino + Generation:Morpholino, data = metaALL)
    ##                        Df SumOfSqs      R2       F Pr(>F)    
    ## Generation              2   1.1709 0.21293 60.9419  0.001 ***
    ## Exposure                1   0.1020 0.01856 10.6216  0.001 ***
    ## Morpholino              1   0.0305 0.00555  3.1788  0.010 ** 
    ## Generation:Exposure     2   0.0900 0.01636  4.6819  0.001 ***
    ## Exposure:Morpholino     1   0.0664 0.01207  6.9086  0.001 ***
    ## Generation:Morpholino   2   0.0811 0.01475  4.2222  0.001 ***
    ## Residual              412   3.9579 0.71978                   
    ## Total                 421   5.4988 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Beta dipsersion results

    # calc dispersion for each generation/treatment group
    metaALL <- metaALL %>% mutate(treatgen = paste0(Treatments, "_", Generation))
    pw.disper <- betadisper(pw_dm, metaALL$treatgen)

    # make file for graphing
    pw.Disper1 <- data.frame(pw.disper$distances)
    colnames(pw.Disper1) <- "dispers"
    pw.Disper2 <- pw.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaALL, 
                 by = "SampleID")

    ## linear model approach
    # step AIC to test model
    testmod_bdisp <- lm(pw.disper$distances ~ Generation + Exposure + Morpholino + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino, data = pw.Disper2)
    AIC_bdispALL <- stepAIC(testmod_bdisp)
    ## Start:  AIC=-2416.51
    ## pw.disper$distances ~ Generation + Exposure + Morpholino + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino    1  0.002433 1.3141 -2417.7
    ## <none>                               1.3116 -2416.5
    ## - Generation:Exposure    2  0.015742 1.3274 -2415.5
    ## - Generation:Morpholino  2  0.032968 1.3446 -2410.0
    ## 
    ## Step:  AIC=-2417.73
    ## pw.disper$distances ~ Generation + Exposure + Morpholino + Generation:Exposure + 
    ##     Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## <none>                               1.3141 -2417.7
    ## - Generation:Exposure    2  0.015684 1.3298 -2416.7
    ## - Generation:Morpholino  2  0.033029 1.3471 -2411.3

    # save in object
    lm_betadisper <- summary(AIC_bdispALL)
    lm_betadisper
    ## 
    ## Call:
    ## lm(formula = pw.disper$distances ~ Generation + Exposure + Morpholino + 
    ##     Generation:Exposure + Generation:Morpholino, data = pw.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.07512 -0.03069 -0.01281  0.01517  0.58193 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    0.083441   0.008142  10.249   <2e-16 ***
    ## GenerationF1                   0.008263   0.011662   0.709   0.4790    
    ## GenerationF2                  -0.016293   0.011844  -1.376   0.1697    
    ## ExposureBaP                   -0.017794   0.009401  -1.893   0.0591 .  
    ## MorpholinoAHR2Mo               0.002962   0.009401   0.315   0.7529    
    ## GenerationF1:ExposureBaP       0.027597   0.013419   2.057   0.0404 *  
    ## GenerationF2:ExposureBaP       0.023332   0.013424   1.738   0.0829 .  
    ## GenerationF1:MorpholinoAHR2Mo -0.030512   0.013415  -2.275   0.0234 *  
    ## GenerationF2:MorpholinoAHR2Mo  0.011703   0.013424   0.872   0.3838    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05641 on 413 degrees of freedom
    ## Multiple R-squared:  0.03911,    Adjusted R-squared:  0.02049 
    ## F-statistic: 2.101 on 8 and 413 DF,  p-value: 0.03462


    # best model was Generation + Exposure + Morpholino + Generation:Exposure + Generation:Morpholino
    # GenerationF1:ExposureBaP and GenerationF1:MorpholinoAHR2Mo were only significant ones

    pw.disper_bars <- ggplot(data = pw.Disper2, aes(x = Generation, y = log(dispers),
                                                    fill = Treatments)) +
      geom_boxplot(outlier.shape = NA) + 
      theme_classic() +
      ylab("log(Dispersion)") + xlab("") +
      geom_point(size = 3, shape = 21, alpha = 0.9,
                 position = position_jitterdodge(jitter.width = 0.07)) +
      labs(fill = "Treatments") +
      scale_fill_brewer(palette = "PuOr", direction = -1,
                        labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 20))

    pw.disper_bars

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/beta dispersion graph-1.png" width="98%" height="98%" />

How does treatment comparison change over generation?

Look at spread of variation of beta diversity distance
(intergenerational beta diversity analysis) using the median distance
between all points in F0, F1, or F2 cluster on ordination and the F0
centroid.

    # read in the script to calculate this metric
    source("/Users/alexieva/Documents/Projects/Writing/Research manuscripts/metagen_zfBaP/Pub_analysis_gitlab/scripts/treatment_beta_div_across_ALLgens.R")

    horiz_treatboxes2 <- ggplot(data = big_data, aes(x = Generation, 
                                                    y = median, 
                                                    color = Treatments,
                                                    group = Treatments)) +
        geom_errorbar(width = 0.1, 
                    aes(ymin = median-sd, ymax = median+sd),
                    position = position_dodge(0.1)) +
      geom_point(size = 3, position = position_dodge(0.1)) +
      geom_line(position = position_dodge(0.1)) +
      labs(y = "Median dissimilarity\nfrom F0 generation",
           x = "") +
      theme_classic() +
      scale_color_brewer(palette = "PuOr", direction = -1,
                        labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 20)) +
      scale_x_discrete(labels = c("F0 (Baseline)", "F1", "F2"))

    horiz_treatboxes2

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/treatment comparison change over generation-1.png" width="98%" height="98%" />

The following sections show beta diversity analysis on each individual
generation (i.e., subsetting out one generation and running
models/ordinations to see trends within each, indpendently of each
other). Much like above, for each generation, we are trying to see if
beta diversity metrics vary with treatment variables.

# F0

Is beta diversity and dispersion of the gut metagenome of F0 zebrafish
associated with exposure, morpholino, or other variables?

    meta_modF0 <- metaF0 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex")) # filter to only important variables we are testing

    F0mod0 <- capscale(pathabund_tab_F0 ~ 1, meta_modF0, distance = "bray")  # Model with intercept only
    F0mod1 <- capscale(pathabund_tab_F0 ~ . + Exposure:Morpholino + Morpholino:Exposure, meta_modF0, distance = "bray")  # Model with all explanatory variables
    pw_ordiF0 <- ordistep(F0mod0, scope = formula(F0mod1)) # this determines what the best model is to run RDA on
    pw_ordiF0

    # best model is pathabund_tab_F0 ~ Morpholino + Exposure + Morpholino:Exposure
    # model inertia was 1.5837; Inertia is scaled Chi-square
    # proportion variance explained is 0.2204

    # calculate RDA (distance based)
    pw_rdaF0 <- capscale(formula = pathabund_tab_F0 ~ Morpholino + Exposure + Morpholino:Exposure, data = meta_modF0, distance = "bray")
    # using bray because euclidean, bray, and jaccard were all correlated via mantel test, and bray is consistent with other analyses in NMDS and other graphs in this analysis
    RsquareAdj(pw_rdaF0)
    ## $r.squared
    ## [1] 0.1946419
    ## 
    ## $adj.r.squared
    ## [1] 0.1773843

    # 0.2036827 of variation explained by this model
    plot(pw_rdaF0)

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/dbRDA run with best model from ordistep F0-1.png" width="98%" height="98%" />

    smry_rdaF0 <- summary(pw_rdaF0)
    pwF0_PC1  <- data.frame(smry_rdaF0$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF0, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    # the biplot scores are the correlations between your environmental variables and axes
    pwF0_PC2  <- data.frame(smry_rdaF0$biplot)     # x and y coordinates for the vectors
    pwF0_PC1$Treatments <- factor(pwF0_PC1$Treatments, 
                                    levels = c("CoMo_DMSO",
                                               "CoMo_BaP", 
                                               "AHR2Mo_DMSO",
                                               "AHR2Mo_BaP"))

    pwrda_plotF0 <- ggplot(pwF0_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = pwF0_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = pwF0_PC2, 
                aes(x = CAP1, y = CAP2, label = c("AHR2 Morpholino", "BaP Exposure",
                                                  "Interaction")), 
                color = "black", size = 4, hjust = c(0.45, 0, 0.45), vjust = 1) +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF0$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF0$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "PuOr", direction = -1,
                        name = "Treatments", 
                        labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 20))

    pwrda_plotF0

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/graph ordination F0-1.png" width="98%" height="98%" />

    pw_dmF0 <- vegdist(pathabund_tab_F0, method = "bray")
    PermExpandMod_pwF0 <- adonis2(pw_dmF0 ~ Morpholino + Exposure + Morpholino:Exposure, data = metaF0)
    PermExpandMod_pwF0
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pw_dmF0 ~ Morpholino + Exposure + Morpholino:Exposure, data = metaF0)
    ##                      Df SumOfSqs      R2      F Pr(>F)    
    ## Morpholino            1  0.14303 0.09031 16.163  0.001 ***
    ## Exposure              1  0.09500 0.05999 10.736  0.001 ***
    ## Morpholino:Exposure   1  0.10680 0.06744 12.070  0.001 ***
    ## Residual            140  1.23883 0.78226                  
    ## Total               143  1.58365 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # morph > exp > morph:exposure

Evaluate beta dispersion at F0 as well.

    # calc dispersion
    pw.disperF0 <- betadisper(pw_dmF0, metaF0$Treatments)

    # make file for graphing and stats
    pwF0.Disper1 <- data.frame(pw.disperF0$distances)
    colnames(pwF0.Disper1) <- "dispers"
    pwF0.Disper2 <- pwF0.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF0, 
                 by = "SampleID")

    ## linear model approach
    # step AIC to test model
    testmod_bdispF0 <- lm(pw.disperF0$distances ~ Exposure + Morpholino + Exposure:Morpholino, data = pwF0.Disper2)
    AIC_bdispF0 <- stepAIC(testmod_bdispF0)
    ## Start:  AIC=-909.32
    ## pw.disperF0$distances ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Morpholino  1 0.0032879 0.24977 -909.41
    ## <none>                             0.24649 -909.32
    ## 
    ## Step:  AIC=-909.41
    ## pw.disperF0$distances ~ Exposure + Morpholino
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - Exposure    1 0.0028539 0.25263 -909.77
    ## <none>                    0.24977 -909.41
    ## - Morpholino  1 0.0068281 0.25660 -907.53
    ## 
    ## Step:  AIC=-909.77
    ## pw.disperF0$distances ~ Morpholino
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## <none>                    0.25263 -909.77
    ## - Morpholino  1 0.0068281 0.25945 -907.93

    # save in object
    lm_betadisperF0 <- summary(AIC_bdispF0)
    lm_betadisperF0
    ## 
    ## Call:
    ## lm(formula = pw.disperF0$distances ~ Morpholino, data = pwF0.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.06061 -0.03025 -0.01081  0.02018  0.22838 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.076272   0.004971  15.344   <2e-16 ***
    ## MorpholinoAHR2Mo 0.013772   0.007030   1.959   0.0521 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04218 on 142 degrees of freedom
    ## Multiple R-squared:  0.02632,    Adjusted R-squared:  0.01946 
    ## F-statistic: 3.838 on 1 and 142 DF,  p-value: 0.05206

    # Morpholino best model and MorpholinoAHR2Mo almost significant (nothing else)

    pw.disper_barsF0 <- ggplot(data = pwF0.Disper2, aes(x = Treatments, 
                                              y = dispers)) +
      geom_boxplot(aes(fill = Treatments), alpha = 0.6, 
                   outlier.shape = NA) + 
      theme_classic() +
      ylab("Dispersion") + xlab("") +
      geom_point(aes(fill = Treatments), size = 3, shape = 21, position = position_jitter(
        width = 0.25
      )) +
      labs(fill = "Treatments")

    pw.disper_barsF0

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F0 beta dispersion graph-1.png" width="98%" height="98%" />

# F1

Is beta diversity and dispersion of the gut metagenome of F1 zebrafish
associated with exposure, morpholino, or other variables?

    meta_modF1 <- metaF1 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex", "BMI")) # filter to only important variables we are testing

    F1mod0 <- capscale(pathabund_tab_F1 ~ 1, meta_modF1, distance = "bray")  # Model with intercept only
    F1mod1 <- capscale(pathabund_tab_F1 ~ . + Exposure:Morpholino + Morpholino:Exposure + Sex:BMI, meta_modF1, distance = "bray")  # Model with all explanatory variables
    pw_ordiF1 <- ordistep(F1mod0, scope = formula(F1mod1)) # this determines what the best model is to run RDA on
    pw_ordiF1

    # best model is pathabund_tab_F1 ~ Morpholino + Exposure
    # model inertia was 2.25210; Inertia is scaled Chi-square
    # proportion variance explained is 0.20189

    # calculate RDA (distance based)
    pw_rdaF1 <- capscale(formula = pathabund_tab_F1 ~ Morpholino + Exposure, data = meta_modF1, distance = "bray")
    # using bray because euclidean, bray, and jaccard were all correlated via mantel test, and bray is consistent with other analyses in NMDS and other graphs in this analysis
    RsquareAdj(pw_rdaF1)
    ## $r.squared
    ## [1] 0.1836533
    ## 
    ## $adj.r.squared
    ## [1] 0.1716482

    # 0.1901491 of variation explained by this model

    smry_rdaF1 <- summary(pw_rdaF1)
    pwF1_PC1  <- data.frame(smry_rdaF1$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF1, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    # the biplot scores are the correlations between your environmental variables and axes
    pwF1_PC2  <- data.frame(smry_rdaF1$biplot)     # x and y coordinates for the vectors
    pwrda_plotF1 <- ggplot(pwF1_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = pwF1_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "blue", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = pwF1_PC2, 
                aes(x = CAP1, y = CAP2, label = c(rownames(pwF1_PC2))), 
                color = "black") +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF1$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF1$cont$importance[2, "CAP2"], digits = 2),"%)"))

    # add centroids to ordination
    metaF1_treat <- dplyr::select(metaF1, c("SampleID", "Treatments")) # select only these for join in next step
    centDistrda_pwF1 <- as.data.frame(smry_rdaF1$sites) %>% # take only sites rda data
      dplyr::select(c("CAP1", "CAP2")) %>% # select axis points for rda ordination
      rownames_to_column("SampleID") %>% 
      inner_join(metaF1_treat, by = "SampleID") %>% # join with Treatment metadata
      group_by(Treatments) %>% # group Treatments
      summarize(CentroidCAP1 = mean(CAP1), # calculate mean of each treatment's axes
                CentroidCAP2 = mean(CAP2)) %>%
      ungroup() %>% # ungroup
      as.data.frame()
    # check if these are in the right place on rda
    pwrda_plot_centrF1 <- pwrda_plotF1 + geom_point(data = centDistrda_pwF1,
                              aes(x = CentroidCAP1, y = CentroidCAP2,
                                  fill = Treatments),
                              shape = 21, 
                              size = 3, 
                              color = "black",
                              stroke = 1)
    pwrda_plot_centrF1

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/graph ordination F1-1.png" width="98%" height="98%" />

    pw_dmF1 <- vegdist(pathabund_tab_F1, method = "bray")
    PermExpandMod_pwF1 <- adonis2(pw_dmF1 ~ Morpholino + Exposure, data = metaF1)
    PermExpandMod_pwF1
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pw_dmF1 ~ Morpholino + Exposure, data = metaF1)
    ##             Df SumOfSqs      R2       F Pr(>F)    
    ## Morpholino   1  0.37431 0.16621 28.2680  0.001 ***
    ## Exposure     1  0.07694 0.03416  5.8104  0.003 ** 
    ## Residual   136  1.80085 0.79963                   
    ## Total      138  2.25210 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # morph > exp

Evaluate beta dispersion for F1 fish.

    # calc dispersion
    pw.disperF1 <- betadisper(pw_dmF1, metaF1$Treatments)

    # make file for graphing
    pwF1.Disper1 <- data.frame(pw.disperF1$distances)
    colnames(pwF1.Disper1) <- "dispers"
    pwF1.Disper2 <- pwF1.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF1, 
                 by = "SampleID")

    ## linear model approach
    # step AIC to test model
    testmod_bdispF1 <- lm(pw.disperF1$distances ~ Exposure + Morpholino + Exposure:Morpholino, data = pwF1.Disper2)
    AIC_bdispF1 <- stepAIC(testmod_bdispF1)
    ## Start:  AIC=-722.3
    ## pw.disperF1$distances ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq     RSS     AIC
    ## - Exposure:Morpholino  1 0.0014433 0.72799 -724.02
    ## <none>                             0.72654 -722.30
    ## 
    ## Step:  AIC=-724.02
    ## pw.disperF1$distances ~ Exposure + Morpholino
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## - Exposure    1 0.0003786 0.72836 -725.95
    ## <none>                    0.72799 -724.02
    ## - Morpholino  1 0.0191851 0.74717 -722.41
    ## 
    ## Step:  AIC=-725.95
    ## pw.disperF1$distances ~ Morpholino
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## <none>                    0.72836 -725.95
    ## - Morpholino  1  0.019146 0.74751 -724.34

    # save in object
    lm_betadisperF1 <- summary(AIC_bdispF1)
    lm_betadisperF1
    ## 
    ## Call:
    ## lm(formula = pw.disperF1$distances ~ Morpholino, data = pwF1.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.06337 -0.03591 -0.01473  0.00674  0.57771 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       0.099326   0.008715  11.397   <2e-16 ***
    ## MorpholinoAHR2Mo -0.023473   0.012369  -1.898   0.0598 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07291 on 137 degrees of freedom
    ## Multiple R-squared:  0.02561,    Adjusted R-squared:  0.0185 
    ## F-statistic: 3.601 on 1 and 137 DF,  p-value: 0.05984

    # Morpholino best model and MorpholinoAHR2Mo almost significant (nothing else)

    pw.disper_barsF1 <- ggplot(data = pwF1.Disper2, aes(x = Treatments, 
                                              y = dispers)) +
      geom_boxplot(aes(fill = Treatments), alpha = 0.6, 
                   outlier.shape = NA) + 
      theme_classic() +
      ylab("Dispersion") + xlab("") +
      geom_point(aes(fill = Treatments), size = 3, shape = 21, position = position_jitter(
        width = 0.25
      )) +
      labs(fill = "Treatments")

    pw.disper_barsF1

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F1 beta dispersion graph-1.png" width="98%" height="98%" />

# F2

Is beta diversity and dispersion of the gut metagenome of F2 zebrafish
associated with exposure, morpholino, or other variables?

    meta_modF2 <- metaF2 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex", "BMI")) # filter to only important variables we are testing

    F2mod0 <- capscale(pathabund_tab_F2 ~ 1, meta_modF2, distance = "bray")  # Model with intercept only
    F2mod1 <- capscale(pathabund_tab_F2 ~ . + Exposure:Morpholino + Morpholino:Exposure + Sex:BMI, meta_modF2, distance = "bray")  # Model with all explanatory variables
    pw_ordiF2 <- ordistep(F2mod0, scope = formula(F2mod1)) # this determines what the best model is to run RDA on
    pw_ordiF2

    # best model is pathabund_tab_F2 ~ Morpholino + Exposure + Morpholino:Exposure
    # model inertia was 1.5020; Inertia is scaled Chi-square
    # proportion variance explained is 0.2722

    # calculate RDA (distance based)
    pw_rdaF2 <- capscale(formula = pathabund_tab_F2 ~ Morpholino + Exposure + Morpholino:Exposure, data = meta_modF2, distance = "bray")
    # using bray because euclidean, bray, and jaccard were all correlated via mantel test, and bray is consistent with other analyses in NMDS and other graphs in this analysis
    RsquareAdj(pw_rdaF2)
    ## $r.squared
    ## [1] 0.2424956
    ## 
    ## $adj.r.squared
    ## [1] 0.2256622

    # 0.25603 of variation explained by this model

    smry_rdaF2 <- summary(pw_rdaF2)
    pwF2_PC1  <- data.frame(smry_rdaF2$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF2, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    # the biplot scores are the correlations between your environmental variables and axes
    pwF2_PC2  <- data.frame(smry_rdaF2$biplot)     # x and y coordinates for the vectors
    pwrda_plotF2 <- ggplot(pwF2_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = pwF2_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "blue", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = pwF2_PC2, 
                aes(x = CAP1, y = CAP2, label = c(rownames(pwF2_PC2))), 
                color = "black") +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF2$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF2$cont$importance[2, "CAP2"], digits = 2),"%)"))

    # add centroids to ordination
    metaF2_treat <- dplyr::select(metaF2, c("SampleID", "Treatments")) # select only these for join in next step
    centDistrda_pwF2 <- as.data.frame(smry_rdaF2$sites) %>% # take only sites rda data
      dplyr::select(c("CAP1", "CAP2")) %>% # select axis points for rda ordination
      rownames_to_column("SampleID") %>% 
      inner_join(metaF2_treat, by = "SampleID") %>% # join with Treatment metadata
      group_by(Treatments) %>% # group Treatments
      summarize(CentroidCAP1 = mean(CAP1), # calculate mean of each treatment's axes
                CentroidCAP2 = mean(CAP2)) %>%
      ungroup() %>% # ungroup
      as.data.frame()
    # check if these are in the right place on rda
    pwrda_plot_centrF2 <- pwrda_plotF2 + geom_point(data = centDistrda_pwF2,
                              aes(x = CentroidCAP1, y = CentroidCAP2,
                                  fill = Treatments),
                              shape = 21, 
                              size = 3, 
                              color = "black",
                              stroke = 1)
    pwrda_plot_centrF2

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/graph ordination F2-1.png" width="98%" height="98%" />

    pw_dmF2 <- vegdist(pathabund_tab_F2, method = "bray")
    PermExpandMod_pwF2 <- adonis2(pw_dmF2 ~ Morpholino + Exposure + Morpholino:Exposure, data = metaF2)
    PermExpandMod_pwF2
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pw_dmF2 ~ Morpholino + Exposure + Morpholino:Exposure, data = metaF2)
    ##                      Df SumOfSqs     R2       F Pr(>F)    
    ## Morpholino            1  0.32008 0.2131 39.3762  0.001 ***
    ## Exposure              1  0.05017 0.0334  6.1716  0.001 ***
    ## Morpholino:Exposure   1  0.03439 0.0229  4.2313  0.001 ***
    ## Residual            135  1.09737 0.7306                   
    ## Total               138  1.50201 1.0000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # morph > exp > interaction

Evaluate beta dispersion for F2 fish.

    # calc dispersion
    pw.disperF2 <- betadisper(pw_dmF2, metaF2$Treatments)

    # make file for graphing
    pwF2.Disper1 <- data.frame(pw.disperF2$distances)
    colnames(pwF2.Disper1) <- "dispers"
    pwF2.Disper2 <- pwF2.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF2, 
                 by = "SampleID")

    ## linear model approach
    # step AIC to test model
    testmod_bdispF2 <- lm(pw.disperF2$distances ~ Exposure + Morpholino + Exposure:Morpholino, data = pwF2.Disper2)
    AIC_bdispF2 <- stepAIC(testmod_bdispF2)
    ## Start:  AIC=-878.94
    ## pw.disperF2$distances ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df  Sum of Sq     RSS     AIC
    ## - Exposure:Morpholino  1 5.2404e-05 0.23547 -880.91
    ## <none>                              0.23541 -878.94
    ## 
    ## Step:  AIC=-880.91
    ## pw.disperF2$distances ~ Exposure + Morpholino
    ## 
    ##              Df Sum of Sq     RSS     AIC
    ## <none>                    0.23547 -880.91
    ## - Exposure    1 0.0098204 0.24529 -877.23
    ## - Morpholino  1 0.0111999 0.24667 -876.45

    # save in object
    lm_betadisperF2 <- summary(AIC_bdispF2)
    lm_betadisperF2
    ## 
    ## Call:
    ## lm(formula = pw.disperF2$distances ~ Exposure + Morpholino, data = pwF2.Disper2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.060313 -0.026429 -0.008092  0.013833  0.192168 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.060822   0.006345   9.585   <2e-16 ***
    ## ExposureBaP      0.016833   0.007068   2.382   0.0186 *  
    ## MorpholinoAHR2Mo 0.017977   0.007068   2.543   0.0121 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04161 on 136 degrees of freedom
    ## Multiple R-squared:  0.07925,    Adjusted R-squared:  0.06571 
    ## F-statistic: 5.853 on 2 and 136 DF,  p-value: 0.003645

    # Exposure + Morpholino best model and MorpholinoAHR2Mo almost significant (nothing else)

    pw.disper_barsF2 <- ggplot(data = pwF2.Disper2, aes(x = Treatments, 
                                              y = dispers)) +
      geom_boxplot(aes(fill = Treatments), alpha = 0.6, 
                   outlier.shape = NA) + 
      theme_classic() +
      ylab("Dispersion") + xlab("") +
      geom_point(aes(fill = Treatments), size = 3, shape = 21, position = position_jitter(
        width = 0.25
      )) +
      labs(fill = "Treatments")
    pw.disper_barsF2

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/F2 beta dispersion graph-1.png" width="98%" height="98%" />

## Various useful figures

Paneled beta dispersion figure of each generation’s dispersion (most
were not significant or were very low impact on trends)

    pw.disper_barsF0 <- pw.disper_barsF0 + 
      scale_y_continuous(limits = c(0, 0.7))
    pw.disper_barsF1 <- pw.disper_barsF1 + 
      scale_y_continuous(limits = c(0, 0.7))
    pw.disper_barsF2 <- pw.disper_barsF2 + 
      scale_y_continuous(limits = c(0, 0.7))
    pw_disperbars_bygen <- ggarrange(pw.disper_barsF0, pw.disper_barsF1, pw.disper_barsF2, ncol = 3, nrow = 1)
    pw_disperbars_bygen

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/paneled beta dispersion fig-1.png" width="98%" height="98%" />

dbRDA of each generation in a paneled supplemental figure, based on the
optimal models for each subset generation.

    ########################### F0 ########################### 
    # pre plot data frame
    smry_rdaF0 <- summary(capscale(formula = pathabund_tab_F0 ~ Morpholino + Exposure + Morpholino:Exposure, 
                                   data = meta_modF0,
                                   distance = "bray"))
    F0_PC1  <- data.frame(smry_rdaF0$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF0, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    F0_PC2  <- data.frame(smry_rdaF0$biplot)     # x and y coordinates for the vectors
    # reorder x-axis
    F0_PC1$Treatments <- factor(F0_PC1$Treatments, 
                            levels = c("CoMo_DMSO",
                                       "CoMo_BaP", 
                                       "AHR2Mo_DMSO",
                                       "AHR2Mo_BaP"))

    rda_plotF0 <- ggplot(F0_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = F0_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = F0_PC2, 
                aes(x = CAP1, y = CAP2, label = c("AHR2 Morpholino", "BaP Exposure", "Interaction")), 
                color = "black", size = 4, hjust = 0.5, vjust = 0.5) +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF0$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF0$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "PuOr", direction = -1,
                         name = "Treatments",
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      theme(text = element_text(size = 20))

    rda_plotF0

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/manuscript suppl figure of each individual generation rda-1.png" width="98%" height="98%" />



    ########################### F1 ########################### 
    # pre plot data frame
    smry_rdaF1 <- summary(capscale(formula = pathabund_tab_F1 ~ Morpholino + Exposure, 
                                   data = meta_modF1,
                                   distance = "bray"))
    F1_PC1  <- data.frame(smry_rdaF1$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF1, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    F1_PC2  <- data.frame(smry_rdaF1$biplot)     # x and y coordinates for the vectors
    # reorder x-axis
    F1_PC1$Treatments <- factor(F1_PC1$Treatments, 
                            levels = c("CoMo_DMSO",
                                       "CoMo_BaP", 
                                       "AHR2Mo_DMSO",
                                       "AHR2Mo_BaP"))

    rda_plotF1 <- ggplot(F1_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = F1_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = F1_PC2, 
                aes(x = CAP1, y = CAP2, label = c("AHR2 Morpholino", "BaP Exposure")), 
                color = "black", size = 4, hjust = 1, vjust = 0.5) +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF1$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF1$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "PuOr", direction = -1,
                         name = "Treatments",
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      scale_shape_discrete(labels = c("Female", "Male")) +
      theme(text = element_text(size = 20))

    rda_plotF1

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/manuscript suppl figure of each individual generation rda-2.png" width="98%" height="98%" />



    ########################### F2 ########################### 
    # pre plot data frame
    smry_rdaF2 <- summary(capscale(formula = pathabund_tab_F2 ~ Morpholino + Exposure + Morpholino:Exposure, 
                                   data = meta_modF2,
                                   distance = "bray"))
    F2_PC1  <- data.frame(smry_rdaF2$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF2, c(SampleID, Treatments, Exposure, Morpholino)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    F2_PC2  <- data.frame(smry_rdaF2$biplot)     # x and y coordinates for the vectors
    # reorder x-axis
    F2_PC1$Treatments <- factor(F2_PC1$Treatments, 
                            levels = c("CoMo_DMSO",
                                       "CoMo_BaP", 
                                       "AHR2Mo_DMSO",
                                       "AHR2Mo_BaP"))

    rda_plotF2 <- ggplot(F2_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = F2_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = F2_PC2, 
                aes(x = CAP1, y = CAP2, label = c("AHR2 Morpholino", "BaP Exposure", "Interaction")), 
                color = "black", size = 4, hjust = c(0.5, 0, 0), vjust = 0.5) +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF2$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF2$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "PuOr", direction = -1,
                         name = "Treatments",
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      scale_shape_discrete(labels = c("Female", "Male")) +
      theme(text = element_text(size = 20))

    rda_plotF2

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/manuscript suppl figure of each individual generation rda-3.png" width="98%" height="98%" />

Heatmap of all stats - behavior, alpha, and beta diversity - per
treatment

    # read in data and clean up for graphing
    sumstats_beh <- read.csv(paste0(input_dir, "SumStatsHeat_behavior.csv"),
                             header = T) %>%
      mutate(Treatment = replace(Treatment, Treatment == "CoMo _BaP", "AHR2Mo - / BaP +"),
             Treatment = replace(Treatment, Treatment == "AHRMo_DMSO", "AHR2Mo + / BaP -"),
             Treatment = replace(Treatment, Treatment == "AHRMo_BaP", "AHR2Mo + / BaP +"),
             Metric = replace(Metric, Metric == "freeswim", "Freeswim"),
             Metric = replace(Metric, Metric == "shoaling_iid", "Shoaling iid"),
             Metric = replace(Metric, Metric == "shoaling_nnd", "Shoaling nnd"),
             Metric = replace(Metric, Metric == "shoaling_speed", "Shoaling speed"))
    sumstats_beh$Metric <- factor(x = sumstats_beh$Metric,
                                 levels = c("Freeswim",
                                            "Shoaling nnd",
                                            "Shoaling iid",
                                            "Shoaling speed"))

    sumstats_alphadiv <- read.csv(paste0(input_dir, "SumStatsHeat_alphadiv.csv"),
                             header = T) %>%
      mutate(Treatment = replace(Treatment, Treatment == "CoMo _BaP", "AHR2Mo - / BaP +"),
             Treatment = replace(Treatment, Treatment == "AHRMo_DMSO", "AHR2Mo + / BaP -"),
             Treatment = replace(Treatment, Treatment == "AHRMo_BaP", "AHR2Mo + / BaP +"),
             Metric = replace(Metric, Metric == "Richness_16S", "Taxonomic Richness"),
             Metric = replace(Metric, Metric == "Shannon_16S", "Taxonomic Shannon Diversity"),
             Metric = replace(Metric, Metric == "Richness_MG", "Pathway Richness"),
             Metric = replace(Metric, Metric == "Shannon_MG", "Pathway Shannon Diversity"))
    sumstats_alphadiv$Metric <- factor(x = sumstats_alphadiv$Metric,
                                 levels = c("Taxonomic Richness",
                                            "Taxonomic Shannon Diversity",
                                            "Pathway Richness",
                                            "Pathway Shannon Diversity"))

    sumstats_betadiv <- read.csv(paste0(input_dir, "SumStatsHeat_betadiv.csv"),
                             header = T) %>%
      mutate(Treatment = replace(Treatment, Treatment == "CoMo _BaP", "AHR2Mo - / BaP +"),
             Treatment = replace(Treatment, Treatment == "AHRMo_DMSO", "AHR2Mo + / BaP -"),
             Treatment = replace(Treatment, Treatment == "AHRMo_BaP", "AHR2Mo + / BaP +"),
             Metric = replace(Metric, Metric == "BC_16S", "Taxonomic Bray-Curtis Dissimilarity"),
             Metric = replace(Metric, Metric == "BC_MG", "Pathway Bray-Curtis Dissimilarity"))
    sumstats_betadiv$Metric <- factor(x = sumstats_betadiv$Metric,
                                 levels = c("Taxonomic Bray-Curtis Dissimilarity",
                                            "Pathway Bray-Curtis Dissimilarity"))

    # make a heatmap per measure
    hm_stats_beh <- ggplot(data = sumstats_beh, aes(x = Metric, y = Treatment)) +
      geom_tile(aes(fill = Statistic), colour = "white", linetype = 1) +
      theme_classic() +
      scale_fill_gradient2(low = "white", high = "darkgreen", 
                           na.value = "white")  + 
      labs(y = "Treatment compared to Control", x = "", fill = "KS test\nStatistic") +
      geom_text(aes(label = Significance)) + 
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y.right = element_text(angle = 0),
            text = element_text(size = 23),
            plot.margin = unit(c(0,0.2,0,1), 'lines')) +
      facet_grid("Generation") +
      ggtitle("Behavioral metrics")

    hm_stats_alpha <- ggplot(data = sumstats_alphadiv, aes(x = Metric, y = Treatment)) +
      geom_tile(aes(fill = Estimate), colour = "white", linetype = 1) +
      theme_classic() +
      scale_fill_gradient2(low = "tan", mid = "white", high = "darkgreen", 
                           midpoint = 0, na.value = "white")  + 
      labs(y = "Treatment compared to Control", x = "", fill = "Linear Model\nEstimate") +
      geom_text(aes(label = Significance)) + 
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y.right = element_text(angle = 0),
            text = element_text(size = 23),
            plot.margin = unit(c(0,0.2,0,1), 'lines')) +
      facet_grid("Generation") +
      ggtitle("Alpha diversity metrics")

    hm_stats_beta <- ggplot(data = sumstats_betadiv, aes(x = Metric, y = Treatment)) +
      geom_tile(aes(fill = Estimate), colour = "white", linetype = 1) +
      theme_classic() +
      scale_fill_gradient2(low = "tan", mid = "white", high = "darkgreen", 
                           midpoint = 0, na.value = "white")  + 
      labs(y = "Treatment compared to Control", x = "", fill = "Linear Model\nEstimate") +
      geom_text(aes(label = Significance)) + 
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y.right = element_text(angle = 0),
            text = element_text(size = 23),
            plot.margin = unit(c(0,0.2,0,1), 'lines')) +
      facet_grid("Generation") +
      ggtitle("Beta diversity metrics")

    # combine into a paneled fig and save
    hm_sumstates <- ggarrange(hm_stats_beh, hm_stats_alpha, hm_stats_beta, 
                              ncol = 3, nrow = 1,
                              align = "hv",
                              labels = c("A", "B", "C"))
    hm_sumstates

<img src="zfBaP_PubAnalysis_files/figure-markdown_strict/paneled heatmap of all stats-1.png" width="98%" height="98%" />
