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
    home_dir <- "/Users/alexieva/Documents/Projects/Writing/Research manuscripts/metagen_zfBaP"
    input_dir <- "/Users/alexieva/Documents/Projects/Writing/Research manuscripts/metagen_zfBaP/Pub_analysis_gitlab/input_files_pub/"
    output_dir <- "/Users/alexieva/Documents/Projects/Analysis/metagen_zfBaP/16S_Ebony/"

## load in 16S data from Ebony

-   this data has been pre-processed via the DADA2 pipeline and Ebony
    quality filtered out singletons, mitochondria, and chloroplasts

<!-- -->

    ## whole dataset (all generations) files

    # read in 16S data object for whole dataset
    # convert to data frame from phyloseq object Ebony gave me
    taxa_psobj <- readRDS(paste0(input_dir, "ps.qc.rds"))
    taxa_otutab <- as.data.frame(as.matrix(taxa_psobj@otu_table))
    # metadata for whole dataset
    metaALL <- as.data.frame(as.matrix(taxa_psobj@sam_data)) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), 
             Morph = relevel(factor(Morph), "CoMo"),
             Treatments = relevel(factor(Treatments), "CoMo+DMSO"),
             Generation = factor(Generation),
             Sex = factor(Sex)) %>% # relevel and factor so these are accurate
      rename_at("Morph", ~"Morpholino")
    metaALL$Treatments <- gsub("\\+", "_", metaALL$Treatments)

    ## F0
    F0_meta <- dplyr::filter(metaALL, Generation == "F0")
    F0_taxa_psobj <- readRDS(paste0(home_dir, "ps.qc.f0.rds"))
    F0_taxa_otutab <- as.data.frame(as.matrix(F0_taxa_psobj@otu_table))
    metaF0 <- as.data.frame(as.matrix(F0_taxa_psobj@sam_data)) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), 
             Morph = relevel(factor(Morph), "CoMo"),
             Treatments = relevel(factor(Treatments), "CoMo+DMSO"),
             Generation = factor(Generation),
             Sex = factor(Sex)) %>% # relevel and factor so these are accurate
      rename_at("Morph", ~"Morpholino")
    metaF0$Treatments <- gsub("\\+", "_", metaF0$Treatments)

    ## F1
    F1_meta <- dplyr::filter(metaALL, Generation == "F1")
    F1_taxa_psobj <- readRDS(paste0(home_dir, "ps.qc.f1.rds"))
    F1_taxa_otutab <- as.data.frame(as.matrix(F1_taxa_psobj@otu_table))
    metaF1 <- as.data.frame(as.matrix(F1_taxa_psobj@sam_data)) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), 
             Morph = relevel(factor(Morph), "CoMo"),
             Treatments = relevel(factor(Treatments), "CoMo+DMSO"),
             Generation = factor(Generation),
             Sex = factor(Sex)) %>% # relevel and factor so these are accurate
      rename_at("Morph", ~"Morpholino")
    metaF1$Treatments <- gsub("\\+", "_", metaF1$Treatments)

    ## F2
    F2_meta <- dplyr::filter(metaALL, Generation == "F2")
    F2_taxa_psobj <- readRDS(paste0(home_dir, "ps.qc.f2.rds"))
    F2_taxa_otutab <- as.data.frame(as.matrix(F2_taxa_psobj@otu_table))
    metaF2 <- as.data.frame(as.matrix(F2_taxa_psobj@sam_data)) %>%
      mutate(Exposure = relevel(factor(Exposure), "DMSO"), 
             Morph = relevel(factor(Morph), "CoMo"),
             Treatments = relevel(factor(Treatments), "CoMo+DMSO"),
             Generation = factor(Generation),
             Sex = factor(Sex)) %>% # relevel and factor so these are accurate
      rename_at("Morph", ~"Morpholino")
    metaF2$Treatments <- gsub("\\+", "_", metaF2$Treatments)

## Alpha diversity

# Whole dataset

Can gut microbiome taxonomic alpha diverity be predicted by covariates,
in particular generation, exposure, and morpholino treatment?

    tax_richness <- specnumber(taxa_otutab)
    tax_shannon <- diversity(taxa_otutab, "shannon", base = exp(1))
    tax_simpson <- diversity(taxa_otutab, "simpson")
    alphadiv <- data.frame(tax_richness, tax_shannon, tax_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(metaALL, by = "SampleID") %>%
      dplyr::select(c("SampleID", "tax_richness", "tax_shannon", "tax_simpson", "Exposure", "Morpholino", "Generation", "Sex", "Treatments")) # filter to only important variables we are testing
    alphadiv
    ##         SampleID tax_richness tax_shannon tax_simpson Exposure Morpholino
    ## 1       BaP_0001           18   1.9043087   0.8095260      BaP     AHR2Mo
    ## 2       BaP_0002           18   1.8046092   0.7757695      BaP     AHR2Mo
    ## 3       BaP_0003           11   1.4751795   0.7105426     DMSO     AHR2Mo
    ## 4       BaP_0004           25   1.2613050   0.5446668     DMSO     AHR2Mo
    ## 5       BaP_0005           17   1.5067301   0.6663175     DMSO     AHR2Mo
    ## 6       BaP_0006           13   0.6188942   0.2594314     DMSO     AHR2Mo
    ## 7       BaP_0007           21   1.5170818   0.6999072     DMSO     AHR2Mo
    ## 8       BaP_0008           19   1.6607160   0.7340606     DMSO     AHR2Mo
    ## 9       BaP_0009           25   1.0537078   0.4455745     DMSO       CoMo
    ## 10      BaP_0010           40   1.8399199   0.7588862     DMSO       CoMo
    ## 11      BaP_0011           22   1.7670211   0.7529178     DMSO       CoMo
    ## 12      BaP_0012           22   1.8079990   0.8021222     DMSO       CoMo
    ## 13      BaP_0013           14   1.2095781   0.5412487     DMSO       CoMo
    ## 14      BaP_0014           86   1.9706742   0.7003931     DMSO       CoMo
    ## 15      BaP_0015           19   1.1628763   0.5890970     DMSO       CoMo
    ## 16      BaP_0016           39   1.8126558   0.7489156     DMSO       CoMo
    ## 17      BaP_0017           83   2.6841609   0.8648838      BaP       CoMo
    ## 18      BaP_0018           21   1.5068673   0.6794985      BaP       CoMo
    ## 19      BaP_0019           10   1.0512185   0.5314333     DMSO     AHR2Mo
    ## 20      BaP_0020           13   1.6128215   0.7261125     DMSO     AHR2Mo
    ## 21      BaP_0021           13   1.3914185   0.6307450     DMSO     AHR2Mo
    ## 22      BaP_0022           11   1.4019605   0.6706061     DMSO     AHR2Mo
    ## 23      BaP_0023           11   1.2529471   0.5496405     DMSO     AHR2Mo
    ## 24      BaP_0024           13   1.4084568   0.6584834     DMSO     AHR2Mo
    ## 25      BaP_0025           17   1.5781985   0.7269485     DMSO       CoMo
    ## 26      BaP_0026           12   1.3812271   0.6463410     DMSO       CoMo
    ## 27      BaP_0027           14   1.1236284   0.5495108     DMSO       CoMo
    ## 28      BaP_0028           11   1.2224680   0.6196388     DMSO       CoMo
    ## 29      BaP_0029           11   0.6422178   0.2744240     DMSO       CoMo
    ## 30      BaP_0030           11   1.3441607   0.6149830     DMSO       CoMo
    ## 31      BaP_0031           11   1.5601041   0.7293716     DMSO       CoMo
    ## 32      BaP_0032           10   1.5125483   0.6921326     DMSO       CoMo
    ## 33      BaP_0033           14   1.5559083   0.7273313     DMSO       CoMo
    ## 34      BaP_0034           11   1.5570319   0.7292421     DMSO       CoMo
    ## 35      BaP_0035           15   1.6449805   0.7551149     DMSO       CoMo
    ## 36      BaP_0036           13   1.4391952   0.6708733     DMSO       CoMo
    ## 37      BaP_0037           16   1.2757459   0.5622978     DMSO     AHR2Mo
    ## 38      BaP_0038           15   1.4314688   0.6539241     DMSO     AHR2Mo
    ## 39      BaP_0039           18   1.4227330   0.6235380     DMSO     AHR2Mo
    ## 40      BaP_0040           15   1.3540076   0.6405306     DMSO     AHR2Mo
    ## 41      BaP_0041           22   1.6076996   0.6927714     DMSO     AHR2Mo
    ## 42      BaP_0042           19   1.6392542   0.6932340     DMSO     AHR2Mo
    ## 43      BaP_0043           23   1.8047307   0.7514470     DMSO     AHR2Mo
    ## 44      BaP_0044           19   1.4213272   0.6450931     DMSO     AHR2Mo
    ## 45      BaP_0045           16   1.2819910   0.5927275     DMSO       CoMo
    ## 46      BaP_0046           25   1.7335938   0.7183714     DMSO       CoMo
    ## 47      BaP_0047           18   1.6367647   0.6794797     DMSO       CoMo
    ## 48      BaP_0048           20   1.5739728   0.6595144     DMSO       CoMo
    ## 49      BaP_0049           24   1.8541737   0.7496501     DMSO       CoMo
    ## 50      BaP_0050           30   1.7450228   0.7304204     DMSO       CoMo
    ## 51      BaP_0051           24   1.8448173   0.8018402     DMSO       CoMo
    ## 52      BaP_0052           23   1.4070826   0.6576872     DMSO       CoMo
    ## 53      BaP_0053           19   1.7323135   0.7368580     DMSO     AHR2Mo
    ## 54      BaP_0054           31   1.7907850   0.7494194     DMSO     AHR2Mo
    ## 55      BaP_0055           11   1.4899875   0.7075503     DMSO     AHR2Mo
    ## 56      BaP_0056           10   0.7240748   0.3372249     DMSO     AHR2Mo
    ## 57      BaP_0057           10   1.3983658   0.7053094     DMSO     AHR2Mo
    ## 58      BaP_0058           12   1.4607582   0.6859999     DMSO     AHR2Mo
    ## 59      BaP_0059           10   1.2163893   0.6325380     DMSO     AHR2Mo
    ## 60      BaP_0060           15   1.1607740   0.4850843     DMSO     AHR2Mo
    ## 61      BaP_0061           11   0.9883705   0.4652649     DMSO       CoMo
    ## 62      BaP_0062           16   1.0280139   0.4284470     DMSO       CoMo
    ## 63      BaP_0063            8   1.3655938   0.6735804     DMSO       CoMo
    ## 64      BaP_0064           14   1.1532284   0.5001923     DMSO       CoMo
    ## 65      BaP_0065           15   1.1261625   0.5135362     DMSO       CoMo
    ## 66      BaP_0066           17   0.9325314   0.3812879     DMSO       CoMo
    ## 67      BaP_0067           10   1.4807195   0.7117886     DMSO       CoMo
    ## 68      BaP_0068           19   1.0989520   0.4611377     DMSO       CoMo
    ## 69      BaP_0069           14   1.3220686   0.5840625     DMSO       CoMo
    ## 70      BaP_0070           13   1.6237753   0.7550168     DMSO       CoMo
    ## 71      BaP_0071            9   1.1621968   0.5928558     DMSO       CoMo
    ## 72      BaP_0072           14   1.5510806   0.6747229      BaP       CoMo
    ## 73      BaP_0073           17   0.9210665   0.4118116     DMSO     AHR2Mo
    ## 74      BaP_0074           16   0.5077373   0.1991810     DMSO     AHR2Mo
    ## 75      BaP_0075           11   0.8843224   0.3796460      BaP     AHR2Mo
    ## 76      BaP_0076           13   1.1102764   0.5234061      BaP     AHR2Mo
    ## 77      BaP_0077           12   1.2423913   0.6245837      BaP     AHR2Mo
    ## 78      BaP_0078           15   1.4341145   0.6314839      BaP     AHR2Mo
    ## 79      BaP_0079           13   1.2787306   0.6323861      BaP     AHR2Mo
    ## 80      BaP_0080           13   0.5900458   0.2553297      BaP     AHR2Mo
    ## 81      BaP_0081           15   0.9098571   0.4128258     DMSO       CoMo
    ## 82      BaP_0082           34   2.2943894   0.8527067     DMSO       CoMo
    ## 83      BaP_0083           17   0.7864446   0.3316524     DMSO       CoMo
    ## 84      BaP_0084            9   0.2282591   0.0790422     DMSO       CoMo
    ## 85      BaP_0085           22   0.6780049   0.2589378     DMSO       CoMo
    ## 86      BaP_0086           22   1.4914075   0.6900570     DMSO       CoMo
    ## 87      BaP_0087           22   0.6362403   0.2464652     DMSO       CoMo
    ## 88      BaP_0088           19   0.5877726   0.2230730     DMSO       CoMo
    ## 89      BaP_0089           16   1.4845507   0.7092555      BaP     AHR2Mo
    ## 90      BaP_0090           26   1.0614801   0.4594206      BaP     AHR2Mo
    ## 91      BaP_0091           14   1.6387108   0.7453882     DMSO     AHR2Mo
    ## 92      BaP_0092           19   1.1598175   0.4685972     DMSO     AHR2Mo
    ## 93      BaP_0093           12   1.5393032   0.6904846     DMSO     AHR2Mo
    ## 94      BaP_0094           15   0.9172201   0.3726186     DMSO     AHR2Mo
    ## 95      BaP_0095           15   1.4679705   0.6229566     DMSO     AHR2Mo
    ## 96      BaP_0096           14   1.2951597   0.5676459     DMSO     AHR2Mo
    ## 97      BaP_0097           11   1.7900889   0.8056914      BaP       CoMo
    ## 98      BaP_0098           19   1.1648889   0.4808637      BaP     AHR2Mo
    ## 99      BaP_0099           14   1.5210748   0.6909720      BaP     AHR2Mo
    ## 100     BaP_0100           21   1.4759416   0.6374041      BaP     AHR2Mo
    ## 101     BaP_0101           17   1.6343917   0.7146165      BaP     AHR2Mo
    ## 102     BaP_0102           15   1.3618935   0.6504646      BaP     AHR2Mo
    ## 103     BaP_0103            7   0.2844285   0.1059817      BaP     AHR2Mo
    ## 104     BaP_0104           15   1.1505353   0.5055484      BaP     AHR2Mo
    ## 105     BaP_0105           17   1.2686579   0.5789640      BaP       CoMo
    ## 106     BaP_0106           22   1.4529187   0.6119164      BaP       CoMo
    ## 107     BaP_0107           19   1.5296770   0.6947166      BaP       CoMo
    ## 108     BaP_0108           22   1.2960939   0.5209255      BaP       CoMo
    ## 109     BaP_0109           19   1.5707789   0.6712854      BaP     AHR2Mo
    ## 110     BaP_0110           18   1.8693331   0.7897186      BaP     AHR2Mo
    ## 111     BaP_0111           17   2.2136365   0.8664748      BaP     AHR2Mo
    ## 112     BaP_0112           16   1.4751848   0.6344341      BaP     AHR2Mo
    ## 113     BaP_0113           20   1.7883009   0.7514739      BaP     AHR2Mo
    ## 114     BaP_0114           19   1.2977889   0.5405955      BaP     AHR2Mo
    ## 115     BaP_0115           17   1.7658332   0.7512136      BaP     AHR2Mo
    ## 116     BaP_0116           20   2.0178751   0.8165626      BaP     AHR2Mo
    ## 117     BaP_0117           17   1.9804538   0.8153927     DMSO       CoMo
    ## 118     BaP_0118           19   2.0619582   0.8316066     DMSO       CoMo
    ## 119     BaP_0119           20   1.7785959   0.7649543     DMSO       CoMo
    ## 120     BaP_0120           40   2.7443870   0.9055234     DMSO       CoMo
    ## 121     BaP_0121           17   1.6836094   0.7142154      BaP       CoMo
    ## 122     BaP_0122           14   1.4623184   0.6786593      BaP       CoMo
    ## 123     BaP_0123           20   1.9103262   0.7953088      BaP       CoMo
    ## 124     BaP_0124           16   1.5024854   0.6395335      BaP       CoMo
    ## 125     BaP_0125           18   1.7705105   0.7454335     DMSO     AHR2Mo
    ## 126     BaP_0126           22   2.0481538   0.8281015     DMSO     AHR2Mo
    ## 127     BaP_0127           12   1.5212087   0.6870889     DMSO     AHR2Mo
    ## 128     BaP_0128           13   1.7106230   0.7770207     DMSO     AHR2Mo
    ## 129     BaP_0129            9   1.3491707   0.6300909     DMSO     AHR2Mo
    ## 130     BaP_0130            9   1.7068443   0.7746385     DMSO     AHR2Mo
    ## 131     BaP_0131           11   1.4665259   0.6641323     DMSO     AHR2Mo
    ## 132     BaP_0132           10   1.3414603   0.6168983     DMSO     AHR2Mo
    ## 133     BaP_0133           15   1.5557641   0.7035133      BaP     AHR2Mo
    ## 134     BaP_0134           10   1.7733486   0.8046963      BaP     AHR2Mo
    ## 135     BaP_0135           17   1.6612666   0.7253661      BaP     AHR2Mo
    ## 136     BaP_0136           14   1.7117134   0.7833875      BaP     AHR2Mo
    ## 137     BaP_0137           13   1.7029994   0.7734793     DMSO       CoMo
    ## 138     BaP_0138           14   1.6222456   0.7569204     DMSO       CoMo
    ## 139     BaP_0139           11   1.3733754   0.6156915     DMSO       CoMo
    ## 140     BaP_0140           11   1.5254011   0.7328212     DMSO       CoMo
    ## 141     BaP_0141           12   1.3943281   0.6546422      BaP       CoMo
    ## 142     BaP_0142           10   1.5592462   0.7320276      BaP     AHR2Mo
    ## 143     BaP_0143           13   1.6226867   0.7373885      BaP       CoMo
    ## 144     BaP_0144           14   1.7932130   0.7779819      BaP       CoMo
    ## 145     BaP_0433           11   1.5400814   0.6786236      BaP     AHR2Mo
    ## 146     BaP_0434           13   1.4016651   0.6213876      BaP     AHR2Mo
    ## 147     BaP_0435           12   1.7230774   0.7631621      BaP     AHR2Mo
    ## 148     BaP_0436           19   1.3042965   0.5231613      BaP     AHR2Mo
    ## 149     BaP_0437           13   1.7245850   0.7585671      BaP     AHR2Mo
    ## 150     BaP_0438           19   1.5794031   0.6510486      BaP     AHR2Mo
    ## 151     BaP_0439           18   1.8483205   0.7443729      BaP     AHR2Mo
    ## 152     BaP_0440           15   1.3145918   0.5503501      BaP     AHR2Mo
    ## 153     BaP_0441           15   1.7411699   0.7473252      BaP     AHR2Mo
    ## 154     BaP_0442           13   1.1033586   0.5466737      BaP     AHR2Mo
    ## 155     BaP_0443           13   1.5628293   0.6611785      BaP       CoMo
    ## 156     BaP_0444           12   1.3807386   0.6552241      BaP       CoMo
    ## 157     BaP_0445           15   1.1741261   0.5816562      BaP       CoMo
    ## 158     BaP_0446           16   1.7321591   0.7566231      BaP       CoMo
    ## 159     BaP_0447           13   1.2836558   0.5385074      BaP       CoMo
    ## 160     BaP_0448           14   1.8982637   0.7981676      BaP       CoMo
    ## 161     BaP_0449           24   0.8019297   0.3515296     DMSO     AHR2Mo
    ## 162     BaP_0450           14   0.8658108   0.3568893     DMSO     AHR2Mo
    ## 163     BaP_0451           17   1.2850795   0.5295749     DMSO     AHR2Mo
    ## 164     BaP_0452           18   1.7698049   0.7464540     DMSO     AHR2Mo
    ## 165     BaP_0453           14   1.5076940   0.6795778     DMSO     AHR2Mo
    ## 166     BaP_0454           13   1.1511536   0.5087171     DMSO     AHR2Mo
    ## 167     BaP_0455           15   1.8584973   0.7948014      BaP     AHR2Mo
    ## 168     BaP_0456           21   1.8036572   0.7700729      BaP     AHR2Mo
    ## 169     BaP_0457           20   1.3619889   0.6399455      BaP       CoMo
    ## 170     BaP_0458           11   1.2945837   0.6110593     DMSO     AHR2Mo
    ## 171     BaP_0459           16   1.8863330   0.7932580     DMSO     AHR2Mo
    ## 172     BaP_0460           16   1.9276707   0.8214444     DMSO     AHR2Mo
    ## 173     BaP_0461           14   1.8685375   0.8079892     DMSO     AHR2Mo
    ## 174     BaP_0462           15   1.8874619   0.7996027     DMSO     AHR2Mo
    ## 175     BaP_0463           17   1.7603842   0.7728881     DMSO     AHR2Mo
    ## 176     BaP_0464           16   1.6262117   0.7112714     DMSO     AHR2Mo
    ## 177     BaP_0465           16   1.4686872   0.6782274     DMSO       CoMo
    ## 178     BaP_0466           15   1.7159484   0.7398164     DMSO       CoMo
    ## 179     BaP_0467           13   1.2291729   0.5096855     DMSO       CoMo
    ## 180     BaP_0468           18   1.7159057   0.7246432     DMSO       CoMo
    ## 181     BaP_0541           12   0.8702130   0.3836692      BaP     AHR2Mo
    ## 182     BaP_0542           14   1.4349457   0.6212981      BaP     AHR2Mo
    ## 183     BaP_0543           14   1.9831535   0.8257172      BaP     AHR2Mo
    ## 184     BaP_0544           18   1.5901155   0.6851979      BaP     AHR2Mo
    ## 185     BaP_0545           13   1.4100273   0.6882498      BaP     AHR2Mo
    ## 186     BaP_0546           20   1.4481003   0.6181108      BaP     AHR2Mo
    ## 187     BaP_0547           15   1.8730731   0.7858866      BaP     AHR2Mo
    ## 188     BaP_0548           16   1.4555081   0.5880449      BaP     AHR2Mo
    ## 189     BaP_0549           15   1.6631217   0.7397746      BaP       CoMo
    ## 190     BaP_0550           18   1.8342936   0.7746328      BaP       CoMo
    ## 191     BaP_0551           13   1.7058252   0.7519537      BaP       CoMo
    ## 192     BaP_0552           14   1.3610035   0.6190912      BaP       CoMo
    ## 193     BaP_0553           15   1.7542992   0.7393483      BaP       CoMo
    ## 194     BaP_0554           31   1.8164412   0.7507869      BaP       CoMo
    ## 195     BaP_0555           14   1.4736822   0.6203809      BaP       CoMo
    ## 196     BaP_0556           15   1.3838688   0.6027038      BaP       CoMo
    ## 197     BaP_0557           14   1.0749063   0.4623406      BaP     AHR2Mo
    ## 198     BaP_0558           12   1.2545657   0.5409020      BaP     AHR2Mo
    ## 199     BaP_0559           13   1.3720506   0.6077584      BaP     AHR2Mo
    ## 200     BaP_0560           15   1.5424006   0.6625499      BaP     AHR2Mo
    ## 201     BaP_0561           12   1.3894178   0.5840512      BaP     AHR2Mo
    ## 202     BaP_0562           16   1.2948532   0.5516323      BaP     AHR2Mo
    ## 203     BaP_0563           12   1.0945083   0.4897972      BaP     AHR2Mo
    ## 204     BaP_0564           13   1.0638495   0.5086387      BaP     AHR2Mo
    ## 205     BaP_0565           14   1.5201629   0.6604762     DMSO     AHR2Mo
    ## 206     BaP_0566           15   0.5924145   0.2281548     DMSO     AHR2Mo
    ## 207     BaP_0567           14   1.2646891   0.5592312     DMSO     AHR2Mo
    ## 208     BaP_0568           25   1.5505481   0.6384222     DMSO     AHR2Mo
    ## 209     BaP_0569           16   1.8921836   0.7918377     DMSO     AHR2Mo
    ## 210     BaP_0570           20   0.9596588   0.3612343     DMSO     AHR2Mo
    ## 211     BaP_0571           21   1.5552347   0.6818153     DMSO     AHR2Mo
    ## 212     BaP_0572           22   1.1347429   0.4590758     DMSO     AHR2Mo
    ## 213     BaP_0573           17   1.8727897   0.7940596      BaP       CoMo
    ## 214     BaP_0574           24   1.5832070   0.6631870      BaP       CoMo
    ## 215     BaP_0575           18   1.7975205   0.7557423      BaP       CoMo
    ## 216     BaP_0576           18   1.5783798   0.6720348      BaP       CoMo
    ## 217     BaP_0649           12   1.3953017   0.6182775      BaP     AHR2Mo
    ## 218     BaP_0650           12   1.7083842   0.7591331      BaP     AHR2Mo
    ## 219     BaP_0651           19   1.6804522   0.7202230      BaP     AHR2Mo
    ## 220     BaP_0652           12   1.7409314   0.7798406      BaP     AHR2Mo
    ## 221     BaP_0653           18   1.4405575   0.6647490      BaP     AHR2Mo
    ## 222     BaP_0654           15   1.8885544   0.8056192      BaP     AHR2Mo
    ## 223     BaP_0655           18   1.5955236   0.7004320      BaP     AHR2Mo
    ## 224     BaP_0656           12   1.2561565   0.5464938      BaP     AHR2Mo
    ## 225     BaP_0657           16   1.7785188   0.7541276      BaP       CoMo
    ## 226     BaP_0658          208   2.8916923   0.8670128      BaP       CoMo
    ## 227     BaP_0659           11   1.2333396   0.6281330      BaP       CoMo
    ## 228     BaP_0660           13   1.3551697   0.5759812      BaP       CoMo
    ## 229     BaP_0661           17   1.3538517   0.5958574      BaP       CoMo
    ## 230     BaP_0662           16   1.6466672   0.7090466      BaP       CoMo
    ## 231     BaP_0663           16   1.1884733   0.5083128      BaP       CoMo
    ## 232     BaP_0664           18   1.9376368   0.8176444      BaP       CoMo
    ## 233     BaP_0665           11   1.0357631   0.5434388      BaP     AHR2Mo
    ## 234     BaP_0666           13   0.8639422   0.3587170      BaP     AHR2Mo
    ## 235     BaP_0667           16   0.4595846   0.2104500      BaP     AHR2Mo
    ## 236     BaP_0668           17   1.9120868   0.7872042      BaP     AHR2Mo
    ## 237     BaP_0669           16   1.2385989   0.5832955      BaP     AHR2Mo
    ## 238     BaP_0670           19   1.4657946   0.6164830      BaP     AHR2Mo
    ## 239     BaP_0671           32   1.0112294   0.4021154      BaP     AHR2Mo
    ## 240     BaP_0672           12   0.9052721   0.4046150      BaP     AHR2Mo
    ## 241     BaP_0673           17   1.5996000   0.7099769      BaP       CoMo
    ## 242     BaP_0674           18   1.5427745   0.7051863      BaP       CoMo
    ## 243     BaP_0675           23   1.8661832   0.7858353      BaP       CoMo
    ## 244     BaP_0676           16   1.8917307   0.8103155      BaP       CoMo
    ## 245     BaP_0677           14   1.1312231   0.5923337      BaP       CoMo
    ## 246     BaP_0678           20   1.9713178   0.8174723      BaP       CoMo
    ## 247     BaP_0679           12   1.2338612   0.5576651      BaP       CoMo
    ## 248     BaP_0680           22   1.7206802   0.7435792      BaP       CoMo
    ## 249     BaP_0683           24   1.9502490   0.8126746      BaP       CoMo
    ## 250     BaP_0757           12   1.4911574   0.6846831     DMSO       CoMo
    ## 251     BaP_0758           20   1.6891783   0.7106677     DMSO       CoMo
    ## 252     BaP_0759           16   1.4918534   0.6563517     DMSO       CoMo
    ## 253     BaP_0760           14   1.7415647   0.7501041     DMSO       CoMo
    ## 254     BaP_0761           19   2.0068434   0.8165190     DMSO       CoMo
    ## 255     BaP_0762           17   1.7746827   0.7540968     DMSO       CoMo
    ## 256     BaP_0763           16   1.4411254   0.6681289     DMSO       CoMo
    ## 257     BaP_0764           13   1.3212162   0.5937944     DMSO       CoMo
    ## 258     BaP_0765           13   1.3578895   0.6536736      BaP       CoMo
    ## 259     BaP_0766           29   1.9911005   0.7774242      BaP       CoMo
    ## 260     BaP_0767           15   1.3441663   0.6440316      BaP       CoMo
    ## 261     BaP_0768           13   1.5376708   0.6997437      BaP       CoMo
    ## 262     BaP_0769           15   1.4746777   0.6503399      BaP       CoMo
    ## 263     BaP_0770           21   1.8536349   0.7678408      BaP       CoMo
    ## 264     BaP_0771           13   1.1551838   0.4954541      BaP       CoMo
    ## 265     BaP_0772           52   1.9194794   0.7473449      BaP       CoMo
    ## 266     BaP_0773           17   1.4243220   0.6144829      BaP     AHR2Mo
    ## 267     BaP_0774           22   1.8422666   0.7410393      BaP     AHR2Mo
    ## 268     BaP_0775           14   1.0849415   0.4592906      BaP     AHR2Mo
    ## 269     BaP_0776           20   1.6597852   0.6836521      BaP     AHR2Mo
    ## 270     BaP_0777           15   1.8507113   0.7795969      BaP     AHR2Mo
    ## 271     BaP_0778           11   1.3267478   0.6457926      BaP     AHR2Mo
    ## 272     BaP_0779           15   1.8333689   0.7809226      BaP     AHR2Mo
    ## 273     BaP_0780           13   1.8065587   0.7689507      BaP     AHR2Mo
    ## 274     BaP_0781           15   1.7750671   0.7841541      BaP       CoMo
    ## 275     BaP_0782           16   1.8136750   0.7948110      BaP       CoMo
    ## 276     BaP_0783           10   0.3862750   0.1393885      BaP       CoMo
    ## 277     BaP_0784           17   1.7700392   0.7464504      BaP       CoMo
    ## 278     BaP_0785           14   1.9834065   0.8162065      BaP       CoMo
    ## 279     BaP_0787           16   2.0014974   0.8050209      BaP       CoMo
    ## 280     BaP_0788           21   1.7542037   0.7235587      BaP       CoMo
    ## 281     BaP_0789           19   1.6760964   0.7509448      BaP       CoMo
    ## 282     BaP_0790           35   1.2382260   0.5491490      BaP       CoMo
    ## 283     BaP_0791           28   2.2014666   0.8522743      BaP       CoMo
    ## 284     BaP_0792           21   1.7274964   0.7184387      BaP       CoMo
    ## 285   BaP_0865_1           28   0.8440744   0.3721503      BaP     AHR2Mo
    ## 286   BaP_0866_2           17   1.4576230   0.6489273      BaP     AHR2Mo
    ## 287   BaP_0867_3           26   1.5905825   0.6858998     DMSO       CoMo
    ## 288   BaP_0868_4           21   1.8156332   0.7559129     DMSO       CoMo
    ## 289   BaP_0869_5           25   1.2702292   0.5107292     DMSO       CoMo
    ## 290   BaP_0870_6           22   1.2830253   0.5810241     DMSO       CoMo
    ## 291   BaP_0871_7           19   1.1182158   0.4605329     DMSO       CoMo
    ## 292   BaP_0872_8           19   1.8238365   0.7984123     DMSO       CoMo
    ## 293   BaP_0873_9           16   1.3204311   0.5839353     DMSO       CoMo
    ## 294  BaP_0874_10           19   1.7265464   0.7123286     DMSO       CoMo
    ## 295  BaP_0875_11           33   2.0411677   0.7855329     DMSO       CoMo
    ## 296  BaP_0876_12           15   1.5176551   0.6915855     DMSO       CoMo
    ## 297  BaP_0877_13           22   1.8002468   0.7363844     DMSO       CoMo
    ## 298  BaP_0878_14           20   1.9266435   0.7952160     DMSO       CoMo
    ## 299  BaP_0879_15           18   1.4656249   0.6180212      BaP       CoMo
    ## 300  BaP_0880_16           18   1.0036599   0.4181813      BaP       CoMo
    ## 301  BaP_0881_17           18   1.3703869   0.6542948      BaP       CoMo
    ## 302  BaP_0882_18           19   1.5674132   0.6648947      BaP       CoMo
    ## 303  BaP_0883_19           21   1.4355634   0.5760808      BaP       CoMo
    ## 304  BaP_0884_20           20   1.4996051   0.6844768      BaP       CoMo
    ## 305  BaP_0885_21           31   1.6325804   0.6886196      BaP       CoMo
    ## 306  BaP_0886_22           20   1.6875573   0.7446084      BaP       CoMo
    ## 307  BaP_0887_23           19   1.1979789   0.4735800      BaP       CoMo
    ## 308  BaP_0888_24           20   1.8373614   0.7972553      BaP       CoMo
    ## 309  BaP_0889_25           16   1.0767715   0.4639642     DMSO     AHR2Mo
    ## 310  BaP_0890_26           32   1.4121768   0.5947607     DMSO     AHR2Mo
    ## 311  BaP_0891_27           32   1.9345137   0.7479137     DMSO     AHR2Mo
    ## 312  BaP_0893_29           19   0.7701844   0.3418286     DMSO     AHR2Mo
    ## 313  BaP_0895_31           27   1.6657049   0.6901205     DMSO     AHR2Mo
    ## 314  BaP_0897_33           16   1.1874874   0.5146849     DMSO     AHR2Mo
    ## 315  BaP_0899_35           23   1.3214634   0.5316480     DMSO     AHR2Mo
    ## 316  BaP_0937_37           14   1.1310148   0.5115819     DMSO       CoMo
    ## 317  BaP_0938_38           20   1.7282664   0.7616901     DMSO       CoMo
    ## 318  BaP_0939_39           18   1.6125481   0.7354112     DMSO       CoMo
    ## 319  BaP_0940_40           26   1.7100552   0.7046484     DMSO       CoMo
    ## 320  BaP_0941_41           13   0.6902676   0.3073755     DMSO       CoMo
    ## 321  BaP_0942_42           23   1.7218268   0.7314220     DMSO       CoMo
    ## 322  BaP_0943_43           15   0.7134299   0.3625318     DMSO       CoMo
    ## 323  BaP_0944_44           40   1.5703546   0.6065130     DMSO       CoMo
    ## 324  BaP_0945_45           14   0.9542451   0.4328225      BaP       CoMo
    ## 325  BaP_0946_46           24   1.4712019   0.6823407      BaP       CoMo
    ## 326  BaP_0947_47           16   1.0479021   0.4553281      BaP       CoMo
    ## 327  BaP_0948_48           20   1.3128642   0.6271243      BaP       CoMo
    ## 328  BaP_0949_49           20   1.2369623   0.5095466      BaP       CoMo
    ## 329  BaP_0950_50           21   1.1560473   0.4602688      BaP       CoMo
    ## 330  BaP_0952_52           35   1.4869488   0.6050766      BaP       CoMo
    ## 331  BaP_0953_53           22   1.5857298   0.7186602      BaP       CoMo
    ## 332  BaP_0954_54           20   1.4906294   0.7091686      BaP       CoMo
    ## 333  BaP_0955_55           16   1.1986457   0.4953309      BaP       CoMo
    ## 334  BaP_0956_56           25   0.9458416   0.3852497      BaP       CoMo
    ## 335  BaP_0957_57           16   1.2275094   0.5922697     DMSO     AHR2Mo
    ## 336  BaP_0958_58           13   0.6270656   0.2902314     DMSO     AHR2Mo
    ## 337  BaP_0959_59           18   1.4064048   0.6645318     DMSO     AHR2Mo
    ## 338  BaP_0961_61           32   1.5646916   0.7034385     DMSO     AHR2Mo
    ## 339  BaP_0962_62           19   0.9369043   0.4073664     DMSO     AHR2Mo
    ## 340  BaP_0963_63           20   1.3137428   0.6313591     DMSO     AHR2Mo
    ## 341  BaP_0964_64           21   1.5084249   0.6612268     DMSO     AHR2Mo
    ## 342  BaP_0965_65           17   1.4332569   0.6450117     DMSO     AHR2Mo
    ## 343  BaP_0966_66           21   1.5028965   0.6286678     DMSO     AHR2Mo
    ## 344  BaP_0967_67           32   1.6730300   0.6719081     DMSO     AHR2Mo
    ## 345  BaP_0968_68           24   1.5761689   0.7214420     DMSO     AHR2Mo
    ## 346  BaP_0969_69           30   1.6438472   0.7246569      BaP     AHR2Mo
    ## 347  BaP_0970_70           24   2.0961639   0.8348115      BaP     AHR2Mo
    ## 348  BaP_0971_71           21   1.6205830   0.6858555      BaP     AHR2Mo
    ## 349  BaP_0972_72           23   1.4823534   0.6859854      BaP     AHR2Mo
    ## 350  BaP_1009_73           20   1.2529448   0.5941495     DMSO       CoMo
    ## 351  BaP_1010_74           18   1.2436749   0.6349144     DMSO       CoMo
    ## 352  BaP_1011_75           16   1.1960174   0.6142970     DMSO       CoMo
    ## 353  BaP_1012_76           15   0.6324907   0.2816366     DMSO       CoMo
    ## 354  BaP_1013_77           14   0.7745159   0.3569255     DMSO       CoMo
    ## 355  BaP_1014_78           22   0.9021513   0.3746210     DMSO       CoMo
    ## 356  BaP_1015_79           18   1.2427808   0.5908974     DMSO       CoMo
    ## 357  BaP_1016_80           24   2.0850612   0.8310963     DMSO       CoMo
    ## 358  BaP_1017_81           15   1.3949389   0.6373362      BaP     AHR2Mo
    ## 359  BaP_1018_82           17   1.9408218   0.7989419      BaP     AHR2Mo
    ## 360  BaP_1019_83           26   1.8720988   0.7814942      BaP       CoMo
    ## 361  BaP_1020_84           21   1.0560262   0.4761014      BaP       CoMo
    ## 362  BaP_1021_85           19   1.1693111   0.5682190      BaP       CoMo
    ## 363  BaP_1022_86           23   1.3713448   0.6485070      BaP       CoMo
    ## 364  BaP_1023_87           21   0.7536806   0.3125816      BaP       CoMo
    ## 365  BaP_1024_88           17   1.9786158   0.8082451      BaP       CoMo
    ## 366  BaP_1025_89           24   1.5874615   0.6804009     DMSO     AHR2Mo
    ## 367  BaP_1026_90           15   0.9376780   0.4329432     DMSO     AHR2Mo
    ## 368  BaP_1027_91           30   1.3806167   0.6538250     DMSO     AHR2Mo
    ## 369  BaP_1028_92           19   1.2579480   0.6124797     DMSO     AHR2Mo
    ## 370  BaP_1029_93           20   1.1975762   0.5593290     DMSO     AHR2Mo
    ## 371  BaP_1030_94           23   1.6215225   0.6782986     DMSO     AHR2Mo
    ## 372  BaP_1031_95           20   1.8544046   0.7865852     DMSO     AHR2Mo
    ## 373  BaP_1032_96           20   1.6017831   0.6923882     DMSO     AHR2Mo
    ## 374  BaP_1033_97           23   1.3247200   0.5760505     DMSO     AHR2Mo
    ## 375  BaP_1034_98           28   0.7665856   0.2927918     DMSO     AHR2Mo
    ## 376  BaP_1035_99           26   1.0973758   0.4570096     DMSO     AHR2Mo
    ## 377 BaP_1036_100           20   1.5453808   0.6464839     DMSO     AHR2Mo
    ## 378 BaP_1037_101           25   0.9123431   0.3755418     DMSO     AHR2Mo
    ## 379 BaP_1038_102           26   1.6657357   0.7051391     DMSO     AHR2Mo
    ## 380 BaP_1039_103           28   0.8000462   0.3135108      BaP     AHR2Mo
    ## 381 BaP_1040_104           21   1.1852220   0.5118023      BaP     AHR2Mo
    ## 382 BaP_1041_105           37   1.4715371   0.6305432      BaP     AHR2Mo
    ## 383 BaP_1042_106           22   1.5813461   0.6877553      BaP     AHR2Mo
    ## 384 BaP_1043_107           24   1.1656256   0.5598563      BaP     AHR2Mo
    ## 385 BaP_1044_108           20   0.8957975   0.3467599      BaP     AHR2Mo
    ## 386 BaP_1081_109           21   1.6209010   0.7316597     DMSO       CoMo
    ## 387 BaP_1082_110           19   1.8602330   0.7468747     DMSO       CoMo
    ## 388 BaP_1083_111           24   1.0740621   0.5709533     DMSO       CoMo
    ## 389 BaP_1084_112           32   1.9425344   0.7715106     DMSO       CoMo
    ## 390 BaP_1085_113            8   1.0891589   0.5826127     DMSO       CoMo
    ## 391 BaP_1086_114           21   1.7411009   0.7123417     DMSO       CoMo
    ## 392 BaP_1087_115           17   1.4740000   0.6875803     DMSO       CoMo
    ## 393 BaP_1088_116           27   2.1075956   0.8265048     DMSO       CoMo
    ## 394 BaP_1089_117           20   1.4636259   0.6164434      BaP       CoMo
    ## 395 BaP_1090_118           46   2.1777498   0.7873784      BaP       CoMo
    ## 396 BaP_1091_119           15   1.1056648   0.5078319      BaP       CoMo
    ## 397 BaP_1092_120           26   1.9116244   0.7886228      BaP       CoMo
    ## 398 BaP_1093_121           19   1.1439063   0.4730745      BaP       CoMo
    ## 399 BaP_1094_122           31   1.3384365   0.6452775      BaP       CoMo
    ## 400 BaP_1095_123           18   1.3375702   0.6528500      BaP       CoMo
    ## 401 BaP_1096_124           24   1.7392803   0.7505465      BaP       CoMo
    ## 402 BaP_1097_125           60   2.1366823   0.7982259     DMSO     AHR2Mo
    ## 403 BaP_1098_126           20   1.5250785   0.6714095     DMSO     AHR2Mo
    ## 404 BaP_1099_127           41   2.0457012   0.8289968     DMSO     AHR2Mo
    ## 405 BaP_1100_128           31   1.4011972   0.5367085     DMSO     AHR2Mo
    ## 406 BaP_1101_129           28   1.5570434   0.6288432     DMSO     AHR2Mo
    ## 407 BaP_1102_130           30   1.9261776   0.7307202     DMSO     AHR2Mo
    ## 408 BaP_1103_131           19   1.6182805   0.6588414     DMSO     AHR2Mo
    ## 409 BaP_1104_132           19   0.9984163   0.5239839     DMSO     AHR2Mo
    ## 410 BaP_1105_133           20   1.4547542   0.6040293      BaP     AHR2Mo
    ## 411 BaP_1106_134           52   1.4894458   0.5544966      BaP     AHR2Mo
    ## 412 BaP_1107_135           19   1.2076802   0.5234009      BaP     AHR2Mo
    ## 413 BaP_1108_136           47   1.6350254   0.6856778      BaP     AHR2Mo
    ## 414 BaP_1109_137           20   1.3695078   0.6584481      BaP     AHR2Mo
    ## 415 BaP_1110_138           22   1.6343283   0.6849217      BaP     AHR2Mo
    ## 416 BaP_1111_139           20   1.1242670   0.4960874      BaP     AHR2Mo
    ## 417 BaP_1112_140           19   0.9384733   0.4328709      BaP     AHR2Mo
    ## 418 BaP_1113_141           22   1.4661090   0.7018937      BaP     AHR2Mo
    ## 419 BaP_1114_142           39   1.4523074   0.5501694      BaP     AHR2Mo
    ## 420 BaP_1115_143           23   1.8146699   0.7471965      BaP     AHR2Mo
    ## 421 BaP_1116_144           33   2.2385886   0.7913963      BaP     AHR2Mo
    ##     Generation Sex  Treatments
    ## 1           F2   F  AHR2Mo_BaP
    ## 2           F2   F  AHR2Mo_BaP
    ## 3           F0   F AHR2Mo_DMSO
    ## 4           F0   F AHR2Mo_DMSO
    ## 5           F0   F AHR2Mo_DMSO
    ## 6           F0   M AHR2Mo_DMSO
    ## 7           F0   F AHR2Mo_DMSO
    ## 8           F0   M AHR2Mo_DMSO
    ## 9           F1   F   CoMo_DMSO
    ## 10          F1   M   CoMo_DMSO
    ## 11          F1   F   CoMo_DMSO
    ## 12          F1   M   CoMo_DMSO
    ## 13          F1   F   CoMo_DMSO
    ## 14          F1   M   CoMo_DMSO
    ## 15          F1   F   CoMo_DMSO
    ## 16          F1   M   CoMo_DMSO
    ## 17          F1   F    CoMo_BaP
    ## 18          F1   M    CoMo_BaP
    ## 19          F1   F AHR2Mo_DMSO
    ## 20          F1   M AHR2Mo_DMSO
    ## 21          F1   F AHR2Mo_DMSO
    ## 22          F1   M AHR2Mo_DMSO
    ## 23          F1   F AHR2Mo_DMSO
    ## 24          F1   M AHR2Mo_DMSO
    ## 25          F2   F   CoMo_DMSO
    ## 26          F2   M   CoMo_DMSO
    ## 27          F2   F   CoMo_DMSO
    ## 28          F2   M   CoMo_DMSO
    ## 29          F2   F   CoMo_DMSO
    ## 30          F2   M   CoMo_DMSO
    ## 31          F2   F   CoMo_DMSO
    ## 32          F2   M   CoMo_DMSO
    ## 33          F2   F   CoMo_DMSO
    ## 34          F2   M   CoMo_DMSO
    ## 35          F2   F   CoMo_DMSO
    ## 36          F2   M   CoMo_DMSO
    ## 37          F0   F AHR2Mo_DMSO
    ## 38          F0   M AHR2Mo_DMSO
    ## 39          F0   F AHR2Mo_DMSO
    ## 40          F0   F AHR2Mo_DMSO
    ## 41          F0   F AHR2Mo_DMSO
    ## 42          F0   M AHR2Mo_DMSO
    ## 43          F0   F AHR2Mo_DMSO
    ## 44          F0   M AHR2Mo_DMSO
    ## 45          F1   F   CoMo_DMSO
    ## 46          F1   M   CoMo_DMSO
    ## 47          F1   F   CoMo_DMSO
    ## 48          F1   M   CoMo_DMSO
    ## 49          F1   F   CoMo_DMSO
    ## 50          F1   M   CoMo_DMSO
    ## 51          F1   F   CoMo_DMSO
    ## 52          F1   M   CoMo_DMSO
    ## 53          F1   F AHR2Mo_DMSO
    ## 54          F1   M AHR2Mo_DMSO
    ## 55          F1   F AHR2Mo_DMSO
    ## 56          F1   M AHR2Mo_DMSO
    ## 57          F1   F AHR2Mo_DMSO
    ## 58          F1   M AHR2Mo_DMSO
    ## 59          F1   F AHR2Mo_DMSO
    ## 60          F1   M AHR2Mo_DMSO
    ## 61          F2   F   CoMo_DMSO
    ## 62          F2   M   CoMo_DMSO
    ## 63          F2   F   CoMo_DMSO
    ## 64          F2   M   CoMo_DMSO
    ## 65          F2   F   CoMo_DMSO
    ## 66          F2   M   CoMo_DMSO
    ## 67          F2   F   CoMo_DMSO
    ## 68          F2   M   CoMo_DMSO
    ## 69          F2   F   CoMo_DMSO
    ## 70          F2   F   CoMo_DMSO
    ## 71          F2   F   CoMo_DMSO
    ## 72          F2   F    CoMo_BaP
    ## 73          F0   F AHR2Mo_DMSO
    ## 74          F0   M AHR2Mo_DMSO
    ## 75          F0   F  AHR2Mo_BaP
    ## 76          F0   M  AHR2Mo_BaP
    ## 77          F0   F  AHR2Mo_BaP
    ## 78          F0   M  AHR2Mo_BaP
    ## 79          F0   F  AHR2Mo_BaP
    ## 80          F0   M  AHR2Mo_BaP
    ## 81          F1   F   CoMo_DMSO
    ## 82          F1   M   CoMo_DMSO
    ## 83          F1   F   CoMo_DMSO
    ## 84          F1   M   CoMo_DMSO
    ## 85          F1   F   CoMo_DMSO
    ## 86          F1   M   CoMo_DMSO
    ## 87          F1   F   CoMo_DMSO
    ## 88          F1   M   CoMo_DMSO
    ## 89          F2   M  AHR2Mo_BaP
    ## 90          F2   M  AHR2Mo_BaP
    ## 91          F1   F AHR2Mo_DMSO
    ## 92          F1   M AHR2Mo_DMSO
    ## 93          F1   F AHR2Mo_DMSO
    ## 94          F1   M AHR2Mo_DMSO
    ## 95          F1   F AHR2Mo_DMSO
    ## 96          F1   M AHR2Mo_DMSO
    ## 97          F2   F    CoMo_BaP
    ## 98          F2   F  AHR2Mo_BaP
    ## 99          F1   F  AHR2Mo_BaP
    ## 100         F1   M  AHR2Mo_BaP
    ## 101         F1   F  AHR2Mo_BaP
    ## 102         F1   M  AHR2Mo_BaP
    ## 103         F1   F  AHR2Mo_BaP
    ## 104         F1   M  AHR2Mo_BaP
    ## 105         F2   F    CoMo_BaP
    ## 106         F2   M    CoMo_BaP
    ## 107         F2   M    CoMo_BaP
    ## 108         F2   F    CoMo_BaP
    ## 109         F0   F  AHR2Mo_BaP
    ## 110         F0   M  AHR2Mo_BaP
    ## 111         F0   F  AHR2Mo_BaP
    ## 112         F0   M  AHR2Mo_BaP
    ## 113         F0   F  AHR2Mo_BaP
    ## 114         F0   M  AHR2Mo_BaP
    ## 115         F0   F  AHR2Mo_BaP
    ## 116         F0   M  AHR2Mo_BaP
    ## 117         F1   F   CoMo_DMSO
    ## 118         F1   M   CoMo_DMSO
    ## 119         F1   F   CoMo_DMSO
    ## 120         F1   M   CoMo_DMSO
    ## 121         F1   F    CoMo_BaP
    ## 122         F1   M    CoMo_BaP
    ## 123         F1   F    CoMo_BaP
    ## 124         F1   M    CoMo_BaP
    ## 125         F1   F AHR2Mo_DMSO
    ## 126         F1   M AHR2Mo_DMSO
    ## 127         F1   F AHR2Mo_DMSO
    ## 128         F1   M AHR2Mo_DMSO
    ## 129         F1   F AHR2Mo_DMSO
    ## 130         F1   M AHR2Mo_DMSO
    ## 131         F1   F AHR2Mo_DMSO
    ## 132         F1   M AHR2Mo_DMSO
    ## 133         F1   F  AHR2Mo_BaP
    ## 134         F1   M  AHR2Mo_BaP
    ## 135         F1   F  AHR2Mo_BaP
    ## 136         F1   M  AHR2Mo_BaP
    ## 137         F2   F   CoMo_DMSO
    ## 138         F2   M   CoMo_DMSO
    ## 139         F2   F   CoMo_DMSO
    ## 140         F2   M   CoMo_DMSO
    ## 141         F2   F    CoMo_BaP
    ## 142         F2   M  AHR2Mo_BaP
    ## 143         F2   M    CoMo_BaP
    ## 144         F2   F    CoMo_BaP
    ## 145         F0   F  AHR2Mo_BaP
    ## 146         F0   M  AHR2Mo_BaP
    ## 147         F0   F  AHR2Mo_BaP
    ## 148         F0   M  AHR2Mo_BaP
    ## 149         F0   F  AHR2Mo_BaP
    ## 150         F0   M  AHR2Mo_BaP
    ## 151         F0   F  AHR2Mo_BaP
    ## 152         F0   M  AHR2Mo_BaP
    ## 153         F2   F  AHR2Mo_BaP
    ## 154         F2   F  AHR2Mo_BaP
    ## 155         F1   F    CoMo_BaP
    ## 156         F1   M    CoMo_BaP
    ## 157         F1   F    CoMo_BaP
    ## 158         F1   M    CoMo_BaP
    ## 159         F1   F    CoMo_BaP
    ## 160         F1   M    CoMo_BaP
    ## 161         F1   F AHR2Mo_DMSO
    ## 162         F1   M AHR2Mo_DMSO
    ## 163         F1   F AHR2Mo_DMSO
    ## 164         F1   M AHR2Mo_DMSO
    ## 165         F1   F AHR2Mo_DMSO
    ## 166         F1   F AHR2Mo_DMSO
    ## 167         F1   F  AHR2Mo_BaP
    ## 168         F1   M  AHR2Mo_BaP
    ## 169         F2   M    CoMo_BaP
    ## 170         F2   F AHR2Mo_DMSO
    ## 171         F2   M AHR2Mo_DMSO
    ## 172         F2   F AHR2Mo_DMSO
    ## 173         F2   M AHR2Mo_DMSO
    ## 174         F2   F AHR2Mo_DMSO
    ## 175         F2   M AHR2Mo_DMSO
    ## 176         F2   F AHR2Mo_DMSO
    ## 177         F2   F   CoMo_DMSO
    ## 178         F2   M   CoMo_DMSO
    ## 179         F2   F   CoMo_DMSO
    ## 180         F2   M   CoMo_DMSO
    ## 181         F0   F  AHR2Mo_BaP
    ## 182         F0   M  AHR2Mo_BaP
    ## 183         F0   F  AHR2Mo_BaP
    ## 184         F0   M  AHR2Mo_BaP
    ## 185         F0   M  AHR2Mo_BaP
    ## 186         F0   F  AHR2Mo_BaP
    ## 187         F0   F  AHR2Mo_BaP
    ## 188         F0   M  AHR2Mo_BaP
    ## 189         F1   F    CoMo_BaP
    ## 190         F1   M    CoMo_BaP
    ## 191         F1   M    CoMo_BaP
    ## 192         F1   M    CoMo_BaP
    ## 193         F1   F    CoMo_BaP
    ## 194         F1   M    CoMo_BaP
    ## 195         F1   F    CoMo_BaP
    ## 196         F1   M    CoMo_BaP
    ## 197         F1   F  AHR2Mo_BaP
    ## 198         F1   M  AHR2Mo_BaP
    ## 199         F1   F  AHR2Mo_BaP
    ## 200         F1   M  AHR2Mo_BaP
    ## 201         F1   F  AHR2Mo_BaP
    ## 202         F1   M  AHR2Mo_BaP
    ## 203         F1   F  AHR2Mo_BaP
    ## 204         F1   M  AHR2Mo_BaP
    ## 205         F2   M AHR2Mo_DMSO
    ## 206         F2   F AHR2Mo_DMSO
    ## 207         F2   M AHR2Mo_DMSO
    ## 208         F2   F AHR2Mo_DMSO
    ## 209         F2   M AHR2Mo_DMSO
    ## 210         F2   F AHR2Mo_DMSO
    ## 211         F2   M AHR2Mo_DMSO
    ## 212         F2   F AHR2Mo_DMSO
    ## 213         F2   M    CoMo_BaP
    ## 214         F2   F    CoMo_BaP
    ## 215         F2   M    CoMo_BaP
    ## 216         F2   F    CoMo_BaP
    ## 217         F2   M  AHR2Mo_BaP
    ## 218         F2   M  AHR2Mo_BaP
    ## 219         F0   F  AHR2Mo_BaP
    ## 220         F0   M  AHR2Mo_BaP
    ## 221         F0   M  AHR2Mo_BaP
    ## 222         F0   M  AHR2Mo_BaP
    ## 223         F0   M  AHR2Mo_BaP
    ## 224         F0   M  AHR2Mo_BaP
    ## 225         F1   F    CoMo_BaP
    ## 226         F1   M    CoMo_BaP
    ## 227         F1   F    CoMo_BaP
    ## 228         F1   M    CoMo_BaP
    ## 229         F1   F    CoMo_BaP
    ## 230         F1   M    CoMo_BaP
    ## 231         F1   F    CoMo_BaP
    ## 232         F1   M    CoMo_BaP
    ## 233         F1   F  AHR2Mo_BaP
    ## 234         F1   M  AHR2Mo_BaP
    ## 235         F1   F  AHR2Mo_BaP
    ## 236         F1   M  AHR2Mo_BaP
    ## 237         F1   F  AHR2Mo_BaP
    ## 238         F1   M  AHR2Mo_BaP
    ## 239         F1   F  AHR2Mo_BaP
    ## 240         F1   M  AHR2Mo_BaP
    ## 241         F2   M    CoMo_BaP
    ## 242         F2   F    CoMo_BaP
    ## 243         F2   M    CoMo_BaP
    ## 244         F2   F    CoMo_BaP
    ## 245         F2   M    CoMo_BaP
    ## 246         F2   F    CoMo_BaP
    ## 247         F2   M    CoMo_BaP
    ## 248         F2   F    CoMo_BaP
    ## 249         F2   M    CoMo_BaP
    ## 250         F1   F   CoMo_DMSO
    ## 251         F1   M   CoMo_DMSO
    ## 252         F1   F   CoMo_DMSO
    ## 253         F1   M   CoMo_DMSO
    ## 254         F1   F   CoMo_DMSO
    ## 255         F1   M   CoMo_DMSO
    ## 256         F1   F   CoMo_DMSO
    ## 257         F1   M   CoMo_DMSO
    ## 258         F1   F    CoMo_BaP
    ## 259         F1   M    CoMo_BaP
    ## 260         F1   F    CoMo_BaP
    ## 261         F1   M    CoMo_BaP
    ## 262         F1   F    CoMo_BaP
    ## 263         F1   M    CoMo_BaP
    ## 264         F1   F    CoMo_BaP
    ## 265         F1   M    CoMo_BaP
    ## 266         F1   F  AHR2Mo_BaP
    ## 267         F1   M  AHR2Mo_BaP
    ## 268         F1   F  AHR2Mo_BaP
    ## 269         F1   M  AHR2Mo_BaP
    ## 270         F1   F  AHR2Mo_BaP
    ## 271         F1   M  AHR2Mo_BaP
    ## 272         F1   F  AHR2Mo_BaP
    ## 273         F1   M  AHR2Mo_BaP
    ## 274         F2   M    CoMo_BaP
    ## 275         F2   F    CoMo_BaP
    ## 276         F2   M    CoMo_BaP
    ## 277         F2   F    CoMo_BaP
    ## 278         F2   M    CoMo_BaP
    ## 279         F2   M    CoMo_BaP
    ## 280         F2   F    CoMo_BaP
    ## 281         F2   M    CoMo_BaP
    ## 282         F2   F    CoMo_BaP
    ## 283         F2   M    CoMo_BaP
    ## 284         F2   F    CoMo_BaP
    ## 285         F2   F  AHR2Mo_BaP
    ## 286         F2   F  AHR2Mo_BaP
    ## 287         F0   F   CoMo_DMSO
    ## 288         F0   M   CoMo_DMSO
    ## 289         F0   F   CoMo_DMSO
    ## 290         F0   M   CoMo_DMSO
    ## 291         F0   F   CoMo_DMSO
    ## 292         F0   M   CoMo_DMSO
    ## 293         F0   F   CoMo_DMSO
    ## 294         F0   M   CoMo_DMSO
    ## 295         F0   F   CoMo_DMSO
    ## 296         F0   M   CoMo_DMSO
    ## 297         F0   F   CoMo_DMSO
    ## 298         F0   M   CoMo_DMSO
    ## 299         F0   F    CoMo_BaP
    ## 300         F0   F    CoMo_BaP
    ## 301         F0   F    CoMo_BaP
    ## 302         F0   M    CoMo_BaP
    ## 303         F0   F    CoMo_BaP
    ## 304         F0   M    CoMo_BaP
    ## 305         F0   F    CoMo_BaP
    ## 306         F0   M    CoMo_BaP
    ## 307         F0   F    CoMo_BaP
    ## 308         F0   M    CoMo_BaP
    ## 309         F2   M AHR2Mo_DMSO
    ## 310         F2   F AHR2Mo_DMSO
    ## 311         F2   M AHR2Mo_DMSO
    ## 312         F2   F AHR2Mo_DMSO
    ## 313         F2   M AHR2Mo_DMSO
    ## 314         F2   F AHR2Mo_DMSO
    ## 315         F2   M AHR2Mo_DMSO
    ## 316         F0   F   CoMo_DMSO
    ## 317         F0   M   CoMo_DMSO
    ## 318         F0   F   CoMo_DMSO
    ## 319         F0   M   CoMo_DMSO
    ## 320         F0   F   CoMo_DMSO
    ## 321         F0   M   CoMo_DMSO
    ## 322         F0   F   CoMo_DMSO
    ## 323         F0   M   CoMo_DMSO
    ## 324         F0   F    CoMo_BaP
    ## 325         F0   M    CoMo_BaP
    ## 326         F0   F    CoMo_BaP
    ## 327         F0   M    CoMo_BaP
    ## 328         F0   F    CoMo_BaP
    ## 329         F0   M    CoMo_BaP
    ## 330         F0   M    CoMo_BaP
    ## 331         F0   F    CoMo_BaP
    ## 332         F0   M    CoMo_BaP
    ## 333         F0   F    CoMo_BaP
    ## 334         F0   M    CoMo_BaP
    ## 335         F0   F AHR2Mo_DMSO
    ## 336         F0   M AHR2Mo_DMSO
    ## 337         F0   F AHR2Mo_DMSO
    ## 338         F2   F AHR2Mo_DMSO
    ## 339         F2   M AHR2Mo_DMSO
    ## 340         F2   F AHR2Mo_DMSO
    ## 341         F2   M AHR2Mo_DMSO
    ## 342         F2   F AHR2Mo_DMSO
    ## 343         F2   M AHR2Mo_DMSO
    ## 344         F2   F AHR2Mo_DMSO
    ## 345         F2   M AHR2Mo_DMSO
    ## 346         F2   F  AHR2Mo_BaP
    ## 347         F2   M  AHR2Mo_BaP
    ## 348         F2   F  AHR2Mo_BaP
    ## 349         F2   M  AHR2Mo_BaP
    ## 350         F0   F   CoMo_DMSO
    ## 351         F0   M   CoMo_DMSO
    ## 352         F0   M   CoMo_DMSO
    ## 353         F0   F   CoMo_DMSO
    ## 354         F0   F   CoMo_DMSO
    ## 355         F0   M   CoMo_DMSO
    ## 356         F0   F   CoMo_DMSO
    ## 357         F0   M   CoMo_DMSO
    ## 358         F2   M  AHR2Mo_BaP
    ## 359         F2   M  AHR2Mo_BaP
    ## 360         F0   F    CoMo_BaP
    ## 361         F0   M    CoMo_BaP
    ## 362         F0   F    CoMo_BaP
    ## 363         F0   M    CoMo_BaP
    ## 364         F0   F    CoMo_BaP
    ## 365         F0   M    CoMo_BaP
    ## 366         F0   F AHR2Mo_DMSO
    ## 367         F0   M AHR2Mo_DMSO
    ## 368         F0   F AHR2Mo_DMSO
    ## 369         F0   M AHR2Mo_DMSO
    ## 370         F0   F AHR2Mo_DMSO
    ## 371         F0   M AHR2Mo_DMSO
    ## 372         F0   F AHR2Mo_DMSO
    ## 373         F0   M AHR2Mo_DMSO
    ## 374         F2   F AHR2Mo_DMSO
    ## 375         F2   M AHR2Mo_DMSO
    ## 376         F2   F AHR2Mo_DMSO
    ## 377         F2   M AHR2Mo_DMSO
    ## 378         F2   F AHR2Mo_DMSO
    ## 379         F2   M AHR2Mo_DMSO
    ## 380         F2   F  AHR2Mo_BaP
    ## 381         F2   M  AHR2Mo_BaP
    ## 382         F2   F  AHR2Mo_BaP
    ## 383         F2   M  AHR2Mo_BaP
    ## 384         F2   F  AHR2Mo_BaP
    ## 385         F2   M  AHR2Mo_BaP
    ## 386         F0   F   CoMo_DMSO
    ## 387         F0   M   CoMo_DMSO
    ## 388         F0   F   CoMo_DMSO
    ## 389         F0   M   CoMo_DMSO
    ## 390         F0   F   CoMo_DMSO
    ## 391         F0   M   CoMo_DMSO
    ## 392         F0   F   CoMo_DMSO
    ## 393         F0   F   CoMo_DMSO
    ## 394         F0   F    CoMo_BaP
    ## 395         F0   M    CoMo_BaP
    ## 396         F0   F    CoMo_BaP
    ## 397         F0   M    CoMo_BaP
    ## 398         F0   F    CoMo_BaP
    ## 399         F0   M    CoMo_BaP
    ## 400         F0   F    CoMo_BaP
    ## 401         F0   M    CoMo_BaP
    ## 402         F0   F AHR2Mo_DMSO
    ## 403         F0   M AHR2Mo_DMSO
    ## 404         F0   F AHR2Mo_DMSO
    ## 405         F0   M AHR2Mo_DMSO
    ## 406         F0   F AHR2Mo_DMSO
    ## 407         F0   M AHR2Mo_DMSO
    ## 408         F0   M AHR2Mo_DMSO
    ## 409         F0   M AHR2Mo_DMSO
    ## 410         F2   F  AHR2Mo_BaP
    ## 411         F2   M  AHR2Mo_BaP
    ## 412         F2   F  AHR2Mo_BaP
    ## 413         F2   M  AHR2Mo_BaP
    ## 414         F2   F  AHR2Mo_BaP
    ## 415         F2   M  AHR2Mo_BaP
    ## 416         F2   F  AHR2Mo_BaP
    ## 417         F2   M  AHR2Mo_BaP
    ## 418         F2   F  AHR2Mo_BaP
    ## 419         F2   M  AHR2Mo_BaP
    ## 420         F2   F  AHR2Mo_BaP
    ## 421         F2   M  AHR2Mo_BaP

    # step AIC to test which model
    testmod_rich <- lm(tax_richness ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino + Sex:Generation + Sex:Exposure + Sex:Morpholino, data = alphadiv)
    AIC_wholedata <- stepAIC(testmod_rich)
    ## Start:  AIC=2101.22
    ## tax_richness ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Sex:Generation + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                         Df Sum of Sq   RSS    AIC
    ## - Exposure:Sex           1      38.3 57702 2099.5
    ## - Generation:Sex         2     465.1 58129 2100.6
    ## <none>                               57664 2101.2
    ## - Exposure:Morpholino    1     291.1 57955 2101.3
    ## - Generation:Exposure    2     594.1 58258 2101.5
    ## - Morpholino:Sex         1     535.2 58199 2103.1
    ## - Generation:Morpholino  2    3901.7 61566 2124.8
    ## 
    ## Step:  AIC=2099.5
    ## tax_richness ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Generation:Sex + 
    ##     Morpholino:Sex
    ## 
    ##                         Df Sum of Sq   RSS    AIC
    ## - Generation:Sex         2     465.2 58168 2098.9
    ## <none>                               57702 2099.5
    ## - Exposure:Morpholino    1     289.5 57992 2099.6
    ## - Generation:Exposure    2     593.2 58296 2099.8
    ## - Morpholino:Sex         1     534.8 58237 2101.4
    ## - Generation:Morpholino  2    3892.4 61595 2123.0
    ## 
    ## Step:  AIC=2098.88
    ## tax_richness ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Morpholino:Sex
    ## 
    ##                         Df Sum of Sq   RSS    AIC
    ## <none>                               58168 2098.9
    ## - Exposure:Morpholino    1     290.1 58458 2099.0
    ## - Generation:Exposure    2     603.0 58771 2099.2
    ## - Morpholino:Sex         1     556.3 58724 2100.9
    ## - Generation:Morpholino  2    3912.9 62081 2122.3

    summary(AIC_wholedata)
    ## 
    ## Call:
    ## lm(formula = tax_richness ~ Generation + Exposure + Morpholino + 
    ##     Sex + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino + 
    ##     Morpholino:Sex, data = alphadiv)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -15.513  -4.507  -1.460   2.058 180.722 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    19.3848     1.9767   9.807  < 2e-16 ***
    ## GenerationF1                    0.5712     2.4408   0.234 0.815073    
    ## GenerationF2                   -8.4352     2.5117  -3.358 0.000858 ***
    ## ExposureBaP                    -0.7502     2.3159  -0.324 0.746152    
    ## MorpholinoAHR2Mo                1.1709     2.5653   0.456 0.648305    
    ## SexM                            4.5570     1.6558   2.752 0.006186 ** 
    ## GenerationF1:ExposureBaP        3.5157     2.8317   1.242 0.215108    
    ## GenerationF2:ExposureBaP        5.8387     2.8592   2.042 0.041783 *  
    ## GenerationF1:MorpholinoAHR2Mo  -6.0585     2.8319  -2.139 0.032997 *  
    ## GenerationF2:MorpholinoAHR2Mo   8.8692     2.8589   3.102 0.002053 ** 
    ## ExposureBaP:MorpholinoAHR2Mo   -3.3254     2.3282  -1.428 0.153969    
    ## MorpholinoAHR2Mo:SexM          -4.6052     2.3286  -1.978 0.048631 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 11.93 on 409 degrees of freedom
    ## Multiple R-squared:  0.09982,    Adjusted R-squared:  0.07561 
    ## F-statistic: 4.123 on 11 and 409 DF,  p-value: 8.967e-06

    # step AIC to test which model
    testmod_shan <- lm(tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + Generation:Morpholino + Exposure:Morpholino + Sex:Generation + Sex:Exposure + Sex:Morpholino, data = alphadiv)
    AIC_wholedata <- stepAIC(testmod_shan)
    ## Start:  AIC=-833.24
    ## tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Sex:Generation + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Generation:Sex         2   0.13865 54.311 -836.17
    ## - Exposure:Sex           1   0.00563 54.178 -835.20
    ## - Generation:Exposure    2   0.26869 54.441 -835.16
    ## - Exposure:Morpholino    1   0.07799 54.250 -834.64
    ## <none>                               54.172 -833.24
    ## - Generation:Morpholino  2   0.87780 55.050 -830.48
    ## - Morpholino:Sex         1   0.61692 54.789 -830.48
    ## 
    ## Step:  AIC=-836.17
    ## tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Exposure:Sex + 
    ##     Morpholino:Sex
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex           1   0.00546 54.316 -838.13
    ## - Generation:Exposure    2   0.27052 54.582 -838.08
    ## - Exposure:Morpholino    1   0.07997 54.391 -837.55
    ## <none>                               54.311 -836.17
    ## - Generation:Morpholino  2   0.88222 55.193 -833.38
    ## - Morpholino:Sex         1   0.62532 54.936 -833.35
    ## 
    ## Step:  AIC=-838.13
    ## tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + 
    ##     Generation:Morpholino + Exposure:Morpholino + Morpholino:Sex
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Generation:Exposure    2   0.27101 54.587 -840.03
    ## - Exposure:Morpholino    1   0.08028 54.397 -839.50
    ## <none>                               54.316 -838.13
    ## - Generation:Morpholino  2   0.88030 55.197 -835.36
    ## - Morpholino:Sex         1   0.62550 54.942 -835.31
    ## 
    ## Step:  AIC=-840.03
    ## tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Morpholino + 
    ##     Exposure:Morpholino + Morpholino:Sex
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino    1   0.07244 54.660 -841.47
    ## <none>                               54.587 -840.03
    ## - Generation:Morpholino  2   0.87779 55.465 -837.31
    ## - Morpholino:Sex         1   0.64602 55.233 -837.08
    ## 
    ## Step:  AIC=-841.47
    ## tax_shannon ~ Generation + Exposure + Morpholino + Sex + Generation:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## <none>                               54.660 -841.47
    ## - Exposure               1   0.49758 55.157 -839.66
    ## - Generation:Morpholino  2   0.88047 55.540 -838.75
    ## - Morpholino:Sex         1   0.66247 55.322 -838.40

    summary(AIC_wholedata)
    ## 
    ## Call:
    ## lm(formula = tax_shannon ~ Generation + Exposure + Morpholino + 
    ##     Sex + Generation:Morpholino + Morpholino:Sex, data = alphadiv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.40548 -0.21086  0.01633  0.23802  1.18914 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    1.30500    0.05239  24.912  < 2e-16 ***
    ## GenerationF1                   0.14370    0.06095   2.358 0.018846 *  
    ## GenerationF2                   0.06806    0.06254   1.088 0.277104    
    ## ExposureBaP                    0.06882    0.03554   1.937 0.053474 .  
    ## MorpholinoAHR2Mo               0.11124    0.07046   1.579 0.115142    
    ## SexM                           0.18503    0.05056   3.660 0.000285 ***
    ## GenerationF1:MorpholinoAHR2Mo -0.22229    0.08648  -2.570 0.010511 *  
    ## GenerationF2:MorpholinoAHR2Mo -0.09791    0.08731  -1.121 0.262743    
    ## MorpholinoAHR2Mo:SexM         -0.15876    0.07105  -2.235 0.025980 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3642 on 412 degrees of freedom
    ## Multiple R-squared:  0.06752,    Adjusted R-squared:  0.04942 
    ## F-statistic: 3.729 on 8 and 412 DF,  p-value: 0.0003157

    alphadiv$Treatments <- factor(alphadiv$Treatments, 
                                    levels = c("CoMo_DMSO",
                                               "CoMo_BaP", 
                                               "AHR2Mo_DMSO",
                                               "AHR2Mo_BaP"))
    alphadiv <- alphadiv[alphadiv$tax_richness < 200,] # remove extreme point
    tax_richplotint <- ggplot(alphadiv, aes(x = Treatments, y = tax_richness,
                                        fill = Generation)) + 
      geom_boxplot() +
      theme_classic() +
      labs(x = NULL, y = "ASV Richness") +
      scale_fill_brewer(palette = "PuBu") +
      scale_x_discrete(labels = c("AhR2Mo - / BaP -",
                                   "AhR2Mo - / BaP +",
                                   "AhR2Mo + / BaP -",
                                   "AhR2Mo + / BaP +")) +
      theme(text = element_text(size = 25))

    tax_richplotint

<img src="zfBaP_16Sdiversity_analysis_files/figure-markdown_strict/manuscript figure with whole dataset alpha diversity across generations-1.png" width="98%" height="98%" />

# F0

Can F0 gut microbiome taxonomic alpha diveristy be predicted by
exposure, morpholino, and other covariates?

    tax_richness <- specnumber(F0_taxa_otutab)
    tax_shannon <- diversity(F0_taxa_otutab, "shannon", base = exp(1))
    tax_simpson <- diversity(F0_taxa_otutab, "simpson")
    F0alphadiv <- data.frame(tax_richness, tax_shannon, tax_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(metaALL, by = "SampleID") %>%
      dplyr::select(c("SampleID", "tax_richness", "tax_shannon", "tax_simpson", "Exposure", "Morpholino", "Generation", "Sex", "Treatments")) # filter to only important variables we are testing
    F0alphadiv
    ##         SampleID tax_richness tax_shannon tax_simpson Exposure Morpholino
    ## 1       BaP_0003           14   1.4702960   0.7070562     DMSO     AHR2Mo
    ## 2       BaP_0004           26   1.2874125   0.5591892     DMSO     AHR2Mo
    ## 3       BaP_0005           16   1.5075086   0.6649067     DMSO     AHR2Mo
    ## 4       BaP_0006           13   0.6520477   0.2725672     DMSO     AHR2Mo
    ## 5       BaP_0007           21   1.5023649   0.6952833     DMSO     AHR2Mo
    ## 6       BaP_0008           18   1.6551501   0.7373428     DMSO     AHR2Mo
    ## 7       BaP_0037           17   1.2778151   0.5661378     DMSO     AHR2Mo
    ## 8       BaP_0038           17   1.4028953   0.6460898     DMSO     AHR2Mo
    ## 9       BaP_0039           17   1.4152523   0.6194654     DMSO     AHR2Mo
    ## 10      BaP_0040           17   1.3627698   0.6397405     DMSO     AHR2Mo
    ## 11      BaP_0041           20   1.6011533   0.6914055     DMSO     AHR2Mo
    ## 12      BaP_0042           20   1.6844358   0.7015727     DMSO     AHR2Mo
    ## 13      BaP_0043           22   1.8229502   0.7557027     DMSO     AHR2Mo
    ## 14      BaP_0044           19   1.4358387   0.6438256     DMSO     AHR2Mo
    ## 15      BaP_0073           14   0.9215374   0.4120802     DMSO     AHR2Mo
    ## 16      BaP_0074           15   0.5066940   0.2013811     DMSO     AHR2Mo
    ## 17      BaP_0075           11   0.8492558   0.3626213      BaP     AHR2Mo
    ## 18      BaP_0076           13   1.0985514   0.5195078      BaP     AHR2Mo
    ## 19      BaP_0077           10   1.2358170   0.6237886      BaP     AHR2Mo
    ## 20      BaP_0078           14   1.4393861   0.6343657      BaP     AHR2Mo
    ## 21      BaP_0079           12   1.2805703   0.6302899      BaP     AHR2Mo
    ## 22      BaP_0080           13   0.5910716   0.2567960      BaP     AHR2Mo
    ## 23      BaP_0109           18   1.5624955   0.6668380      BaP     AHR2Mo
    ## 24      BaP_0110           17   1.8751160   0.7903943      BaP     AHR2Mo
    ## 25      BaP_0111           17   2.2319408   0.8680638      BaP     AHR2Mo
    ## 26      BaP_0112           16   1.4544187   0.6240423      BaP     AHR2Mo
    ## 27      BaP_0113           21   1.8301305   0.7642390      BaP     AHR2Mo
    ## 28      BaP_0114           19   1.2773132   0.5351615      BaP     AHR2Mo
    ## 29      BaP_0115           18   1.7594538   0.7457759      BaP     AHR2Mo
    ## 30      BaP_0116           20   1.9929124   0.8125506      BaP     AHR2Mo
    ## 31      BaP_0433           11   1.5575403   0.6828301      BaP     AHR2Mo
    ## 32      BaP_0434           14   1.3945754   0.6255032      BaP     AHR2Mo
    ## 33      BaP_0435           11   1.7073390   0.7561776      BaP     AHR2Mo
    ## 34      BaP_0436           20   1.2839119   0.5159371      BaP     AHR2Mo
    ## 35      BaP_0437           13   1.7483138   0.7649595      BaP     AHR2Mo
    ## 36      BaP_0438           20   1.5791932   0.6508409      BaP     AHR2Mo
    ## 37      BaP_0439           18   1.8572981   0.7490220      BaP     AHR2Mo
    ## 38      BaP_0440           15   1.3164538   0.5500543      BaP     AHR2Mo
    ## 39      BaP_0541           12   0.8470690   0.3682266      BaP     AHR2Mo
    ## 40      BaP_0542           14   1.4197915   0.6128809      BaP     AHR2Mo
    ## 41      BaP_0543           14   1.9499746   0.8182923      BaP     AHR2Mo
    ## 42      BaP_0544           19   1.5998551   0.6883619      BaP     AHR2Mo
    ## 43      BaP_0545           13   1.4340580   0.6922337      BaP     AHR2Mo
    ## 44      BaP_0546           19   1.4475621   0.6199393      BaP     AHR2Mo
    ## 45      BaP_0547           14   1.8720482   0.7851371      BaP     AHR2Mo
    ## 46      BaP_0548           15   1.4380125   0.5836530      BaP     AHR2Mo
    ## 47      BaP_0651           20   1.6918286   0.7238284      BaP     AHR2Mo
    ## 48      BaP_0652           13   1.7435292   0.7808425      BaP     AHR2Mo
    ## 49      BaP_0653           18   1.4167367   0.6582153      BaP     AHR2Mo
    ## 50      BaP_0654           15   1.8935820   0.8077640      BaP     AHR2Mo
    ## 51      BaP_0655           17   1.5990087   0.7019642      BaP     AHR2Mo
    ## 52      BaP_0656           12   1.2389436   0.5439267      BaP     AHR2Mo
    ## 53    BaP_0867_3           26   1.6016488   0.6902775     DMSO       CoMo
    ## 54    BaP_0868_4           20   1.8278036   0.7560452     DMSO       CoMo
    ## 55    BaP_0869_5           24   1.2710792   0.5073932     DMSO       CoMo
    ## 56    BaP_0870_6           22   1.2780600   0.5847440     DMSO       CoMo
    ## 57    BaP_0871_7           18   1.1141393   0.4612583     DMSO       CoMo
    ## 58    BaP_0872_8           18   1.8126536   0.7964482     DMSO       CoMo
    ## 59    BaP_0873_9           16   1.2954591   0.5708421     DMSO       CoMo
    ## 60   BaP_0874_10           18   1.7220749   0.7150318     DMSO       CoMo
    ## 61   BaP_0875_11           36   2.0461820   0.7887246     DMSO       CoMo
    ## 62   BaP_0876_12           16   1.5036209   0.6826519     DMSO       CoMo
    ## 63   BaP_0877_13           23   1.8005267   0.7358851     DMSO       CoMo
    ## 64   BaP_0878_14           21   1.9496799   0.7982555     DMSO       CoMo
    ## 65   BaP_0879_15           18   1.4389469   0.6104028      BaP       CoMo
    ## 66   BaP_0880_16           18   1.0286000   0.4279351      BaP       CoMo
    ## 67   BaP_0881_17           18   1.3848331   0.6598513      BaP       CoMo
    ## 68   BaP_0882_18           19   1.5801275   0.6688335      BaP       CoMo
    ## 69   BaP_0883_19           21   1.4462549   0.5806136      BaP       CoMo
    ## 70   BaP_0884_20           21   1.5205938   0.6899993      BaP       CoMo
    ## 71   BaP_0885_21           30   1.6273398   0.6889371      BaP       CoMo
    ## 72   BaP_0886_22           20   1.6950093   0.7464219      BaP       CoMo
    ## 73   BaP_0887_23           19   1.1855460   0.4725967      BaP       CoMo
    ## 74   BaP_0888_24           20   1.8296810   0.7951587      BaP       CoMo
    ## 75   BaP_0937_37           14   1.1163170   0.5083295     DMSO       CoMo
    ## 76   BaP_0938_38           19   1.7465840   0.7657325     DMSO       CoMo
    ## 77   BaP_0939_39           18   1.6046709   0.7339528     DMSO       CoMo
    ## 78   BaP_0940_40           26   1.7293748   0.7140150     DMSO       CoMo
    ## 79   BaP_0941_41           13   0.6928288   0.3123534     DMSO       CoMo
    ## 80   BaP_0942_42           23   1.7138286   0.7263249     DMSO       CoMo
    ## 81   BaP_0943_43           16   0.6945077   0.3573019     DMSO       CoMo
    ## 82   BaP_0944_44           37   1.5741116   0.6071224     DMSO       CoMo
    ## 83   BaP_0945_45           15   0.9663698   0.4352287      BaP       CoMo
    ## 84   BaP_0946_46           22   1.4999376   0.6882600      BaP       CoMo
    ## 85   BaP_0947_47           16   1.0485485   0.4559376      BaP       CoMo
    ## 86   BaP_0948_48           20   1.3381241   0.6361579      BaP       CoMo
    ## 87   BaP_0949_49           23   1.2342674   0.5086389      BaP       CoMo
    ## 88   BaP_0950_50           21   1.1551340   0.4618511      BaP       CoMo
    ## 89   BaP_0952_52           31   1.4736219   0.6047398      BaP       CoMo
    ## 90   BaP_0953_53           19   1.5840769   0.7190812      BaP       CoMo
    ## 91   BaP_0954_54           24   1.4786693   0.7067847      BaP       CoMo
    ## 92   BaP_0955_55           17   1.1713751   0.4824578      BaP       CoMo
    ## 93   BaP_0956_56           24   0.9566900   0.3903520      BaP       CoMo
    ## 94   BaP_0957_57           16   1.2117269   0.5827917     DMSO     AHR2Mo
    ## 95   BaP_0958_58           13   0.6406406   0.2955753     DMSO     AHR2Mo
    ## 96   BaP_0959_59           17   1.4194704   0.6631430     DMSO     AHR2Mo
    ## 97   BaP_1009_73           21   1.2746465   0.5999351     DMSO       CoMo
    ## 98   BaP_1010_74           17   1.2371015   0.6347883     DMSO       CoMo
    ## 99   BaP_1011_75           16   1.2231518   0.6226317     DMSO       CoMo
    ## 100  BaP_1012_76           16   0.6432861   0.2855730     DMSO       CoMo
    ## 101  BaP_1013_77           14   0.7814755   0.3640987     DMSO       CoMo
    ## 102  BaP_1014_78           24   0.8806131   0.3650752     DMSO       CoMo
    ## 103  BaP_1015_79           19   1.2397495   0.5931710     DMSO       CoMo
    ## 104  BaP_1016_80           23   2.0834412   0.8305337     DMSO       CoMo
    ## 105  BaP_1019_83           24   1.8759874   0.7838554      BaP       CoMo
    ## 106  BaP_1020_84           21   1.0972252   0.4870963      BaP       CoMo
    ## 107  BaP_1021_85           19   1.1508315   0.5578235      BaP       CoMo
    ## 108  BaP_1022_86           23   1.3550011   0.6418281      BaP       CoMo
    ## 109  BaP_1023_87           19   0.7911481   0.3278321      BaP       CoMo
    ## 110  BaP_1024_88           18   1.9929967   0.8115773      BaP       CoMo
    ## 111  BaP_1025_89           23   1.6036216   0.6830050     DMSO     AHR2Mo
    ## 112  BaP_1026_90           14   0.9204839   0.4320783     DMSO     AHR2Mo
    ## 113  BaP_1027_91           30   1.3594108   0.6474358     DMSO     AHR2Mo
    ## 114  BaP_1028_92           18   1.2704234   0.6188625     DMSO     AHR2Mo
    ## 115  BaP_1029_93           20   1.1823856   0.5525838     DMSO     AHR2Mo
    ## 116  BaP_1030_94           23   1.6168636   0.6736201     DMSO     AHR2Mo
    ## 117  BaP_1031_95           20   1.8476469   0.7837364     DMSO     AHR2Mo
    ## 118  BaP_1032_96           19   1.5874008   0.6836337     DMSO     AHR2Mo
    ## 119 BaP_1081_109           23   1.6139935   0.7289081     DMSO       CoMo
    ## 120 BaP_1082_110           19   1.8810233   0.7557052     DMSO       CoMo
    ## 121 BaP_1083_111           26   1.0790460   0.5723466     DMSO       CoMo
    ## 122 BaP_1084_112           34   1.9608982   0.7750959     DMSO       CoMo
    ## 123 BaP_1085_113            8   1.0740533   0.5782383     DMSO       CoMo
    ## 124 BaP_1086_114           22   1.7421617   0.7065174     DMSO       CoMo
    ## 125 BaP_1087_115           17   1.4617438   0.6837158     DMSO       CoMo
    ## 126 BaP_1088_116           27   2.1007223   0.8250435     DMSO       CoMo
    ## 127 BaP_1089_117           20   1.4376723   0.6117348      BaP       CoMo
    ## 128 BaP_1090_118           47   2.1839469   0.7866647      BaP       CoMo
    ## 129 BaP_1091_119           15   1.1103716   0.5064895      BaP       CoMo
    ## 130 BaP_1092_120           27   1.9013622   0.7872869      BaP       CoMo
    ## 131 BaP_1093_121           19   1.1619930   0.4802705      BaP       CoMo
    ## 132 BaP_1094_122           33   1.3587696   0.6478143      BaP       CoMo
    ## 133 BaP_1095_123           18   1.3469446   0.6548077      BaP       CoMo
    ## 134 BaP_1096_124           23   1.7177060   0.7456826      BaP       CoMo
    ## 135 BaP_1097_125           61   2.0996829   0.7925030     DMSO     AHR2Mo
    ## 136 BaP_1098_126           21   1.5106550   0.6633081     DMSO     AHR2Mo
    ## 137 BaP_1099_127           43   2.0224400   0.8267781     DMSO     AHR2Mo
    ## 138 BaP_1100_128           29   1.4187523   0.5433279     DMSO     AHR2Mo
    ## 139 BaP_1101_129           29   1.5362398   0.6238639     DMSO     AHR2Mo
    ## 140 BaP_1102_130           30   1.8841590   0.7220295     DMSO     AHR2Mo
    ## 141 BaP_1103_131           18   1.6234920   0.6625203     DMSO     AHR2Mo
    ## 142 BaP_1104_132           19   0.9935530   0.5177490     DMSO     AHR2Mo
    ##     Generation Sex  Treatments
    ## 1           F0   F AHR2Mo_DMSO
    ## 2           F0   F AHR2Mo_DMSO
    ## 3           F0   F AHR2Mo_DMSO
    ## 4           F0   M AHR2Mo_DMSO
    ## 5           F0   F AHR2Mo_DMSO
    ## 6           F0   M AHR2Mo_DMSO
    ## 7           F0   F AHR2Mo_DMSO
    ## 8           F0   M AHR2Mo_DMSO
    ## 9           F0   F AHR2Mo_DMSO
    ## 10          F0   F AHR2Mo_DMSO
    ## 11          F0   F AHR2Mo_DMSO
    ## 12          F0   M AHR2Mo_DMSO
    ## 13          F0   F AHR2Mo_DMSO
    ## 14          F0   M AHR2Mo_DMSO
    ## 15          F0   F AHR2Mo_DMSO
    ## 16          F0   M AHR2Mo_DMSO
    ## 17          F0   F  AHR2Mo_BaP
    ## 18          F0   M  AHR2Mo_BaP
    ## 19          F0   F  AHR2Mo_BaP
    ## 20          F0   M  AHR2Mo_BaP
    ## 21          F0   F  AHR2Mo_BaP
    ## 22          F0   M  AHR2Mo_BaP
    ## 23          F0   F  AHR2Mo_BaP
    ## 24          F0   M  AHR2Mo_BaP
    ## 25          F0   F  AHR2Mo_BaP
    ## 26          F0   M  AHR2Mo_BaP
    ## 27          F0   F  AHR2Mo_BaP
    ## 28          F0   M  AHR2Mo_BaP
    ## 29          F0   F  AHR2Mo_BaP
    ## 30          F0   M  AHR2Mo_BaP
    ## 31          F0   F  AHR2Mo_BaP
    ## 32          F0   M  AHR2Mo_BaP
    ## 33          F0   F  AHR2Mo_BaP
    ## 34          F0   M  AHR2Mo_BaP
    ## 35          F0   F  AHR2Mo_BaP
    ## 36          F0   M  AHR2Mo_BaP
    ## 37          F0   F  AHR2Mo_BaP
    ## 38          F0   M  AHR2Mo_BaP
    ## 39          F0   F  AHR2Mo_BaP
    ## 40          F0   M  AHR2Mo_BaP
    ## 41          F0   F  AHR2Mo_BaP
    ## 42          F0   M  AHR2Mo_BaP
    ## 43          F0   M  AHR2Mo_BaP
    ## 44          F0   F  AHR2Mo_BaP
    ## 45          F0   F  AHR2Mo_BaP
    ## 46          F0   M  AHR2Mo_BaP
    ## 47          F0   F  AHR2Mo_BaP
    ## 48          F0   M  AHR2Mo_BaP
    ## 49          F0   M  AHR2Mo_BaP
    ## 50          F0   M  AHR2Mo_BaP
    ## 51          F0   M  AHR2Mo_BaP
    ## 52          F0   M  AHR2Mo_BaP
    ## 53          F0   F   CoMo_DMSO
    ## 54          F0   M   CoMo_DMSO
    ## 55          F0   F   CoMo_DMSO
    ## 56          F0   M   CoMo_DMSO
    ## 57          F0   F   CoMo_DMSO
    ## 58          F0   M   CoMo_DMSO
    ## 59          F0   F   CoMo_DMSO
    ## 60          F0   M   CoMo_DMSO
    ## 61          F0   F   CoMo_DMSO
    ## 62          F0   M   CoMo_DMSO
    ## 63          F0   F   CoMo_DMSO
    ## 64          F0   M   CoMo_DMSO
    ## 65          F0   F    CoMo_BaP
    ## 66          F0   F    CoMo_BaP
    ## 67          F0   F    CoMo_BaP
    ## 68          F0   M    CoMo_BaP
    ## 69          F0   F    CoMo_BaP
    ## 70          F0   M    CoMo_BaP
    ## 71          F0   F    CoMo_BaP
    ## 72          F0   M    CoMo_BaP
    ## 73          F0   F    CoMo_BaP
    ## 74          F0   M    CoMo_BaP
    ## 75          F0   F   CoMo_DMSO
    ## 76          F0   M   CoMo_DMSO
    ## 77          F0   F   CoMo_DMSO
    ## 78          F0   M   CoMo_DMSO
    ## 79          F0   F   CoMo_DMSO
    ## 80          F0   M   CoMo_DMSO
    ## 81          F0   F   CoMo_DMSO
    ## 82          F0   M   CoMo_DMSO
    ## 83          F0   F    CoMo_BaP
    ## 84          F0   M    CoMo_BaP
    ## 85          F0   F    CoMo_BaP
    ## 86          F0   M    CoMo_BaP
    ## 87          F0   F    CoMo_BaP
    ## 88          F0   M    CoMo_BaP
    ## 89          F0   M    CoMo_BaP
    ## 90          F0   F    CoMo_BaP
    ## 91          F0   M    CoMo_BaP
    ## 92          F0   F    CoMo_BaP
    ## 93          F0   M    CoMo_BaP
    ## 94          F0   F AHR2Mo_DMSO
    ## 95          F0   M AHR2Mo_DMSO
    ## 96          F0   F AHR2Mo_DMSO
    ## 97          F0   F   CoMo_DMSO
    ## 98          F0   M   CoMo_DMSO
    ## 99          F0   M   CoMo_DMSO
    ## 100         F0   F   CoMo_DMSO
    ## 101         F0   F   CoMo_DMSO
    ## 102         F0   M   CoMo_DMSO
    ## 103         F0   F   CoMo_DMSO
    ## 104         F0   M   CoMo_DMSO
    ## 105         F0   F    CoMo_BaP
    ## 106         F0   M    CoMo_BaP
    ## 107         F0   F    CoMo_BaP
    ## 108         F0   M    CoMo_BaP
    ## 109         F0   F    CoMo_BaP
    ## 110         F0   M    CoMo_BaP
    ## 111         F0   F AHR2Mo_DMSO
    ## 112         F0   M AHR2Mo_DMSO
    ## 113         F0   F AHR2Mo_DMSO
    ## 114         F0   M AHR2Mo_DMSO
    ## 115         F0   F AHR2Mo_DMSO
    ## 116         F0   M AHR2Mo_DMSO
    ## 117         F0   F AHR2Mo_DMSO
    ## 118         F0   M AHR2Mo_DMSO
    ## 119         F0   F   CoMo_DMSO
    ## 120         F0   M   CoMo_DMSO
    ## 121         F0   F   CoMo_DMSO
    ## 122         F0   M   CoMo_DMSO
    ## 123         F0   F   CoMo_DMSO
    ## 124         F0   M   CoMo_DMSO
    ## 125         F0   F   CoMo_DMSO
    ## 126         F0   F   CoMo_DMSO
    ## 127         F0   F    CoMo_BaP
    ## 128         F0   M    CoMo_BaP
    ## 129         F0   F    CoMo_BaP
    ## 130         F0   M    CoMo_BaP
    ## 131         F0   F    CoMo_BaP
    ## 132         F0   M    CoMo_BaP
    ## 133         F0   F    CoMo_BaP
    ## 134         F0   M    CoMo_BaP
    ## 135         F0   F AHR2Mo_DMSO
    ## 136         F0   M AHR2Mo_DMSO
    ## 137         F0   F AHR2Mo_DMSO
    ## 138         F0   M AHR2Mo_DMSO
    ## 139         F0   F AHR2Mo_DMSO
    ## 140         F0   M AHR2Mo_DMSO
    ## 141         F0   M AHR2Mo_DMSO
    ## 142         F0   M AHR2Mo_DMSO

    # step AIC to test which model
    testmod_richF0 <- lm(tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F0alphadiv)
    AIC_wholedata_richF0 <- stepAIC(testmod_richF0)
    ## Start:  AIC=530.55
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## <none>                             5396.6 530.55
    ## - Exposure:Sex         1    133.87 5530.5 532.03
    ## - Morpholino:Sex       1    248.78 5645.4 534.95
    ## - Exposure:Morpholino  1    406.78 5803.4 538.87

    summary(AIC_wholedata_richF0)
    ## 
    ## Call:
    ## lm(formula = tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino, data = F0alphadiv)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.016  -3.665  -1.039   2.395  37.963 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   20.0161     1.3637  14.677  < 2e-16 ***
    ## ExposureBaP                   -0.9774     1.8147  -0.539  0.59102    
    ## MorpholinoAHR2Mo               3.0205     1.7975   1.680  0.09520 .  
    ## SexM                           1.7307     1.8331   0.944  0.34679    
    ## ExposureBaP:MorpholinoAHR2Mo  -6.7901     2.1286  -3.190  0.00177 ** 
    ## ExposureBaP:SexM               3.8957     2.1288   1.830  0.06946 .  
    ## MorpholinoAHR2Mo:SexM         -5.3107     2.1288  -2.495  0.01381 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.323 on 135 degrees of freedom
    ## Multiple R-squared:  0.2016, Adjusted R-squared:  0.1661 
    ## F-statistic: 5.681 on 6 and 135 DF,  p-value: 2.727e-05

    # step AIC to test which model
    testmod_shanF0 <- lm(tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F0alphadiv)
    AIC_wholedata_shanF0 <- stepAIC(testmod_shanF0)
    ## Start:  AIC=-294.51
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex         1   0.00169 16.173 -296.50
    ## <none>                             16.171 -294.51
    ## - Exposure:Morpholino  1   0.28273 16.454 -294.05
    ## - Morpholino:Sex       1   1.95579 18.127 -280.30
    ## 
    ## Step:  AIC=-296.5
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## <none>                             16.173 -296.50
    ## - Exposure:Morpholino  1   0.28176 16.454 -296.04
    ## - Morpholino:Sex       1   1.95749 18.130 -282.27

    summary(AIC_wholedata_shanF0)
    ## 
    ## Call:
    ## lm(formula = tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex, data = F0alphadiv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.84957 -0.20039 -0.00461  0.22078  0.79011 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   1.31062    0.06928  18.917  < 2e-16 ***
    ## ExposureBaP                  -0.05531    0.08187  -0.676 0.500413    
    ## MorpholinoAHR2Mo              0.17244    0.09804   1.759 0.080842 .  
    ## SexM                          0.30530    0.08193   3.726 0.000284 ***
    ## ExposureBaP:MorpholinoAHR2Mo  0.17864    0.11605   1.539 0.126057    
    ## MorpholinoAHR2Mo:SexM        -0.47104    0.11610  -4.057 8.33e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3448 on 136 degrees of freedom
    ## Multiple R-squared:  0.1302, Adjusted R-squared:  0.09818 
    ## F-statistic:  4.07 on 5 and 136 DF,  p-value: 0.001786

    pal_shan <- c("#5F3C98", "#E76100") # pallete for morpholino

    taxF0_richplot <- ggplot(F0alphadiv) +
      geom_boxplot(aes(x = factor(Exposure), y = tax_richness, fill = Morpholino),
                   position = position_dodge(0.7), width = 0.5, alpha = 0.5,
                   outlier.shape = NA) +
      geom_point(aes(x = factor(Exposure), y = tax_richness, color = Morpholino),
                 position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.7, seed = 3)) +
      stat_smooth(aes(x = as.integer(Exposure), y = tax_richness, color = Morpholino, fill = Morpholino),
                  method = lm, se = T,
                  position = position_dodge(0.7)) +
      stat_regline_equation(aes(x = as.integer(Exposure), y = tax_richness, color = Morpholino),
                            label.x.npc = "left") +
      labs(x = "", y = "Taxonomic Richness") +
      theme_classic() +
      scale_fill_manual(values = pal_shan) +
      scale_color_manual(values = pal_shan, 
                         labels = c("Control morpholino",
                                    "AhR2 morpholino")) +
      theme(text = element_text(size = 20),
            legend.title = element_blank()) +
      guides(fill = "none") +
      ylim(0, range(F0alphadiv$tax_richness)[2]+10) +
      facet_wrap("Sex")

    taxF0_richplot

<img src="zfBaP_16Sdiversity_analysis_files/figure-markdown_strict/make manuscript fig of F0 16S richness-1.png" width="98%" height="98%" />

# F1

Can F1 gut microbiome taxonomic alpha diveristy be predicted by
exposure, morpholino, and other covariates?

    tax_richness <- specnumber(F1_taxa_otutab)
    tax_shannon <- diversity(F1_taxa_otutab, "shannon", base = exp(1))
    tax_simpson <- diversity(F1_taxa_otutab, "simpson")
    F1alphadiv <- data.frame(tax_richness, tax_shannon, tax_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(metaALL, by = "SampleID") %>%
      dplyr::select(c("SampleID", "tax_richness", "tax_shannon", "tax_simpson", "Exposure", "Morpholino", "Generation", "Sex", "Treatments")) # filter to only important variables we are testing
    F1alphadiv
    ##     SampleID tax_richness tax_shannon tax_simpson Exposure Morpholino
    ## 1   BaP_0009           26   1.0550000  0.44580036     DMSO       CoMo
    ## 2   BaP_0010           40   1.8451294  0.75681460     DMSO       CoMo
    ## 3   BaP_0011           22   1.7775087  0.75395438     DMSO       CoMo
    ## 4   BaP_0012           25   1.8395631  0.80548266     DMSO       CoMo
    ## 5   BaP_0013           14   1.2274427  0.54726444     DMSO       CoMo
    ## 6   BaP_0014           88   1.9489449  0.69153414     DMSO       CoMo
    ## 7   BaP_0015           19   1.1601736  0.58720282     DMSO       CoMo
    ## 8   BaP_0016           38   1.8273867  0.74688796     DMSO       CoMo
    ## 9   BaP_0017           81   2.6688424  0.86346058      BaP       CoMo
    ## 10  BaP_0018           22   1.4891461  0.67375562      BaP       CoMo
    ## 11  BaP_0019           11   1.0397329  0.52664836     DMSO     AHR2Mo
    ## 12  BaP_0020           13   1.6036382  0.72185798     DMSO     AHR2Mo
    ## 13  BaP_0021           13   1.3752092  0.61815616     DMSO     AHR2Mo
    ## 14  BaP_0022           11   1.4045409  0.67108270     DMSO     AHR2Mo
    ## 15  BaP_0023           11   1.2280969  0.53416564     DMSO     AHR2Mo
    ## 16  BaP_0024           13   1.3990105  0.65710200     DMSO     AHR2Mo
    ## 17  BaP_0045           16   1.2837731  0.59229858     DMSO       CoMo
    ## 18  BaP_0046           24   1.7118275  0.71145026     DMSO       CoMo
    ## 19  BaP_0047           18   1.6166866  0.67324536     DMSO       CoMo
    ## 20  BaP_0048           20   1.5827749  0.66564594     DMSO       CoMo
    ## 21  BaP_0049           23   1.8530236  0.75010548     DMSO       CoMo
    ## 22  BaP_0050           29   1.7314137  0.72549000     DMSO       CoMo
    ## 23  BaP_0051           24   1.8663432  0.80464032     DMSO       CoMo
    ## 24  BaP_0052           21   1.4246801  0.66462826     DMSO       CoMo
    ## 25  BaP_0053           19   1.7211848  0.73814090     DMSO     AHR2Mo
    ## 26  BaP_0054           31   1.7824508  0.74736260     DMSO     AHR2Mo
    ## 27  BaP_0055           11   1.4789403  0.70526496     DMSO     AHR2Mo
    ## 28  BaP_0056           10   0.7300687  0.34187032     DMSO     AHR2Mo
    ## 29  BaP_0057           10   1.3826247  0.70375550     DMSO     AHR2Mo
    ## 30  BaP_0058           13   1.4447184  0.68003122     DMSO     AHR2Mo
    ## 31  BaP_0059           10   1.2321743  0.63872098     DMSO     AHR2Mo
    ## 32  BaP_0060           15   1.1440807  0.47777496     DMSO     AHR2Mo
    ## 33  BaP_0081           14   0.9327751  0.42545580     DMSO       CoMo
    ## 34  BaP_0082           33   2.2960205  0.85363422     DMSO       CoMo
    ## 35  BaP_0083           18   0.7726131  0.32815926     DMSO       CoMo
    ## 36  BaP_0084            9   0.2217718  0.07713258     DMSO       CoMo
    ## 37  BaP_0085           24   0.6489038  0.24284646     DMSO       CoMo
    ## 38  BaP_0086           20   1.4842003  0.69202694     DMSO       CoMo
    ## 39  BaP_0087           21   0.6222945  0.24494974     DMSO       CoMo
    ## 40  BaP_0088           19   0.5948992  0.22673664     DMSO       CoMo
    ## 41  BaP_0091           15   1.6293541  0.74038680     DMSO     AHR2Mo
    ## 42  BaP_0092           19   1.1455033  0.46563098     DMSO     AHR2Mo
    ## 43  BaP_0093           12   1.5605224  0.70102114     DMSO     AHR2Mo
    ## 44  BaP_0094           15   0.8796767  0.35777704     DMSO     AHR2Mo
    ## 45  BaP_0095           15   1.4540506  0.61654080     DMSO     AHR2Mo
    ## 46  BaP_0096           14   1.2521365  0.54694464     DMSO     AHR2Mo
    ## 47  BaP_0099           14   1.5355727  0.69497208      BaP     AHR2Mo
    ## 48  BaP_0100           22   1.4721796  0.63156740      BaP     AHR2Mo
    ## 49  BaP_0101           15   1.6550513  0.72283624      BaP     AHR2Mo
    ## 50  BaP_0102           14   1.3638806  0.64795984      BaP     AHR2Mo
    ## 51  BaP_0103            9   0.2955045  0.11006620      BaP     AHR2Mo
    ## 52  BaP_0104           15   1.1582457  0.50705410      BaP     AHR2Mo
    ## 53  BaP_0117           17   2.0001880  0.81786700     DMSO       CoMo
    ## 54  BaP_0118           19   2.0618074  0.83156476     DMSO       CoMo
    ## 55  BaP_0119           20   1.7847182  0.76583526     DMSO       CoMo
    ## 56  BaP_0120           40   2.7363946  0.90564940     DMSO       CoMo
    ## 57  BaP_0121           17   1.6876949  0.71693862      BaP       CoMo
    ## 58  BaP_0122           14   1.4536910  0.67662792      BaP       CoMo
    ## 59  BaP_0123           23   1.9244152  0.79702102      BaP       CoMo
    ## 60  BaP_0124           15   1.4995987  0.63822224      BaP       CoMo
    ## 61  BaP_0125           18   1.7734295  0.74552458     DMSO     AHR2Mo
    ## 62  BaP_0126           21   2.0539329  0.82999658     DMSO     AHR2Mo
    ## 63  BaP_0127           12   1.5199397  0.68549254     DMSO     AHR2Mo
    ## 64  BaP_0128           14   1.7079675  0.77312174     DMSO     AHR2Mo
    ## 65  BaP_0129           10   1.3721468  0.63337654     DMSO     AHR2Mo
    ## 66  BaP_0130            9   1.7056329  0.77346718     DMSO     AHR2Mo
    ## 67  BaP_0131           11   1.4547375  0.66388460     DMSO     AHR2Mo
    ## 68  BaP_0132           11   1.3384380  0.61798558     DMSO     AHR2Mo
    ## 69  BaP_0133           15   1.5558135  0.70374590      BaP     AHR2Mo
    ## 70  BaP_0134           10   1.7727817  0.80577684      BaP     AHR2Mo
    ## 71  BaP_0135           17   1.6594704  0.71973032      BaP     AHR2Mo
    ## 72  BaP_0136           13   1.7134776  0.78245524      BaP     AHR2Mo
    ## 73  BaP_0443           12   1.5300453  0.65001746      BaP       CoMo
    ## 74  BaP_0444           13   1.3757184  0.65290140      BaP       CoMo
    ## 75  BaP_0445           14   1.1814603  0.58293076      BaP       CoMo
    ## 76  BaP_0446           15   1.7575532  0.76201578      BaP       CoMo
    ## 77  BaP_0447           13   1.2925654  0.54425876      BaP       CoMo
    ## 78  BaP_0448           14   1.8854291  0.79692338      BaP       CoMo
    ## 79  BaP_0449           25   0.7796233  0.33810144     DMSO     AHR2Mo
    ## 80  BaP_0450           14   0.8514230  0.35106038     DMSO     AHR2Mo
    ## 81  BaP_0451           18   1.2966828  0.53130260     DMSO     AHR2Mo
    ## 82  BaP_0452           18   1.7722585  0.74707008     DMSO     AHR2Mo
    ## 83  BaP_0453           14   1.4821911  0.67318604     DMSO     AHR2Mo
    ## 84  BaP_0454           13   1.1560871  0.50834240     DMSO     AHR2Mo
    ## 85  BaP_0455           15   1.8600101  0.79424888      BaP     AHR2Mo
    ## 86  BaP_0456           20   1.8114307  0.77150948      BaP     AHR2Mo
    ## 87  BaP_0549           15   1.6392488  0.73341560      BaP       CoMo
    ## 88  BaP_0550           16   1.8471078  0.78044812      BaP       CoMo
    ## 89  BaP_0551           14   1.7002237  0.74904136      BaP       CoMo
    ## 90  BaP_0552           13   1.3573809  0.62100126      BaP       CoMo
    ## 91  BaP_0553           17   1.7760802  0.74532272      BaP       CoMo
    ## 92  BaP_0554           31   1.8026709  0.74625802      BaP       CoMo
    ## 93  BaP_0555           13   1.4664865  0.62153678      BaP       CoMo
    ## 94  BaP_0556           15   1.3980040  0.60987172      BaP       CoMo
    ## 95  BaP_0557           14   1.0731049  0.45756680      BaP     AHR2Mo
    ## 96  BaP_0558           12   1.2680555  0.54777022      BaP     AHR2Mo
    ## 97  BaP_0559           13   1.3759648  0.61208432      BaP     AHR2Mo
    ## 98  BaP_0560           17   1.5159011  0.65138650      BaP     AHR2Mo
    ## 99  BaP_0561           14   1.4162226  0.59187968      BaP     AHR2Mo
    ## 100 BaP_0562           17   1.3054448  0.55924794      BaP     AHR2Mo
    ## 101 BaP_0563           12   1.1049245  0.49125908      BaP     AHR2Mo
    ## 102 BaP_0564           13   1.0601210  0.51206330      BaP     AHR2Mo
    ## 103 BaP_0657           16   1.8141153  0.76634080      BaP       CoMo
    ## 104 BaP_0658          203   2.9147982  0.87040306      BaP       CoMo
    ## 105 BaP_0659           11   1.2155143  0.62327674      BaP       CoMo
    ## 106 BaP_0660           13   1.3411320  0.56967466      BaP       CoMo
    ## 107 BaP_0661           16   1.3776567  0.60863176      BaP       CoMo
    ## 108 BaP_0662           16   1.6655793  0.71321830      BaP       CoMo
    ## 109 BaP_0663           16   1.1652533  0.50125680      BaP       CoMo
    ## 110 BaP_0664           17   1.9403434  0.81700230      BaP       CoMo
    ## 111 BaP_0665           11   1.0370641  0.54701554      BaP     AHR2Mo
    ## 112 BaP_0666           13   0.8657569  0.35727638      BaP     AHR2Mo
    ## 113 BaP_0667           13   0.4462159  0.20232516      BaP     AHR2Mo
    ## 114 BaP_0668           17   1.9206588  0.78634382      BaP     AHR2Mo
    ## 115 BaP_0669           16   1.2221205  0.57974126      BaP     AHR2Mo
    ## 116 BaP_0670           19   1.4562172  0.61107550      BaP     AHR2Mo
    ## 117 BaP_0671           28   1.0241325  0.40369958      BaP     AHR2Mo
    ## 118 BaP_0672           12   0.9007960  0.40286256      BaP     AHR2Mo
    ## 119 BaP_0757           12   1.4769542  0.68008148     DMSO       CoMo
    ## 120 BaP_0758           19   1.7049445  0.71408212     DMSO       CoMo
    ## 121 BaP_0759           16   1.4836858  0.65114226     DMSO       CoMo
    ## 122 BaP_0760           14   1.7361974  0.74997588     DMSO       CoMo
    ## 123 BaP_0761           18   2.0092002  0.81690188     DMSO       CoMo
    ## 124 BaP_0762           16   1.7601725  0.75141852     DMSO       CoMo
    ## 125 BaP_0763           16   1.4268191  0.66385444     DMSO       CoMo
    ## 126 BaP_0764           13   1.3397004  0.60181034     DMSO       CoMo
    ## 127 BaP_0765           13   1.3414715  0.65060314      BaP       CoMo
    ## 128 BaP_0766           31   1.9933107  0.77709058      BaP       CoMo
    ## 129 BaP_0767           15   1.3458115  0.64848714      BaP       CoMo
    ## 130 BaP_0768           13   1.5450093  0.70006290      BaP       CoMo
    ## 131 BaP_0769           15   1.4872148  0.65206434      BaP       CoMo
    ## 132 BaP_0770           19   1.8596406  0.77284898      BaP       CoMo
    ## 133 BaP_0771           13   1.1533317  0.49433218      BaP       CoMo
    ## 134 BaP_0772           52   1.9371315  0.74690198      BaP       CoMo
    ## 135 BaP_0773           16   1.4599104  0.62565726      BaP     AHR2Mo
    ## 136 BaP_0774           22   1.8513868  0.74285336      BaP     AHR2Mo
    ## 137 BaP_0775           14   1.0572848  0.44616986      BaP     AHR2Mo
    ## 138 BaP_0776           19   1.6377811  0.67833482      BaP     AHR2Mo
    ## 139 BaP_0777           15   1.8377131  0.77481692      BaP     AHR2Mo
    ## 140 BaP_0778           11   1.3310334  0.64588096      BaP     AHR2Mo
    ## 141 BaP_0779           15   1.8357555  0.77939650      BaP     AHR2Mo
    ## 142 BaP_0780           13   1.8066978  0.76762930      BaP     AHR2Mo
    ##     Generation Sex  Treatments
    ## 1           F1   F   CoMo_DMSO
    ## 2           F1   M   CoMo_DMSO
    ## 3           F1   F   CoMo_DMSO
    ## 4           F1   M   CoMo_DMSO
    ## 5           F1   F   CoMo_DMSO
    ## 6           F1   M   CoMo_DMSO
    ## 7           F1   F   CoMo_DMSO
    ## 8           F1   M   CoMo_DMSO
    ## 9           F1   F    CoMo_BaP
    ## 10          F1   M    CoMo_BaP
    ## 11          F1   F AHR2Mo_DMSO
    ## 12          F1   M AHR2Mo_DMSO
    ## 13          F1   F AHR2Mo_DMSO
    ## 14          F1   M AHR2Mo_DMSO
    ## 15          F1   F AHR2Mo_DMSO
    ## 16          F1   M AHR2Mo_DMSO
    ## 17          F1   F   CoMo_DMSO
    ## 18          F1   M   CoMo_DMSO
    ## 19          F1   F   CoMo_DMSO
    ## 20          F1   M   CoMo_DMSO
    ## 21          F1   F   CoMo_DMSO
    ## 22          F1   M   CoMo_DMSO
    ## 23          F1   F   CoMo_DMSO
    ## 24          F1   M   CoMo_DMSO
    ## 25          F1   F AHR2Mo_DMSO
    ## 26          F1   M AHR2Mo_DMSO
    ## 27          F1   F AHR2Mo_DMSO
    ## 28          F1   M AHR2Mo_DMSO
    ## 29          F1   F AHR2Mo_DMSO
    ## 30          F1   M AHR2Mo_DMSO
    ## 31          F1   F AHR2Mo_DMSO
    ## 32          F1   M AHR2Mo_DMSO
    ## 33          F1   F   CoMo_DMSO
    ## 34          F1   M   CoMo_DMSO
    ## 35          F1   F   CoMo_DMSO
    ## 36          F1   M   CoMo_DMSO
    ## 37          F1   F   CoMo_DMSO
    ## 38          F1   M   CoMo_DMSO
    ## 39          F1   F   CoMo_DMSO
    ## 40          F1   M   CoMo_DMSO
    ## 41          F1   F AHR2Mo_DMSO
    ## 42          F1   M AHR2Mo_DMSO
    ## 43          F1   F AHR2Mo_DMSO
    ## 44          F1   M AHR2Mo_DMSO
    ## 45          F1   F AHR2Mo_DMSO
    ## 46          F1   M AHR2Mo_DMSO
    ## 47          F1   F  AHR2Mo_BaP
    ## 48          F1   M  AHR2Mo_BaP
    ## 49          F1   F  AHR2Mo_BaP
    ## 50          F1   M  AHR2Mo_BaP
    ## 51          F1   F  AHR2Mo_BaP
    ## 52          F1   M  AHR2Mo_BaP
    ## 53          F1   F   CoMo_DMSO
    ## 54          F1   M   CoMo_DMSO
    ## 55          F1   F   CoMo_DMSO
    ## 56          F1   M   CoMo_DMSO
    ## 57          F1   F    CoMo_BaP
    ## 58          F1   M    CoMo_BaP
    ## 59          F1   F    CoMo_BaP
    ## 60          F1   M    CoMo_BaP
    ## 61          F1   F AHR2Mo_DMSO
    ## 62          F1   M AHR2Mo_DMSO
    ## 63          F1   F AHR2Mo_DMSO
    ## 64          F1   M AHR2Mo_DMSO
    ## 65          F1   F AHR2Mo_DMSO
    ## 66          F1   M AHR2Mo_DMSO
    ## 67          F1   F AHR2Mo_DMSO
    ## 68          F1   M AHR2Mo_DMSO
    ## 69          F1   F  AHR2Mo_BaP
    ## 70          F1   M  AHR2Mo_BaP
    ## 71          F1   F  AHR2Mo_BaP
    ## 72          F1   M  AHR2Mo_BaP
    ## 73          F1   F    CoMo_BaP
    ## 74          F1   M    CoMo_BaP
    ## 75          F1   F    CoMo_BaP
    ## 76          F1   M    CoMo_BaP
    ## 77          F1   F    CoMo_BaP
    ## 78          F1   M    CoMo_BaP
    ## 79          F1   F AHR2Mo_DMSO
    ## 80          F1   M AHR2Mo_DMSO
    ## 81          F1   F AHR2Mo_DMSO
    ## 82          F1   M AHR2Mo_DMSO
    ## 83          F1   F AHR2Mo_DMSO
    ## 84          F1   F AHR2Mo_DMSO
    ## 85          F1   F  AHR2Mo_BaP
    ## 86          F1   M  AHR2Mo_BaP
    ## 87          F1   F    CoMo_BaP
    ## 88          F1   M    CoMo_BaP
    ## 89          F1   M    CoMo_BaP
    ## 90          F1   M    CoMo_BaP
    ## 91          F1   F    CoMo_BaP
    ## 92          F1   M    CoMo_BaP
    ## 93          F1   F    CoMo_BaP
    ## 94          F1   M    CoMo_BaP
    ## 95          F1   F  AHR2Mo_BaP
    ## 96          F1   M  AHR2Mo_BaP
    ## 97          F1   F  AHR2Mo_BaP
    ## 98          F1   M  AHR2Mo_BaP
    ## 99          F1   F  AHR2Mo_BaP
    ## 100         F1   M  AHR2Mo_BaP
    ## 101         F1   F  AHR2Mo_BaP
    ## 102         F1   M  AHR2Mo_BaP
    ## 103         F1   F    CoMo_BaP
    ## 104         F1   M    CoMo_BaP
    ## 105         F1   F    CoMo_BaP
    ## 106         F1   M    CoMo_BaP
    ## 107         F1   F    CoMo_BaP
    ## 108         F1   M    CoMo_BaP
    ## 109         F1   F    CoMo_BaP
    ## 110         F1   M    CoMo_BaP
    ## 111         F1   F  AHR2Mo_BaP
    ## 112         F1   M  AHR2Mo_BaP
    ## 113         F1   F  AHR2Mo_BaP
    ## 114         F1   M  AHR2Mo_BaP
    ## 115         F1   F  AHR2Mo_BaP
    ## 116         F1   M  AHR2Mo_BaP
    ## 117         F1   F  AHR2Mo_BaP
    ## 118         F1   M  AHR2Mo_BaP
    ## 119         F1   F   CoMo_DMSO
    ## 120         F1   M   CoMo_DMSO
    ## 121         F1   F   CoMo_DMSO
    ## 122         F1   M   CoMo_DMSO
    ## 123         F1   F   CoMo_DMSO
    ## 124         F1   M   CoMo_DMSO
    ## 125         F1   F   CoMo_DMSO
    ## 126         F1   M   CoMo_DMSO
    ## 127         F1   F    CoMo_BaP
    ## 128         F1   M    CoMo_BaP
    ## 129         F1   F    CoMo_BaP
    ## 130         F1   M    CoMo_BaP
    ## 131         F1   F    CoMo_BaP
    ## 132         F1   M    CoMo_BaP
    ## 133         F1   F    CoMo_BaP
    ## 134         F1   M    CoMo_BaP
    ## 135         F1   F  AHR2Mo_BaP
    ## 136         F1   M  AHR2Mo_BaP
    ## 137         F1   F  AHR2Mo_BaP
    ## 138         F1   M  AHR2Mo_BaP
    ## 139         F1   F  AHR2Mo_BaP
    ## 140         F1   M  AHR2Mo_BaP
    ## 141         F1   F  AHR2Mo_BaP
    ## 142         F1   M  AHR2Mo_BaP

    # step AIC to test which model
    testmod_richF1 <- lm(tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F1alphadiv)
    AIC_wholedata_richF1 <- stepAIC(testmod_richF1)
    ## Start:  AIC=831.02
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq   RSS    AIC
    ## - Exposure:Morpholino  1      0.18 44779 829.02
    ## - Exposure:Sex         1      2.70 44781 829.03
    ## - Morpholino:Sex       1    582.25 45361 830.85
    ## <none>                             44779 831.02
    ## 
    ## Step:  AIC=829.02
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq   RSS    AIC
    ## - Exposure:Sex    1      2.75 44781 827.03
    ## - Morpholino:Sex  1    583.34 45362 828.86
    ## <none>                        44779 829.02
    ## 
    ## Step:  AIC=827.03
    ## tax_richness ~ Exposure + Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq   RSS    AIC
    ## - Exposure        1     23.14 44805 825.10
    ## - Morpholino:Sex  1    582.27 45364 826.86
    ## <none>                        44781 827.03
    ## 
    ## Step:  AIC=825.1
    ## tax_richness ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq   RSS    AIC
    ## - Morpholino:Sex  1    582.08 45387 824.93
    ## <none>                        44805 825.10
    ## 
    ## Step:  AIC=824.93
    ## tax_richness ~ Morpholino + Sex
    ## 
    ##              Df Sum of Sq   RSS    AIC
    ## <none>                    45387 824.93
    ## - Sex         1    931.61 46318 825.82
    ## - Morpholino  1   2604.89 47992 830.86

    summary(AIC_wholedata_richF1)
    ## 
    ## Call:
    ## lm(formula = tax_richness ~ Morpholino + Sex, data = F1alphadiv)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -16.977  -6.407  -2.407   2.040 177.023 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        20.853      2.639   7.901 7.61e-13 ***
    ## MorpholinoAHR2Mo   -8.570      3.034  -2.824  0.00543 ** 
    ## SexM                5.125      3.034   1.689  0.09344 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 18.07 on 139 degrees of freedom
    ## Multiple R-squared:  0.074,  Adjusted R-squared:  0.06068 
    ## F-statistic: 5.554 on 2 and 139 DF,  p-value: 0.00478

    # step AIC to test which model
    testmod_shanF1 <- lm(tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F1alphadiv)
    AIC_wholedata_shanF1 <- stepAIC(testmod_shanF1)
    ## Start:  AIC=-251.64
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex         1  0.010389 21.880 -253.58
    ## - Exposure:Morpholino  1  0.112958 21.983 -252.91
    ## - Morpholino:Sex       1  0.200192 22.070 -252.35
    ## <none>                             21.870 -251.64
    ## 
    ## Step:  AIC=-253.58
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino  1   0.11505 21.995 -254.83
    ## - Morpholino:Sex       1   0.19879 22.079 -254.29
    ## <none>                             21.880 -253.58
    ## 
    ## Step:  AIC=-254.83
    ## tax_shannon ~ Exposure + Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Exposure        1  0.082198 22.077 -256.30
    ## - Morpholino:Sex  1  0.207711 22.203 -255.50
    ## <none>                        21.995 -254.83
    ## 
    ## Step:  AIC=-256.3
    ## tax_shannon ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Morpholino:Sex  1    0.2075 22.285 -256.97
    ## <none>                        22.077 -256.30
    ## 
    ## Step:  AIC=-256.97
    ## tax_shannon ~ Morpholino + Sex
    ## 
    ##              Df Sum of Sq    RSS     AIC
    ## <none>                    22.285 -256.97
    ## - Sex         1   0.88301 23.168 -253.46
    ## - Morpholino  1   1.29803 23.583 -250.93

    summary(AIC_wholedata_shanF1)
    ## 
    ## Call:
    ## lm(formula = tax_shannon ~ Morpholino + Sex, data = F1alphadiv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.43377 -0.22664  0.03848  0.24934  1.25926 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       1.49776    0.05848  25.610  < 2e-16 ***
    ## MorpholinoAHR2Mo -0.19131    0.06724  -2.845  0.00511 ** 
    ## SexM              0.15778    0.06723   2.347  0.02034 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4004 on 139 degrees of freedom
    ## Multiple R-squared:  0.09145,    Adjusted R-squared:  0.07838 
    ## F-statistic: 6.996 on 2 and 139 DF,  p-value: 0.001274

# F2

Can F2 gut microbiome taxonomic alpha diveristy be predicted by
exposure, morpholino, and other covariates?

    tax_richness <- specnumber(F2_taxa_otutab)
    tax_shannon <- diversity(F2_taxa_otutab, "shannon", base = exp(1))
    tax_simpson <- diversity(F2_taxa_otutab, "simpson")
    F2alphadiv <- data.frame(tax_richness, tax_shannon, tax_simpson) %>%
      rownames_to_column("SampleID") %>%
      inner_join(metaALL, by = "SampleID") %>%
      dplyr::select(c("SampleID", "tax_richness", "tax_shannon", "tax_simpson", "Exposure", "Morpholino", "Generation", "Sex", "Treatments")) # filter to only important variables we are testing
    F2alphadiv
    ##         SampleID tax_richness tax_shannon tax_simpson Exposure Morpholino
    ## 1       BaP_0001           18   1.8842464   0.8075514      BaP     AHR2Mo
    ## 2       BaP_0002           18   1.8002448   0.7764163      BaP     AHR2Mo
    ## 3       BaP_0025           18   1.5655002   0.7205704     DMSO       CoMo
    ## 4       BaP_0026           14   1.3802800   0.6474444     DMSO       CoMo
    ## 5       BaP_0027           14   1.1409894   0.5557470     DMSO       CoMo
    ## 6       BaP_0028           11   1.1898860   0.6110916     DMSO       CoMo
    ## 7       BaP_0029           13   0.6495473   0.2803063     DMSO       CoMo
    ## 8       BaP_0030           12   1.3477925   0.6170500     DMSO       CoMo
    ## 9       BaP_0031           11   1.5556001   0.7289101     DMSO       CoMo
    ## 10      BaP_0032           12   1.5507609   0.7042284     DMSO       CoMo
    ## 11      BaP_0033           14   1.5563702   0.7260304     DMSO       CoMo
    ## 12      BaP_0034           11   1.5701407   0.7330302     DMSO       CoMo
    ## 13      BaP_0035           16   1.6407887   0.7544659     DMSO       CoMo
    ## 14      BaP_0036           13   1.4443289   0.6730159     DMSO       CoMo
    ## 15      BaP_0061           10   0.9981632   0.4733440     DMSO       CoMo
    ## 16      BaP_0062           16   0.9981814   0.4127875     DMSO       CoMo
    ## 17      BaP_0063            9   1.3654900   0.6727401     DMSO       CoMo
    ## 18      BaP_0064           13   1.1546729   0.4989884     DMSO       CoMo
    ## 19      BaP_0065           14   1.1321478   0.5126835     DMSO       CoMo
    ## 20      BaP_0066           17   0.9270646   0.3812343     DMSO       CoMo
    ## 21      BaP_0067           10   1.4897255   0.7159194     DMSO       CoMo
    ## 22      BaP_0068           18   1.0766016   0.4546126     DMSO       CoMo
    ## 23      BaP_0069           14   1.3341302   0.5927315     DMSO       CoMo
    ## 24      BaP_0070           14   1.6492183   0.7592375     DMSO       CoMo
    ## 25      BaP_0071            9   1.1400649   0.5861441     DMSO       CoMo
    ## 26      BaP_0072           14   1.5567159   0.6801061      BaP       CoMo
    ## 27      BaP_0089           16   1.4822967   0.7098267      BaP     AHR2Mo
    ## 28      BaP_0090           27   1.0888159   0.4600446      BaP     AHR2Mo
    ## 29      BaP_0097           12   1.7816676   0.8024716      BaP       CoMo
    ## 30      BaP_0098           19   1.1707273   0.4844849      BaP     AHR2Mo
    ## 31      BaP_0105           17   1.3074853   0.5955495      BaP       CoMo
    ## 32      BaP_0106           19   1.4666155   0.6206107      BaP       CoMo
    ## 33      BaP_0107           19   1.5213766   0.6935917      BaP       CoMo
    ## 34      BaP_0108           22   1.2831230   0.5178601      BaP       CoMo
    ## 35      BaP_0137           14   1.7080262   0.7740097     DMSO       CoMo
    ## 36      BaP_0138           14   1.6340009   0.7627410     DMSO       CoMo
    ## 37      BaP_0139           11   1.3776681   0.6184465     DMSO       CoMo
    ## 38      BaP_0140           12   1.5548633   0.7400038     DMSO       CoMo
    ## 39      BaP_0141           11   1.3951352   0.6523612      BaP       CoMo
    ## 40      BaP_0142           10   1.5347047   0.7258042      BaP     AHR2Mo
    ## 41      BaP_0143           12   1.6224437   0.7385849      BaP       CoMo
    ## 42      BaP_0144           15   1.7657048   0.7714085      BaP       CoMo
    ## 43      BaP_0441           15   1.7424743   0.7476064      BaP     AHR2Mo
    ## 44      BaP_0442           13   1.1012776   0.5473905      BaP     AHR2Mo
    ## 45      BaP_0457           20   1.3625515   0.6361144      BaP       CoMo
    ## 46      BaP_0458           11   1.3175854   0.6193976     DMSO     AHR2Mo
    ## 47      BaP_0459           14   1.8994412   0.7977018     DMSO     AHR2Mo
    ## 48      BaP_0460           15   1.9257784   0.8210840     DMSO     AHR2Mo
    ## 49      BaP_0461           14   1.8677108   0.8062589     DMSO     AHR2Mo
    ## 50      BaP_0462           14   1.8895180   0.7995876     DMSO     AHR2Mo
    ## 51      BaP_0463           16   1.7669817   0.7734431     DMSO     AHR2Mo
    ## 52      BaP_0464           16   1.6136572   0.7076914     DMSO     AHR2Mo
    ## 53      BaP_0465           16   1.4354739   0.6693849     DMSO       CoMo
    ## 54      BaP_0466           16   1.7284371   0.7423759     DMSO       CoMo
    ## 55      BaP_0467           13   1.2279367   0.5088657     DMSO       CoMo
    ## 56      BaP_0468           19   1.7022848   0.7212924     DMSO       CoMo
    ## 57      BaP_0565           14   1.5365619   0.6651360     DMSO     AHR2Mo
    ## 58      BaP_0566           15   0.6034781   0.2324267     DMSO     AHR2Mo
    ## 59      BaP_0567           14   1.3027313   0.5714217     DMSO     AHR2Mo
    ## 60      BaP_0568           24   1.5510862   0.6411151     DMSO     AHR2Mo
    ## 61      BaP_0569           18   1.8743919   0.7894480     DMSO     AHR2Mo
    ## 62      BaP_0570           25   0.9681408   0.3629287     DMSO     AHR2Mo
    ## 63      BaP_0571           21   1.5594541   0.6860283     DMSO     AHR2Mo
    ## 64      BaP_0572           24   1.1389270   0.4604448     DMSO     AHR2Mo
    ## 65      BaP_0573           18   1.8792881   0.7949445      BaP       CoMo
    ## 66      BaP_0574           24   1.5767703   0.6591339      BaP       CoMo
    ## 67      BaP_0575           20   1.8133528   0.7559692      BaP       CoMo
    ## 68      BaP_0576           17   1.5572708   0.6697547      BaP       CoMo
    ## 69      BaP_0649           12   1.3959252   0.6114784      BaP     AHR2Mo
    ## 70      BaP_0650           12   1.7099016   0.7607734      BaP     AHR2Mo
    ## 71      BaP_0673           17   1.6179174   0.7158344      BaP       CoMo
    ## 72      BaP_0674           17   1.5439392   0.7041643      BaP       CoMo
    ## 73      BaP_0675           22   1.8587988   0.7858640      BaP       CoMo
    ## 74      BaP_0676           15   1.8751114   0.8072020      BaP       CoMo
    ## 75      BaP_0677           14   1.1320326   0.5923136      BaP       CoMo
    ## 76      BaP_0678           20   1.9762084   0.8191111      BaP       CoMo
    ## 77      BaP_0679           12   1.2320345   0.5585303      BaP       CoMo
    ## 78      BaP_0680           21   1.7064263   0.7408483      BaP       CoMo
    ## 79      BaP_0683           23   1.9529767   0.8143508      BaP       CoMo
    ## 80      BaP_0781           15   1.7777810   0.7832698      BaP       CoMo
    ## 81      BaP_0782           17   1.8185205   0.7960746      BaP       CoMo
    ## 82      BaP_0783           10   0.3862890   0.1397118      BaP       CoMo
    ## 83      BaP_0784           18   1.7791698   0.7469604      BaP       CoMo
    ## 84      BaP_0785           14   1.9754552   0.8164374      BaP       CoMo
    ## 85      BaP_0787           16   2.0122591   0.8105631      BaP       CoMo
    ## 86      BaP_0788           20   1.7506400   0.7214839      BaP       CoMo
    ## 87      BaP_0789           20   1.6571081   0.7473505      BaP       CoMo
    ## 88      BaP_0790           34   1.2402774   0.5564423      BaP       CoMo
    ## 89      BaP_0791           29   2.2000657   0.8524004      BaP       CoMo
    ## 90      BaP_0792           22   1.7286617   0.7170064      BaP       CoMo
    ## 91    BaP_0865_1           24   0.8203734   0.3555586      BaP     AHR2Mo
    ## 92    BaP_0866_2           17   1.4670094   0.6498193      BaP     AHR2Mo
    ## 93   BaP_0889_25           17   1.0552701   0.4552095     DMSO     AHR2Mo
    ## 94   BaP_0890_26           30   1.4007176   0.5890197     DMSO     AHR2Mo
    ## 95   BaP_0891_27           32   1.9782598   0.7565867     DMSO     AHR2Mo
    ## 96   BaP_0893_29           18   0.7800389   0.3514297     DMSO     AHR2Mo
    ## 97   BaP_0895_31           27   1.6628248   0.6890247     DMSO     AHR2Mo
    ## 98   BaP_0897_33           17   1.1586724   0.5067864     DMSO     AHR2Mo
    ## 99   BaP_0899_35           22   1.2873921   0.5222688     DMSO     AHR2Mo
    ## 100  BaP_0961_61           31   1.5890008   0.7102563     DMSO     AHR2Mo
    ## 101  BaP_0962_62           17   0.9291130   0.4017830     DMSO     AHR2Mo
    ## 102  BaP_0963_63           19   1.3014286   0.6293603     DMSO     AHR2Mo
    ## 103  BaP_0964_64           21   1.5071439   0.6579072     DMSO     AHR2Mo
    ## 104  BaP_0965_65           17   1.4221702   0.6450400     DMSO     AHR2Mo
    ## 105  BaP_0966_66           21   1.5128708   0.6294222     DMSO     AHR2Mo
    ## 106  BaP_0967_67           29   1.6640119   0.6732072     DMSO     AHR2Mo
    ## 107  BaP_0968_68           24   1.5525491   0.7151785     DMSO     AHR2Mo
    ## 108  BaP_0969_69           30   1.6262657   0.7215035      BaP     AHR2Mo
    ## 109  BaP_0970_70           24   2.0794836   0.8296562      BaP     AHR2Mo
    ## 110  BaP_0971_71           21   1.6091983   0.6796824      BaP     AHR2Mo
    ## 111  BaP_0972_72           23   1.5040702   0.6896538      BaP     AHR2Mo
    ## 112  BaP_1017_81           15   1.4195093   0.6427196      BaP     AHR2Mo
    ## 113  BaP_1018_82           17   1.9647775   0.8034947      BaP     AHR2Mo
    ## 114  BaP_1033_97           22   1.3365637   0.5752839     DMSO     AHR2Mo
    ## 115  BaP_1034_98           26   0.7779672   0.2928568     DMSO     AHR2Mo
    ## 116  BaP_1035_99           26   1.1223753   0.4695991     DMSO     AHR2Mo
    ## 117 BaP_1036_100           22   1.5683346   0.6533996     DMSO     AHR2Mo
    ## 118 BaP_1037_101           24   0.9354576   0.3830119     DMSO     AHR2Mo
    ## 119 BaP_1038_102           26   1.6539237   0.7045263     DMSO     AHR2Mo
    ## 120 BaP_1039_103           29   0.7729422   0.3018056      BaP     AHR2Mo
    ## 121 BaP_1040_104           21   1.1900909   0.5202002      BaP     AHR2Mo
    ## 122 BaP_1041_105           37   1.4698499   0.6284842      BaP     AHR2Mo
    ## 123 BaP_1042_106           22   1.5864927   0.6927114      BaP     AHR2Mo
    ## 124 BaP_1043_107           22   1.1468418   0.5549435      BaP     AHR2Mo
    ## 125 BaP_1044_108           20   0.8643444   0.3350710      BaP     AHR2Mo
    ## 126 BaP_1105_133           20   1.4584559   0.6068385      BaP     AHR2Mo
    ## 127 BaP_1106_134           50   1.4535804   0.5468717      BaP     AHR2Mo
    ## 128 BaP_1107_135           18   1.2057342   0.5192935      BaP     AHR2Mo
    ## 129 BaP_1108_136           44   1.6132923   0.6825551      BaP     AHR2Mo
    ## 130 BaP_1109_137           21   1.3671175   0.6545204      BaP     AHR2Mo
    ## 131 BaP_1110_138           23   1.6646842   0.6975124      BaP     AHR2Mo
    ## 132 BaP_1111_139           21   1.0955069   0.4804426      BaP     AHR2Mo
    ## 133 BaP_1112_140           18   0.9699140   0.4528078      BaP     AHR2Mo
    ## 134 BaP_1113_141           22   1.4602832   0.6990529      BaP     AHR2Mo
    ## 135 BaP_1114_142           41   1.4422312   0.5413321      BaP     AHR2Mo
    ## 136 BaP_1115_143           23   1.7834792   0.7389577      BaP     AHR2Mo
    ## 137 BaP_1116_144           34   2.2274082   0.7912829      BaP     AHR2Mo
    ##     Generation Sex  Treatments
    ## 1           F2   F  AHR2Mo_BaP
    ## 2           F2   F  AHR2Mo_BaP
    ## 3           F2   F   CoMo_DMSO
    ## 4           F2   M   CoMo_DMSO
    ## 5           F2   F   CoMo_DMSO
    ## 6           F2   M   CoMo_DMSO
    ## 7           F2   F   CoMo_DMSO
    ## 8           F2   M   CoMo_DMSO
    ## 9           F2   F   CoMo_DMSO
    ## 10          F2   M   CoMo_DMSO
    ## 11          F2   F   CoMo_DMSO
    ## 12          F2   M   CoMo_DMSO
    ## 13          F2   F   CoMo_DMSO
    ## 14          F2   M   CoMo_DMSO
    ## 15          F2   F   CoMo_DMSO
    ## 16          F2   M   CoMo_DMSO
    ## 17          F2   F   CoMo_DMSO
    ## 18          F2   M   CoMo_DMSO
    ## 19          F2   F   CoMo_DMSO
    ## 20          F2   M   CoMo_DMSO
    ## 21          F2   F   CoMo_DMSO
    ## 22          F2   M   CoMo_DMSO
    ## 23          F2   F   CoMo_DMSO
    ## 24          F2   F   CoMo_DMSO
    ## 25          F2   F   CoMo_DMSO
    ## 26          F2   F    CoMo_BaP
    ## 27          F2   M  AHR2Mo_BaP
    ## 28          F2   M  AHR2Mo_BaP
    ## 29          F2   F    CoMo_BaP
    ## 30          F2   F  AHR2Mo_BaP
    ## 31          F2   F    CoMo_BaP
    ## 32          F2   M    CoMo_BaP
    ## 33          F2   M    CoMo_BaP
    ## 34          F2   F    CoMo_BaP
    ## 35          F2   F   CoMo_DMSO
    ## 36          F2   M   CoMo_DMSO
    ## 37          F2   F   CoMo_DMSO
    ## 38          F2   M   CoMo_DMSO
    ## 39          F2   F    CoMo_BaP
    ## 40          F2   M  AHR2Mo_BaP
    ## 41          F2   M    CoMo_BaP
    ## 42          F2   F    CoMo_BaP
    ## 43          F2   F  AHR2Mo_BaP
    ## 44          F2   F  AHR2Mo_BaP
    ## 45          F2   M    CoMo_BaP
    ## 46          F2   F AHR2Mo_DMSO
    ## 47          F2   M AHR2Mo_DMSO
    ## 48          F2   F AHR2Mo_DMSO
    ## 49          F2   M AHR2Mo_DMSO
    ## 50          F2   F AHR2Mo_DMSO
    ## 51          F2   M AHR2Mo_DMSO
    ## 52          F2   F AHR2Mo_DMSO
    ## 53          F2   F   CoMo_DMSO
    ## 54          F2   M   CoMo_DMSO
    ## 55          F2   F   CoMo_DMSO
    ## 56          F2   M   CoMo_DMSO
    ## 57          F2   M AHR2Mo_DMSO
    ## 58          F2   F AHR2Mo_DMSO
    ## 59          F2   M AHR2Mo_DMSO
    ## 60          F2   F AHR2Mo_DMSO
    ## 61          F2   M AHR2Mo_DMSO
    ## 62          F2   F AHR2Mo_DMSO
    ## 63          F2   M AHR2Mo_DMSO
    ## 64          F2   F AHR2Mo_DMSO
    ## 65          F2   M    CoMo_BaP
    ## 66          F2   F    CoMo_BaP
    ## 67          F2   M    CoMo_BaP
    ## 68          F2   F    CoMo_BaP
    ## 69          F2   M  AHR2Mo_BaP
    ## 70          F2   M  AHR2Mo_BaP
    ## 71          F2   M    CoMo_BaP
    ## 72          F2   F    CoMo_BaP
    ## 73          F2   M    CoMo_BaP
    ## 74          F2   F    CoMo_BaP
    ## 75          F2   M    CoMo_BaP
    ## 76          F2   F    CoMo_BaP
    ## 77          F2   M    CoMo_BaP
    ## 78          F2   F    CoMo_BaP
    ## 79          F2   M    CoMo_BaP
    ## 80          F2   M    CoMo_BaP
    ## 81          F2   F    CoMo_BaP
    ## 82          F2   M    CoMo_BaP
    ## 83          F2   F    CoMo_BaP
    ## 84          F2   M    CoMo_BaP
    ## 85          F2   M    CoMo_BaP
    ## 86          F2   F    CoMo_BaP
    ## 87          F2   M    CoMo_BaP
    ## 88          F2   F    CoMo_BaP
    ## 89          F2   M    CoMo_BaP
    ## 90          F2   F    CoMo_BaP
    ## 91          F2   F  AHR2Mo_BaP
    ## 92          F2   F  AHR2Mo_BaP
    ## 93          F2   M AHR2Mo_DMSO
    ## 94          F2   F AHR2Mo_DMSO
    ## 95          F2   M AHR2Mo_DMSO
    ## 96          F2   F AHR2Mo_DMSO
    ## 97          F2   M AHR2Mo_DMSO
    ## 98          F2   F AHR2Mo_DMSO
    ## 99          F2   M AHR2Mo_DMSO
    ## 100         F2   F AHR2Mo_DMSO
    ## 101         F2   M AHR2Mo_DMSO
    ## 102         F2   F AHR2Mo_DMSO
    ## 103         F2   M AHR2Mo_DMSO
    ## 104         F2   F AHR2Mo_DMSO
    ## 105         F2   M AHR2Mo_DMSO
    ## 106         F2   F AHR2Mo_DMSO
    ## 107         F2   M AHR2Mo_DMSO
    ## 108         F2   F  AHR2Mo_BaP
    ## 109         F2   M  AHR2Mo_BaP
    ## 110         F2   F  AHR2Mo_BaP
    ## 111         F2   M  AHR2Mo_BaP
    ## 112         F2   M  AHR2Mo_BaP
    ## 113         F2   M  AHR2Mo_BaP
    ## 114         F2   F AHR2Mo_DMSO
    ## 115         F2   M AHR2Mo_DMSO
    ## 116         F2   F AHR2Mo_DMSO
    ## 117         F2   M AHR2Mo_DMSO
    ## 118         F2   F AHR2Mo_DMSO
    ## 119         F2   M AHR2Mo_DMSO
    ## 120         F2   F  AHR2Mo_BaP
    ## 121         F2   M  AHR2Mo_BaP
    ## 122         F2   F  AHR2Mo_BaP
    ## 123         F2   M  AHR2Mo_BaP
    ## 124         F2   F  AHR2Mo_BaP
    ## 125         F2   M  AHR2Mo_BaP
    ## 126         F2   F  AHR2Mo_BaP
    ## 127         F2   M  AHR2Mo_BaP
    ## 128         F2   F  AHR2Mo_BaP
    ## 129         F2   M  AHR2Mo_BaP
    ## 130         F2   F  AHR2Mo_BaP
    ## 131         F2   M  AHR2Mo_BaP
    ## 132         F2   F  AHR2Mo_BaP
    ## 133         F2   M  AHR2Mo_BaP
    ## 134         F2   F  AHR2Mo_BaP
    ## 135         F2   M  AHR2Mo_BaP
    ## 136         F2   F  AHR2Mo_BaP
    ## 137         F2   M  AHR2Mo_BaP

    # step AIC to test which model
    testmod_richF2 <- lm(tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F2alphadiv)
    AIC_wholedata_richF2 <- stepAIC(testmod_richF2)
    ## Start:  AIC=504.56
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Exposure:Sex         1     2.223 4920.6 502.63
    ## - Morpholino:Sex       1     5.055 4923.5 502.70
    ## - Exposure:Morpholino  1    57.008 4975.4 504.14
    ## <none>                             4918.4 504.56
    ## 
    ## Step:  AIC=502.63
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Morpholino:Sex       1     4.890 4925.5 500.76
    ## - Exposure:Morpholino  1    56.474 4977.1 502.19
    ## <none>                             4920.6 502.63
    ## 
    ## Step:  AIC=500.76
    ## tax_richness ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Sex                  1     7.703 4933.2 498.98
    ## - Exposure:Morpholino  1    55.666 4981.2 500.30
    ## <none>                             4925.5 500.76
    ## 
    ## Step:  AIC=498.98
    ## tax_richness ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS    AIC
    ## - Exposure:Morpholino  1    56.705 4989.9 498.54
    ## <none>                             4933.2 498.98
    ## 
    ## Step:  AIC=498.54
    ## tax_richness ~ Exposure + Morpholino
    ## 
    ##              Df Sum of Sq    RSS    AIC
    ## <none>                    4989.9 498.54
    ## - Exposure    1    367.53 5357.4 506.28
    ## - Morpholino  1   1162.30 6152.2 525.23

    summary(AIC_wholedata_richF2)
    ## 
    ## Call:
    ## lm(formula = tax_richness ~ Exposure + Morpholino, data = F2alphadiv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.3054  -3.4708  -0.4708   2.5292  26.6946 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       14.1933     0.9331  15.211  < 2e-16 ***
    ## ExposureBaP        3.2774     1.0432   3.142  0.00207 ** 
    ## MorpholinoAHR2Mo   5.8346     1.0444   5.587 1.24e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.102 on 134 degrees of freedom
    ## Multiple R-squared:  0.2312, Adjusted R-squared:  0.2197 
    ## F-statistic: 20.15 on 2 and 134 DF,  p-value: 2.239e-08

    # step AIC to test which model
    testmod_shanF2 <- lm(tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                       data = F2alphadiv)
    AIC_wholedata_shanF2 <- stepAIC(testmod_shanF2)
    ## Start:  AIC=-298.96
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex         1   0.02724 13.979 -300.69
    ## - Morpholino:Sex       1   0.19908 14.151 -299.02
    ## <none>                             13.952 -298.96
    ## - Exposure:Morpholino  1   0.43014 14.382 -296.80
    ## 
    ## Step:  AIC=-300.69
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Morpholino:Sex
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Morpholino:Sex       1   0.20297 14.182 -300.72
    ## <none>                             13.979 -300.69
    ## - Exposure:Morpholino  1   0.43594 14.415 -298.48
    ## 
    ## Step:  AIC=-300.72
    ## tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## <none>                             14.182 -300.72
    ## - Sex                  1   0.26230 14.444 -300.20
    ## - Exposure:Morpholino  1   0.42118 14.603 -298.71

    summary(AIC_wholedata_shanF2)
    ## 
    ## Call:
    ## lm(formula = tax_shannon ~ Exposure + Morpholino + Sex + Exposure:Morpholino, 
    ##     data = F2alphadiv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.27843 -0.19973  0.04635  0.21582  0.73434 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   1.32257    0.06409  20.638  < 2e-16 ***
    ## ExposureBaP                   0.25454    0.08144   3.125  0.00218 ** 
    ## MorpholinoAHR2Mo              0.05061    0.08036   0.630  0.52988    
    ## SexM                          0.08760    0.05607   1.562  0.12057    
    ## ExposureBaP:MorpholinoAHR2Mo -0.22227    0.11226  -1.980  0.04979 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3278 on 132 degrees of freedom
    ## Multiple R-squared:  0.09614,    Adjusted R-squared:  0.06875 
    ## F-statistic:  3.51 on 4 and 132 DF,  p-value: 0.009303

## Beta diversity and dispersion

# Whole dataset

    ## check model with all possible options (have checked sex with ordistep in the past and doesn't come out important to model)
    meta_modALL <- metaALL %>%
      dplyr::select(c("Exposure", "Morpholino", "Generation", "Sex"))
    ALLmod0 <- capscale(taxa_otutab ~ 1, meta_modALL, distance = "bray")  # Model with intercept only
    ALLmod1 <- capscale(taxa_otutab ~ . + Exposure:Morpholino + Morpholino:Exposure + Generation:Exposure + Generation:Morpholino + Sex: Generation + Sex:Exposure + Sex:Morpholino, meta_modALL, distance = "bray")  # Model with all explanatory variables
    tax_ordiALL <- ordistep(ALLmod0, scope = formula(ALLmod1)) # this determines what the best model is to run RDA on
    tax_ordiALL

    tax_dm <- vegdist(taxa_otutab, method = "bray")
    PermExpandMod_tax <- adonis2(tax_dm ~ Generation + Exposure + Sex + Generation:Exposure, 
                                 data = metaALL)
    PermExpandMod_tax
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = tax_dm ~ Generation + Exposure + Sex + Generation:Exposure, data = metaALL)
    ##                      Df SumOfSqs      R2      F Pr(>F)    
    ## Generation            2    2.304 0.03233 7.1571  0.001 ***
    ## Exposure              1    0.835 0.01171 5.1851  0.002 ** 
    ## Sex                   1    0.484 0.00680 3.0084  0.026 *  
    ## Generation:Exposure   2    1.001 0.01405 3.1092  0.008 ** 
    ## Residual            414   66.635 0.93512                  
    ## Total               420   71.258 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # calc dispersion
    metaALL$group <- as.factor(paste0(metaALL$Treatments, "_", metaALL$Generation))
    tax.disper <- betadisper(tax_dm, metaALL$group)

    # make file for graphing
    tax.Disper1 <- data.frame(tax.disper$distances)
    colnames(tax.Disper1) <- "dispers"
    tax.Disper2 <- tax.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaALL, 
                 by = "SampleID")

    ## NOT WORKING BECAUSE variable lengths differ (found for 'Generation')??
    ## linear model
    # step AIC to test model
    testmod_bdisp <- lm(tax.disper$distances ~ Generation + Exposure + Morpholino + Sex + Generation:Exposure + Exposure:Morpholino + Generation:Morpholino,
                        data = tax.Disper2)
    AIC_bdispALL <- stepAIC(testmod_bdisp)
    ## Start:  AIC=-1630.86
    ## tax.disper$distances ~ Generation + Exposure + Morpholino + Sex + 
    ##     Generation:Exposure + Exposure:Morpholino + Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Generation:Exposure    2  0.026324 8.3291 -1633.5
    ## - Exposure:Morpholino    1  0.008741 8.3116 -1632.4
    ## - Sex                    1  0.035025 8.3378 -1631.1
    ## - Generation:Morpholino  2  0.078126 8.3809 -1630.9
    ## <none>                               8.3028 -1630.9
    ## 
    ## Step:  AIC=-1633.53
    ## tax.disper$distances ~ Generation + Exposure + Morpholino + Sex + 
    ##     Exposure:Morpholino + Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino    1  0.009635 8.3388 -1635.0
    ## - Sex                    1  0.035446 8.3646 -1633.7
    ## - Generation:Morpholino  2  0.077807 8.4069 -1633.6
    ## <none>                               8.3291 -1633.5
    ## 
    ## Step:  AIC=-1635.04
    ## tax.disper$distances ~ Generation + Exposure + Morpholino + Sex + 
    ##     Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Exposure               1  0.004376 8.3431 -1636.8
    ## - Sex                    1  0.035205 8.3740 -1635.3
    ## - Generation:Morpholino  2  0.077780 8.4165 -1635.1
    ## <none>                               8.3388 -1635.0
    ## 
    ## Step:  AIC=-1636.82
    ## tax.disper$distances ~ Generation + Morpholino + Sex + Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Sex                    1  0.034362 8.3775 -1637.1
    ## - Generation:Morpholino  2  0.077381 8.4205 -1636.9
    ## <none>                               8.3431 -1636.8
    ## 
    ## Step:  AIC=-1637.09
    ## tax.disper$distances ~ Generation + Morpholino + Generation:Morpholino
    ## 
    ##                         Df Sum of Sq    RSS     AIC
    ## - Generation:Morpholino  2  0.074859 8.4524 -1637.3
    ## <none>                               8.3775 -1637.1
    ## 
    ## Step:  AIC=-1637.35
    ## tax.disper$distances ~ Generation + Morpholino
    ## 
    ##              Df Sum of Sq    RSS     AIC
    ## - Generation  2  0.046343 8.4987 -1639.0
    ## <none>                    8.4524 -1637.3
    ## - Morpholino  1  0.139969 8.5923 -1632.4
    ## 
    ## Step:  AIC=-1639.04
    ## tax.disper$distances ~ Morpholino
    ## 
    ##              Df Sum of Sq    RSS     AIC
    ## <none>                    8.4987 -1639.0
    ## - Morpholino  1   0.13649 8.6352 -1634.3

    # save in object
    lm_betadisper <- summary(AIC_bdispALL)
    lm_betadisper
    ## 
    ## Call:
    ## lm(formula = tax.disper$distances ~ Morpholino, data = tax.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.26920 -0.10580 -0.01735  0.08660  0.44792 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.343350   0.009875  34.770  < 2e-16 ***
    ## MorpholinoAHR2Mo 0.036013   0.013883   2.594  0.00982 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1424 on 419 degrees of freedom
    ## Multiple R-squared:  0.01581,    Adjusted R-squared:  0.01346 
    ## F-statistic: 6.729 on 1 and 419 DF,  p-value: 0.009818

# F0

    ## check model with all possible options (have checked sex with ordistep in the past and doesn't come out important to model)
    meta_modF0 <- metaF0 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex"))
    F0mod0 <- capscale(F0_taxa_otutab ~ 1, meta_modF0, distance = "bray")  # Model with intercept only
    F0mod1 <- capscale(F0_taxa_otutab ~ . + Exposure:Morpholino + Morpholino:Exposure + Sex:Exposure + Sex:Morpholino, meta_modF0, distance = "bray")  # Model with all explanatory variables
    tax_ordiF0 <- ordistep(F0mod0, scope = formula(F0mod1)) # this determines what the best model is to run RDA on
    tax_ordiF0

    F0tax_dm <- vegdist(F0_taxa_otutab, method = "bray")
    PermExpandMod_taxF0 <- adonis2(F0tax_dm ~ Morpholino + Sex, 
                                 data = metaF0)
    PermExpandMod_taxF0
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = F0tax_dm ~ Morpholino + Sex, data = metaF0)
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## Morpholino   1   1.3780 0.06135 9.2937  0.001 ***
    ## Sex          1   0.4732 0.02107 3.1912  0.017 *  
    ## Residual   139  20.6102 0.91758                  
    ## Total      141  22.4614 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # calc dispersion
    F0tax.disper <- betadisper(F0tax_dm, metaF0$Treatments)

    # make file for graphing
    F0tax.Disper1 <- data.frame(F0tax.disper$distances)
    colnames(F0tax.Disper1) <- "dispers"
    F0tax.Disper2 <- F0tax.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF0, 
                 by = "SampleID")

    ## linear model
    # step AIC to test model
    testmod_bdispF0 <- lm(F0tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                        data = F0tax.Disper2)
    AIC_bdispF0 <- stepAIC(testmod_bdispF0)
    ## Start:  AIC=-556.47
    ## F0tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino  1  0.008619 2.5646 -557.99
    ## - Exposure:Sex         1  0.018097 2.5740 -557.47
    ## <none>                             2.5560 -556.47
    ## - Morpholino:Sex       1  0.057770 2.6137 -555.30
    ## 
    ## Step:  AIC=-557.99
    ## F0tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Sex + 
    ##     Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex    1  0.017448 2.5820 -559.03
    ## <none>                        2.5646 -557.99
    ## - Morpholino:Sex  1  0.055479 2.6200 -556.95
    ## 
    ## Step:  AIC=-559.03
    ## F0tax.disper$distances ~ Exposure + Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Exposure        1  0.000657 2.5827 -560.99
    ## <none>                        2.5820 -559.03
    ## - Morpholino:Sex  1  0.056350 2.6384 -557.96
    ## 
    ## Step:  AIC=-560.99
    ## F0tax.disper$distances ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## <none>                        2.5827 -560.99
    ## - Morpholino:Sex  1  0.055935 2.6386 -559.95

    # save in object
    lm_betadisperF0 <- summary(AIC_bdispF0)
    lm_betadisperF0
    ## 
    ## Call:
    ## lm(formula = F0tax.disper$distances ~ Morpholino + Sex + Morpholino:Sex, 
    ##     data = F0tax.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.20375 -0.10309 -0.02493  0.08446  0.33215 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.33874    0.02249  15.062   <2e-16 ***
    ## MorpholinoAHR2Mo       0.01301    0.03226   0.403   0.6874    
    ## SexM                  -0.01798    0.03250  -0.553   0.5809    
    ## MorpholinoAHR2Mo:SexM  0.07943    0.04594   1.729   0.0861 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1368 on 138 degrees of freedom
    ## Multiple R-squared:  0.0623, Adjusted R-squared:  0.04192 
    ## F-statistic: 3.056 on 3 and 138 DF,  p-value: 0.03052

    # pre plot data frame
    smry_rdaF0 <- summary(capscale(formula = F0_taxa_otutab ~ Morpholino + Sex, 
                                   data = metaF0,
                                   distance = "bray"))
    F0_PC1  <- data.frame(smry_rdaF0$sites[,1:2]) %>%      # these are the x, y coordinates for the sample points (e.g. BaP_####)
      rownames_to_column("SampleID") %>%
      inner_join(dplyr::select(metaF0, c(SampleID, Treatments, Exposure, Morpholino, Sex)), by = "SampleID") %>%
      column_to_rownames("SampleID")
    F0_PC2  <- data.frame(smry_rdaF0$biplot)     # x and y coordinates for the vectors
    # reorder x-axis
    F0_PC1$Treatments <- factor(F0_PC1$Treatments, 
                            levels = c("CoMo_DMSO",
                                       "CoMo_BaP", 
                                       "AHR2Mo_DMSO",
                                       "AHR2Mo_BaP"))

    rda_plotF0 <- ggplot(F0_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = Treatments, shape = Sex), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = F0_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = F0_PC2, 
                aes(x = CAP1, y = CAP2, label = c("AHR2 Morpholino", "Sex (Male)")), 
                color = "black", size = 4, hjust = c(0.55, 0.5), vjust = c(-0.5, -0.45)) +
      labs(x = paste0("CAP1 (",round(100*smry_rdaF0$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rdaF0$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "PuOr", direction = -1,
                         name = "Treatments",
                         labels = c("AHR2Mo - / BaP -",
                                   "AHR2Mo - / BaP +",
                                   "AHR2Mo + / BaP -",
                                   "AHR2Mo + / BaP +")) +
      scale_shape_discrete(labels = c("Female", "Male")) +
      theme(text = element_text(size = 20))

    rda_plotF0

<img src="zfBaP_16Sdiversity_analysis_files/figure-markdown_strict/figure of F0 beta diversity bray-curtis rda-1.png" width="98%" height="98%" />

# F1

    ## check model with all possible options (have checked sex with ordistep in the past and doesn't come out important to model)
    meta_modF1 <- metaF1 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex"))
    F1mod0 <- capscale(F1_taxa_otutab ~ 1, meta_modF1, distance = "bray")  # Model with intercept only
    F1mod1 <- capscale(F1_taxa_otutab ~ . + Exposure:Morpholino + Morpholino:Exposure + Sex:Exposure + Sex:Morpholino, meta_modF1, distance = "bray")  # Model with all explanatory variables
    tax_ordiF1 <- ordistep(F1mod0, scope = formula(F1mod1)) # this determines what the best model is to run RDA on
    tax_ordiF1

    F1tax_dm <- vegdist(F1_taxa_otutab, method = "bray")
    PermExpandMod_taxF1 <- adonis2(F1tax_dm ~ Exposure + Sex, 
                                 data = metaF1)
    PermExpandMod_taxF1
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = F1tax_dm ~ Exposure + Sex, data = metaF1)
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## Exposure   1   1.0658 0.04367 6.4650  0.001 ***
    ## Sex        1   0.4256 0.01744 2.5815  0.034 *  
    ## Residual 139  22.9149 0.93889                  
    ## Total    141  24.4063 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # calc dispersion
    F1tax.disper <- betadisper(F1tax_dm, metaF1$Treatments)

    # make file for graphing
    F1tax.Disper1 <- data.frame(F1tax.disper$distances)
    colnames(F1tax.Disper1) <- "dispers"
    F1tax.Disper2 <- F1tax.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF1, 
                 by = "SampleID")

    ## linear model
    # step AIC to test model
    testmod_bdispF1 <- lm(F1tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                        data = F1tax.Disper2)
    AIC_bdispF1 <- stepAIC(testmod_bdispF1)
    ## Start:  AIC=-536.34
    ## F1tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Morpholino:Sex       1  0.001011 2.9463 -538.29
    ## - Exposure:Sex         1  0.019969 2.9652 -537.38
    ## <none>                             2.9453 -536.34
    ## - Exposure:Morpholino  1  0.047608 2.9929 -536.06
    ## 
    ## Step:  AIC=-538.29
    ## F1tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Exposure:Sex
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex         1  0.019831 2.9661 -539.34
    ## <none>                             2.9463 -538.29
    ## - Exposure:Morpholino  1  0.048052 2.9943 -537.99
    ## 
    ## Step:  AIC=-539.34
    ## F1tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Sex                  1  0.014992 2.9811 -540.62
    ## <none>                             2.9661 -539.34
    ## - Exposure:Morpholino  1  0.049876 3.0160 -538.97
    ## 
    ## Step:  AIC=-540.62
    ## F1tax.disper$distances ~ Exposure + Morpholino + Exposure:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## <none>                             2.9811 -540.62
    ## - Exposure:Morpholino  1  0.049831 3.0309 -540.27

    # save in object
    lm_betadisperF1 <- summary(AIC_bdispF1)
    lm_betadisperF1
    ## 
    ## Call:
    ## lm(formula = F1tax.disper$distances ~ Exposure + Morpholino + 
    ##     Exposure:Morpholino, data = F1tax.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.25896 -0.09227 -0.02793  0.08959  0.45769 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   0.39873    0.02450  16.277   <2e-16 ***
    ## ExposureBaP                  -0.04656    0.03464  -1.344    0.181    
    ## MorpholinoAHR2Mo             -0.03777    0.03515  -1.075    0.284    
    ## ExposureBaP:MorpholinoAHR2Mo  0.07495    0.04935   1.519    0.131    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.147 on 138 degrees of freedom
    ## Multiple R-squared:  0.01751,    Adjusted R-squared:  -0.003852 
    ## F-statistic: 0.8196 on 3 and 138 DF,  p-value: 0.4851

# F2

    ## check model with all possible options (have checked sex with ordistep in the past and doesn't come out important to model)
    meta_modF2 <- metaF2 %>%
      dplyr::select(c("Exposure", "Morpholino", "Sex"))
    F2mod0 <- capscale(F2_taxa_otutab ~ 1, meta_modF2, distance = "bray")  # Model with intercept only
    F2mod1 <- capscale(F2_taxa_otutab ~ . + Exposure:Morpholino + Morpholino:Exposure + Sex:Exposure + Sex:Morpholino, meta_modF2, distance = "bray")  # Model with all explanatory variables
    tax_ordiF2 <- ordistep(F2mod0, scope = formula(F2mod1)) # this determines what the best model is to run RDA on
    tax_ordiF2

    F2tax_dm <- vegdist(F2_taxa_otutab, method = "bray")
    PermExpandMod_taxF2 <- adonis2(F2tax_dm ~ Morpholino + Exposure + Morpholino:Exposure, 
                                 data = metaF2)
    PermExpandMod_taxF2
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = F2tax_dm ~ Morpholino + Exposure + Morpholino:Exposure, data = metaF2)
    ##                      Df SumOfSqs      R2      F Pr(>F)    
    ## Morpholino            1   1.3338 0.06051 9.0837  0.001 ***
    ## Exposure              1   0.5963 0.02705 4.0611  0.004 ** 
    ## Morpholino:Exposure   1   0.5847 0.02652 3.9821  0.012 *  
    ## Residual            133  19.5289 0.88592                  
    ## Total               136  22.0437 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # calc dispersion
    F2tax.disper <- betadisper(F2tax_dm, metaF2$Treatments)

    # make file for graphing
    F2tax.Disper1 <- data.frame(F2tax.disper$distances)
    colnames(F2tax.Disper1) <- "dispers"
    F2tax.Disper2 <- F2tax.Disper1 %>%
      rownames_to_column("SampleID") %>%
      inner_join(y = metaF2, 
                 by = "SampleID")

    ## linear model
    # step AIC to test model
    testmod_bdispF2 <- lm(F2tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + Sex:Exposure + Sex:Morpholino, 
                        data = F2tax.Disper2)
    AIC_bdispF2 <- stepAIC(testmod_bdispF2)
    ## Start:  AIC=-530.94
    ## F2tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Morpholino + 
    ##     Sex:Exposure + Sex:Morpholino
    ## 
    ##                       Df Sum of Sq    RSS     AIC
    ## - Exposure:Morpholino  1  0.000025 2.5660 -532.93
    ## - Exposure:Sex         1  0.023688 2.5897 -531.68
    ## <none>                             2.5660 -530.94
    ## - Morpholino:Sex       1  0.039150 2.6052 -530.86
    ## 
    ## Step:  AIC=-532.93
    ## F2tax.disper$distances ~ Exposure + Morpholino + Sex + Exposure:Sex + 
    ##     Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Exposure:Sex    1  0.023665 2.5897 -533.68
    ## <none>                        2.5660 -532.93
    ## - Morpholino:Sex  1  0.039228 2.6053 -532.86
    ## 
    ## Step:  AIC=-533.68
    ## F2tax.disper$distances ~ Exposure + Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## - Exposure        1  0.031565 2.6213 -534.02
    ## <none>                        2.5897 -533.68
    ## - Morpholino:Sex  1  0.040776 2.6305 -533.54
    ## 
    ## Step:  AIC=-534.02
    ## F2tax.disper$distances ~ Morpholino + Sex + Morpholino:Sex
    ## 
    ##                  Df Sum of Sq    RSS     AIC
    ## <none>                        2.6213 -534.02
    ## - Morpholino:Sex  1  0.042553 2.6638 -533.81

    # save in object
    lm_betadisperF2 <- summary(AIC_bdispF2)
    lm_betadisperF2
    ## 
    ## Call:
    ## lm(formula = F2tax.disper$distances ~ Morpholino + Sex + Morpholino:Sex, 
    ##     data = F2tax.Disper2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.25480 -0.10884 -0.00589  0.08368  0.37313 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.33132    0.02408  13.761   <2e-16 ***
    ## MorpholinoAHR2Mo       0.09495    0.03357   2.828   0.0054 ** 
    ## SexM                  -0.02255    0.03486  -0.647   0.5188    
    ## MorpholinoAHR2Mo:SexM -0.07063    0.04807  -1.469   0.1441    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1404 on 133 degrees of freedom
    ## Multiple R-squared:  0.09773,    Adjusted R-squared:  0.07738 
    ## F-statistic: 4.802 on 3 and 133 DF,  p-value: 0.003298
