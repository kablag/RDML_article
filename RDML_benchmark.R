library(RDML)
library(chipPCR)
library(MBmca)
library(ggplot2)
library(qpcR)
data(C54)
library(microbenchmark)

do_rdml <- function(n_exps, exp_data) {
  # Add method to RDML class for VideoScan data preprocessing and Cq calculation
  # Argument 'last.cycle' - cycle limit
  # Argument 'bg.range' - bg.range of CPP function
  
  # Create a data frame of metadata 
  descr <- data.frame(
    fdata.name = paste0("D", 1L:n_exps),
    exp.id = rep("exp1", n_exps),
    run.id = paste0("run", 1L:n_exps),
    react.id = rep(1, n_exps),
    sample = paste0("D", 1L:n_exps),
    target = rep("MLC-2v", n_exps),
    target.dyeId = rep("Cy5", n_exps),
    stringsAsFactors = FALSE
  )
  
  # Create an empty RDML object
  video.scan <- RDML$new()
  
  # Add fluorescence data and metadata to the RDML object from a given source
  # Fare the sake of easyness we use the C54 dataset from the chipPCR package.
  dat <- cbind(exp_data[, 1], exp_data[, sample(2L:ncol(exp_data), n_exps, replace = TRUE)])
  colnames(dat) <- c("cyc", paste0("D", 1L:n_exps))

  video.scan$SetFData(dat, descr)
  
  # Add experimentator information
  video.scan$experimenter <- 
    list(
      experimenterType$new(
        idType$new("SR"),
        "Stefan",
        "Roediger",
        "stefan.roediger@b-tu.de"
      ),
      experimenterType$new(
        idType$new("CD"),
        "Claudi",
        "Deutschmann"
      )
    )
  
  # Add a reference to documentation
  video.scan$documentation <- list(
    documentationType$new(
      idType$new("Roediger et al. 2013"),
      paste("A Highly Versatile Microscope Imaging Technology Platform for the Multiplex",
            "Real-Time Detection of Biomolecules and Autoimmune Antibodies. S. Roediger,",
            "P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M.",
            "Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder. Advances in",
            "Biochemical Bioengineering/Biotechnology. 133:33-74, 2013.",
            "http://www.ncbi.nlm.nih.gov/pubmed/22437246")
    )
  )
  
  cdna <- 
    cdnaSynthesisMethodType$new(
      enzyme = "SuperScript II",
      primingMethod =
        primingMethodType$new("oligo-dt"),
      dnaseTreatment = TRUE
    )
  
  for(i in paste0("D", 1L:n_exps)) {
    video.scan$sample[[i]]$description <- "Input stock cDNA was used undiluted (D1)"
    video.scan$sample[[i]]$cdnaSynthesisMethod <- cdna
  }
  
  
  
  
  video.scan$target$`MLC-2v`$xRef <- list(
    xRefType$new("uniprot",
                 "P10916")
  )
  video.scan$target$`MLC-2v`$sequences <- 
    sequencesType$new(
      forwardPrimer <- oligoType$new(
        sequence = "ACAGGGATGGCTTCATTGAC"),
      reversePrimer <- oligoType$new(
        sequence = "ATGCGTTGAGAATGGTTTCC"),
      probe1 <- oligoType$new(
        threePrimeTag = "Atto647N",
        sequence = "CAGGGTCCGCTCCCTTAAGTTTCTCC",
        fivePrimeTag = "BHQ2")
    )
  
  tcc <- 
    thermalCyclingConditionsType$new(
      idType$new("Amplification"),
      experimenter = list(
        idReferencesType$new("SR"),
        idReferencesType$new("CD")
      ),
      step = 
        list(
          stepType$new(
            nr = 1,
            temperature = temperatureType$new(95,
                                              600)
          ),
          stepType$new(
            nr = 2,
            temperature = temperatureType$new(95,
                                              40)
          ),
          stepType$new(
            nr = 3,
            temperature = temperatureType$new(58.5,
                                              90)
          ),
          stepType$new(
            nr = 4,
            temperature = temperatureType$new(68.5,
                                              90)
          ),
          stepType$new(
            nr = 5,
            loop = loopType$new(goto = 2,
                                repeat.n = 49)
          )
        )
    )
  video.scan$thermalCyclingConditions <- list(
    tcc
  )
  
  #add description of the experiment
  video.scan$experiment$exp1$description <- 
    paste("The aim was to amplify MLC-2v in the VideoScan and to monitor with a",
          "hydrolysis probe for MLC-2v. The primer sequences for MLC-2v were taken",
          "from Roediger et al. (2013). The amplification was detected in solution of",
          "the 1 HCU (see Roediger et al. 2013 for details). A 20 micro L PCR reaction",
          "was composed of 250 nM primer (forward and reverse), 1x Maxima Probe qPCR",
          "Master Mix (Fermentas), 1 micro L template (MLC-2v amplification product in",
          "different dilutions), 50 nM hydrolysis probe probe for MLC-2v and A.",
          "bidest. During the amplification, fluorescence was measured at 59.5 degree",
          "Celsius. The Cy5 channel was used to monitor the MLC-2v specific hydrolysis",
          "probe. Input stock cDNA was used undiluted (D1). D2 was 1/1000 and D3",
          "1/1000000 diluted in A. bidest. The D1, D2, and D3 have different numbers",
          "measure points and D2 contains a missing value at cycle 37.")
  
  for(i in 1L:n_exps)
    video.scan$experiment$exp1$run[[paste0("run", i)]]$thermalCyclingConditions
  
  video.scan
}

benchmark <- microbenchmark(do_rdml(1, reps384), 
                            do_rdml(8, reps384),
                            do_rdml(20, reps384), 
                            do_rdml(30, reps384),
                            do_rdml(32, reps384),
                            do_rdml(50, reps384),
                            do_rdml(96, reps384),
                            do_rdml(100, reps384),
                            do_rdml(300, reps384),
                            do_rdml(384, reps384),
                            do_rdml(765, reps384),
                            do_rdml(1000, reps384), 
                            times = 100)
load("benchmark_win.RData")
swin <- summary(benchmark_win, unit = "s")
#number of reactions
nr <- as.vector(na.omit(as.numeric(unlist(strsplit(unlist(strsplit(levels(swin[["expr"]]), "(", fixed = TRUE)), ",")))))
bench_df <- data.frame(nr = nr, swin[, c("min", "lq", "mean", "median", "uq", "max")], 
                       os = rep("Windows", length(nr)))

ggplot(bench_df, aes(x = nr, y = mean, colour = os)) +
  geom_point(size = 3) +
  scale_x_continuous("Number of experiments") +
  scale_y_continuous("Time [s]") +
  scale_color_discrete("Operating\nsystem") +
  #stat_smooth(method = "lm", se = FALSE) + #not really linear
  stat_smooth(method = "loess")
