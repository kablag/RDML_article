library(RDML)
library(chipPCR)
data(C54)

# Create a data frame of metadata
descr <- data.frame(
  fdata.name = c("D1", "D2", "D3"),
  exp.id = c("exp1", "exp1", "exp1"),
  run.id = c("run1", "run2", "run3"),
  react.id = c(1, 2, 3),
  sample = c("Stock cDNA", "1/10 cDNA", "1/100 cDNA"),
  target = c("MLC-2v", "MLC-2v", "MLC-2v"),
  target.dyeId = c("Cy5", "Cy5", "Cy5"),
  stringsAsFactors = FALSE
)

# Create the RDML object
video.scan <- RDML$new()

# Add metadata to the RDML object
video.scan$SetFData(C54, descr)

# Add experimentator information
video.scan$experimenter <- 
  list(
    experimenterType$new(
      idType$new("SR"),
      "Stefan",
      "Roediger",
      "stefan.roediger@hs-lausitz.de"
    ),
    experimenterType$new(
      idType$new("CD"),
      "Claudia",
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
          "Biochemical Bioengineering/Biotechnology. 133:33–74, 2013.",
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
video.scan$sample$`Stock cDNA`$description <- "Input stock cDNA was used undiluted (D1)"
video.scan$sample$`Stock cDNA`$cdnaSynthesisMethod <- cdna
video.scan$sample$`1/10 cDNA`$description <- "1/1000 diluted in A. bidest"
video.scan$sample$`1/10 cDNA`$cdnaSynthesisMethod <- cdna
video.scan$sample$`1/100 cDNA`$description <- "1/1000000 diluted in A. bidest"
video.scan$sample$`1/100 cDNA`$cdnaSynthesisMethod <- cdna

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
video.scan$experiment$exp1$run$run1$thermalCyclingConditions <- idReferencesType$new("Amplification")
video.scan$experiment$exp1$run$run2$thermalCyclingConditions <- idReferencesType$new("Amplification")
video.scan$experiment$exp1$run$run3$thermalCyclingConditions <- idReferencesType$new("Amplification")

#visualise RDML object
video.scan$AsDendrogram()

