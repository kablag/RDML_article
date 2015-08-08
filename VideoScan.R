data(C54)
descr <- data.frame(
    fdata.name = c("D1", "D2", "D3"),
    exp.id = c("exp1", "exp1", "exp1"),
    run.id = c("run1", "run1", "run1"),
    react.id = c(1, 2, 3),
    sample = c("Stock cDNA", "1/10 cDNA", "1/100 cDNA"),
    target = c("MLC-2v", "MLC-2v", "MLC-2v"),
    target.dyeId = c("Cy5", "Cy5", "Cy5"),
    stringsAsFactors = FALSE
)
video.scan <- RDML$new()
video.scan$SetFData(C54, descr)

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
video.scan$documentation <- list(
  documentationType$new(
    idType$new("Roediger et al. 2013"),
    "A Highly Versatile Microscope Imaging Technology Platform for the Multiplex Real-Time Detection of Biomolecules and Autoimmune Antibodies. S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder. Advances in Biochemical Bioengineering/Biotechnology. 133:33–74, 2013. http://www.ncbi.nlm.nih.gov/pubmed/22437246"
  )
)
video.scan$sample$`Stock cDNA`$description <- "Input stock cDNA was used undiluted (D1)"
video.scan$sample$`1/10 cDNA`$description <- "1/1000 diluted in A. bidest"
video.scan$sample$`1/100 cDNA`$description <- "1/1000000 diluted in A. bidest"
video.scan$target$`MLC-2v`$xRef <- list(
  xRefType$new("uniprot",
               "P10916")
)
video.scan$experiment$exp1$description <- 
  "qPCR Experiment for the amplification of MLC-2v using the VideoScan heating/cooling-unit. The aim was to amplify MLC-2v in the VideoScan and to monitor with a hydrolysis probe for MLC-2v. The primer sequences for MLC-2v were taken from Roediger et al. (2013). The amplification was detected in solution of the 1 HCU (see Roediger et al. 2013 for details). A 20 micro L PCR reaction was composed of 250 nM primer (forward and reverse), 1x Maxima Probe qPCR Master Mix (Fermentas), 1 micro L template (MLC-2v amplification product in different dilutions), 50 nM hydrolysis probe probe for MLC-2v and A. bidest. During the amplification, fluorescence was measured at 59.5 degree Celsius. The Cy5 channel was used to monitor the MLC-2v specific hydrolysis probe. Input stock cDNA was used undiluted (D1). D2 was 1/1000 and D3 1/1000000 diluted in A. bidest. The D1, D2, and D3 have different numbers measure points and D2 contains a missing value at cycle 37."
