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
