
Miao.collectGenotype <- function() {
    patient <- Miao.loadPatientDetail()
    neoDetail <- Miao.loadNeoDetail()

    patient <- patient[,c("patient_id", "Tumor_Sample_Barcode")]
    neoDetail <- neoDetail[,c("Tumor_Sample_Barcode", "HLA")]
    neoDetail <- merge(patient, neoDetail, by = "Tumor_Sample_Barcode")

    byHLA <-
        by(neoDetail,
           neoDetail$patient_id,
           function(x) paste(sort(unique(x$HLA)), collapse = " "))

    genoFrame <-
        data.frame(patient_id = names(byHLA), Genotype = as.character(byHLA))

    genoFrame$Genotype <-
        gsub("HLA-", "", genoFrame$Genotype)
    
    genoFrame$Genotype <-
        gsub(":", "", genoFrame$Genotype)

    genoFrame
}

Miao.compileRawData <- function() {
    dataDir <- Miao.dataDir()

    table5  <- Miao.loadTSV(file.path(dataDir, "Miao_SupTable5.tsv"))
    table10 <- Miao.loadTSV(file.path(dataDir, "Miao_SupTable10.tsv"))

    JamIO.save(table5,  file.path(dataDir, "Miao_SupTable5.RData"))
    JamIO.save(table10, file.path(dataDir, "Miao_SupTable10.RData"))
}

Miao.compileMaster <- function() {
    master <- Miao.loadPatientDetail()

    genoFrame <- Miao.collectGenotype()
    presentFrame <- Miao.loadPatientPresentation()

    master <- merge(master, genoFrame, by = "patient_id")
    master <- merge(master, presentFrame, by = "patient_id")
    
    mutDetail <- Miao.loadMutDetail()
    mutBurden <- Miao.computeTMB(mutDetail)

    master <- merge(master, mutBurden, by = "pair_id")

    ## Remove these cancer types with only one observation...
    master <- subset(master, !(cancer_type %in% c("Anal", "Sarcoma")))

    ## Compute z-scores for HLA presentation and tumor mutational
    ## burden...
    master$zAGE <- Miao.zscore(master$age_start_io)
    master$zHLA <- Miao.zscore(master$presentRate.mean)

    master$logTMB <- log(master$nonSilentCount)
    master$zTMB   <- Miao.zscore(master$logTMB)

    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "logTMB", "zTMB.By"))

    ## Events have occurred if the observation is censored...
    master$os_event <- master$os_censor
    master$pfs_event <- master$pfs_censor

    ## Melanoma will be the reference...
    master$Bladder <- as.numeric(master$cancer_type == "Bladder")
    master$HNSCC   <- as.numeric(master$cancer_type == "HNSCC")
    master$Lung    <- as.numeric(master$cancer_type == "Lung")

    ## "Anti-CTLA-4" will be the reference...
    master$PD1  <- as.numeric(master$drug_type == "anti-PD-1/anti-PD-L1")
    master$Both <- as.numeric(master$drug_type == "anti-CTLA-4 + anti-PD-1/PD-L1")

    master
}

Miao.computeTMB <- function(mutDetail) {
    mutDetail <- subset(mutDetail, clonal_dm == 1)

    aggFunc <- function(slice) {
        data.frame(pair_id = slice$pair_id[1],
                   missenseCount = sum(slice$Variant_Classification == "Missense_Mutation"),
                   nonSilentCount = sum(slice$Variant_Classification %in%
                                        c("Nonstop_Mutation",
                                          "In_Frame_Ins",
                                          "In_Frame_Del",
                                          "Frame_Shift_Ins",
                                          "Frame_Shift_Del",
                                          "Missense_Mutation")))
    }
    
    result <- do.call(rbind, by(mutDetail, mutDetail$pair_id, aggFunc))
    rownames(result) <- NULL

    result
}

Miao.cox <- function() {
    require(survival)
    master <- Miao.loadMaster()

    coxOS <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zHLA + Bladder + HNSCC + Lung + PD1 + Both, data = master)

    coxPFS <-
        coxph(Surv(pfs_days, pfs_event) ~
                  zTMB.By + zHLA + Bladder + HNSCC + Lung + PD1 + Both, data = master)

    coxOS  <- Miao.coxFrame(coxOS)
    coxPFS <- Miao.coxFrame(coxPFS)

    merged <- merge(coxOS, coxPFS, all = TRUE, by = "Covariate", suffixes = c(".OS", ".PFS"))

    rownames(merged) <- merged$Covariate
    merged$Covariate <- NULL

    ##merged <- merged[c("zHLA", "zTMB", "Bladder", "HNSCC", "Lung"),]
    merged <- merged[c("zHLA", "zTMB", "Bladder", "HNSCC", "Lung", "PD1", "Both"),]
    merged
}

Miao.coxFrame <- function(coxModel) {
    coeffMat <- summary(coxModel)$coefficients
    coxFrame <- 
        data.frame(Covariate   = rownames(coeffMat),
                   HazardRatio = coeffMat[,2],
                   HR_CI95_Lo  = exp(coeffMat[,1] - 2.0 * coeffMat[,3]),
                   HR_CI95_Up  = exp(coeffMat[,1] + 2.0 * coeffMat[,3]),
                   PValue      = coeffMat[,5])

    rownames(coxFrame) <- NULL
    coxFrame
}

Miao.coxPlot <- function(survType = "OS") {
    par(las = 1)
    par(fig = c(0.2, 1.0, 0.0, 1.0))

    plot(c(0.2, 8.0), c(1.5, 8.5),
         log  = "x",
         type = "n",
         axes = FALSE,
         xlab = "Hazard ratio",
         ylab = "")
    axis(1, at = c(0.2, 0.5, 1.0, 2.0, 4.0, 8.0))
    lines(c(1, 1), c(-2, 12), lty = 3)

    cox1 <- Miao.cox()
    errY <- 0.075

    plotCoeff <- function(rowName, height) {
        x1 <- cox1[rowName, sprintf("HR_CI95_Lo.%s", survType)]
        x2 <- cox1[rowName, sprintf("HR_CI95_Up.%s", survType)]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox1[rowName, sprintf("HazardRatio.%s", survType)]
        points(c(x, x), c(height, height), cex = 1.25, col = "black", pch = 16)

        text(x, height + 0.2, Miao.signifCode(cox1[rowName, sprintf("PValue.%s", survType)]), cex = 1.0)
    }

    plotCoeff("zHLA",    8.0)
    plotCoeff("zTMB.By", 7.0)
    plotCoeff("Bladder", 6.0)
    plotCoeff("HNSCC",   5.0)
    plotCoeff("Lung",    4.0)
    plotCoeff("PD1",     3.0)
    plotCoeff("Both",    2.0)
    box()

    text(2.0, 8.0, sprintf("p = %.2f", cox1["zHLA", sprintf("PValue.%s", survType)]), adj = 0, font = 3)

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(1.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    tx <- 0.85

    text(tx, 8.0, "HLA Score", adj = 1)
    text(tx, 7.1, "Mutation Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Bladder", adj = 1)
    text(tx, 5.0, "HNSCC", adj = 1)
    text(tx, 4.0, "Lung", adj = 1)
    text(tx, 3.0, "PD-1", adj = 1)
    text(tx, 2.0, "PD-1 + CTLA-4", adj = 1, cex = 0.8)
}

Miao.dataDir <- function() {
    homeDir <- Sys.getenv("CSB_DATA_VAULT", unset = NA)

    if (is.na(homeDir))
        JamLog.error("Environment variable CSB_DATA_VAULT is not set.");

    file.path(homeDir, "Miao")
}

Miao.loadAllelePresentation <- function() {
    read.csv(file.path(Miao.dataDir(), "Miao_Allele_Present.csv"))
}

Miao.loadMaster <- function() {
    read.csv(file.path(Miao.dataDir(), "Miao_Master.csv"))
}

Miao.loadMutDetail <- function() {
    JamIO.load(file.path(Miao.dataDir(), "Miao_SupTable5.RData"))
}

Miao.loadNeoDetail <- function() {
    JamIO.load(file.path(Miao.dataDir(), "Miao_SupTable10.RData"))
}

Miao.loadPatientDetail <- function() {
    read.csv(file.path(Miao.dataDir(), "Miao_SupTable2.csv"))
}

Miao.loadPatientPresentation <- function() {
    read.csv(file.path(Miao.dataDir(), "Miao_Patient_Present.csv"))
}

Miao.loadTSV <- function(fileName) {
    JamLog.info("Reading table [%s]...", fileName)
    read.table(fileName,  sep = "\t", header = TRUE, strip.white = TRUE, comment.char = "#")
}

Miao.signifCode <- function(p) {
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*",
    ifelse(p < 0.1,   ".", " "))))
}

Miao.writeGenotype <- function() {
    fileName <- file.path(Miao.dataDir(), "Miao_Genotype_Input.csv")
    genoFrame <- Miao.collectGenotype()

    write.csv(genoFrame, fileName, quote = FALSE, row.names = FALSE)
}

Miao.writeMaster <- function() {
    fileName <- file.path(Miao.dataDir(), "Miao_Master.csv")
    masterFrame <- Miao.compileMaster()

    write.csv(masterFrame, fileName, quote = FALSE, row.names = FALSE)
}

Miao.zscore <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

Miao.zBy <- function(dframe, keyCol, byCol, targetCol, scoreName) {

    aggFunc <- function(slice) {
        data.frame(key = slice[,keyCol], score = Miao.zscore(slice[,targetCol]), count = nrow(slice))
    }

    result <- do.call(rbind, by(dframe, dframe[,byCol], aggFunc))

    names(result) <- c(keyCol, scoreName, sprintf("%s.Count", byCol))
    rownames(result) <- NULL

    result
}
