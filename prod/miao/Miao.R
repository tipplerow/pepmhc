
Miao.barplotA <- function(skip = 2) {
    dframe <- Miao.loadAllelePresentation()
    dframe <- dframe[grep("HLA-A", dframe$allele),]
    dframe <- dframe[order(-dframe$presentRate),]
    dframe <- dframe[seq(1, nrow(dframe), skip),]

    par(las = 2)
    par(fig = c(0.05, 1.0, 0.15, 0.85))

    barplot(dframe$presentRate,
            names.arg = substr(dframe$allele, 5, 11),
            cex.names = 0.9,
            ylab = "Presentation rate",
            ylim = c(0.0, 0.10),
            srt = 45)
}

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

    genoFrame$AlleleCount <-
        unlist(lapply(strsplit(genoFrame$Genotype, " "), function(x) length(unique(x))))

    genoFrame$Homozygous <-
        as.numeric(genoFrame$AlleleCount < 6)

    genoFrame
}

Miao.compileRawData <- function() {
    dataDir <- Miao.dataDir()

    table5  <- Miao.loadTSV(file.path(dataDir, "Miao_SupTable5.tsv"))
    table10 <- Miao.loadTSV(file.path(dataDir, "Miao_SupTable10.tsv"))

    JamIO.save(table5,  file.path(dataDir, "Miao_SupTable5.RData"))
    JamIO.save(table10, file.path(dataDir, "Miao_SupTable10.RData"))
}

Miao.compileMaster <- function(threshold = 100) {
    master <- Miao.loadPatientDetail()

    genoFrame <- Miao.collectGenotype()
    presentFrame <- Miao.loadPatientPresentation()

    master <- merge(master, genoFrame, by = "patient_id")
    master <- merge(master, presentFrame, by = "patient_id")
    
    mutDetail <- Miao.loadMutDetail()
    mutBurden <- Miao.computeTMB(mutDetail)

    master <- merge(master, mutBurden, by = "pair_id")

    neoDetail <- Miao.loadNeoDetail()
    neoBurden <- Miao.computeNAB(neoDetail, threshold)

    master <- merge(master, neoBurden, by = "Tumor_Sample_Barcode")

    ## Remove these cancer types with only one observation...
    master <- subset(master, !(cancer_type %in% c("Anal", "Sarcoma")))

    ## Compute z-scores for HLA presentation, tumor mutational burden,
    ## and neoantigen burden...
    master$zAGE <- Miao.zscore(master$age_start_io)
    master$zHLA <- Miao.zscore(master$actualRate)
    master$zNPR <- Miao.zscore(master$neoAgBindingRate)

    master$logTMB <- log(master$nonSilentCount)
    master$zTMB   <- Miao.zscore(master$logTMB)

    master$logNAB <- log(master$neoPeptideTotal)
    master$zNAB   <- Miao.zscore(master$logNAB)

    ## Compute the linear and log-ratio "excess neoantigen binding
    ## rates" (relative to the overall binding of the HLA genotype)...
    master$xNPR.Linear <- master$neoAgBindingRate - master$actualRate
    master$xNPR.LogRat <- log(master$neoAgBindingRate / master$actualRate)

    master$zxNPR.Linear <- Miao.zscore(master$xNPR.Linear)
    master$zxNPR.LogRat <- Miao.zscore(master$xNPR.LogRat)

    ## Compute another measure of excess neoantigen binding as the
    ## residual from a regression model...
    resid <- lm(neoAgBindingRate ~ actualRate, data = master)$resid
    stopifnot(length(resid) == nrow(master))

    master$xNPR.Resid <- resid
    master$zxNPR.Resid <- Miao.zscore(master$xNPR.Resid)

    ## Neoantigen presentation model...
    lmobj <- lm(log(neoAgBindingCount) ~ log(nonSilentCount) + log(actualRate), data = master)

    stopifnot(length(lmobj$fitted.values) == nrow(master))
    stopifnot(length(lmobj$residuals) == nrow(master))

    master$NAP.Actual <- log(master$neoAgBindingCount)
    master$NAP.Fitted <- lmobj$fitted.values
    master$NAP.Excess <- lmobj$residuals

    master$NAP.Actual.z <- Miao.zscore(master$NAP.Actual)
    master$NAP.Fitted.z <- Miao.zscore(master$NAP.Fitted)
    master$NAP.Excess.z <- Miao.zscore(master$NAP.Excess)

    ## "Neutralize" by cancer type...
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "logTMB", "zTMB.By"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "logNAB", "zNAB.By"))

    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "actualRate", "zHLA.By"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "neoAgBindingRate", "zNPR.By"))

    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "xNPR.Linear", "zxNPR.Linear.By"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "xNPR.LogRat", "zxNPR.LogRat.By"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "xNPR.Resid",  "zxNPR.Resid.By"))

    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "NAP.Actual",  "NAP.Actual.zBy"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "NAP.Fitted",  "NAP.Fitted.zBy"))
    master <- merge(master, Miao.zBy(master, "pair_id", "cancer_type", "NAP.Excess",  "NAP.Excess.zBy"))

    ## Death occurred if the overall survival is censored...
    master$os_event <- 1 - master$os_censor

    ## Progression occurred if the PFS duration is shorter than the overall duration...
    master$pfs_event <- as.numeric(master$pfs_days < master$os_days)

    ## Melanoma will be the reference...
    master$Bladder <- as.numeric(master$cancer_type == "Bladder")
    master$HNSCC   <- as.numeric(master$cancer_type == "HNSCC")
    master$Lung    <- as.numeric(master$cancer_type == "Lung")

    ## "Anti-CTLA-4" will be the reference...
    master$PD1  <- as.numeric(master$drug_type == "anti-PD-1/anti-PD-L1")
    master$Both <- as.numeric(master$drug_type == "anti-CTLA-4 + anti-PD-1/PD-L1")

    master
}

Miao.computeNAB <- function(neoDetail, threshold = 100) {
    ##
    ## Find the HLA allele with the strongest binding for each
    ## neopeptide...
    ##
    neoDetail <-
        neoDetail[order(neoDetail$Tumor_Sample_Barcode,
                        neoDetail$pep_mut,
                        neoDetail$affinity_mut),]

    neoDetail <-
        neoDetail[!duplicated(neoDetail[,c("Tumor_Sample_Barcode", "pep_mut")]),]

    aggFunc <- function(slice) {
        data.frame(Tumor_Sample_Barcode = slice$Tumor_Sample_Barcode[1],
                   neoPeptideTotal      = nrow(slice),
                   neoAgBindingCount    = sum(slice$affinity_mut <= threshold),
                   neoAgBindingRate     = mean(slice$affinity_mut <= threshold))
    }

    result <- do.call(rbind, by(neoDetail, neoDetail$Tumor_Sample_Barcode, aggFunc))
    rownames(result) <- NULL

    result
}

Miao.computeTMB <- function(mutDetail) {
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

Miao.cox <- function(master = NULL) {
    require(survival)

    if (is.null(master))
        master <- Miao.loadMaster()

    coxOS <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zHLA.By + Bladder + HNSCC + Lung + PD1, data = master)

    coxPFS <-
        coxph(Surv(pfs_days, pfs_event) ~
                  zTMB.By + zHLA.By + Bladder + HNSCC + Lung + PD1 + Both, data = master)

    coxOS  <- Miao.coxFrame(coxOS)
    coxPFS <- Miao.coxFrame(coxPFS)

    merged <- merge(coxOS, coxPFS, all = TRUE, by = "Covariate", suffixes = c(".OS", ".PFS"))

    rownames(merged) <- merged$Covariate
    merged$Covariate <- NULL

    ##merged <- merged[c("zHLA", "zTMB", "Bladder", "HNSCC", "Lung"),]
    merged <- merged[c("zHLA", "zTMB", "Bladder", "HNSCC", "Lung", "PD1", "Both"),]
    merged
}

Miao.coxNAP <- function(master = NULL, threshold = 100) {
    require(survival)

    if (is.null(master))
        master <- Miao.loadMaster(threshold)

    coxOS <-
        coxph(Surv(os_days, os_event) ~
                  NAP.Actual.zBy + zHLA.By + Bladder + HNSCC + Lung + PD1, data = master)

    coxPFS <-
        coxph(Surv(pfs_days, pfs_event) ~
                  NAP.Actual.zBy + zHLA.By + Bladder + HNSCC + Lung + PD1 + Both, data = master)

    coxOS  <- Miao.coxFrame(coxOS)
    coxPFS <- Miao.coxFrame(coxPFS)

    merged <- merge(coxOS, coxPFS, all = TRUE, by = "Covariate", suffixes = c(".OS", ".PFS"))

    rownames(merged) <- merged$Covariate
    merged$Covariate <- NULL

    ##merged <- merged[c("zHLA", "zTMB", "Bladder", "HNSCC", "Lung"),]
    merged <- merged[c("zHLA", "NAP", "Bladder", "HNSCC", "Lung", "PD1", "Both"),]
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

Miao.coxMatrix <- function(coxModel) {
    coeffMat <- summary(coxModel)$coefficients
    coxFrame <- 
        data.frame(HazardRatio = coeffMat[,2],
                   HR_CI95_Lo  = exp(coeffMat[,1] - 2.0 * coeffMat[,3]),
                   HR_CI95_Up  = exp(coeffMat[,1] + 2.0 * coeffMat[,3]),
                   PValue      = coeffMat[,5])

    rownames(coxFrame) <- rownames(coeffMat)
    as.matrix(coxFrame)
}

Miao.coxModel2 <- function(threshold = 100) {
    require(survival)
    master <- Miao.loadMaster(threshold)

    actual <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zNPR.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    excess <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zxNPR.LogRat.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    actual <- Miao.coxMatrix(actual)
    excess <- Miao.coxMatrix(excess)

    list(actual = actual, excess = excess)
}

Miao.coxModel4 <- function(master = NULL, threshold = 100) {
    require(survival)

    if (is.null(master))
        master <- Miao.loadMaster(threshold)

    cox1 <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    cox2 <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zHLA.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    cox3 <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zNPR.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    ## cox4 <-
    ##     coxph(Surv(os_days, os_event) ~
    ##               zTMB.By + zHLA.By + zxNPR.LogRat.By +
    ##               Bladder + HNSCC + Lung + PD1 + Both, data = master)

    cox4 <-
        coxph(Surv(os_days, os_event) ~
                  zTMB.By + zxNPR.LogRat.By +
                  Bladder + HNSCC + Lung + PD1 + Both, data = master)

    cox1 <- Miao.coxMatrix(cox1)
    cox2 <- Miao.coxMatrix(cox2)
    cox3 <- Miao.coxMatrix(cox3)
    cox4 <- Miao.coxMatrix(cox4)

    list(cox1 = cox1, cox2 = cox2, cox3 = cox3, cox4 = cox4)
}

Miao.coxPlot <- function(master = NULL, survType = "PFS") {
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

    cox1 <- Miao.cox(master)
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

    text(4.5, 8.0, sprintf("p = %.2f", cox1["zHLA", sprintf("PValue.%s", survType)]), adj = 0, font = 3, cex = 0.85)
    text(4.5, 7.0, sprintf("p = %.3f", cox1["zTMB", sprintf("PValue.%s", survType)]), adj = 0, font = 3, cex = 0.85)

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(1.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    tx <- 0.85

    par(xpd = TRUE)
    text(tx, 8.1, "HLA Score", adj = 1)
    text(tx, 7.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 7.1, "Mutation Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Bladder", adj = 1)
    text(tx, 5.0, "HNSCC", adj = 1)
    text(tx, 4.0, "Lung", adj = 1)
    text(tx, 3.0, "PD-1", adj = 1)
    text(tx, 2.0, "PD-1 + CTLA-4", adj = 1, cex = 1.0)
    par(xpd = FALSE)
}

Miao.coxPlotNAP <- function(master = NULL, survType = "PFS") {
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

    cox1 <- Miao.coxNAP(master)
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
    plotCoeff("NAP",     7.0)
    plotCoeff("Bladder", 6.0)
    plotCoeff("HNSCC",   5.0)
    plotCoeff("Lung",    4.0)
    plotCoeff("PD1",     3.0)
    plotCoeff("Both",    2.0)
    box()

    text(4.5, 8.0, sprintf("p = %.3f", cox1["zHLA", sprintf("PValue.%s", survType)]), adj = 0, font = 3, cex = 0.85)
    text(4.5, 7.0, sprintf("p = %.3f", cox1["NAP",  sprintf("PValue.%s", survType)]), adj = 0, font = 3, cex = 0.85)

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(1.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    tx <- 0.85

    par(xpd = TRUE)
    text(tx, 8.1, "HLA Score", adj = 1)
    text(tx, 7.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 7.1, "Neoantigen Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Bladder", adj = 1)
    text(tx, 5.0, "HNSCC", adj = 1)
    text(tx, 4.0, "Lung", adj = 1)
    text(tx, 3.0, "PD-1", adj = 1)
    text(tx, 2.0, "PD-1 + CTLA-4", adj = 1)
    par(xpd = FALSE)
}

Miao.coxPlot2 <- function(threshold = 100) {
    par(las = 1)
    par(fig = c(0.2, 1.0, 0.0, 1.0))

    ##plot(c(0.2, 8.0), c(0.5, 7.5),
    plot(c(0.46, 8.0), c(0.5, 7.5),
         log  = "x",
         type = "n",
         axes = FALSE,
         xlab = "Hazard ratio",
         ylab = "")
    ##axis(1, at = c(0.2, 0.5, 1.0, 2.0, 4.0, 8.0))
    axis(1, at = c(0.5, 1.0, 2.0, 4.0, 8.0))
    lines(c(1, 1), c(-2, 12), lty = 3)

    cox2 <- Miao.coxModel2(threshold)
    errY <- 0.075

    model.cex <- c(actual = 1.25, excess = 1.25)
    model.col <- c(actual = "red", excess = "black")
    model.pch <- c(actual = 15, excess = 16)
    model.dpv <- c(actual = 0.2, excess = -0.2)

    plotCoeff <- function(modelKey, covariate, height, showP = FALSE) {
        x  <- cox2[[modelKey]][covariate, "HazardRatio"]
        x1 <- cox2[[modelKey]][covariate, "HR_CI95_Lo"]
        x2 <- cox2[[modelKey]][covariate, "HR_CI95_Up"]
        pv <- cox2[[modelKey]][covariate, "PValue"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        points(c(x, x), c(height, height),
               cex = model.cex[modelKey],
               col = model.col[modelKey],
               pch = model.pch[modelKey])

        text(x, height + model.dpv[modelKey], Miao.signifCode(pv))

        if (showP)
            text(1.4, height, sprintf("p = %s", Miao.pValueString(pv)), adj = 0, cex = 0.75, font = 3)
    }

    plotCoeffs <- function(covariate, height, showP = FALSE) {
        dy <- 0.1

        plotCoeff("actual", covariate, height + dy, showP)
        plotCoeff("excess", covariate, height - dy, showP)
    }

    plotCoeff("actual", "zNPR.By", 7.1, TRUE)
    plotCoeff("excess", "zxNPR.LogRat.By", 6.9, TRUE)

    ##pActual <- cox2[["actual"]]["zNPR.By", "PValue"]
    ##pExcess <- cox2[["excess"]]["zxNPR.LogRat.By", "PValue"]

    ##text(1.4, 7.1, sprintf("p = %.2f", pActual), adj = 0, cex = 0.75, font = 3)
    ##text(1.4, 6.9, sprintf("p = %.3f", pExcess), adj = 0, cex = 0.75, font = 3)

    plotCoeffs("zTMB.By", 6.0, TRUE)
    plotCoeffs("Bladder", 5.0)
    plotCoeffs("HNSCC",   4.0)
    plotCoeffs("Lung",    3.0)
    plotCoeffs("PD1",     2.0)
    plotCoeffs("Both",    1.0)
    box()

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(0.5, 7.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    px <- 0.99
    tx <- 0.85
    dy <- 0.13

    points(px, 7.0 + dy, cex = model.cex["actual"], col = model.col["actual"], pch = model.pch["actual"])
    points(px, 7.0 - dy, cex = model.cex["excess"], col = model.col["excess"], pch = model.pch["excess"])

    text(tx, 7.0 + dy, "Actual NPR", adj = 1, font = 2)
    text(tx, 7.0 - dy, "Excess NPR", adj = 1, font = 2)

    text(tx, 6.1, "Mutation Load", adj = 1)
    text(tx, 5.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 5.0, "Bladder", adj = 1)
    text(tx, 4.0, "HNSCC", adj = 1)
    text(tx, 3.0, "Lung", adj = 1)
    text(tx, 2.0, "PD-1", adj = 1)
    text(tx, 1.0, "PD-1 + CTLA-4", adj = 1, cex = 0.95)
}

Miao.coxPlot4 <- function(threshold = 100) {
    par(las = 1)
    par(fig = c(0.2, 1.0, 0.0, 1.0))

    plot(c(0.2, 8.0), c(0.5, 9.5),
         log  = "x",
         type = "n",
         axes = FALSE,
         xlab = "Hazard ratio",
         ylab = "")
    axis(1, at = c(0.2, 0.5, 1.0, 2.0, 4.0, 8.0))
    lines(c(1, 1), c(-2, 12), lty = 3)

    cox4 <- Miao.coxModel4(threshold)
    errY <- 0.075

    model.cex <- c(1.2, 1.2, 1.2, 1.5)
    model.col <- c("black", "red", "purple", "blue")
    model.pch <- c(15, 16, 17, 18)

    plotCoeff <- function(modelIndex, Covariate, height, p.digits) {
        x  <- cox4[[modelIndex]][Covariate, "HazardRatio"]
        x1 <- cox4[[modelIndex]][Covariate, "HR_CI95_Lo"]
        x2 <- cox4[[modelIndex]][Covariate, "HR_CI95_Up"]
        pv <- cox4[[modelIndex]][Covariate, "PValue"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        points(c(x, x), c(height, height),
               cex = model.cex[modelIndex],
               col = model.col[modelIndex],
               pch = model.pch[modelIndex])

        text(x1 * 0.9, height, Miao.signifCode(pv), adj = 1, cex = 1.0)

        if (p.digits > 0) {
            fmt <- sprintf("p = %%.%df", p.digits)
            text(x2 * 1.1, height, sprintf(fmt, pv), adj = 0, cex = 0.8, font = 3)
        }
    }

    plotCoeffs <- function(Covariate, height) {
        dy <- 0.15

        plotCoeff(1, Covariate, height + 1.5 * dy, 0)
        plotCoeff(2, Covariate, height + 0.5 * dy, 0)
        plotCoeff(3, Covariate, height - 0.5 * dy, 0)
        plotCoeff(4, Covariate, height - 1.5 * dy, 0)
    }


    ## pHLA <- cox2["zHLA.By",         "PValue.HLA"]
    ## pNPR <- cox2["zxNPR.LogRat.By", "PValue.NPR"]

    ## text(2.2,  8.0 + dy, sprintf("p = %4.2f", pHLA), adj = 0, cex = 0.8, font = 3)
    ## text(0.85, 8.0 - dy, sprintf("p = %4.2f", pNPR), adj = 1, cex = 0.8, font = 3)

    plotCoeff(2, "zHLA.By", 9.0, 2)
    plotCoeff(3, "zNPR.By", 8.0, 2)
    plotCoeff(4, "zxNPR.LogRat.By", 7.0, 3)

    plotCoeffs("zTMB.By", 6.0)
    plotCoeffs("Bladder", 5.0)
    plotCoeffs("HNSCC", 4.0)
    plotCoeffs("Lung", 3.0)
    plotCoeffs("PD1", 2.0)
    plotCoeffs("Both", 1.0)
    box()

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(0.5, 9.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    px <- 0.99
    tx <- 0.85
    dy <- 0.13

    ## points(px, 8.0 + dy, cex = cex.HLA, col = col.HLA, pch = pch.HLA)
    ## points(px, 8.0 - dy, cex = cex.NPR, col = col.NPR, pch = pch.NPR)

    ## text(tx, 8.0 + dy, "HLA Score", adj = 1, font = 2)
    ## text(tx, 8.0 - dy, "HLA Score",  adj = 1, font = 2)

    par(xpd = TRUE)
    text(tx, 9.1, "Self-Ag Binding", adj = 1)
    text(tx, 8.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)

    text(tx, 8.1, "Neo-Ag Binding", adj = 1)
    text(tx, 7.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)

    text(tx, 7.1, "Excess Neo-Ag", adj = 1, font = 2)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)

    text(tx, 6.1, "Mutation Load", adj = 1)
    text(tx, 5.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)

    text(tx, 5.0, "Bladder", adj = 1)
    text(tx, 4.0, "HNSCC", adj = 1)
    text(tx, 3.0, "Lung", adj = 1)
    text(tx, 2.0, "PD-1", adj = 1)
    text(tx, 1.0, "PD-1 + CTLA-4", adj = 1)
    par(xpd = FALSE)
}

Miao.coxHomo <- function(master = NULL, threshold = 100) {
    require(survival)

    if (is.null(master))
        master <- Miao.loadMaster(threshold)

    model1 <- coxph(Surv(pfs_days, pfs_event) ~ Homozygous + zTMB.By + Bladder + HNSCC + Lung + PD1 + Both, data = master)
    model2 <- coxph(Surv(pfs_days, pfs_event) ~ zHLA.By    + zTMB.By + Bladder + HNSCC + Lung + PD1 + Both, data = master)

    frame1 <- Miao.coxFrame(model1)
    frame2 <- Miao.coxFrame(model2)

    master <- merge(frame1, frame2, all = TRUE, by = "Covariate", suffixes = c(".Binary", ".ZScore"))

    rownames(master) <- master$Covariate
    master$Covariate <- NULL

    master <- master[c("Homozygous", "zHLA", "zTMB", "Bladder", "HNSCC", "Lung", "PD1", "Both"),]
    master
}

Miao.coxHomoPlot <- function(master = NULL, threshold = 100) {
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

    cox <- Miao.coxHomo(master, threshold)
    errY <- 0.075

    cex.binary <- 1.25
    col.binary <- "red"
    pch.binary <- 15

    cex.zscore <- 1.25
    col.zscore <- "black"
    pch.zscore <- 16

    plotBinary <- function(rowName, height) {
        x1 <- cox[rowName, "HR_CI95_Lo.Binary"]
        x2 <- cox[rowName, "HR_CI95_Up.Binary"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox[rowName, "HazardRatio.Binary"]
        points(c(x, x), c(height, height), cex = cex.binary, col = col.binary, pch = pch.binary)

        text(x, height + 0.2, Miao.signifCode(cox[rowName, "PValue.Binary"]), cex = 1.0)
    }

    plotZScore <- function(rowName, height) {
        x1 <- cox[rowName, "HR_CI95_Lo.ZScore"]
        x2 <- cox[rowName, "HR_CI95_Up.ZScore"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox[rowName, "HazardRatio.ZScore"]
        points(c(x, x), c(height, height), cex = cex.zscore, col = col.zscore, pch = pch.zscore)

        text(x, height - 0.2, Miao.signifCode(cox[rowName, "PValue.ZScore"]), cex = 1.0)
    }

    dy <- 0.1

    plotBinary("Homozygous", 8.0 + dy)
    plotZScore("zHLA",       8.0 - dy)

    pHomozygous <- cox["Homozygous", "PValue.Binary"]
    pHLAScore   <- cox["zHLA", "PValue.ZScore"]

    xp <- 4.8

    text(xp, 8.13, sprintf("p = %.2f", pHomozygous), adj = 0, cex = 0.8, font = 3)
    text(xp, 7.87, sprintf("p = %.2f", pHLAScore),   adj = 0, cex = 0.8, font = 3)

    plotBinary("zTMB", 7.0 + dy)
    plotZScore("zTMB", 7.0 - dy)

    text(xp, 7.13, sprintf("p = %.3f", cox["zTMB", "PValue.Binary"]), adj = 0, cex = 0.8, font = 3)
    text(xp, 6.87, sprintf("p = %.3f", cox["zTMB", "PValue.ZScore"]), adj = 0, cex = 0.8, font = 3)

    plotBinary("Bladder", 6.0 + dy)
    plotZScore("Bladder", 6.0 - dy)

    plotBinary("HNSCC", 5.0 + dy)
    plotZScore("HNSCC", 5.0 - dy)

    plotBinary("Lung", 4.0 + dy)
    plotZScore("Lung", 4.0 - dy)

    plotBinary("PD1", 3.0 + dy)
    plotZScore("PD1", 3.0 - dy)

    plotBinary("Both", 2.0 + dy)
    plotZScore("Both", 2.0 - dy)

    box()

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(1.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    px <- 0.99
    tx <- 0.85
    dy <- 0.13

    points(px, 8.0 + dy, cex = cex.binary, col = col.binary, pch = pch.binary)
    points(px, 8.0 - dy, cex = cex.zscore, col = col.zscore, pch = pch.zscore)

    par(xpd = TRUE)
    text(tx, 8.0 + dy, "Homozygous", adj = 1, font = 2)
    text(tx, 8.0 - dy, "HLA Score",  adj = 1, font = 2)
    text(tx, 7.1, "Mutation Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Bladder", adj = 1)
    text(tx, 5.0, "HNSCC", adj = 1)
    text(tx, 4.0, "Lung", adj = 1)
    text(tx, 3.0, "PD-1", adj = 1)
    text(tx, 2.0, "PD-1 + CTLA-4", adj = 1)
    par(xpd = FALSE)
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

Miao.loadMaster <- function(threshold = 100) {
    read.csv(Miao.masterFileName(threshold))
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
    read.table(fileName,  sep = "\t", quote = "", header = TRUE, strip.white = TRUE, comment.char = "#")
}

Miao.masterFileName <- function(threshold = 100) {
    file.path(Miao.dataDir(), sprintf("Miao_Master_%03d.csv", threshold))
}

Miao.plotExcessNeo <- function(master = NULL) {
    if (is.null(master))
        master <- Miao.loadMaster()

    par(las = 1)
    par(fig = c(0.05, 1.0, 0.15, 0.85))

    JamPlot.logXY(xlim = c(1, 10000),
                  ylim = c(1, 10000),
                  xlab = "TMB * HLA",
                  ylab = "Presented neoantigen count")

    x <- master$nonSilentCount * master$actualRate
    y <- master$neoAgBindingCount

    points(x, y)
    JamPlot.loglogline(x, y, 10 ^ seq(-1, 4, 0.1))
}

Miao.pValueString <- function(p) {
    fmt <- sprintf("%%.%df", 1 + ceiling(-log10(p)))
    sprintf(fmt, p)
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

Miao.writeMaster <- function(threshold = 100) {
    fileName <- Miao.masterFileName(threshold)
    masterFrame <- Miao.compileMaster(threshold)

    write.csv(masterFrame, fileName, quote = FALSE, row.names = FALSE)
}

Miao.zscore <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

Miao.zBy <- function(dframe, keyCol, byCol, targetCol, scoreName) {

    aggFunc <- function(slice) {
        data.frame(key = slice[,keyCol], score = Miao.zscore(slice[,targetCol]))
    }

    result <- do.call(rbind, by(dframe, dframe[,byCol], aggFunc))

    names(result) <- c(keyCol, scoreName)
    rownames(result) <- NULL

    result
}
