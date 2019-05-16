
Chowell.countHomozygous <- function() {
    genotype <- read.csv(file.path(Chowell.dataDir(), "Chowell_Genotype_Input.csv"))
    alleleList <- strsplit(genotype$HLA_Class_I_Alleles, " ")

    homozygous <- character(0)
    heterozygous <- character(0)

    findHomozygous <- function(x) {
        if (x[1] == x[2])
            homozygous <<- c(homozygous, x[1])
        else
            heterozygous <<- c(heterozygous, x[1], x[2])

        if (x[3] == x[4])
            homozygous <<- c(homozygous, x[3])
        else
            heterozygous <<- c(heterozygous, x[3], x[4])

        if (x[5] == x[6])
            homozygous <<- c(homozygous, x[5])
        else
            heterozygous <<- c(heterozygous, x[5], x[6])
    }

    lapply(alleleList, findHomozygous)

    homoCount   <- tapply(homozygous,   homozygous,   length)
    heteroCount <- tapply(heterozygous, heterozygous, length)

    homoFrame   <- data.frame(allele = names(homoCount),   homoCount   = homoCount)
    heteroFrame <- data.frame(allele = names(heteroCount), heteroCount = heteroCount)

    result <- merge(homoFrame, heteroFrame, all = TRUE)

    result$homoCount <- Filter.replaceNA(result$homoCount, 0)
    result$heteroCount <- Filter.replaceNA(result$heteroCount, 0)

    result$homoRatio <- result$homoCount / result$heteroCount
    result
}


Chowell.cox1 <- function(zByType = TRUE) {
    require(survival)
    cohort <- Chowell.loadCohort1()

    if (zByType) {
        model1 <- coxph(Surv(OS_Months, OS_Event) ~ Homozygous + zTMB.By + zAGE + M1a + M1b + M1c + NSCLC + PD1, data = cohort)
        model2 <- coxph(Surv(OS_Months, OS_Event) ~ zHLA       + zTMB.By + zAGE + M1a + M1b + M1c + NSCLC + PD1, data = cohort)
    }
    else {
        model1 <- coxph(Surv(OS_Months, OS_Event) ~ Homozygous + zTMB + zAGE + M1a + M1b + M1c + NSCLC + PD1, data = cohort)
        model2 <- coxph(Surv(OS_Months, OS_Event) ~ zHLA       + zTMB + zAGE + M1a + M1b + M1c + NSCLC + PD1, data = cohort)
    }

    frame1 <- Chowell.coxFrame(model1)
    frame2 <- Chowell.coxFrame(model2)

    master <- merge(frame1, frame2, all = TRUE, by = "Covariate", suffixes = c(".Binary", ".ZScore"))

    rownames(master) <- master$Covariate
    master$Covariate <- NULL

    master <- master[c("Homozygous", "zHLA", "zTMB", "zAGE", "M1a", "M1b", "M1c", "NSCLC", "PD1"),]
    master
}

Chowell.cox1A <- function() {
    require(survival)

    cohort <- Chowell.loadCohort1()
    model1 <- coxph(Surv(OS_Months, OS_Event) ~ Homozygous + zTMB + zAGE + M1a + M1b + M1c + PD1, data = cohort)
    model2 <- coxph(Surv(OS_Months, OS_Event) ~ zHLA       + zTMB + zAGE + M1a + M1b + M1c + PD1, data = cohort)

    frame1 <- Chowell.coxFrame(model1)
    frame2 <- Chowell.coxFrame(model2)

    master <- merge(frame1, frame2, all = TRUE, by = "Covariate", suffixes = c(".Binary", ".ZScore"))

    rownames(master) <- master$Covariate
    master$Covariate <- NULL

    master <- master[c("Homozygous", "zHLA", "zTMB", "zAGE", "M1a", "M1b", "M1c", "PD1"),]
    master
}

Chowell.coxPlot1 <- function(zByType = TRUE) {
    par(las = 1)
    par(fig = c(0.2, 1.0, 0.0, 1.0))

    plot(c(0.5, 4.0), c(0.5, 8.5),
         log  = "x",
         type = "n",
         axes = FALSE,
         xlab = "Hazard ratio",
         ylab = "")
    axis(1, at = c(0.5, 1.0, 2.0, 4.0))
    lines(c(1, 1), c(-2, 12), lty = 3)

    cox1 <- Chowell.cox1(zByType)
    errY <- 0.075

    cex.binary <- 1.25
    col.binary <- "red"
    pch.binary <- 15

    cex.zscore <- 1.25
    col.zscore <- "black"
    pch.zscore <- 16

    plotBinary <- function(rowName, height) {
        x1 <- cox1[rowName, "HR_CI95_Lo.Binary"]
        x2 <- cox1[rowName, "HR_CI95_Up.Binary"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox1[rowName, "HazardRatio.Binary"]
        points(c(x, x), c(height, height), cex = cex.binary, col = col.binary, pch = pch.binary)

        text(x, height + 0.2, Chowell.signifCode(cox1[rowName, "PValue.Binary"]), cex = 1.0)
    }

    plotZScore <- function(rowName, height) {
        x1 <- cox1[rowName, "HR_CI95_Lo.ZScore"]
        x2 <- cox1[rowName, "HR_CI95_Up.ZScore"]

        lines(c(x1, x2), c(height, height), col = 1)
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox1[rowName, "HazardRatio.ZScore"]
        points(c(x, x), c(height, height), cex = cex.zscore, col = col.zscore, pch = pch.zscore)

        text(x, height - 0.2, Chowell.signifCode(cox1[rowName, "PValue.ZScore"]), cex = 1.0)
    }

    dy <- 0.1

    plotBinary("Homozygous", 8.0 + dy)
    plotZScore("zHLA",       8.0 - dy)

    pHomozygous <- cox1["Homozygous", "PValue.Binary"]
    pHLAScore   <- cox1["zHLA", "PValue.ZScore"]

    text(2.9, 8.13, sprintf("p = %5.3f", pHomozygous), adj = 0, cex = 0.8, font = 3)
    text(2.9, 7.87, sprintf("p = %4.2f", pHLAScore),   adj = 0, cex = 0.8, font = 3)

    plotBinary("zTMB", 7.0 + dy)
    plotZScore("zTMB", 7.0 - dy)

    text(2.9, 7.13, sprintf("p = %7.5f", cox1["zTMB", "PValue.Binary"]), adj = 0, cex = 0.8, font = 3)
    text(2.9, 6.87, sprintf("p = %7.5f", cox1["zTMB", "PValue.ZScore"]), adj = 0, cex = 0.8, font = 3)

    plotBinary("zAGE", 6.0 + dy)
    plotZScore("zAGE", 6.0 - dy)

    plotBinary("M1a", 5.0 + dy)
    plotZScore("M1a", 5.0 - dy)

    plotBinary("M1b", 4.0 + dy)
    plotZScore("M1b", 4.0 - dy)

    plotBinary("M1c", 3.0 + dy)
    plotZScore("M1c", 3.0 - dy)

    plotBinary("NSCLC", 2.0 + dy)
    plotZScore("NSCLC", 2.0 - dy)

    plotBinary("PD1", 1.0 + dy)
    plotZScore("PD1", 1.0 - dy)

    box()

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(0.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    px <- 0.99
    tx <- 0.85
    dy <- 0.13

    points(px, 8.0 + dy, cex = cex.binary, col = col.binary, pch = pch.binary)
    points(px, 8.0 - dy, cex = cex.zscore, col = col.zscore, pch = pch.zscore)

    text(tx, 8.0 + dy, "Homozygous", adj = 1, font = 2)
    text(tx, 8.0 - dy, "HLA Score",  adj = 1, font = 2)

    text(tx, 7.1, "Mutation Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Age", adj = 1)
    text(tx, 5.0, "Stage M1a", adj = 1)
    text(tx, 4.0, "Stage M1b", adj = 1)
    text(tx, 3.0, "Stage M1c", adj = 1)
    text(tx, 2.0, "NSCLC", adj = 1)
    text(tx, 1.0, "PD-1", adj = 1)
}

Chowell.cox2 <- function(zByType = TRUE, minCount = 5) {
    require(survival)
    cohort <- Chowell.loadCohort2()

    if (zByType) {
        model1 <- coxph(Surv(OS_Months, OS_Event) ~ Homozygous
                        + zTMB.By + Age_31_50 + Age_50_60 + Age_61_70 + Age_GT_71 + CTLA4 + PD1,
                        data = cohort, subset = Cancer_Type.Count >= minCount)

        model2 <- coxph(Surv(OS_Months, OS_Event) ~ zHLA
                        + zTMB.By + Age_31_50 + Age_50_60 + Age_61_70 + Age_GT_71 + CTLA4 + PD1,
                        data = cohort, subset = Cancer_Type.Count >= minCount)
    }
    else {
        model1 <- coxph(Surv(OS_Months, OS_Event) ~ Homozygous
                        + zTMB + Age_31_50 + Age_50_60 + Age_61_70 + Age_GT_71 + CTLA4 + PD1, data = cohort)

        model2 <- coxph(Surv(OS_Months, OS_Event) ~ zHLA
                        + zTMB + Age_31_50 + Age_50_60 + Age_61_70 + Age_GT_71 + CTLA4 + PD1, data = cohort)
    }

    frame1 <- Chowell.coxFrame(model1)
    frame2 <- Chowell.coxFrame(model2)

    master <- merge(frame1, frame2, all = TRUE, by = "Covariate", suffixes = c(".Binary", ".ZScore"))

    rownames(master) <- master$Covariate
    master$Covariate <- NULL

    master <- master[c("Homozygous", "zHLA", "zTMB", "Age_31_50", "Age_50_60", "Age_61_70", "Age_GT_71", "CTLA4", "PD1"),]
    master
}

Chowell.coxPlot2 <- function(zByType = TRUE, minCount = 5) {
    par(las = 1)
    par(fig = c(0.2, 1.0, 0.0, 1.0))

    plot(c(0.4, 2.3), c(0.5, 8.5),
         log  = "x",
         type = "n",
         axes = FALSE,
         xlab = "Hazard ratio",
         ylab = "")
    axis(1, at = c(0.5, 1.0, 2.0))
    lines(c(1, 1), c(-2, 12), lty = 3)

    cox2 <- Chowell.cox2(zByType, minCount)
    errY <- 0.075

    cex.binary <- 1.25
    col.binary <- "red"
    pch.binary <- 15

    cex.zscore <- 1.25
    col.zscore <- "black"
    pch.zscore <- 16

    plotBinary <- function(rowName, height) {
        x1 <- cox2[rowName, "HR_CI95_Lo.Binary"]
        x2 <- cox2[rowName, "HR_CI95_Up.Binary"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox2[rowName, "HazardRatio.Binary"]
        points(c(x, x), c(height, height), cex = cex.binary, col = col.binary, pch = pch.binary)

        text(x, height + 0.2, Chowell.signifCode(cox2[rowName, "PValue.Binary"]), cex = 1.0)
    }

    plotZScore <- function(rowName, height) {
        x1 <- cox2[rowName, "HR_CI95_Lo.ZScore"]
        x2 <- cox2[rowName, "HR_CI95_Up.ZScore"]

        lines(c(x1, x2), c(height, height))
        lines(c(x1, x1), c(height - errY, height + errY))
        lines(c(x2, x2), c(height - errY, height + errY))

        x <- cox2[rowName, "HazardRatio.ZScore"]
        points(c(x, x), c(height, height), cex = cex.zscore, col = col.zscore, pch = pch.zscore)

        text(x, height - 0.2, Chowell.signifCode(cox2[rowName, "PValue.ZScore"]), cex = 1.0)
    }

    dy <- 0.1

    plotBinary("Homozygous", 8.0 + dy)
    plotZScore("zHLA",       8.0 - dy)

    pHomozygous <- cox2["Homozygous", "PValue.Binary"]
    pHLAScore   <- cox2["zHLA", "PValue.ZScore"]

    text(1.8, 8.13, sprintf("p = %5.3f", pHomozygous), adj = 0, cex = 0.8, font = 3)
    text(1.8, 7.87, sprintf("p = %4.2f", pHLAScore),   adj = 0, cex = 0.8, font = 3)

    plotBinary("zTMB", 7.0 + dy)
    plotZScore("zTMB", 7.0 - dy)

    text(1.8, 7.13, sprintf("p = %6.4f", cox2["zTMB", "PValue.Binary"]), adj = 0, cex = 0.8, font = 3)
    text(1.8, 6.87, sprintf("p = %6.4f", cox2["zTMB", "PValue.ZScore"]), adj = 0, cex = 0.8, font = 3)

    plotBinary("Age_31_50", 6.0 + dy)
    plotZScore("Age_31_50", 6.0 - dy)

    plotBinary("Age_50_60", 5.0 + dy)
    plotZScore("Age_50_60", 5.0 - dy)

    plotBinary("Age_61_70", 4.0 + dy)
    plotZScore("Age_61_70", 4.0 - dy)

    plotBinary("Age_GT_71", 3.0 + dy)
    plotZScore("Age_GT_71", 3.0 - dy)

    plotBinary("CTLA4", 2.0 + dy)
    plotZScore("CTLA4", 2.0 - dy)

    plotBinary("PD1", 1.0 + dy)
    plotZScore("PD1", 1.0 - dy)

    box()

    par(fig = c(0.0, 0.36, 0.0, 1.0), new = TRUE)
    plot(c(0.0, 1.0), c(0.5, 8.5),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "")

    px <- 0.99
    tx <- 0.85
    dy <- 0.13

    points(px, 8.0 + dy, cex = cex.binary, col = col.binary, pch = pch.binary)
    points(px, 8.0 - dy, cex = cex.zscore, col = col.zscore, pch = pch.zscore)

    text(tx, 8.0 + dy, "Homozygous", adj = 1, font = 2)
    text(tx, 8.0 - dy, "HLA Score",  adj = 1, font = 2)
    text(tx, 7.1, "Mutation Load", adj = 1)
    text(tx, 6.8, "(z-score by cancer type)", adj = 1, font = 3, cex = 0.59)
    text(tx, 6.0, "Age 31-50", adj = 1)
    text(tx, 5.0, "Age 50-60", adj = 1)
    text(tx, 4.0, "Age 61-70", adj = 1)
    text(tx, 3.0, "Age > 71", adj = 1)
    text(tx, 2.0, "CTLA-4", adj = 1)
    text(tx, 1.0, "PD-1/PDL-1", adj = 1)
}

Chowell.coxFrame <- function(coxModel) {
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

Chowell.dataDir <- function() {
    homeDir <- Sys.getenv("CSB_DATA_VAULT", unset = NA)

    if (is.na(homeDir))
        JamLog.error("Environment variable CSB_DATA_VAULT is not set.");

    file.path(homeDir, "Chowell")
}

Chowell.histGenotype <- function() {
    dframe <- Chowell.loadGenotypePresentation()
    dframe <- subset(dframe, actualRate < 0.3)

    par(las = 1)
    par(fig = c(0.05, 1.0, 0.15, 0.85))

    truehist(dframe$actualRate,
             h = 0.01,
             xlim = c(0.0, 0.3),
             xlab = "Presentation rate",
             ylab = "Density")
    box()
}

Chowell.homoPresent <- function() {
    allelePresent <- Chowell.loadAllelePresentation()

    allelePresent$allele <- gsub("HLA-", "", allelePresent$allele)
    allelePresent$allele <- gsub("\\*",  "", allelePresent$allele)
    allelePresent$allele <- gsub(":",    "", allelePresent$allele)

    homozygousCount <- Chowell.countHomozygous();

    workspace <- merge(allelePresent, homozygousCount)
    workspace
}

Chowell.loadAllelePresentation <- function() {
    read.csv(file.path(Chowell.dataDir(), "Chowell_Allele_Present.csv"))
}

Chowell.loadCohort <- function(fileName) {
    cohort <- read.table(fileName, sep = "\t", header = TRUE, strip.white = TRUE)
    cohort
}

Chowell.loadCohort1 <- function() {
    present <- Chowell.loadGenotypePresentation()

    cohort <- Chowell.loadCohort(file.path(Chowell.dataDir(), "Chowell_Cohort1.tsv"))
    cohort <- merge(cohort, present, by = "Sample")

    cohort <- subset(cohort, is.finite(OS_Months))
    cohort <- subset(cohort, is.finite(MutCnt))
    cohort <- subset(cohort, !is.na(Stage_M))

    cohort$M1a <- as.numeric(cohort$Stage_M == "M1a")
    cohort$M1b <- as.numeric(cohort$Stage_M == "M1b")
    cohort$M1c <- as.numeric(cohort$Stage_M == "M1c")

    cohort$NSCLC <- as.numeric(cohort$Cancer_Type == "Non-Small Cell Lung Cancer")
    cohort$PD1   <- as.numeric(cohort$Drug_Class == "PD-1")

    cohort$zTMB <- Chowell.zscore(log(cohort$MutCnt))
    cohort$zAGE <- Chowell.zscore(cohort$Age)
    cohort$zHLA <- Chowell.zscore(cohort$actualRate)

    cohort$logTMB <- log(cohort$MutCnt)

    cohort <- merge(cohort, Chowell.zBy(cohort, "Sample", "Cancer_Type", "logTMB", "zTMB.By"))
    cohort
}

Chowell.loadCohort2 <- function() {
    present <- Chowell.loadGenotypePresentation()

    cohort <- Chowell.loadCohort(file.path(Chowell.dataDir(), "Chowell_Cohort2.tsv"))
    cohort <- merge(cohort, present, by = "Sample")

    cohort$Age_31_50 <- as.numeric(cohort$Age_Group == "31-50")
    cohort$Age_50_60 <- as.numeric(cohort$Age_Group == "50-60")
    cohort$Age_61_70 <- as.numeric(cohort$Age_Group == "61-70")
    cohort$Age_GT_71 <- as.numeric(cohort$Age_Group == ">71")

    cohort$CTLA4 <- as.numeric(cohort$Drug_Class == "CTLA4")
    cohort$PD1   <- as.numeric(cohort$Drug_Class == "PD-1/PDL-1")

    ## Add 1 because one tumor has zero mutations...
    cohort$logTMB <- log(1.0 + cohort$IMPACT_MutCnt)

    cohort$zTMB <- Chowell.zscore(cohort$logTMB)
    cohort$zHLA <- Chowell.zscore(cohort$actualRate)

    cohort <- merge(cohort, Chowell.zBy(cohort, "Sample", "Cancer_Type", "logTMB", "zTMB.By"), all.x = TRUE)
    cohort
}

Chowell.loadGenotypePresentation <- function() {
    read.csv(file.path(Chowell.dataDir(), "Chowell_Patient_Present.csv"))
}

Chowell.plotOverlap <- function(size = 500) {
    dframe <- Chowell.loadGenotypePresentation()

    XWD <- 0.52
    YHT <- 0.60

    par(las = 1)
    par(fig = c(0.0, XWD, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)))

    k <- sample.int(nrow(dframe), size)
    x <- 100 * dframe$idealRate[k]
    y <- 100 * dframe$actualRate[k]

    plot(x, y,
         cex  = 0.50,
         pch  = 16,
         axes = FALSE,
         type = "p",
         xlab = "",
         ylab = "",
         xlim = c(0, 40),
         ylim = c(0, 40))
    lines(c(-10, 50), c(-10, 50), lty = 2)

    axis(1, cex.axis = 0.9)
    axis(2, cex.axis = 0.9)
    box()

    text(2, 38, "A", font = 2)
    mtext("Non-overlapping presentation [%]", side = 1, line = 2.5, cex = 0.9)

    par(xpd = TRUE)
    text(-12, 20, "Actual genotype presentation [%]", cex = 0.9, srt = 90)
    par(xpd = FALSE)

    overlap <- 100 * (dframe$idealRate - dframe$actualRate)
    overlap <- overlap[which(overlap < 20)]

    par(fig = c(1.0 - XWD, 1.0, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)), new = TRUE)
    truehist(overlap, h = 0.5, xlab = "", ylab = "", axes = FALSE)
    axis(1, cex.axis = 0.9)
    axis(2, cex.axis = 0.9)
    box()

    text(18.5, 0.20, "B", font = 2)
    mtext("Presentation overlap loss [%]", side = 1, line = 2.5, cex = 0.9)

    par(xpd = TRUE)
    text(-6.3, 0.11, "Density", cex = 0.9, srt = 90)
    par(xpd = FALSE)
}

Chowell.plotOverlapHomo <- function() {
    chow1 <- Chowell.loadCohort1()
    chow2 <- Chowell.loadCohort2()

    par(las = 1)
    par(fig = c(0.05, 1.0, 0.15, 0.85))

    plot(chow1$idealRate,
         chow1$actualRate,
         xlim = c(0.0, 0.50),
         ylim = c(0.0, 0.30),
         xlab = "Assumed non-overlapping presentation rate",
         ylab = "True genotype presentation rate",
         pch  = 16,
         col  = ifelse(chow1$Homozygous, 2, 1))
}

Chowell.violinHLARate <- function() {
    require(vioplot)
    par(las = 1)

    XWD  <- 0.575
    YHT  <- 0.60
    y.at <- c(0.0, 0.1, 0.2, 0.3)
    ylim <- c(0.0, 0.3)

    cohort1 <- Chowell.loadCohort1()
    cohort2 <- Chowell.loadCohort2()

    plotViolin <- function(cohort) {
        plot(c(-1.0, 1.0),
             ylim,
             axes = FALSE,
             type = "n",
             xlab = "",
             ylab = "")
        vioplot(cohort$actualRate[cohort$Homozygous == 1], at = -0.5, add = TRUE, col = "lightgreen")
        vioplot(cohort$actualRate[cohort$Homozygous == 0], at =  0.5, add = TRUE, col = "lightblue")
        return(0)
        
        vioplot(cohort$zHLA[cohort$Homozygous == 1],
                cohort$zHLA[cohort$Homozygous == 0],
                col = color,
                ylim = c(-3.0, 3.0),
                names = c("Homo", "Hetero"), add = TRUE)
    }

    par(fig = c(0.0, XWD, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)))
    plotViolin(cohort1)

    axis(1, at = c(-0.5, 0.5), labels = c("Homozygous", "Heterozygous"))
    axis(2, at = y.at, labels = TRUE)

    mtext("HLA Presentation Rate", side = 2, line = 2.4, las = 0)
    text(0.90, 0.28, "1", font = 2, cex = 1)
    
    par(fig = c(1.0 - XWD, 1.0, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)), new = TRUE)
    plotViolin(cohort2)

    axis(1, at = c(-0.5, 0.5), labels = c("Homozygous", "Heterozygous"))
    axis(2, at = y.at, labels = FALSE)

    text(0.90, 0.28, "2", font = 2, cex = 1)
}

Chowell.violinHLAZ <- function() {
    require(vioplot)
    par(las = 1)

    XWD <- 0.575
    YHT <- 0.60

    cohort1 <- Chowell.loadCohort1()
    cohort2 <- Chowell.loadCohort2()

    plotViolin <- function(cohort) {
        plot(c(-1.0, 1.0),
             c(-3.0, 3.0),
             axes = FALSE,
             type = "n",
             xlab = "",
             ylab = "")
        vioplot(cohort$zHLA[cohort$Homozygous == 1], at = -0.5, add = TRUE, col = "lightgreen")
        vioplot(cohort$zHLA[cohort$Homozygous == 0], at =  0.5, add = TRUE, col = "lightblue")
        return(0)
        
        vioplot(cohort$zHLA[cohort$Homozygous == 1],
                cohort$zHLA[cohort$Homozygous == 0],
                col = color,
                ylim = c(-3.0, 3.0),
                names = c("Homo", "Hetero"), add = TRUE)
    }

    par(fig = c(0.0, XWD, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)))
    plotViolin(cohort1)

    axis(1, at = c(-0.5, 0.5), labels = c("Homozygous", "Heterozygous"))
    axis(2, at = -3:3, labels = TRUE)

    mtext("HLA Presentation (z-score)", side = 2, line = 2.4, las = 0)
    text(0.90, 2.6, "1", font = 2, cex = 1)
    
    par(fig = c(1.0 - XWD, 1.0, 0.5 * (1.0 - YHT), 0.5 * (1.0 + YHT)), new = TRUE)
    plotViolin(cohort2)

    axis(1, at = c(-0.5, 0.5), labels = c("Homozygous", "Heterozygous"))
    axis(2, at = -3:3, labels = FALSE)

    text(0.90, 2.6, "2", font = 2, cex = 1)
}

Chowell.writeGenotypeFlat <- function(fileName) {
    cohort1 <- Chowell.loadCohort1()
    cohort2 <- Chowell.loadCohort2()

    alleles <- c(cohort1$HLA_Class_I_Alleles, cohort2$HLA_Class_I_Alleles)
    alleles <- strsplit(alleles, ",")
    alleles <- lapply(alleles, function(x) paste(sort(unique(x)), collapse = ","))
    alleles <- unlist(alleles)
    alleles <- sort(unique(alleles))

    writeLines(alleles, fileName)
}

Chowell.writeGenotypeInput <- function(fileName) {
    cohort1 <- Chowell.loadCohort1()
    cohort2 <- Chowell.loadCohort2()

    cohort1 <- cohort1[,c("Sample", "HLA_Class_I_Alleles")]
    cohort2 <- cohort2[,c("Sample", "HLA_Class_I_Alleles")]

    master <- rbind(cohort1, cohort2)
    master$HLA_Class_I_Alleles <- gsub(",", " ", master$HLA_Class_I_Alleles)

    write.csv(master, fileName, quote = FALSE, row.names = FALSE)
}

Chowell.signifCode <- function(p) {
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*",
    ifelse(p < 0.1,   ".", " "))))
}

Chowell.zscore <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

Chowell.zBy <- function(dframe, keyCol, byCol, targetCol, scoreName) {

    aggFunc <- function(slice) {
        data.frame(key = slice[,keyCol], score = Chowell.zscore(slice[,targetCol]), count = nrow(slice))
    }

    result <- do.call(rbind, by(dframe, dframe[,byCol], aggFunc))

    names(result) <- c(keyCol, scoreName, sprintf("%s.Count", byCol))
    rownames(result) <- NULL

    result
}
