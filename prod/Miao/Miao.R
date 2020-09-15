
JENE_MIAO <-
    file.path(JamEnv.getRequired("JENE_HOME"), "prod", "Miao", "Miao.R")

if (!file.exists(JENE_MIAO))
    JamLog.error("Cannot find [%s].", JENE_MIAO)

jamr.loadFile(JENE_MIAO)

## ---------------------------------------------------------------------

Miao.aggregateAffinity <- function(affinityFrame) {
    aggFunc <- function(slice) {
        data.frame(Tumor_Barcode  = slice$Tumor_Barcode[1],
                   Missense_Count = length(unique(slice$Protein_Change)),
                   Novel_Count    = sum(slice$Novel_Neo_Peptide),
                   Visible_Count  = sum(slice$Visible_Mutation))
    }

    aggList  <- by(affinityFrame, affinityFrame$Tumor_Barcode, aggFunc)
    aggFrame <- do.call(rbind, aggList)

    rownames(aggFrame) <- NULL
    aggFrame
}

## ---------------------------------------------------------------------

Miao.buildAffinityRN <- function(aggFrame) {
    pairList <- list()
    pairFrame <- Miao.pairRN()

    for (k in 1:nrow(pairFrame)) {
        resIndex <- which(aggFrame$Tumor_Barcode == pairFrame$Tumor_Barcode.Responder[k])
        nonIndex <- which(aggFrame$Tumor_Barcode == pairFrame$Tumor_Barcode.NonResponder[k])

        if (length(resIndex) != 1)
            next

        if (length(nonIndex) != 1)
            next

        pairRow <-
            data.frame(Tumor_Barcode.Responder    = aggFrame$Tumor_Barcode[resIndex],
                       Tumor_Barcode.NonResponder = aggFrame$Tumor_Barcode[nonIndex],
                       Novel_Count.Responder      = aggFrame$Novel_Count[resIndex],
                       Novel_Count.NonResponder   = aggFrame$Novel_Count[nonIndex],
                       Visible_Count.Responder    = aggFrame$Visible_Count[resIndex],
                       Visible_Count.NonResponder = aggFrame$Visible_Count[nonIndex])

        pairList[[length(pairList) + 1]] <- pairRow
    }

    pairFrame <- do.call(rbind, pairList)
    rownames(pairFrame) <- NULL

    pairFrame
}

Miao.buildChangeFromTo <- function(probFrame, novelThreshold = 0.5) {
    novelFrame <- subset(probFrame, Novel_Neo_Prob > novelThreshold)
    changeFrom <- substr(novelFrame$Protein_Change, 1, 1)
    changeTo   <- substr(novelFrame$Protein_Change, nchar(novelFrame$Protein_Change), nchar(novelFrame$Protein_Change))

    changeFrom <- by(changeFrom, changeFrom, length)
    changeTo   <- by(changeTo,   changeTo,   length)

    changeFrom <- changeFrom / sum(changeFrom)
    changeTo   <- changeTo /   sum(changeTo)

    residues <- sort(union(names(changeFrom), names(changeTo)))
    chgFrame <- data.frame(residue    = residues,
                           changeFrom = changeFrom[residues],
                           changeTo   = changeTo[residues])

    rownames(chgFrame) <- NULL
    chgFrame
}

Miao.buildImGenMaster <- function(affinityFrame, rankPctCenter = 2.0, rankPctWidth = 0.50) {
    aggFunc <- function(slice) {
        data.frame(Tumor_Barcode  = slice$Tumor_Barcode[1],
                   Missense_Count = length(unique(slice$Protein_Change)),
                   Novel_Count    = sum(slice$Novel_Neo_Peptide),
                   Visible_Count  = sum(slice$Visible_Mutation))
    }

    imGenList  <- by(affinityFrame, affinityFrame$Tumor_Barcode, aggFunc)
    imGenFrame <- do.call(rbind, imGenList)

    rownames(imGenFrame) <- NULL

    imGenFrame$zTMB <- Filter.zscore(log(1 + imGenFrame$Missense_Count))
    imGenFrame$zNov <- Filter.zscore(log(1 + imGenFrame$Novel_Count))
    imGenFrame$zVis <- Filter.zscore(log(1 + imGenFrame$Visible_Count))

    coxMod <- Miao.buildCoxModelFrame()


    master$zrTMB <- lm(zTMB ~ Cancer_Type - 1, data = master)$resid
    master$zrNov <- lm(zNov ~ Cancer_Type - 1, data = master)$resid
    master$zrVis <- lm(zVis ~ Cancer_Type - 1, data = master)$resid
    master

    ## coxModelFrame <- Miao.buildCoxModelFrame()
    ## affinityFrame <- Miao.loadMissAffinity(fragLen)

    ## affinityFrame$Neo_Present_Prob <-
    ##     Miao.mhcPresentationProb(affinityFrame$Neo_Affinity_Rank, rankPctCenter, rankPctWidth)

    ## affinityFrame$Self_Present_Prob <-
    ##     Miao.mhcPresentationProb(affinityFrame$Self_Affinity_Rank, rankPctCenter, rankPctWidth)

    ## affinityFrame <- merge(coxModelFrame, affinityFrame, by = "Tumor_Barcode")
    ## affinityFrame
}

## ---------------------------------------------------------------------

Miao.buildProteinChangeNovelProb <- function(affinityFrame) {
    probFunc <- function(slice) {
        slice$Novel_Neo_Peptide <-
            Filter.replaceNA(slice$Novel_Neo_Peptide, 0.0)

        data.frame(Tumor_Barcode  = slice$Tumor_Barcode[1],
                   Hugo_Symbol    = slice$Hugo_Symbol[1],
                   Protein_Change = slice$Protein_Change[1],
                   Novel_Neo_Prob = 1.0 - prod(1.0 - slice$Novel_Neo_Peptide))

    }

    index <- sprintf("%s|%s|%s",
                     affinityFrame$Tumor_Barcode,
                     affinityFrame$Hugo_Symbol,
                     affinityFrame$Protein_Change)

    probList <- by(affinityFrame, index, probFunc)
    probFrame <- do.call(rbind, probList)

    rownames(probFrame) <- NULL
    probFrame
}

## ---------------------------------------------------------------------

Miao.loadMissAffinity <- function(fragLen = 9, rankPctCenter = 2.0, rankPctWidth = 0.50) {
    baseName <- sprintf("Miao_MissAffinity_%d.txt", fragLen)
    fileName <- file.path(Miao.homeDir(), "Miss", baseName)

    affinityFrame <- read.table(JamIO.file(fileName), header = TRUE, sep = "\t")

    affinityFrame$Neo_Present_Prob <-
        Miao.mhcPresentationProb(affinityFrame$Neo_Peptide_Cleave_Prob,
                                 affinityFrame$Neo_Affinity_Rank,
                                 rankPctCenter, rankPctWidth)

    affinityFrame$Self_Present_Prob <-
        Miao.mhcPresentationProb(affinityFrame$Self_Peptide_Cleave_Prob,
                                 affinityFrame$Self_Affinity_Rank,
                                 rankPctCenter, rankPctWidth)

    affinityFrame$Novel_Neo_Peptide <-
        affinityFrame$Neo_Present_Prob * (1.0 - affinityFrame$Self_Present_Prob)

    affinityFrame$Visible_Mutation <-
        (affinityFrame$Neo_Present_Prob *
         affinityFrame$Self_Present_Prob *
         (affinityFrame$Neo_Peptide_Missense_Pos %in% 4:7))

    affinityFrame
}

Miao.loadMissCleavage <- function(fragLen = 9) {
    baseName <- sprintf("Miao_MissCleavage_%d.txt", fragLen)
    fileName <- file.path(Miao.homeDir(), "Miss", baseName)

    read.table(JamIO.file(fileName), header = TRUE, sep = "\t")
}

## ---------------------------------------------------------------------

Miao.maxNovel <- function(affinityFrame) {
    maxFunc <- function(slice) {
        slice[which.max(slice$Novel_Neo_Peptide),]
    }

    index <- sprintf("%s|%s|%s",
                     affinityFrame$Tumor_Barcode,
                     affinityFrame$Hugo_Symbol,
                     affinityFrame$Protein_Change)

    maxList <- by(affinityFrame, index, maxFunc)
    maxFrame <- do.call(rbind, maxList)

    rownames(maxFrame) <- NULL
    maxFrame
}

## ---------------------------------------------------------------------

Miao.mhcPresentationProb <- function(cleaveProb, affinityRank, rankPctCenter = 2.0, rankPctWidth = 0.50) {
    0.5 * cleaveProb * (1.0 - tanh((affinityRank - rankPctCenter) / rankPctWidth))
}

## ---------------------------------------------------------------------

Miao.plotNovelProteinChange <- function(chgFrame) {
    xwd <- 0.90
    yht <- 0.55

    par(las = 1)
    par(fig = c(0.5 * (1.0 - xwd), 0.5 * (1.0 + xwd), 1.0 - yht, 1.0))
    barplot(chgFrame$changeFrom, names.arg = chgFrame$residue,
            ylab = "Frequency", ylim = c(0.0, 0.2),
            col = Miao.residueColor(), density = NULL,
            cex.axis = 0.8, cex.lab = 0.8, cex.names = 0.8)
    legend("topleft", bty = "n", legend = c("FROM"), text.font = 2)

    par(fig = c(0.5 * (1.0 - xwd), 0.5 * (1.0 + xwd), 0.0, yht), new = TRUE)
    barplot(chgFrame$changeTo, names.arg = chgFrame$residue,
            ylab = "Frequency", ylim = c(0.0, 0.2),
            col = Miao.residueColor(), density = NULL,
            cex.axis = 0.8, cex.lab = 0.8, cex.names = 0.8)
    legend("topleft", bty = "n", legend = c("TO"), text.font = 2)

    par(xpd = TRUE)
    legend(17.5, 0.25, bty = "n",
           legend = c("Hydrophobic",
                      "Polar Neutral",
                      "Polar Charged",
                      "Unique"),
           pch = c(15, 15, 15, 15),
           col = c("gray50", "mediumpurple1", "lightgreen", "deepskyblue1"))
    par(xpd = FALSE)
}

## ---------------------------------------------------------------------

Miao.plotCohortTMBHist <- function(probFrame) {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    tmb <- by(probFrame, probFrame$Tumor_Barcode, nrow)

    truehist(log10(tmb),
             prob = FALSE,
             nbin = 16,
             xlab = expression(log[10](TMB)),
             ylab = "Count",
             xlim = c(0, 4),
             ylim = c(0, 40))
    box()
}

## ---------------------------------------------------------------------

Miao.plotMissensePosition <- function(novelFrame, novelThreshold = 0.5) {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    novelFrame <- subset(novelFrame, Novel_Neo_Peptide > novelThreshold)

    freq <- by(novelFrame, novelFrame$Neo_Peptide_Missense_Pos, nrow)
    freq <- freq / nrow(novelFrame)

    barplot(freq, xlab = "Mutation Position", ylab = "Frequency", ylim = c(0.0, 0.35))
}

## ---------------------------------------------------------------------

Miao.plotNovelNeoPeptideCDF <- function(probFrame) {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    probFrame <- probFrame[order(probFrame$Novel_Neo_Prob),]

    probFrame$pdf <- 1 / nrow(probFrame)
    probFrame$cdf <- cumsum(probFrame$pdf)
    probFrame <- subset(probFrame, Novel_Neo_Prob > 0)

    plot(probFrame$Novel_Neo_Prob, probFrame$cdf,
         lwd  = 2,
         type = "l",
         xlab = "Novel missense mutation probability",
         ylab = "Cumulative distribution",
         xlim = c(0, 1),
         ylim = c(0, 1))
}

Miao.plotNovelNeoPeptideHist <- function(probFrame) {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    truehist(probFrame$Novel_Neo_Prob,
             xlab = "Novel missense mutation probability",
             ylab = "Density",
             xlim = c(0, 1))
    box()
}

## ---------------------------------------------------------------------

Miao.plotTumorImmunogenicity <- function(simFrame, meanProb) {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    plot(simFrame$TMB, simFrame$MEDIAN,
         log  = "x",
         pch  = 16,
         xlab = "TMB",
         ylab = "Response probability",
         ylim = c(0.0, 1.0))

    for (k in 1:nrow(simFrame))
        lines(c(simFrame$TMB[k], simFrame$TMB[k]),
              c(simFrame$Q1[k],  simFrame$Q3[k]))

    lines(simFrame$TMB, 1.0 - exp(-simFrame$TMB * meanProb), lty = 3)

    legend("bottomright", bty = "n",
           legend = c("Median response probability",
                      "Interquartile range",
                      "Exponential approximation"), ## "1 - exp[-TMB * mean(novel mut prob)]"
           col = c(1, 1, 1),
           lty = c(NA, 1, 3),
           pch = c(16, NA, NA))
}

Miao.plotTumorImmunogenicityConst <- function() {
    par(las = 1)
    par(fig = c(0.05, 1, 0.125, 0.875))

    TMB <- 1:10000

    JamPlot.logX(xlim = c(1, 10000),
                 xlab = "TMB",
                 ylab = "Response probability",
                 ylim = c(0.0, 1.0))

    lines(TMB, 1.0 - (1.0 - 0.1) ^ TMB, lwd = 2, col = 1)
    lines(TMB, 1.0 - (1.0 - 0.01) ^ TMB, lwd = 2, col = 2)
    lines(TMB, 1.0 - (1.0 - 0.001) ^ TMB, lwd = 2, col = 4)

    legend("topleft", bty = "n",
           ##legend = c("0.001", "0.01", "0.1"),
           legend = c(expression(10 ^ -1),
                      expression(10 ^ -2),
                      expression(10 ^ -3)),
           col = c(1, 2, 4),
           lty = c(1, 1, 1),
           lwd = c(2, 2, 2))
}

## ---------------------------------------------------------------------

Miao.residueColor <- function() {
    c(A = "gray50",
      C = "mediumpurple1",
      D = "lightgreen",
      E = "lightgreen",
      F = "gray50",
      G = "deepskyblue1",
      H = "lightgreen",
      I = "gray50",
      K = "lightgreen",
      L = "gray50",
      M = "gray50",
      N = "mediumpurple1",
      P = "deepskyblue1",
      Q = "mediumpurple1",
      R = "lightgreen",
      S = "mediumpurple1",
      T = "mediumpurple1",
      V = "gray50",
      W = "gray50",
      Y = "gray50")
}

## ---------------------------------------------------------------------

Miao.simulate <- function(novelProbVec, TMB, sampleSize = 100) {
    rownames <- as.character(TMB)
    colnames <- c("TMB", "LOW", "Q1", "MEAN", "MEDIAN", "Q3", "HIGH")

    result <- matrix(nrow = length(rownames),
                     ncol = length(colnames),
                     dimnames = list(rownames, colnames))

    for (indexTMB in seq_along(TMB)) {
        JamLog.info("TMB = [%d]...", TMB[indexTMB])
        tumorSampleProb <- numeric(sampleSize)

        for (indexSample in seq_along(tumorSampleProb)) {
            novelProbSample <- sample(novelProbVec, TMB[indexTMB], replace = TRUE)
            tumorSampleProb[indexSample] <- 1.0 - prod(1.0 - novelProbSample)
        }

        result[indexTMB,] <- c(TMB[indexTMB], as.numeric(summary(tumorSampleProb)))
    }

    result <- as.data.frame(result)
    result
}


