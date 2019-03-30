
Balachandran.dataDir <- function() {
    homeDir <- Sys.getenv("CSB_DATA_VAULT", unset = NA)

    if (is.na(homeDir))
        JamLog.error("Environment variable CSB_DATA_VAULT is not set.");

    file.path(homeDir, "Balachandran")
}

Balachandran.loadSup <- function() {
    read.table(file.path(Balachandran.dataDir(), "Sup1Post.txt"), sep = "\t", header = TRUE)
}

Balachandran.writeAlleleFlat <- function(fileName) {
    master <- Balachandran.loadSup()
    master <- master[!duplicated(master$Sample),]

    alleles <-
        sort(unique(unlist(strsplit(master$HLA, ","))))

    writeLines(alleles, fileName)
}

Balachandran.writeGenotypeInput <- function(fileName) {
    master <- Balachandran.loadSup()
    master <- master[!duplicated(master$Sample),]
    master <- master[,c("Sample", "HLA")]

    master$HLA <-
        gsub(",", " ", master$HLA)

    write.csv(master, fileName, quote = FALSE, row.names = FALSE)
}

Balachandran.signifCode <- function(p) {
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*",
    ifelse(p < 0.1,   ".", " "))))
}

Balachandran.zscore <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
