
package pepmhc.app;

import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import jam.app.JamLogger;
import jam.io.IOUtil;
import jam.fasta.FastaReader;
import jam.fasta.FastaRecord;
import jam.math.DoubleUtil;
import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class AllelePresentCalc {
    private final int peptideLength;

    private final String fastaFile;
    private final String alleleCode;

    private final PredictionMethod predMethod;

    private Predictor predictor;
    private FastaReader fastaReader;
    private PrintWriter reportWriter;
    private Set<Peptide> pepFragments;
    private int binderCount;

    private AllelePresentCalc(String fastaFile, PredictionMethod predMethod, String alleleCode, int peptideLength) {
        this.fastaFile = fastaFile;
        this.predMethod = predMethod;
        this.alleleCode = alleleCode;
        this.peptideLength = peptideLength;
    }

    private static final int LOG_INTERVAL = 1000;

    /**
     * Prefix for the name of the output file, written to the working
     * directory.
     */
    public static final String REPORT_PREFIX = "presentation-rate";

    /**
     * Suffix for the name of the output file, written to the working
     * directory.
     */
    public static final String REPORT_SUFFIX = ".csv";

    /**
     * Threshold IC50 value to be considered "bound".
     */
    public static final double BINDING_THRESHOLD = 500.0;

    private void run() {
        predictor    = Predictor.instance(predMethod);
        fastaReader  = FastaReader.open(fastaFile);
        pepFragments = new HashSet<Peptide>();
        binderCount  = 0;

        for (FastaRecord record : fastaReader)
            processPeptide(record.getPeptide());

        reportWriter = IOUtil.openWriter(outputFile());
        reportWriter.println("predMethod,alleleCode,peptideLength,fragmentCount,binderFrac");
        reportWriter.println(String.format("%s,%s,%d,%d,%8.6f",
                                           predMethod,
                                           alleleCode,
                                           peptideLength,
                                           pepFragments.size(),
                                           DoubleUtil.ratio(binderCount, pepFragments.size())));

        fastaReader.close();
        reportWriter.close();
    }

    private void processPeptide(Peptide peptide) {
        for (Peptide fragment : peptide.nativeFragments(peptideLength))
            processFragment(fragment);
    }

    private void processFragment(Peptide fragment) {
        //
        // Do not double-count any duplicate fragments...
        //
        if (pepFragments.contains(fragment))
            return;

        double ic50 = predictor.predict(alleleCode, fragment).getAffinity();

        if (ic50 < BINDING_THRESHOLD)
            ++binderCount;

        pepFragments.add(fragment);

        if (pepFragments.size() % LOG_INTERVAL == 0)
            JamLogger.info("Processed [%d] peptide fragments...", pepFragments.size());
    }

    private String outputFile() {
        return String.format("%s_%s_%s_%d%s",
                             REPORT_PREFIX,
                             predMethod,
                             alleleCode.replace("*", ""),
                             peptideLength,
                             REPORT_SUFFIX);
    }

    private static void usage() {
        System.err.println("Usage: java pepmhc.app.AllelePresentCalc "
                           + "FASTA_FILE PREDICTION_METHOD ALLELE_CODE PEPTIDE_LENGTH");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 4)
            usage();

        int peptideLength;
        String alleleCode, fastaFile;
        PredictionMethod predMethod;

        fastaFile     = args[0];
        predMethod    = PredictionMethod.valueOf(args[1].toUpperCase());
        alleleCode    = args[2];
        peptideLength = Integer.parseInt(args[3]);

        AllelePresentCalc calculator =
            new AllelePresentCalc(fastaFile, predMethod, alleleCode, peptideLength);

        calculator.run();
    }
}
