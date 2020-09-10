
package pepmhc.app;

import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import jam.app.JamLogger;
import jam.io.IOUtil;
import jam.math.DoubleUtil;

import jene.fasta.FastaPeptideReader;
import jene.fasta.FastaPeptideRecord;
import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityPredictor;

public final class AllelePresentCalc {
    private final int pepLen;
    private final String fastIn;
    private final Allele allele;
    private final AffinityMethod method;

    private AffinityPredictor predictor;
    private FastaPeptideReader fastaReader;
    private PrintWriter reportWriter;
    private Set<Peptide> pepFragments;
    private int binderCount;

    private AllelePresentCalc(String fastIn, AffinityMethod method, Allele allele, int pepLen) {
        this.fastIn = fastIn;
        this.method = method;
        this.allele = allele;
        this.pepLen = pepLen;
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
        predictor    = method.getPredictor();
        fastaReader  = FastaPeptideReader.open(fastIn);
        pepFragments = new HashSet<Peptide>();
        binderCount  = 0;

        for (FastaPeptideRecord record : fastaReader)
            processPeptide(record.getPeptide());

        reportWriter = IOUtil.openWriter(outputFile());
        reportWriter.println("method,allele,pepLen,fragmentCount,binderFrac");
        reportWriter.println(String.format("%s,%s,%d,%d,%8.6f",
                                           method,
                                           allele,
                                           pepLen,
                                           pepFragments.size(),
                                           DoubleUtil.ratio(binderCount, pepFragments.size())));

        fastaReader.close();
        reportWriter.close();
    }

    private void processPeptide(Peptide peptide) {
        for (Peptide fragment : peptide.nativeFragments(pepLen))
            processFragment(fragment);
    }

    private void processFragment(Peptide fragment) {
        //
        // Do not double-count any duplicate fragments...
        //
        if (pepFragments.contains(fragment))
            return;

        double ic50 = predictor.predict(allele, fragment).getAffinity().doubleValue();

        if (ic50 < BINDING_THRESHOLD)
            ++binderCount;

        pepFragments.add(fragment);

        if (pepFragments.size() % LOG_INTERVAL == 0)
            JamLogger.info("Processed [%d] peptide fragments...", pepFragments.size());
    }

    private String outputFile() {
        return String.format("%s_%s_%s_%d%s",
                             REPORT_PREFIX,
                             method,
                             allele.longKey().replace("*", ""),
                             pepLen,
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

        int pepLen;
        String fastIn;
        Allele allele;
        AffinityMethod method;

        fastIn = args[0];
        method = AffinityMethod.valueOf(args[1].toUpperCase());
        allele = Allele.instance(args[2]);
        pepLen = Integer.parseInt(args[3]);

        AllelePresentCalc calculator =
            new AllelePresentCalc(fastIn, method, allele, pepLen);

        calculator.run();
    }
}
