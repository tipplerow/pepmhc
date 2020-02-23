
package pepmhc.app;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import jam.app.JamLogger;
import jam.hla.Allele;
import jam.hugo.HugoPeptideTable;
import jam.io.IOUtil;
import jam.peptide.Peptide;
import jam.util.ListUtil;

import pepmhc.affy.AffinityCache;
import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityRecord;
import pepmhc.stab.StabilityCache;
import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityRecord;

public final class AffinityStabilityCorr {
    private final Allele allele;
    private final String outputFile;
    private final String peptideFile;
    private final int    sampleSize;

    private PrintWriter writer;
    private List<Peptide> peptides;

    private static final AffinityMethod AFFINITY_METHOD = AffinityMethod.NET_MHC_PAN;
    private static final StabilityMethod STABILITY_METHOD = StabilityMethod.NET_MHC_STAB_PAN;

    private AffinityStabilityCorr(String[] args) {
        validate(args);

        this.allele      = Allele.instance(args[0]);
        this.peptideFile = args[1];
        this.outputFile  = args[2];
        this.sampleSize  = resolveSampleSize(args);
    }

    private static void validate(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: pepmhc.app.AffinityStabilityCorr ALLELE PEPTIDE_FILE OUTPUT_FILE [SAMPLE_SIZE]");
            System.exit(1);
        }
    }

    private static int resolveSampleSize(String[] args) {
        if (args.length < 4)
            return Integer.MAX_VALUE;
        else
            return Integer.parseInt(args[3]);
    }

    private void run() {
        openWriter();

        try {
            writeHeader();
            loadPeptides();
            processPeptides();
            JamLogger.info("DONE!");
        }
        finally {
            IOUtil.close(writer);
        }
    }

    private void openWriter() {
        writer = IOUtil.openWriter(outputFile);
    }

    private void writeHeader() {
        writer.println("peptide,affinity,halfLife,affinityPct,stabilityPct");
    }

    private void loadPeptides() {
        HugoPeptideTable table = HugoPeptideTable.load(peptideFile);
        peptides = new ArrayList<Peptide>(table.viewPeptides());

        if (sampleSize < peptides.size()) {
            ListUtil.shuffle(peptides);
            peptides = peptides.subList(0, sampleSize);
        }
    }

    private void processPeptides() {
        List<AffinityRecord> affinityRecords = AffinityCache.instance(AFFINITY_METHOD, allele).get(peptides);
        List<StabilityRecord> stabilityRecords = StabilityCache.instance(STABILITY_METHOD, allele).get(peptides);

        for (int index = 0; index < peptides.size(); ++index)
            writeLine(peptides.get(index), affinityRecords.get(index), stabilityRecords.get(index));
    }

    private void writeLine(Peptide peptide, AffinityRecord affinityRecord, StabilityRecord stabilityRecord) {
        writer.println(String.format("%s,%.3f,%.3f,%.2f,%.2f",
                                     peptide.formatString(),
                                     affinityRecord.getAffinity(),
                                     stabilityRecord.getHalfLife(),
                                     affinityRecord.getPercentile(),
                                     stabilityRecord.getPercentile()));
    }

    public static void main(String[] args) {
        AffinityStabilityCorr app = new AffinityStabilityCorr(args);
        app.run();
    }
}
