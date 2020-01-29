
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

import pepmhc.binder.BindingRecord;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;
import pepmhc.stab.StabilityCache;
import pepmhc.stab.StabilityRecord;

public final class AffinityStabilityCorr {
    private final Allele allele;
    private final String outputFile;
    private final String peptideFile;
    private final int    sampleSize;

    private PrintWriter writer;
    private List<Peptide> peptides;

    private static final PredictionMethod METHOD = PredictionMethod.NET_MHC_PAN;

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
        writer.println("peptide,affinity,halfLife");
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
        List<BindingRecord>   bindingRecords   = AffinityCache.get(METHOD, allele, peptides);
        List<StabilityRecord> stabilityRecords = StabilityCache.get(allele, peptides);

        for (int index = 0; index < peptides.size(); ++index)
            writeLine(peptides.get(index),
                      bindingRecords.get(index).getAffinity(),
                      stabilityRecords.get(index).getHalfLife());
    }

    private void writeLine(Peptide peptide, double affinity, double halfLife) {
        writer.println(String.format("%s,%.3f,%.3f", peptide.formatString(), affinity, halfLife));
    }

    public static void main(String[] args) {
        AffinityStabilityCorr app = new AffinityStabilityCorr(args);
        app.run();
    }
}
