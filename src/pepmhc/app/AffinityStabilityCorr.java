
package pepmhc.app;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import jam.hla.Allele;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;
import pepmhc.stab.StabilityCache;
import pepmhc.stab.StabilityRecord;

public final class AffinityStabilityCorr {
    private final Allele allele;
    private final String outputFile;
    private final String peptideFile;

    private PrintWriter writer;
    private final List<Peptide> peptides = new ArrayList<Peptide>();

    private static final PredictionMethod METHOD = PredictionMethod.NET_MHC_PAN;

    private AffinityStabilityCorr(String[] args) {
        validate(args);

        this.allele      = Allele.instance(args[0]);
        this.peptideFile = args[1];
        this.outputFile  = args[2];
    }

    private static void validate(String[] args) {
        if (args.length != 3) {
            System.err.println("Usage: pepmhc.app.AffinityStabilityCorr ALLELE PEPTIDE_FILE OUTPUT_FILE");
            System.exit(1);
        }
    }

    private void run() {
        openWriter();

        try {
            writeHeader();
            loadPeptides();
            processPeptides();
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
        LineReader reader = LineReader.open(peptideFile);

        for (String line : reader)
            peptides.add(Peptide.parse(line));

        reader.close();
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
