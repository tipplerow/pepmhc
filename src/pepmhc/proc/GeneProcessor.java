
package pepmhc.proc;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import jam.app.JamLogger;
import jam.ensembl.EnsemblDb;
import jam.ensembl.EnsemblRecord;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.HugoSymbol;
import jam.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

/**
 * Reads an input file containing HUGO symbols and writes an output
 * file containing all unique 9-mer and 10-mer peptides predicted to
 * be cleaved by a proteasome and transported into the ER by TAP and
 * thus presented for possible MHC binding.
 */
public final class GeneProcessor {
    private final String inputFile;
    private final String outputFile;
    private final Set<String> fragments;

    private GeneProcessor(String inputFile, String outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.fragments = new TreeSet<String>();
    }

    private static final int LOG_INTERVAL = 1;
    private static final int[] PEPTIDE_LENGTHS = new int[] { 9, 10 };

    private void run() {
        generateFragments();
        writeFragments();
    }

    private void generateFragments() {
        int processed = 0;
        LineReader reader = LineReader.open(inputFile);

        for (String hugo : reader) {
            ++processed;
            processGene(HugoSymbol.instance(hugo));

            if (processed % LOG_INTERVAL == 0)
                JamLogger.info("Processed [%d] genes...", processed);
        }

        reader.close();
    }

    private void processGene(HugoSymbol hugo) {
        Collection<EnsemblRecord> records = EnsemblDb.global().get(hugo);

        for (EnsemblRecord record : records)
            processRecord(record);
    }

    private void processRecord(EnsemblRecord record) {
        Peptide protein = record.getPeptide();

        if (protein.isNative())
            processProtein(protein);
    }

    private void processProtein(Peptide protein) {
        List<Peptide> cleaved = NetChop.chop(protein, PEPTIDE_LENGTHS);
        List<Peptide> transported = TAP.INSTANCE.transport(cleaved);

        for (Peptide fragment : transported)
            fragments.add(fragment.formatString());
    }

    private void writeFragments() {
        JamLogger.info("Writing fragments...");
        PrintWriter writer = IOUtil.openWriter(outputFile);

        for (String fragment : fragments)
            writer.println(fragment);

        writer.close();
    }

    private static void usage() {
        System.err.println("Usage: java pepmhc.proc.GeneProcessor INPUT_FILE OUTPUT_FILE");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 2)
            usage();

        String inputFile = args[0];
        String outputFile = args[1];

        GeneProcessor processor =
            new GeneProcessor(inputFile, outputFile);

        processor.run();
    }
}
