
package pepmhc.proc;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import jam.app.JamLogger;
import jam.ensembl.EnsemblGene;
import jam.ensembl.EnsemblMap;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

/**
 * Reads an input file containing Ensembl gene identifiers and writes
 * an output file containing all unique 9-mer peptides predicted to be
 * cleaved by a proteasome and transported into the ER by TAP and thus
 * presented for possible MHC binding.
 */
public final class GeneProcessor {
    private final String inputFile;
    private final String outputFile;
    private final boolean mutate;
    private final Set<String> fragments;

    private GeneProcessor(String inputFile, String outputFile, boolean mutate) {
        this.mutate = mutate;
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.fragments = new TreeSet<String>();
    }

    private static final int LOG_INTERVAL = 1;
    private static final int PEPTIDE_LENGTH = 9;

    private void run() {
        generateFragments();
        writeFragments();
    }

    private void generateFragments() {
        int processed = 0;
        LineReader reader = LineReader.open(inputFile);

        for (String gene : reader) {
            ++processed;
            processGene(EnsemblGene.instance(gene));

            if (processed % LOG_INTERVAL == 0)
                JamLogger.info("Processed [%d] genes...", processed);
        }

        reader.close();
    }

    private void processGene(EnsemblGene gene) {
        Collection<Peptide> proteins = EnsemblMap.instance().get(gene);

        for (Peptide protein : proteins)
            if (protein.isNative())
                processProtein(protein);
    }

    private void processProtein(Peptide protein) {
        if (mutate)
            protein = protein.mutate();

        List<Peptide> cleaved = NetChop.chop(protein, PEPTIDE_LENGTH);
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
        System.err.println("Usage: java pepmhc.proc.GeneProcessor INPUT_FILE OUTPUT_FILE [MUTATE]");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length < 2 || args.length > 3)
            usage();

        String inputFile = args[0];
        String outputFile = args[1];
        boolean mutate = false;

        if (args.length == 3)
            mutate = Boolean.parseBoolean(args[2]);

        GeneProcessor processor =
            new GeneProcessor(inputFile, outputFile, mutate);

        processor.run();
    }
}
