
package pepmhc.agpro;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.TreeSet;

import jam.app.JamLogger;
import jam.ensembl.EnsemblDb;
import jam.ensembl.EnsemblRecord;
import jam.ensembl.TranscriptBiotype;
import jam.hugo.HugoSymbol;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.Peptide;

/**
 * Simulates antigen processing for a collection of genes and writes
 * the processed antigen peptide fragments to an output file indexed
 * by HUGO symbol.
 */
public final class AgProcDriver {
    private final String configFile;
    private final String symbolFile;
    private final String antigenFile;

    private LineReader reader;
    private PrintWriter writer;
    private AntigenProcessor processor;

    private static final String DELIM = "\t";

    private AgProcDriver(String configFile, String symbolFile, String antigenFile) {
        this.configFile = configFile;
        this.symbolFile = symbolFile;
        this.antigenFile = antigenFile;
    }

    /**
     * Simulates antigen processing for a collection of genes.
     *
     * @param configFile the name of a configuration file that defines
     * the antigen processor.
     *
     * @param symbolFile the name of a file containing HUGO symbols of
     * the genes to process.
     *
     * @param antigenFile the name of a file where the antigen peptide
     * fragments will be written.
     *
     * @throws RuntimeException if any errors occur.
     */
    public static void run(String configFile, String symbolFile, String antigenFile) {
        AgProcDriver driver = new AgProcDriver(configFile, symbolFile, antigenFile);
        driver.run();
    }

    private void run() {
        try {
            reader = LineReader.open(symbolFile);
            writer = IOUtil.openWriter(antigenFile);

            processor = AntigenProcessor.resolve(configFile);

            writeHeader();
            processGenes();
        }
        finally {
            IOUtil.close(reader);
            IOUtil.close(writer);
        }
    }

    private void writeHeader() {
        writer.println("Hugo_Symbol" + DELIM + "Peptide");
    }

    private void processGenes() {
        for (String symbol : reader)
            processGene(symbol);

        JamLogger.info("DONE!");
    }

    private void processGene(String symbol) {
        JamLogger.info("Processing gene [%s]...", symbol);
        writePeptides(symbol, generatePeptides(symbol));
    }

    private void writePeptides(String symbol, Collection<String> peptides) {
        for (String peptide : peptides)
            writer.println(symbol + DELIM + peptide);

        writer.flush();
    }

    private Collection<String> generatePeptides(String symbol) {
        Collection<String> peptides = new TreeSet<String>();
        Collection<EnsemblRecord> records = EnsemblDb.reference().get(HugoSymbol.instance(symbol));

        for (EnsemblRecord record : records)
            peptides.addAll(generatePeptides(record));

        return peptides;
    }

    private Collection<String> generatePeptides(EnsemblRecord record) {
        if (!record.getTranscriptBiotype().equals(TranscriptBiotype.PROTEIN_CODING))
            return Collections.emptyList();

        Peptide peptide = record.getPeptide();

        if (!peptide.isNative())
            return Collections.emptyList();

        Collection<Peptide> processed = processor.process(peptide);
        return Peptide.formatString(processed);
    }

    private static void usage() {
        System.err.println("Usage: java [JVMOPTIONS] jam.agproc.AgProcDriver <PROC_CONFIG> <GENE_INPUT> <ANTIGEN_OUTPUT>");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 3)
            usage();

        run(args[0], args[1], args[2]);
    }
}
