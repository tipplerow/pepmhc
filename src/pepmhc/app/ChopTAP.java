
package pepmhc.app;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.TreeSet;

import jam.io.IOUtil;

import jene.ensembl.EnsemblProteinDb;
import jene.ensembl.EnsemblProteinRecord;
import jene.ensembl.TranscriptBiotype;
import jene.hugo.HugoSymbol;
import jene.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

/**
 * Predicts proteasomal cleavage and TAP transport for all genes in an
 * Ensembl database and writes the cleaved and transported peptides to
 * standard output.
 */
public final class ChopTAP {
    private final PrintWriter writer;
    private final EnsemblProteinDb db;

    private static final int[] PEPTIDE_LENGTHS = new int[] { 9, 10 };

    private ChopTAP(EnsemblProteinDb db, PrintWriter writer) {
        this.db = db;
        this.writer = writer;
    }

    public static void run(EnsemblProteinDb db, PrintWriter writer) {
        ChopTAP chopTAP = new ChopTAP(db, writer);
        chopTAP.run();
    }

    public static void run(String fastaFile, String outputFile) {
        run(EnsemblProteinDb.load(fastaFile), IOUtil.openWriter(outputFile));
    }

    private void run() {
        writeHeader();
        processSymbols();
    }

    private void writeHeader() {
        writer.println("Hugo_Symbol\tPeptide");
    }

    private void processSymbols() {
        Collection<HugoSymbol> symbols = new TreeSet<HugoSymbol>(db.hugoSet());

        for (HugoSymbol symbol : symbols)
            processSymbol(symbol);
    }

    private void processSymbol(HugoSymbol symbol) {
        writePeptides(symbol, generatePeptides(symbol));
    }

    private void writePeptides(HugoSymbol symbol, Collection<String> peptides) {
        for (String peptide : peptides)
            writer.println(symbol.getKey() + "\t" + peptide);

        writer.flush();
    }

    private Collection<String> generatePeptides(HugoSymbol symbol) {
        Collection<String> peptides = new TreeSet<String>();
        Collection<EnsemblProteinRecord> records = db.get(symbol);

        for (EnsemblProteinRecord record : records)
            peptides.addAll(generatePeptides(record));

        return peptides;
    }

    private Collection<String> generatePeptides(EnsemblProteinRecord record) {
        if (!record.getTranscriptBiotype().equals(TranscriptBiotype.PROTEIN_CODING))
            return Collections.emptyList();

        Peptide peptide = record.getPeptide();

        if (!peptide.isNative())
            return Collections.emptyList();

        Collection<Peptide> chopped = NetChop.chop(peptide, PEPTIDE_LENGTHS);
        Collection<Peptide> transported = TAP.consensus().transport(chopped);

        return Peptide.formatString(transported);
    }

    private static void usage() {
        System.err.println("Usage: java [JVMOPTIONS] pepmhc.app.ChopTAP <ENSEMBL_FILE> <OUTPUT_FILE>");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 2)
            usage();

        run(args[0], args[1]);
    }
}
