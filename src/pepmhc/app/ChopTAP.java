
package pepmhc.app;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.TreeSet;

import jam.ensembl.EnsemblDb;
import jam.ensembl.EnsemblRecord;
import jam.peptide.HugoSymbol;
import jam.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

/**
 * Predicts proteasomal cleavage and TAP transport for all genes in an
 * Ensembl database and writes the cleaved and transported peptides to
 * standard output.
 */
public final class ChopTAP {
    private final EnsemblDb db;

    private static final int[] PEPTIDE_LENGTHS = new int[] { 9, 10 };

    private ChopTAP(EnsemblDb db) {
        this.db = db;
    }

    public static void run(EnsemblDb db) {
        ChopTAP chopTAP = new ChopTAP(db);
        chopTAP.run();
    }

    public static void run(File fastaFile) {
        run(EnsemblDb.create(fastaFile));
    }

    public static void run(String fastaFile) {
        run(new File(fastaFile));
    }

    private void run() {
        writeHeader();
        processSymbols();
    }

    private void writeHeader() {
        System.out.println("Hugo_Symbol\tPeptide");
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
            System.out.println(symbol.getKey() + "\t" + peptide);
    }

    private Collection<String> generatePeptides(HugoSymbol symbol) {
        Collection<String> peptides = new TreeSet<String>();
        Collection<EnsemblRecord> records = db.get(symbol);

        for (EnsemblRecord record : records)
            peptides.addAll(generatePeptides(record));

        return peptides;
    }

    private Collection<String> generatePeptides(EnsemblRecord record) {
        Peptide peptide = record.getPeptide();

        if (!peptide.isNative())
            return Collections.emptyList();

        Collection<Peptide> chopped = NetChop.chop(peptide, PEPTIDE_LENGTHS);
        Collection<Peptide> transported = TAP.INSTANCE.transport(chopped);

        return Peptide.formatString(transported);
    }

    private static void usage() {
        System.err.println("Usage: java [JVMOPTIONS] pepmhc.app.ChopTAP <ENSEMBL_FILE>");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 1)
            usage();

        run(args[0]);
    }
}
