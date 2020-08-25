
package pepmhc.app;

import java.util.ArrayList;
import java.util.List;

import jam.app.JamLogger;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.util.ListUtil;

import jean.hla.Allele;
import jean.hugo.HugoPeptideTable;
import jean.peptide.Peptide;

import pepmhc.affy.AffinityCache;
import pepmhc.affy.AffinityMethod;

public final class AffinityBank {
    private final String alleleFile;
    private final String peptideFile;
    private final AffinityMethod predMethod;

    private final List<Allele> alleles = new ArrayList<Allele>();
    private final List<Peptide> peptides = new ArrayList<Peptide>();

    private final static int BATCH_SIZE = 100000;

    private AffinityBank(String[] args) {
        validate(args);

        this.alleleFile  = args[0];
        this.peptideFile = args[1];
        this.predMethod  = AffinityMethod.valueOf(args[2]);
    }

    private static void validate(String[] args) {
        if (args.length != 3) {
            System.err.println("Usage: pepmhc.app.AffinityBank ALLELE_FILE PEPTIDE_FILE PREDICTION_METHOD");
            System.exit(1);
        }
    }

    private void run() {
        loadAlleles();
        loadPeptides();
        processAlleles();
        JamLogger.info("DONE!");
    }

    private void loadAlleles() {
        LineReader reader = LineReader.open(alleleFile);

        for (String line : reader)
            alleles.add(Allele.instance(line));

        reader.close();
    }

    private void loadPeptides() {
        HugoPeptideTable table = HugoPeptideTable.load(peptideFile);
        peptides.addAll(table.viewPeptides());
    }

    private void processAlleles() {
        for (Allele allele : alleles)
            processAllele(allele);
    }

    private void processAllele(Allele allele) {
        //
        // Just get the results from the affinity cache, to enforce
        // calculation on demand, but ignore the returned records...
        //
        List<List<Peptide>> subLists = ListUtil.split(peptides, BATCH_SIZE);

        for (List<Peptide> subList : subLists)
            AffinityCache.instance(predMethod, allele).require(subList);
    }

    public static void main(String[] args) {
        AffinityBank bank = new AffinityBank(args);
        bank.run();
    }
}
