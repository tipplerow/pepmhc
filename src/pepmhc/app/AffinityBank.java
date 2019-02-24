
package pepmhc.app;

import java.util.ArrayList;
import java.util.List;

import jam.hla.Allele;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.Peptide;

import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;

public final class AffinityBank {
    private final String alleleFile;
    private final String peptideFile;
    private final PredictionMethod predMethod;

    private final List<Allele> alleles = new ArrayList<Allele>();
    private final List<Peptide> peptides = new ArrayList<Peptide>();

    private AffinityBank(String[] args) {
        validate(args);

        this.alleleFile  = args[0];
        this.peptideFile = args[1];
        this.predMethod  = PredictionMethod.valueOf(args[2]);
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
    }

    private void loadAlleles() {
        LineReader reader = LineReader.open(alleleFile);

        for (String line : reader)
            alleles.add(Allele.instance(line));

        reader.close();
    }

    private void loadPeptides() {
        LineReader reader = LineReader.open(peptideFile);

        for (String line : reader)
            peptides.add(Peptide.parse(line));

        reader.close();
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
        AffinityCache.get(predMethod, allele, peptides);
    }

    public static void main(String[] args) {
        AffinityBank bank = new AffinityBank(args);
        bank.run();
    }
}
