
package pepmhc.app;

import java.util.List;

import jam.app.JamLogger;
import jam.util.ListUtil;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityStore;

/**
 * Populates the stability store with a collection of alleles and
 * peptides.
 */
public final class StabilityBank {
    private final String alleleFile;
    private final String peptideFile;

    private final StabilityMethod method;

    private List<Allele> alleles;
    private List<Peptide> peptides;

    private final static int BATCH_SIZE = 200000;

    private StabilityBank(String[] args) {
        validate(args);

        this.method = StabilityMethod.valueOf(args[0]);
        this.alleleFile = args[1];
        this.peptideFile = args[2];
    }

    private static void validate(String[] args) {
        if (args.length != 2) {
            System.err.println("Usage: pepmhc.app.StabilityBank STABILITY_METHOD ALLELE_FILE PEPTIDE_FILE");
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
        alleles = Allele.load(alleleFile);
    }

    private void loadPeptides() {
        peptides = Peptide.load(peptideFile);
    }

    private void processAlleles() {
        for (Allele allele : alleles)
            processAllele(allele);
    }

    private void processAllele(Allele allele) {
        //
        // Just get the results from the stability store, to enforce
        // calculation on demand, but ignore the returned records...
        //
        try {
            List<List<Peptide>> subLists = ListUtil.split(peptides, BATCH_SIZE);

            for (List<Peptide> subList : subLists)
                StabilityStore.instance(method, allele).get(subList);
        }
        catch (Exception ex) {
            JamLogger.error("Stability calculation failed for allele [%s].", allele);
            JamLogger.error(ex.getMessage());
        }
    }

    public static void main(String[] args) {
        StabilityBank bank = new StabilityBank(args);
        bank.run();
    }
}
