
package pepmhc.stab;

import java.util.ArrayList;
import java.util.List;

import jam.app.JamLogger;
import jam.hla.Allele;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.peptide.Peptide;
import jam.util.ListUtil;

/**
 * Populates the stability cache with a collection of alleles and
 * peptides.
 */
public final class StabilityBank {
    private final String alleleFile;
    private final String peptideFile;

    private final List<Allele> alleles = new ArrayList<Allele>();
    private final List<Peptide> peptides = new ArrayList<Peptide>();

    private final static int BATCH_SIZE = 200000;

    private StabilityBank(String[] args) {
        validate(args);

        this.alleleFile  = args[0];
        this.peptideFile = args[1];
    }

    private static void validate(String[] args) {
        if (args.length != 2) {
            System.err.println("Usage: jam.stab.StabilityBank ALLELE_FILE PEPTIDE_FILE");
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
        // Just get the results from the stability cache, to enforce
        // calculation on demand, but ignore the returned records...
        //
        try {
            List<List<Peptide>> subLists = ListUtil.split(peptides, BATCH_SIZE);

            for (List<Peptide> subList : subLists)
                StabilityCache.get(allele, subList);
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
