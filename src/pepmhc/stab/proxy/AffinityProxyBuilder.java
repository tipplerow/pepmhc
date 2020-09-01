
package pepmhc.stab.proxy;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.io.LineReader;
import jam.math.JamRandom;
import jam.util.ListUtil;

import jene.hla.Allele;
import jene.hugo.HugoPeptideTable;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityRecord;
import pepmhc.affy.AffinityStore;

import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityRecord;
import pepmhc.stab.StabilityStore;

/**
 * Estimates the parameters in affinity-proxy models.
 */
public final class AffinityProxyBuilder {
    private final Allele allele;

    private final AffinityMethod affinityMethod;
    private final StabilityMethod stabilityMethod;

    private final AffinityStore affinityStore;
    private final StabilityStore stabilityStore;

    private final int sampleSize;
    private final String peptideFile;

    private final SimpleRegression regression;
    
    private static List<Peptide> peptideList = null;

    private static final long RANDOM_SEED = 20200131;

    private AffinityProxyBuilder(Allele allele, AffinityMethod affinityMethod) {
        this.allele = allele;

        this.affinityMethod = affinityMethod;
        this.stabilityMethod = StabilityMethod.NET_MHC_STAB_PAN;

        this.affinityStore = AffinityStore.instance(affinityMethod, allele);
        this.stabilityStore = StabilityStore.instance(stabilityMethod, allele);

        this.peptideFile = JamProperties.getRequired(PROXY_PEPTIDE_FILE_PROPERTY);
        this.sampleSize  = JamProperties.getRequiredInt(PROXY_SAMPLE_SIZE_PROPERTY);

        this.regression = new SimpleRegression();
        assert regression.hasIntercept();
    }

    /**
     * System property naming the file that contains a HUGO-peptide
     * table that will provide peptide fragments for model building.
     */
    public static final String PROXY_PEPTIDE_FILE_PROPERTY = "pepmhc.stab.proxyPeptideFile";

    /**
     * System property defining the maximum number of peptides to use
     * in model building.
     */
    public static final String PROXY_SAMPLE_SIZE_PROPERTY = "pepmhc.stab.proxySampleSize";

    /**
     * Builds and stores affinity-proxy models for HLA alleles
     * contained in a file.
     *
     * @param fileName the name of a file containing the HLA alleles
     * of interest.
     *
     * @param affinityMethod the affinity prediction method to use.
     */
    public static void build(String fileName, AffinityMethod affinityMethod) {
        /*
        LineReader reader = LineReader.open(fileName);

        try {
            for (String line : reader) {
                Allele allele = Allele.instance(line);
                JamLogger.info("Processing [%s]...", allele);

                AffinityProxyModel model = build(allele, affinityMethod);
                AffinityProxyDb.instance().add(model);
            }

            JamLogger.info("DONE!");
        }
        finally {
            reader.close();
        }
        */
    }

    /**
     * Builds an affinity-proxy model for a given HLA allele and
     * affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @param affinityMethod the affinity prediction method to use.
     *
     * @return an affinity-proxy model with parameters determined for
     * the the specified HLA allele and affinity prediction method.
     */
    public static AffinityProxyModel build(Allele allele, AffinityMethod affinityMethod) {
        AffinityProxyBuilder builder = new AffinityProxyBuilder(allele, affinityMethod);
        return builder.build();
    }

    private AffinityProxyModel build() {
        selectPeptides();
        addSampleData();

        double intercept = regression.getIntercept();
        double coefficient = regression.getSlope();

        return new AffinityProxyModel(allele, affinityMethod, intercept, coefficient);
    }

    private void selectPeptides() {
        if (peptideList != null)
            return;

        HugoPeptideTable peptideTable = HugoPeptideTable.load(peptideFile);
        peptideList = new ArrayList<Peptide>(peptideTable.viewPeptides());

        if (sampleSize < peptideList.size()) {
            JamRandom generator = JamRandom.generator(RANDOM_SEED);
            ListUtil.shuffle(peptideList, generator);
            peptideList = peptideList.subList(0, sampleSize);
        }
    }

    private void addSampleData() {
        List<AffinityRecord> affinityRecords = affinityStore.get(peptideList);
        List<StabilityRecord> stabilityRecords = stabilityStore.get(peptideList);

        for (int k = 0; k < peptideList.size(); ++k) {
            //
            // The linear regression model is: log(H) = intercept + coeff * log(A)
            //
            double affinity = affinityRecords.get(k).getAffinity();
            double halfLife = stabilityRecords.get(k).getHalfLife();

            if (affinity <= 1.0E-08)
                continue;

            if (halfLife <= 1.0E-08)
                continue;

            double logA = Math.log(affinity);
            double logH = Math.log(halfLife);

            regression.addData(logA, logH);
        }
    }

    private static void usage() {
        System.err.print("Usage: pepmhc.stab.AffinityProxyBuilder ");
        System.err.print("ALLELE_FILE AFFINITY_METHOD PEPTIDE_FILE SAMPLE_SIZE");
        System.err.println();
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length != 4)
            usage();
            
        String alleleFile = args[0];
        AffinityMethod affinityMethod = AffinityMethod.valueOf(args[1]);

        System.setProperty(PROXY_PEPTIDE_FILE_PROPERTY, args[2]);
        System.setProperty(PROXY_SAMPLE_SIZE_PROPERTY,  args[3]);

        build(alleleFile, affinityMethod);
    }
}
