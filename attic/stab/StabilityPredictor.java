
package pepmhc.stab;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import jam.app.JamEnv;
import jam.hla.Allele;
import jam.io.LineReader;
import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class StabilityPredictor {
    private final Allele allele;
    private final Collection<Peptide> peptides;
    private final List<StabilityRecord> stabRecords;

    private List<Peptide> pep8  = new ArrayList<Peptide>();
    private List<Peptide> pep9  = new ArrayList<Peptide>();
    private List<Peptide> pep10 = new ArrayList<Peptide>();
    private List<Peptide> pep11 = new ArrayList<Peptide>();

    private List<StabilityRecord> stab8;
    private List<StabilityRecord> stab9;
    private List<StabilityRecord> stab10;
    private List<StabilityRecord> stab11;

    private static Set<Allele> coverage = null;

    private static final String COVERAGE_FILE = "data/allele/coverage_netMHCstabpan.txt";

    private static final PredictionMethod METHOD = PredictionMethod.NET_MHC_STAB_PAN;

    private StabilityPredictor(Allele allele, Collection<Peptide> peptides) {
        this.allele = allele;
        this.peptides = peptides;
        this.stabRecords = new ArrayList<StabilityRecord>(peptides.size());
    }

    /**
     * Identifies alleles covered by this stability predictor.
     *
     * @param allele an allele of interest.
     *
     * @return {@code true} iff the specified allele is covered by
     * this stability predictor.
     */
    public static boolean isCovered(Allele allele) {
        return viewCoverage().contains(allele);
    }

    /**
     * Returns the set of alleles covered by this stability predictor.
     *
     * @return the set of alleles covered by this stability predictor.
     */
    public static Set<Allele> viewCoverage() {
        if (coverage == null)
            coverage = loadCoverage();

        return coverage;
    }

    private static Set<Allele> loadCoverage() {
        Set<Allele> alleles = new HashSet<Allele>();
        LineReader  reader  = LineReader.open(resolveCoverageFile());

        for (String allele : reader)
            alleles.add(Allele.instance(allele));

        return Collections.unmodifiableSet(alleles);
    }

    private static File resolveCoverageFile() {
        return new File(JamEnv.getRequired("PEPMHC_HOME"), COVERAGE_FILE);
    }

    /**
     * Determines whether the underlying executable program or
     * prediction engine is available to the JVM.
     *
     * @return {@code true} iff the underlying executable program
     * or prediction engine is available to the JVM.
     */
    public static boolean isInstalled() {
        return Predictor.isInstalled(METHOD);
    }

    /**
     * Predicts the stability of a peptide-MHC complex for a given
     * allele.
     *
     * @param allele the code of the MHC allele presenting the
     * peptide.
     *
     * @param peptide the peptide being presented.
     *
     * @return the stability record for the allele and peptide.
     */
    public static StabilityRecord predict(Allele allele, Peptide peptide) {
        return StabilityRecord.convert(Predictor.instance(METHOD).predict(allele, peptide));
    }

    /**
     * Predicts the stability of peptide-MHC complexes for a given
     * allele.
     *
     * @param allele the code of the MHC allele presenting the
     * peptides.
     *
     * @param peptides the peptides being presented.
     *
     * @return a list of stability records for the allele and
     * peptides.
     */
    public static List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        StabilityPredictor predictor =
            new StabilityPredictor(allele, peptides);

        return predictor.predict();
    }

    private List<StabilityRecord> predict() {
        //
        // Peptides passed to the netMHCstabpan executable must have
        // the same length, so we need to separate the peptide input
        // by length and then recompose the results in the original
        // peptide order...
        //
        mapPeptides();
        predictUniform();
        reduceStabRecords();

        return stabRecords;
    }

    private void mapPeptides() {
        for (Peptide peptide : peptides) {
            switch (peptide.length()) {
            case 8:
                pep8.add(peptide);
                break;

            case 9:
                pep9.add(peptide);
                break;

            case 10:
                pep10.add(peptide);
                break;

            case 11:
                pep11.add(peptide);
                break;

            default:
                throw new IllegalArgumentException("Invalid peptide length.");
            }
        }
    }

    private void predictUniform() {
        stab8  = predictUniform(pep8);
        stab9  = predictUniform(pep9);
        stab10 = predictUniform(pep10);
        stab11 = predictUniform(pep11);
    }

    private List<StabilityRecord> predictUniform(List<Peptide> uniform) {
        if (uniform.isEmpty())
            return Collections.emptyList();
        else
            return StabilityBatchProcess.run(allele, uniform);
    }

    private void reduceStabRecords() {
        Iterator<StabilityRecord> iter8  = stab8.iterator();
        Iterator<StabilityRecord> iter9  = stab9.iterator();
        Iterator<StabilityRecord> iter10 = stab10.iterator();
        Iterator<StabilityRecord> iter11 = stab11.iterator();

        for (Peptide peptide : peptides) {
            switch (peptide.length()) {
            case 8:
                stabRecords.add(iter8.next());
                break;

            case 9:
                stabRecords.add(iter9.next());
                break;

            case 10:
                stabRecords.add(iter10.next());
                break;

            case 11:
                stabRecords.add(iter11.next());
                break;

            default:
                throw new IllegalArgumentException("Invalid peptide length.");
            }
        }

        // Verify that all uniform-length results have been processed...
        assert !iter8.hasNext();
        assert !iter9.hasNext();
        assert !iter10.hasNext();
        assert !iter11.hasNext();
        assert stabRecords.size() == peptides.size();
    }
}
