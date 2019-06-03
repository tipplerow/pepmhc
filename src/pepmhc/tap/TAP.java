
package pepmhc.tap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.io.FileUtil;
import jam.peptide.Peptide;

import pepmhc.engine.smm.StabilizedMatrix;

/**
 * Predicts peptides that will be transported into the endoplasmic
 * reticulum by the transporter associated with antigen processing
 * (TAP) using a method reported by Peters et al., J. Immunol. 171,
 * 1741--1749 (2003).
 */
public final class TAP {
    private final double threshold;
    private final StabilizedMatrix matrix;

    private TAP() {
        this.matrix = StabilizedMatrix.load(resolveMatrixFile());
        this.threshold = resolveThreshold();
    }

    private static String resolveMatrixFile() {
        return FileUtil.join(JamEnv.getRequired("PEPMHC_HOME"), "data", "tap", "consensus.txt");
    }

    /**
     * The optimal shrinkage factor for the contribution of N-terminal
     * residues to the binding score as reported by Peters.
     */
    public static final double SCORE_ALPHA = 0.2;

    /**
     * Name of the system property that defines the threshold score
     * for transport: peptides must score BELOW this threshold to be
     * transported.
     */
    public static final String THRESHOLD_SCORE_PROPERTY = "pepmhc.tap.thresholdScore";

    /**
     * Default value for the threshold score for transport: 
     * Peters reports that 98% of known epitopes fall below this
     * threshold.
     */
    public static final double THRESHOLD_SCORE_DEFAULT = 1.0;

    private static double resolveThreshold() {
        return JamProperties.getOptionalDouble(THRESHOLD_SCORE_PROPERTY,
                                               THRESHOLD_SCORE_DEFAULT);
    }

    /**
     * The single global instance.
     */
    public static final TAP INSTANCE = new TAP();

    /**
     * Identifies peptides that will be transported by TAP.
     *
     * @param peptide the peptide of interest.
     *
     * @return {@code true} iff the input peptide will be transported
     * by TAP.
     *
     * @throws IllegalArgumentException unless the peptide has length
     * nine or greater.
     */
    public boolean isTransported(Peptide peptide) {
        return score(peptide) <= threshold;
    }

    /**
     * Computes the TAP binding score for a peptide.
     *
     * @param peptide the peptide being processed.
     *
     * @return the TAP binding score for the given peptide.
     *
     * @throws IllegalArgumentException unless the peptide has length
     * nine or greater.
     */
    public double score(Peptide peptide) {
        //
        // This code implements Equation 3 from Peters et al.,
        // J. Immunol. 171, 1741--1749 (2003) with L = 9 and
        // alpha = 0.2.
        //
        if (peptide.length() < 9)
            throw new IllegalArgumentException("Peptide must have at least nine residues.");

        // Start with the contribution from the C terminus...
        double result = matrix.getElement(peptide.at(peptide.length() - 1), 8);

        // Add the rescaled contributions from the first three
        // N-terminal residues...
        for (int nterm = 0; nterm < 3; ++nterm)
            result += SCORE_ALPHA * matrix.getElement(peptide.at(nterm), nterm);

        return result;
    }

    /**
     * Identifies peptides that will be transported by TAP.
     *
     * @param peptides the peptide of interests.
     *
     * @return a list of peptides that are transported by TAP.
     *
     * @throws IllegalArgumentException if any peptide has length less
     * than nine.
     */
    public List<Peptide> transport(Collection<Peptide> peptides) {
        List<Peptide> transported = new ArrayList<Peptide>(peptides.size());

        for (Peptide peptide : peptides)
            if (isTransported(peptide))
                transported.add(peptide);

        return transported;
    }
}
