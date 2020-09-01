
package pepmhc.tap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.io.FileUtil;
import jam.math.DoubleRange;

import jene.peptide.Peptide;
import jene.peptide.Residue;

/**
 * Predicts peptides that will be transported into the endoplasmic
 * reticulum by the transporter associated with antigen processing
 * (TAP).  The method was developed by Peters et al. 
 * <em>J. Immunol.</em> <b>171</b>, 1741-1749 (2003).
 */
public final class TAP {
    private final double alpha;
    private final double threshold;
    private final TAPMatrix matrix;

    private static TAP consensus = null;

    private TAP(TAPMatrix matrix, double alpha, double threshold) {
        this.alpha = alpha;
        this.matrix = matrix;
        this.threshold = threshold;

        ALPHA_RANGE.validate("Alpha factor", alpha);
        THRESHOLD_SCORE_RANGE.validate("Score threshold", threshold);
    }

    /**
     * Name of the system property that defines the shrinkage factor
     * for the contribution of N-terminal residues.
     */
    public static final String ALPHA_PROPERTY = "pepmhc.tap.alpha";

    /**
     * The optimal shrinkage factor for the contribution of N-terminal
     * residues to the binding score as reported by Peters.
     */
    public static final double ALPHA_DEFAULT = 0.2;

    /**
     * Valid range of alpha shrinkage factors.
     */
    public static final DoubleRange ALPHA_RANGE = DoubleRange.FRACTIONAL;

    /**
     * Name of the system property that defines the threshold score
     * for transport. Peptides must score BELOW this threshold to be
     * transported.
     */
    public static final String THRESHOLD_SCORE_PROPERTY = "pepmhc.tap.thresholdScore";

    /**
     * Default value for the threshold score for transport.  Peters
     * reports that 98% of known epitopes fall below 1.0, while the
     * value -1.0 provides a reasonable trade-off between the false
     * positive and false negative rate.
     */
    public static final double THRESHOLD_SCORE_DEFAULT = 1.0;

    /**
     * Valid range of score thresholds.
     */
    public static final DoubleRange THRESHOLD_SCORE_RANGE = DoubleRange.closed(-2.0, 2.0);

    /**
     * Creates a new TAP scorer.
     *
     * @param alpha the alpha shrinkage factor for N-terminal
     * residues.
     *
     * @param threshold the score threshold for transport.
     */
    public TAP(double alpha, double threshold) {
        this(TAPMatrix.consensus(), alpha, threshold);
    }

    /**
     * Returns the consensus TAP scorer.
     *
     * @return the consensus TAP scorer.
     */
    public static TAP consensus() {
        if (consensus == null)
            consensus = createConsensus();

        return consensus;
    }

    private static TAP createConsensus() {
        return new TAP(TAPMatrix.consensus(), resolveAlpha(), resolveThreshold());
    }

    private static double resolveAlpha() {
        return JamProperties.getOptionalDouble(ALPHA_PROPERTY,
                                               ALPHA_RANGE,
                                               ALPHA_DEFAULT);
    }

    private static double resolveThreshold() {
        return JamProperties.getOptionalDouble(THRESHOLD_SCORE_PROPERTY,
                                               THRESHOLD_SCORE_RANGE,
                                               THRESHOLD_SCORE_DEFAULT);
    }

    /**
     * Returns the alpha shrinkage factor for N-terminal residues.
     *
     * @return the alpha shrinkage factor for N-terminal residues.
     */
    public double getAlpha() {
        return alpha;
    }

    /**
     * Returns the score threshold required for successful TAP
     * transport.
     *
     * @return the score threshold required for successful TAP
     * transport.
     */
    public double getThreshold() {
        return threshold;
    }

    /**
     * Returns the underlying scoring matrix.
     *
     * @return the underlying scoring matrix.
     */
    public TAPMatrix getMatrix() {
        return matrix;
    }

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
        return peptide.isNative() && score(peptide) <= threshold;
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
        int L = peptide.length();

        if (L < 9)
            throw new IllegalArgumentException("Peptide must have at least nine residues.");

        // Start with the contribution from the C terminus...
        int    cterm = L - 1;
        double score = matrix.get(peptide.get(cterm), TAPPosition.CTerm);

        // Add the alpha-scaled contributions from the N-terminal residues...
        score += alpha * matrix.get(peptide.get(0), TAPPosition.NTerm1);
        score += alpha * matrix.get(peptide.get(1), TAPPosition.NTerm2);
        score += alpha * matrix.get(peptide.get(2), TAPPosition.NTerm3);

        return score;
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
