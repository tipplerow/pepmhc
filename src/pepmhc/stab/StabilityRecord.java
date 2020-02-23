
package pepmhc.stab;

import jam.math.DoubleRange;
import jam.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC stability measurement or
 * prediction.
 */
public final class StabilityRecord extends BindRecord {
    private final double halfLife;
    private final double percentile;

    /**
     * Creates a new stability record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     */
    public StabilityRecord(Peptide peptide, double halfLife) {
        this(peptide, halfLife, Double.NaN);
    }

    /**
     * Creates a new stability record.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     *
     * @param percentile the percentile rank of the half-life
     * (relative to other peptides of the same length bound to
     * the same allele).
     */
    public StabilityRecord(Peptide peptide, double halfLife, double percentile) {
        super(peptide);

        this.halfLife = halfLife;
        this.percentile = percentile;

        validate();
    }

    private void validate() {
        DoubleRange.NON_NEGATIVE.validate("half-life", halfLife);
        DoubleRange.PERCENTILE.validate("percentile rank", percentile);
    }

    /**
     * Returns the half-life of the MHC-bound state (hours).
     *
     * @return the half-life of the MHC-bound state (hours).
     */
    public double getHalfLife() {
        return halfLife;
    }

    /**
     * Returns the percentile rank of the half-life (relative to other
     * peptides of the same length bound to the same allele).
     *
     * @return the percentile rank of the half-life (relative to other
     * peptides of the same length bound to the same allele).
     */
    public double getPercentile() {
        return percentile;
    }

    @Override public String toString() {
        return String.format("StabilityRecord(%s, %.2f, %.2f)", peptide.formatString(), halfLife, percentile);
    }
}
