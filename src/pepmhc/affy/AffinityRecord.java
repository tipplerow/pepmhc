
package pepmhc.affy;

import jam.math.DoubleRange;

import jean.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC affinity measurement or
 * prediction.
 */
public final class AffinityRecord extends BindRecord {
    private final double affinity;
    private final double percentile;

    /**
     * Creates a new affinity record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     */
    public AffinityRecord(Peptide peptide, double affinity) {
        this(peptide, affinity, Double.NaN);
    }

    /**
     * Creates a new affinity record.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     *
     * @param percentile the percentile rank of the binding affinity
     * (relative to other peptides of the same length binding to the
     * same allele).
     */
    public AffinityRecord(Peptide peptide, double affinity, double percentile) {
        super(peptide);

        this.affinity = affinity;
        this.percentile = percentile;

        validate();
    }

    private void validate() {
        DoubleRange.NON_NEGATIVE.validate("binding affinity", affinity);
        DoubleRange.PERCENTILE.validate("percentile rank", percentile);
    }

    /**
     * Returns the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     *
     * @return the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     */
    public double getAffinity() {
        return affinity;
    }

    /**
     * Returns the percentile rank of the affinity (relative to other
     * peptides of the same length binding to the same allele).
     *
     * @return the percentile rank of the affinity.
     */
    public double getPercentile() {
        return percentile;
    }

    /**
     * Identifies records with a valid affinity.
     *
     * @return {@code false} iff the affinity is {@code Double.NaN}.
     */
    public boolean isAffinitySet() {
        return !Double.isNaN(affinity);
    }

    /**
     * Identifies records with a valid percentile.
     *
     * @return {@code false} iff the percentile is {@code Double.NaN}.
     */
    public boolean isPercentileSet() {
        return !Double.isNaN(percentile);
    }

    @Override public String toString() {
        return String.format("AffinityRecord(%s, %.1f, %.2f)", peptide.formatString(), affinity, percentile);
    }
}
